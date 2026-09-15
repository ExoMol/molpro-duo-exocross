import copy
import json
import math
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

from pecfit.__main__ import main
from pecfit.duo import Runner, assess_energies, converged_checkpoint, parse_en, parse_output, write_json
from pecfit.model import FitError, HARTREE_TO_CM, Observation, convert_pec, emo, import_template, number, observations, render, shape_check
from pecfit.pipeline import DEFAULTS, Pipeline, load, merge_config, pec_diagnostics

ROOT = Path(__file__).resolve().parents[1]
CO = ROOT/"CO/sample_outputs/CO_duo_fit_01.inp"


class ModelTests(unittest.TestCase):
    def setUp(self):
        self.model = import_template(CO)

    def test_sample_coverage_and_basis_convention(self):
        m = self.model
        self.assertEqual(len(m.observations), 2293)
        self.assertEqual(max(o.j for o in m.observations), 123)
        self.assertEqual(max(o.v for o in m.observations), 41)
        self.assertEqual(len(m.pec), 301)

    def test_round_trip_does_not_modify_data(self):
        m = self.model
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)/"input.inp"
            p.write_text(render(m,m.parameters,jmax=123,vmax=80,iterations=0))
            other = import_template(p)
        self.assertEqual(m.parameters, other.parameters)
        self.assertEqual(len(m.observations), len(other.observations))
        self.assertAlmostEqual(m.observations[-1].energy, other.observations[-1].energy)
        self.assertEqual(m.pec, other.pec)

    def test_hartree_conversion_and_minimum_shift(self):
        result = convert_pec([(1.,-2.), (2.,-3.), (3.,-1.)], "hartree", True)
        self.assertEqual(result[1][1], 0)
        self.assertEqual(result[2][1], 2*HARTREE_TO_CM)

    def test_already_shifted_grid_is_not_shifted_again(self):
        rows = [(1.,12.),(2.,3.),(3.,20.)]
        self.assertEqual(convert_pec(rows,"cm-1",False), rows)

    def test_invalid_grid_is_rejected(self):
        for rows in ([(1.,2.),(1.,3.),(2.,4.)],[(2.,3.),(1.,2.),(3.,5.)]):
            with self.assertRaises(FitError):
                convert_pec(rows,"cm-1",False)

    def test_emo_uses_absolute_asymptote(self):
        p = {**self.model.parameters, "V0":123.,"DE":1123.,"RE":1.,"B0":2.}
        self.assertEqual(emo(1.,p),123.)
        self.assertAlmostEqual(emo(2.,p),123+1000*(1-math.exp(-2))**2)

    def test_emo_branch_changes_at_rref(self):
        p = {**self.model.parameters, "RREF":1.4,"RE":1.1,"NL":0.,"NR":1.,"B0":2.,"B1":3.}
        for r, beta in ((1.3,2.),(1.5,2+3*(1.5**6-1.4**6)/(1.5**6+1.4**6))):
            self.assertAlmostEqual(emo(r,p),p["DE"]*(1-math.exp(-beta*(r-1.1)))**2)

    def test_shape_guard_rejects_extra_well(self):
        p = {**self.model.parameters,"RE":1.1,"B0":2.,"B1":-10.}
        self.assertFalse(shape_check(p,self.model.grid)["ok"])

    def test_bad_assignments_fail_before_duo(self):
        valid = ["0 e 1 0 X 0 0 0 0 1", "1 e 1 4 X 0 0 0 0 1"]
        for bad in ("1 f 1 4 X 0 0 0 0 1", "1 e 2 4 X 0 0 0 0 1", "1 e 1 4 X 0 1 0 0 1", "1 e 1 4 X 0 0 0 0 -1"):
            with self.assertRaises(FitError):
                observations([valid[0],bad])
        with self.assertRaises(FitError):
            observations(valid+[valid[-1]])
        self.assertEqual(observations([valid[0],"1 - 1 4 X 0 0 0 0 1"])[1].parity,"e")

    def test_nan_is_rejected(self):
        with self.assertRaises(FitError):
            observations(["0 e 1 0 X 0 0 0 0 1", "1 e 1 NaN X 0 0 0 0 1"])

    def test_fortran_exponent_and_overflow_fields(self):
        self.assertEqual(number("1.234-100"),1.234e-100)
        with self.assertRaises(FitError):
            number("***********")

    def test_abinitio_stage_contains_only_computable_levels(self):
        text = render(self.model,self.model.parameters,jmax=0,vmax=1,iterations=4)
        self.assertNotIn("values\n\n",text.lower())
        self.assertEqual(len(text.split("energies\n")[1].strip().splitlines()),2)


EN = """Iteration = 0
    1    1      0.0 +           0.0000        0.0000        0.0000   0.96E-17  ( X 0 0 0.0 0.0 )( X 0 0 0.0 0.0 )
    2    2      0.0 +        2143.2711     2143.1711        0.1000   0.00E+00  ( X 1 0 0.0 0.0 )( X 1 0 0.0 0.0 )
"""


class OutputTests(unittest.TestCase):
    def parse(self, text=EN):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"fit.en"
            p.write_text(text)
            return parse_en(p)

    def test_last_iteration_only(self):
        first=EN.replace("2143.1711","2100.0000")
        result=self.parse(first+EN.replace("Iteration = 0","Iteration = 2"))
        self.assertEqual(result[(0,1)]["calculated"],2143.1711)

    def test_all_observations_count_even_if_duo_weight_is_zero(self):
        obs=[Observation(0,"e",1,0.,"X",0,1.),Observation(0,"e",2,2143.271107,"X",1,4.)]
        stats, rows=assess_energies(self.parse(),obs,0)
        self.assertEqual(stats["count"],2)
        self.assertEqual(stats["excluded_by_duo_count"],1)
        self.assertAlmostEqual(stats["weighted_rms_cm"],math.sqrt(4*(.100007**2)/5))
        self.assertAlmostEqual(rows[1]["residual_cm"],.100007)

    def test_missing_and_misassigned_levels_fail(self):
        obs=[Observation(0,"e",1,0.,"X",0,1.),Observation(0,"e",2,2143.271107,"X",1,4.)]
        for text in (EN.splitlines()[0]+"\n"+EN.splitlines()[1],EN.replace(")( X 1", ")( X 2")):
            with self.assertRaises(FitError):
                assess_energies(self.parse(text),obs,0)

    def test_ignore_echo_and_rounded_parameters(self):
        output="""POTEN X
Values
RE 999
end
Iteration = 1

Robust fitting ...
Watson parameter = 0.00100000
... done!

Parameters:

POTEN X X
name X1Sigma+
type EMO
Values
RE 1.12832173137337D+00 fit
DE 9.27101841221794E+04 fit
end

Fitted parameters (rounded):
POTEN X X
Values
RE 1.0 1.13
end
-->| 1 | 42 | 2 | 1.2E-4 | 0.14 | 10.0 | 1.0E-8 |
"""
        iterations,history=parse_output(output)
        self.assertEqual(iterations[-1]["parameters"]["RE"],1.12832173137337)
        self.assertEqual(history[-1]["stability"],1e-8)


class ConfigTests(unittest.TestCase):
    def test_unknown_options_and_bad_ranges_fail(self):
        for config in ({"typo":1},{"experiment":{"typo":1}}, {"vmax":0},
                       {"experiment":{"fit_scale":2}}, {"jmax_schedule":[1.5]},
                       {"experiment":{"robust":"0"}}, {"experimental_de_cm":-1},
                       {"pec_shift_minimum":"false"}, {"vmax":True},
                       {"validation":{"range_extension_angstrom":[0]}},
                       {"experiment":{"target_metric":"duo_weighted_rms"}}):
            with self.assertRaises(FitError):
                merge_config(config)

    def test_schedule_always_includes_full_data(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"config.json"
            write_json(path,{"template":str(CO),"jmax_schedule":[0,10]})
            self.assertEqual(load(path)[1]["jmax_schedule"],[0,10,123])

    def test_resume_rejects_changed_configuration(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"config.json"
            write_json(path,{"template":str(CO)})
            output=Path(d)/"run"
            Pipeline(path,executable=os.sys.executable,output=output)
            Pipeline(path,executable=os.sys.executable,output=output,resume=True)
            write_json(path,{"template":str(CO),"vmax":90})
            with self.assertRaises(FitError):
                Pipeline(path,executable=os.sys.executable,output=output,resume=True)

    def test_existing_output_not_overwritten(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"config.json"
            write_json(path,{"template":str(CO)})
            output=Path(d)/"run"
            output.mkdir()
            (output/"precious.txt").write_text("preserve")
            with self.assertRaises(FitError):
                Pipeline(path,executable=os.sys.executable,output=output)
            self.assertEqual((output/"precious.txt").read_text(),"preserve")


class RecoveryTests(unittest.TestCase):
    def setUp(self):
        self.pipeline = object.__new__(Pipeline)
        p = self.pipeline
        p.config = merge_config({"experiment":{"max_restarts":1,"max_retries":1,"target_metric":"weighted_rms_cm"}})
        p.p = {"RE": 1.1}
        p.root = Path("run")
        p.physical = lambda _: {"ok": True}
        p.record = lambda **kwargs: None
        self.stats = {"weighted_rms_cm":1.,"rms_cm":2.,"max_abs_cm":20.,"count":41}

    def test_worsening_candidate_rolls_back_and_reduces_step(self):
        p = self.pipeline
        scales=[]
        def run(*args,**kwargs):
            scales.append(kwargs["scale"])
            return SimpleNamespace(parameters={"RE":1.2},history=[{"stability":1.0}])
        p.run_duo=run
        first=True
        def evaluate(*args,**kwargs):
            nonlocal first
            stats=self.stats if first else {**self.stats,"weighted_rms_cm":10.}
            first=False
            return stats,None,SimpleNamespace(directory=Path("run/check"))
        p.evaluate=evaluate
        p.fit_experiment(10,"trial")
        self.assertEqual(p.p,{"RE":1.1})
        self.assertEqual(scales,[.1,.05])

    def test_line_search_recovers_a_physical_improving_step(self):
        p=self.pipeline
        p.config["experiment"]["line_search_steps"]=4
        p.physical=lambda q: None if q["RE"]<=1.16 else (_ for _ in ()).throw(FitError("bad shape"))
        p.evaluate=lambda label,q,j:(dict(self.stats,weighted_rms_cm=.5),None,SimpleNamespace(directory=Path("check")))
        candidate,stats,_,fraction=p.validated_step({"RE":1.2},self.stats,10,"trial")
        self.assertAlmostEqual(candidate["RE"],1.15)
        self.assertEqual(fraction,.5)
        self.assertEqual(stats["weighted_rms_cm"],.5)

    def test_cutoff_waits_until_all_residuals_are_small(self):
        for maximum,expected in ((20.,None),(9.,10)):
            p=self.pipeline
            stats={**self.stats,"max_abs_cm":maximum}
            p.evaluate=lambda *a,**k:(stats,None,SimpleNamespace(directory=Path("run/check")))
            thresholds=[]
            def run(*args,**kwargs):
                thresholds.append(kwargs["threshold"])
                return SimpleNamespace(parameters={"RE":1.1},history=[{"stability":1e-9}])
            p.run_duo=run
            p.fit_experiment(10,"trial")
            self.assertEqual(thresholds,[expected])

    def test_stationary_worse_fit_retains_model_without_step_retries(self):
        p=self.pipeline
        calls=[]
        def run(*args,**kwargs):
            calls.append(kwargs["scale"])
            return SimpleNamespace(parameters={"RE":1.2},history=[{"stability":1e-10}])
        p.run_duo=run
        sequence=iter([self.stats,{**self.stats,"weighted_rms_cm":2.}])
        p.evaluate=lambda *a,**k:(next(sequence),None,SimpleNamespace(directory=Path("run/check")))
        p.fit_experiment(10,"trial")
        self.assertEqual(p.p,{"RE":1.1})
        self.assertEqual(calls,[.1])


@unittest.skipUnless(os.environ.get("PEC_FIT_DUO"),"Set PEC_FIT_DUO to run real Duo smoke tests")
class RealDuoTests(unittest.TestCase):
    def test_monitored_convergence_and_fresh_validation(self):
        model=import_template(CO)
        text=render(model,model.parameters,jmax=0,vmax=1,iterations=5000,
                    abi_factor=1e7,energy_factor=1e-10)
        with tempfile.TemporaryDirectory() as d:
            runner=Runner(os.environ["PEC_FIT_DUO"],d,poll_seconds=0.1,resume=True)
            result=runner.run("fit",text,fitting=True,stability_tolerance=1.0)
            self.assertIsNotNone(result.controlled_stop)
            self.assertLess(result.history[-1]["iteration"],5000)
            self.assertIsNone(converged_checkpoint((result.directory/"run.out").read_text(),1.0,5000))
            cached=runner.run("fit",text,fitting=True,stability_tolerance=1.0)
            self.assertEqual(cached.parameters,result.parameters)
            check=runner.run("check",render(model,result.parameters,jmax=0,vmax=1,iterations=0))
            self.assertEqual(parse_en(check.directory/"fit.en")[(0,0)]["calculated"],0.)

    def test_timeout_is_a_failure(self):
        model=import_template(CO)
        text=render(model,model.parameters,jmax=0,vmax=1,iterations=5000,
                    abi_factor=1e7,energy_factor=1e-10)
        with tempfile.TemporaryDirectory() as d:
            runner=Runner(os.environ["PEC_FIT_DUO"],d,timeout=0.05)
            with self.assertRaisesRegex(FitError,"timed out"):
                runner.run("timeout",text,fitting=True)

    def test_straight_through_and_cache_integrity(self):
        model=import_template(ROOT/"CO/sample_outputs/CO_duo_fit_06.inp")
        text=render(model,model.parameters,jmax=0,vmax=80,iterations=0)
        with tempfile.TemporaryDirectory() as d:
            runner=Runner(os.environ["PEC_FIT_DUO"],d,resume=True)
            result=runner.run("smoke",text)
            levels=parse_en(result.directory/"fit.en")
            stats,_=assess_energies(levels,model.observations,0)
            self.assertEqual(stats["count"],41)
            grid_rows = [list(map(float,line.split())) for line in
                         (result.directory/"Potential_functions.dat").read_text().splitlines() if line.strip()]
            # The radii in this file are rounded to 10 decimal places.
            self.assertLess(max(abs(emo(row[0],model.parameters)-row[1]) for row in grid_rows),0.001)
            stamp=(result.directory/"completed.json").read_text()
            runner.run("smoke",text)
            self.assertEqual((result.directory/"completed.json").read_text(),stamp)
            (result.directory/"fit.en").write_text("corrupt")
            runner.run("smoke",text)
            self.assertTrue(parse_en(result.directory/"fit.en"))

    def test_fortran_stop_without_outputs_is_not_success(self):
        with tempfile.TemporaryDirectory() as d:
            runner=Runner(os.environ["PEC_FIT_DUO"],d)
            with self.assertRaises(FitError):
                runner.run("invalid","this_is_not_a_duo_keyword\n")


if __name__ == "__main__":
    unittest.main()
