import math
import os
from pathlib import Path
import tempfile
import unittest
from dataclasses import replace

from pecfit.duo import Runner, parse_en, parse_output
from pecfit.forms import bob, bob_shape, kind, potential, validate_parameters, centrifugal_shape
from pecfit.model import FitError, import_template, render, shape_check

ROOT=Path(__file__).resolve().parents[1]
MLJ={"TE":0.,"RE":1.3,"AE":38000.,"R12":4.,"DELTA":2.25,"PHIINF":0.5,"N":6.,"PHI0":-3.6,"PHI1":0.5,"PHI2":0.}
BR={"BR_RE":1.3,"BR_BETA":0.8,"BR_GAMMA":0.02,"BR_P":6.,"BR_B0":0.001,"BR_B1":-0.002,"BR_B2":0.,"BR_BINF":0.}


class FormTests(unittest.TestCase):
    def test_monotonic_wall_can_still_be_rejected_as_excessive(self):
        m=import_template(ROOT/'CO/sample_outputs/CO_duo_fit_01.inp')
        p={k:(0. if k.startswith('B') else v) for k,v in m.parameters.items()};p['B0']=20.
        self.assertTrue(shape_check(p,(.6,6))['ok'])
        self.assertFalse(shape_check(p,(.6,6),max_potential=1e8)['ok'])
    def test_centrifugal_shape_rejects_outer_well(self):
        p={"TE":0.,"RE":1.27455640639563,"AE":38255.1524340919,"R12":4.,"DELTA":2.25,"PHIINF":.530187972651004,"N":6.,
           **dict(zip([f"PHI{i}" for i in range(11)],[-3.65995834582118,1.13743430496408,.491611780988227,.745856971414173,.823847849551519,
                       2.18508703298665,-4.99575624592288,-2.98669221752415,29.7813885452912,-44.0793555363314,16.011170555391]))}
        self.assertTrue(shape_check(p,(.5,12))["ok"])
        self.assertFalse(centrifugal_shape(p,"1.00782503223 34.968852682",41,(.5,12))["ok"])
        self.assertTrue(centrifugal_shape(MLJ,"1.00782503223 34.968852682",41,(.5,12))["ok"])
    def test_mlj_equilibrium_and_asymptotic_convention(self):
        self.assertEqual(potential(MLJ["RE"],MLJ),0.)
        c6=2*MLJ["AE"]*MLJ["RE"]**6*math.exp(-2*MLJ["PHIINF"])
        r=50.
        self.assertAlmostEqual((MLJ["AE"]-potential(r,MLJ))*r**6/c6,math.exp(2*MLJ["PHIINF"]-MLJ["PHIINF"]*2*(r-1.3)/(r+1.3)),places=4)
        self.assertTrue(shape_check(MLJ,(.6,12))["ok"])

    def test_bob_decay_and_bounds(self):
        self.assertEqual(bob(1.3,BR),.001)
        self.assertAlmostEqual(bob(100,BR),0.,places=12)
        self.assertTrue(bob_shape(BR,(.6,12))["ok"])
        self.assertFalse(bob_shape({**BR,"BR_B0":-2.},(.6,12))["ok"])

    def test_object_roundtrip_does_not_mix_coefficients(self):
        m=import_template(ROOT/"CO/sample_outputs/CO_duo_fit_01.inp")
        p={**MLJ,**BR}
        fitted=["RE","AE","PHI0","PHI1","PHI2","BR_B0","BR_B1","BR_B2"]
        m=replace(m,parameters=p,fitted=fitted)
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"input.inp"
            path.write_text(render(m,p,jmax=0,vmax=80,iterations=1).replace('bob-rot X','\n\nbob-rot X'))
            restored=import_template(path)
        self.assertEqual(restored.parameters,p)
        self.assertEqual(restored.fitted,fitted)
        text="Iteration = 1\nParameters:\n"+render(m,p,jmax=0,vmax=80,iterations=1).split("poten X\n",1)[1]
        text="Iteration = 1\nParameters:\nPOTEN X\n"+text.split("Parameters:\n",1)[1].split("abinitio poten",1)[0]+"Fitted parameters (rounded):\n"
        self.assertEqual(parse_output(text)[0][0]["parameters"],p)

    def test_disallowed_structure_and_missing_coefficients(self):
        for p in ({**MLJ,"TE":1.},{**MLJ,"N":1.5},{**MLJ,"PHI4":0.},{**MLJ,**BR,"BR_BINF":1.}):
            with self.assertRaises(FitError): validate_parameters(p,["RE"])
        with self.assertRaises(FitError):validate_parameters(MLJ,["PHIINF"])

    @unittest.skipUnless(os.environ.get("PEC_FIT_DUO"),"Set PEC_FIT_DUO for native function checks")
    def test_native_mlj_bob_values_and_checkpoint(self):
        base=import_template(ROOT/"CO/sample_outputs/CO_duo_fit_01.inp")
        p={**MLJ,**BR}
        model=replace(base,parameters=p,fitted=["RE","PHI0","BR_B0"],grid=(.7,5.,301),
                      observations=[o for o in base.observations if o.v<3 and o.j<=1])
        with tempfile.TemporaryDirectory() as d:
            runner=Runner(os.environ["PEC_FIT_DUO"],d)
            run=runner.run("mlj_bob",render(model,p,jmax=1,vmax=30,iterations=1,scale=.001),fitting=True)
            self.assertIn("BR_B0",run.parameters)
            self.assertIn("PHI0",run.parameters)
            self.assertEqual(run.parameters["BR_RE"],1.3)
            check=runner.run("check",render(model,p,jmax=1,vmax=30,iterations=0))
            rows=[list(map(float,line.split())) for line in (check.directory/"Potential_functions.dat").read_text().splitlines() if line.strip()]
            self.assertLess(max(abs(potential(row[0],p)-row[1]) for row in rows),.001)
            # A constant dimensionless BOB rescales the J(J+1)/r^2 operator.
            # J=0 eigenvalues must be unchanged when this field is removed.
            without={k:v for k,v in p.items() if not k.startswith("BR_")}
            clean=replace(model,parameters=without,fitted=["RE","PHI0"])
            no=runner.run("no_bob",render(clean,without,jmax=1,vmax=30,iterations=0))
            a,b=parse_en(check.directory/"fit.en"),parse_en(no.directory/"fit.en")
            self.assertEqual(a[(0,1)]["calculated"],b[(0,1)]["calculated"])
            self.assertNotEqual(a[(1,0)]["calculated"],b[(1,0)]["calculated"])
