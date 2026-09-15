import importlib.util
import math
import os
from pathlib import Path
import tempfile
import unittest
from dataclasses import replace

AVAILABLE=importlib.util.find_spec('numpy') is not None and importlib.util.find_spec('scipy') is not None
if AVAILABLE:
    import numpy as np
    from pecfit.localized import Spectrum,curves,native_check,precise_duo_energies,select_localized,optimize_penalized,physical
from pecfit.duo import Runner
from pecfit.forms import potential,bob
from pecfit.model import FitError,Observation,import_template
from pecfit.pipeline import merge_config

ROOT=Path(__file__).resolve().parents[1]
P={'TE':0.,'RE':1.27454677,'AE':37231.89594,'R12':2.5,'DELTA':2.,'PHIINF':.5200237764,'N':6.,
   'PHI0':-3.9755931298,'PHI1':.4132556507,'PHI2':-.6684514359,'PHI3':-.9047392949,
   'PHI4':-1.2947826268,'PHI5':2.9887171627,'PHI6':-.119832624,'PHI7':-19.5673843459,
   'PHI8':9.9504524501,'PHI9':17.6092973464,'PHI10':-10.1890184622,
   'BR_RE':1.27454677,'BR_BETA':.8,'BR_GAMMA':.02,'BR_P':6.,'BR_B0':.001,'BR_B1':-.002,'BR_BINF':0.}


@unittest.skipUnless(AVAILABLE,'optional NumPy/SciPy dependencies')
class LocalizedTests(unittest.TestCase):
    def model(self):
        model=import_template(ROOT/'CO/sample_outputs/CO_duo_fit_01.inp')
        observations=[Observation(j,'e',v+1,0. if (j,v)==(0,0) else 2886.*v+10.*j*(j+1),'X',v,100/math.sqrt((j+1)*(v+1)))
                      for j in (0,1,2) for v in (0,1,2)]
        return replace(model,atoms='H Cl',masses='1.00782503223 34.968852682',parameters=dict(P),
                       fitted=['AE','PHI0','BR_B0'],observations=observations,grid=(.65,5.8,251))

    def test_vector_forms_match_scalar(self):
        emo=import_template(ROOT/'CO/sample_outputs/CO_duo_fit_01.inp').parameters
        for p in (P,emo,{**emo,'RREF':1.5,'PL':2.,'NL':1.},{**emo,'RREF':0.}):
            r=np.linspace(.8,8,99);v,q=curves(r,p)
            np.testing.assert_allclose(v,[potential(float(x),p) for x in r],rtol=1e-12,atol=1e-8)
            np.testing.assert_allclose(q,[bob(float(x),p) for x in r],rtol=1e-12,atol=1e-14)

    def test_hellmann_feynman_derivative_includes_zero_reference(self):
        spec=Spectrum(self.model(),(.65,5.8,251));h=1e-5
        plus={**P,'PHI0':P['PHI0']+h};minus={**P,'PHI0':P['PHI0']-h}
        dv=(curves(spec.r,plus)[0]-curves(spec.r,minus)[0])/(2*h)
        dq=np.zeros_like(dv)
        _,_,jac,_=spec.evaluate(P,(dv[:,None],dq[:,None]))
        numeric=(spec.evaluate(plus)[0]-spec.evaluate(minus)[0])/(2*h)
        np.testing.assert_allclose(jac[:,0],numeric,rtol=2e-7,atol=2e-5)
        self.assertAlmostEqual(jac[0,0],0.,places=10)

    def test_missing_localized_level_fails_instead_of_matching_by_energy(self):
        r=np.arange(1.,5.);eff=np.array([0.,3.,2.,1.]);vec=np.eye(4)
        with self.assertRaises(FitError):
            select_localized(np.array([0.,4.,5.,6.]),vec,eff,r,1.,[0,1])

    def test_penalized_refinement_recovers_synthetic_levels_and_physical_shape(self):
        model=self.model();model.fitted=['PHI0','BR_B0']
        target={**P,'PHI0':P['PHI0']+.001,'BR_B0':P['BR_B0']+.00001}
        exact=Spectrum(model,model.grid).evaluate(target)[0]
        model.observations=[replace(o,energy=float(e)) for o,e in zip(model.observations,exact)]
        config=merge_config({'validation':{'shape_range_angstrom':[.65,5.8],'potential_ceiling_cm':1e8}})
        before=Spectrum(model,model.grid).evaluate(model.parameters)[0]-exact
        with tempfile.TemporaryDirectory() as d:
            p,history=optimize_penalized(model,config,Path(d),max_iterations=10)
        after=Spectrum(model,model.grid).evaluate(p)[0]-exact
        self.assertLess(np.linalg.norm(after),np.linalg.norm(before)/100)
        self.assertTrue(physical(model,p,config)['ok'])
        self.assertTrue(history)

    def test_precise_energy_duplicate_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'energy.dat';p.write_text('0.0 1 0.000000000000 1 0 0 0.0 0.0 0.0 + ||X\n'*2)
            with self.assertRaises(FitError):precise_duo_energies(p)

    @unittest.skipUnless(os.environ.get('PEC_FIT_DUO'),'set PEC_FIT_DUO for native validation')
    def test_native_sinc_energies_with_rotational_bob(self):
        with tempfile.TemporaryDirectory() as d:
            model=self.model();runner=Runner(os.environ['PEC_FIT_DUO'],Path(d)/'duo',threads=1)
            config=merge_config({'vmax':80})
            rows,_,_=native_check(model,P,model.grid,config,runner,'verify')
            self.assertEqual(len(rows),9)
            self.assertLess(max(abs(r['python_duo_delta_cm']) for r in rows),1e-5)


if __name__=='__main__':unittest.main()
