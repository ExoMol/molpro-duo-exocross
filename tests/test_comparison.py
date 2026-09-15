import json
from pathlib import Path
import tempfile
import unittest

from compare_pec_models import compare
from pecfit.duo import write_csv,write_json
from pecfit.model import import_template,render
from pecfit.pipeline import load

ROOT=Path(__file__).resolve().parents[1]


class ComparisonTests(unittest.TestCase):
    def test_subset_selection_does_not_certify_failed_full_dataset(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);m=import_template(ROOT/'CO/sample_outputs/CO_duo_fit_01.inp')
            p=dict(m.parameters)
            for k in p:
                if k.startswith('B'):p[k]=2. if k=='B0' else 0.
            inp=root/'seed.inp';inp.write_text(render(m,p,jmax=max(o.j for o in m.observations),vmax=80,iterations=0))
            cfg=root/'one.json';cfg.write_text(json.dumps({'template':'seed.inp','output_directory':'run'}))
            _,c,m,sources=load(cfg);run=root/'run';run.mkdir()
            write_json(run/'manifest.json',{'config':c,'sources':sources})
            write_json(run/'result.json',{'status':'needs_review','parameters':p,'shape':{'ok':True},'Jmax':max(o.j for o in m.observations),
                       'grid':m.grid,'experiment':{'count':len(m.observations),'rms_cm':.1,'weighted_rms_cm':.1,'max_abs_cm':.1},
                       'basis_max_delta_cm':1.,'basis_converged':False,'basis_tolerance_cm':.001})
            write_csv(run/'final_residuals.csv',[{'J':o.j,'v':o.v,'residual_cm':.1,'input_weight':o.weight} for o in m.observations])
            write_csv(run/'basis_convergence.csv',[{'J':o.j,'v':o.v,'difference_cm':0. if o.energy<=30000 else 1.} for o in m.observations])
            result=compare([cfg],root/'comparison',calibration_max_energy=30000)
            self.assertIsNone(result['best'])
            self.assertEqual(result['best_calibration'],'one')
            self.assertFalse(result['rows'][0]['eligible'])
            self.assertTrue(result['rows'][0]['calibration_eligible'])


if __name__=='__main__':unittest.main()
