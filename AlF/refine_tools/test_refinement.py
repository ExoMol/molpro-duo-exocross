import math
from pathlib import Path
import unittest
from types import SimpleNamespace
import numpy as np
import coupled
from configure import configure
from refine import irls_weights
from barrier import repulsive, components

class ScientificChecks(unittest.TestCase):
    def test_repulsive_inverse_power_sign_and_units(self):
        # Atomic-unit conversion: potential at exactly one bohr must be -C6 Eh.
        bohr=.529177210903; eh=219474.63136320
        p=[0,0,.3,8,0,0,0,0,0,-75.23*eh*bohr**6,0,0]
        self.assertAlmostEqual(float(repulsive(np.array([bohr]),p)[0])/eh,-75.23,places=10)

    def test_coupled_lower_branch_is_matrix_eigenvalue(self):
        a=[43950,1.65,80000,-1,6,6,0,0,1.7]
        b=[55564,2e7,.3,8,0,0,0,0,-73495,-362562,0,0]
        c=[3000,2.55,0,-1,6,6,4,4,.8,0,0,0,0]
        p=a+b+c+[1]
        v1,v2,w=components(np.array([1.65,2.54,8]),p)
        for x,y,z in zip(v1,v2,w):
            exact=np.linalg.eigvalsh([[x,z],[z,y]])[0]
            self.assertAlmostEqual(.5*(x+y-math.hypot(x-y,2*z)),exact,places=8)

    def test_irls_uses_original_weights_and_keeps_inactive_zero(self):
        rows=[dict(input_weight=w,residual_cm=0) for w in (1,2,0)]
        self.assertEqual(irls_weights(rows,.00001),[1/3,2/3,0])
        rows[1]['residual_cm']=10
        weights=irls_weights(rows,.00001)
        self.assertLess(weights[1],weights[0]*1e-6)
        self.assertEqual(weights[2],0)

    def test_duplicate_named_parameters_transfer_by_position(self):
        text='poten 2\nname "A"\nlambda 1\nmult 1\ntype coupled-pec\nvalues\nB6 -2\nB6 -3 fit\nend\n'
        changed=configure(text,values=['POTEN:2:2:1=-4'])
        output=coupled.transfer(text,coupled.fields(changed))
        self.assertEqual(coupled.fields(output)['POTEN',2,2].values,[-2,-4])

    def test_ordinary_weighted_duplicate_merge_is_equivalent(self):
        observations=[(10.,1.),(10.54,1.),(9.9,.5)]
        weight=sum(w for _,w in observations)
        mean=sum(e*w for e,w in observations)/weight
        constant=sum(w*(e-mean)**2 for e,w in observations)
        for predicted in (8.,10.,10.3,12.):
            self.assertAlmostEqual(sum(w*(e-predicted)**2 for e,w in observations),weight*(mean-predicted)**2+constant)

if __name__=='__main__': unittest.main()
