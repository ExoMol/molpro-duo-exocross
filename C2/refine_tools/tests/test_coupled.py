import copy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import coupled
import curves
import prepare


MODEL = '''nstates 2
poten 1
name "X (test)"
type TWO_COUPLED_EMOS
values
V0 0
RE 1 fit (prior 2)
B0 3
V0 100
RE 2 fit
B0 4
end
(unterminated single-line comment
spin-orbit-x 1 2
name "<X|SO|a>"
type BOBLEROY
morphing
values
RE 1.3
RREF -1
P 2
NT 0
B0 1 fit
BINF 0
end
poten 2
name "a"
type EMO
values
V0 100
end
FITTING
JLIST 0 - 1
itmax 0
energies
0 + 1 0.00001 1 0 0 0 0 100
1 - 1 100 2 0 -1 -1 -2 2
end
'''


class CoupledTests(unittest.TestCase):
    def test_repeated_parameters_are_positional_and_metadata_survives(self):
        src = coupled.fields(MODEL)
        target = copy.deepcopy(src)
        target[('POTEN',1,1)].values[4] = 2.05
        out = coupled.transfer(MODEL, target)
        self.assertEqual(coupled.fields(out)[('POTEN',1,1)].values[1], 1)
        self.assertEqual(coupled.fields(out)[('POTEN',1,1)].values[4], 2.05)
        self.assertIn('morphing', out)
        self.assertIn('fit (prior 2)', out)
        self.assertEqual(out.split('FITTING')[1], MODEL.split('FITTING')[1])

    def test_duo_comment_scope_and_class_alias(self):
        parsed = coupled.fields(MODEL)
        self.assertEqual(len(parsed), 3)
        self.assertIn(('SPINORBIT',1,2), parsed)
        self.assertEqual(parsed[('POTEN',1,1)].name, 'X (test)')

    def test_parameter_layout_changes_fail(self):
        cp = coupled.fields(MODEL)
        cp[('POTEN',1,1)].labels.reverse()
        with self.assertRaisesRegex(ValueError, 'layout'):
            coupled.transfer(MODEL, cp)

    def test_checkpoint_ignores_echo_and_rounded_section(self):
        block = 'Iteration = 1\nParameters:\n'+MODEL.split('FITTING')[0]+'Fitted parameters (rounded):\n'
        self.assertEqual(len(coupled.checkpoints(MODEL+'\n'+block)), 1)

    def test_full_observation_coverage_required(self):
        with self.assertRaisesRegex(ValueError, 'coverage'):
            coupled.assess(MODEL, '', 1)

    def test_original_weights_retained_when_duo_zeroes_weight(self):
        output = '1 1 0.0 + 0.0000 0.0000 0.0000 1E-1 (1 0 0 0 0)(1 0 0 0 0)\n2 1 1.0 - 100.0000 101.0000 -1.0000 0E0 (2 0 1 1 2)(2 0 -1 -1 -2)\n'
        stats, rows = coupled.assess(MODEL, output, 1)
        self.assertEqual(stats['positive_weight_count'], 2)
        self.assertEqual(stats['zero_duo_weight_count'], 1)
        self.assertEqual(stats['assignment_mismatches'], 0)
        self.assertAlmostEqual(stats['weighted_rms_cm'], ((100*1e-10+2)/102)**.5)

    def test_potential_and_coupled_components(self):
        p = [0,1.2,40000,-1,2,2,0,0,2]
        q = [10000,1.4,50000,-1,2,2,0,0,2]
        self.assertEqual(curves.value('EMO',p,1.2), 0)
        low = curves.value('TWO_COUPLED_EMOS',p+q+[0,1,1.3,100,1],1.3)
        high = curves.value('TWO_COUPLED_EMOS',p+q+[0,1,1.3,100,2],1.3)
        self.assertAlmostEqual(low+high,curves.value('EMO',p,1.3)+curves.value('EMO',q,1.3))
        self.assertLess(low, high)

    def test_negative_energies_and_nonfinite_values(self):
        self.assertEqual(coupled.number('-2D-3'), -.002)
        with self.assertRaises(ValueError):
            coupled.number('nan')

    def test_duplicate_preparation_is_explicit_and_weighted(self):
        model = MODEL.replace('0 + 1 0.00001 1 0 0 0 0 100', '0 + 1 0.00001 1 0 0 0 0 100\n0 + 1 0.00041 1 0 0 0 0 1')
        unchanged, audit = prepare.prepare(model)
        self.assertEqual(len(coupled.observations(unchanged)), 3)
        self.assertFalse(audit['duplicate_groups'])
        result, audit = prepare.prepare(model, merge_duplicates=True)
        rows = coupled.observations(result)
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]['input_weight'],101)
        self.assertAlmostEqual(rows[0]['observed_cm'],(.00001*100+.00041)/101, places=9)
        self.assertEqual(len(audit['duplicate_groups'][0]['original_rows']),2)

    def test_inconsistent_duplicate_preparation_fails(self):
        model = MODEL.replace('0 + 1 0.00001 1 0 0 0 0 100', '0 + 1 0.00001 1 0 0 0 0 100\n0 + 1 1.5 1 0 0 0 0 1')
        with self.assertRaisesRegex(ValueError, 'Inconsistent'):
            prepare.prepare(model, merge_duplicates=True)


if __name__ == '__main__':
    unittest.main()
