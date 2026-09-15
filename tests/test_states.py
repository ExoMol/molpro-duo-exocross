import bz2
import json
import math
from pathlib import Path
import tempfile
import unittest

from pecfit.__main__ import main
from pecfit.model import FitError, observations
from pecfit.pipeline import load, merge_config
from pecfit.states import read_states


SAMPLE = "1 0.000000 8 0 0\n42 2885.976400 8 0 1\n710 36515.850700 664 41 11\n"


class StatesTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.path = self.root/"test.states"
        self.path.write_text(SAMPLE)

    def test_mapping_weights_and_level_number(self):
        obs = read_states(self.path)
        self.assertEqual([(o.j, o.v, o.n, o.state, o.parity) for o in obs],
                         [(0,0,1,"X","e"),(0,1,2,"X","e"),(41,11,12,"X","e")])
        self.assertEqual(obs[1].energy, 2885.9764)
        self.assertEqual(obs[0].weight, 100.)
        self.assertAlmostEqual(obs[2].weight, 100/math.sqrt(12)/math.sqrt(42))
        self.assertEqual(observations([o.line() for o in obs])[-1].n, 12)

    def test_compressed_and_explicit_columns(self):
        path = self.root/"test.states.bz2"
        with bz2.open(path, "wt") as stream:
            stream.write("1 0 8 0 label 0\n42 2885.976400 8 0 label 1\n")
        obs = read_states(path, v_column=6, weight_scale=10)
        self.assertAlmostEqual(obs[1].weight, 10/math.sqrt(2))

    def test_bad_rows_are_not_silently_skipped(self):
        for row in ("2 2 8 0", "2 NaN 8 0 1", "2 -1 8 0 1", "1 2 8 0 1",
                    "2 2 8 0.5 1", "2 2 8 0 -1", "2 2 0 0 1", "2 2 8 0 0"):
            with self.subTest(row=row):
                self.path.write_text("1 0 8 0 0\n"+row)
                with self.assertRaises(FitError):
                    read_states(self.path)

    def test_invalid_mapping_or_weight(self):
        for kwargs in ({"v_column":4},{"j_column":3},{"v_column":5.5},{"weight_scale":0},{"weight_scale":float("nan")}):
            with self.subTest(kwargs=kwargs), self.assertRaises(FitError):
                read_states(self.path, **kwargs)
        for config in ({"states_j_column":3},{"states_v_column":4},{"states_weight_scale":0}):
            with self.assertRaises(FitError):
                merge_config(config)

    def test_cli_conversion_preserves_existing_output(self):
        output = self.root/"energies.txt"
        self.assertEqual(main(["convert-states",str(self.path),str(output)]),0)
        saved = output.read_bytes()
        self.assertEqual(main(["convert-states",str(self.path),str(output)]),1)
        self.assertEqual(saved,output.read_bytes())

    def test_init_reads_live_states_and_hashes_both_sources(self):
        pec = self.root/"pec.txt"
        pec.write_text("0.8 -2.0\n1.3 -3.0\n3.5 -2.5\n")
        config = self.root/"fit.json"
        self.assertEqual(main(["init","--config",str(config),"--pec",str(pec),"--states",str(self.path),
                               "--atoms","H","Cl"]),0)
        _, c, m, sources = load(config)
        self.assertEqual(len(m.observations),3)
        self.assertEqual(c["states_file"],"test.states")
        self.assertIn(str(pec),sources)
        self.assertIn(str(self.path),sources)
        previous = sources[str(self.path)]
        self.path.write_text(SAMPLE.replace("2885.976400","2885.970000"))
        _, _, m, sources = load(config)
        self.assertEqual(m.observations[1].energy,2885.97)
        self.assertNotEqual(previous,sources[str(self.path)])


if __name__ == "__main__":
    unittest.main()
