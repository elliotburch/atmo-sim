"""Launcher integration checks; run with python3 -m unittest -v."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import run_planets as launcher


class LauncherTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.sim = launcher.load_library(launcher.HERE / "atmo-sim.py")
        cls.planets = launcher.load_planets(launcher.HERE / "example_planets.json", cls.sim)

    def test_named_run_matches_library(self):
        with tempfile.TemporaryDirectory() as tmp:
            process = subprocess.run([
                sys.executable, str(launcher.HERE / "run_planets.py"),
                "--planet", "eArTh", "--epochs", "4", "--output", tmp,
            ], cwd=tmp, text=True, capture_output=True)
            self.assertEqual(process.returncode, 0, process.stderr)
            with (Path(tmp) / "summary.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 1)
            p = self.planets[0]
            state, bulk, redox, _ = self.sim.evolve(p, epochs=4)
            expected = self.sim.result_row(p, state, bulk, redox)
            for key, value in expected.items():
                if type(value) is float:
                    self.assertEqual(float(rows[0][key]), value, key)
                else:
                    self.assertEqual(rows[0][key], str(value), key)
            for filename in ("history", "columns", "layers", "reservoirs", "vertical_profiles"):
                with (Path(tmp) / f"{filename}.csv").open() as stream:
                    self.assertTrue(list(csv.DictReader(stream)), filename)

    def test_invalid_cli_does_not_write_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "should-not-exist"
            for args in (["--planet", "missing"], ["--all", "--epochs", "0"],
                         ["--all", "--planet", "Earth"]):
                result = subprocess.run([
                    sys.executable, str(launcher.HERE / "run_planets.py"),
                    *args, "--output", str(output),
                ], text=True, capture_output=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(output.exists())

    def test_invalid_json_inputs(self):
        original = json.loads((launcher.HERE / "example_planets.json").read_text())
        bad_entries = [
            {**original["planets"][0], "mass_e": -1},
            {**original["planets"][0], "flux": "1361"},
            {**original["planets"][0], "crust": {"typo": 1}},
            {**original["planets"][0], "tidally_locked": "false"},
        ]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.json"
            for entry in bad_entries:
                path.write_text(json.dumps(dict(schema_version=1, planets=[entry])))
                with self.assertRaises(ValueError):
                    launcher.load_planets(path, self.sim)
            path.write_text(json.dumps(dict(schema_version=1, planets=[original["planets"][0]] * 2)))
            with self.assertRaisesRegex(ValueError, "Duplicate"):
                launcher.load_planets(path, self.sim)


if __name__ == "__main__":
    unittest.main()
