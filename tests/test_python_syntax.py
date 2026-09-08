from __future__ import annotations

import ast
from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]


class TestPythonSyntax(unittest.TestCase):
    def test_scripts_parse_as_python(self):
        scripts = [
            "scripts/average_chisq.py",
            "scripts/dust_amp_plot.py",
            "scripts/likelihood_plot.py",
            "scripts/make_mean_maps.py",
            "scripts/make_synch_diff_maps.py",
            "scripts/metrop_test.py",
            "scripts/parameter_plotter.py",
            "scripts/plot_dang.py",
        ]
        for script_relpath in scripts:
            with self.subTest(script=script_relpath):
                source = (REPO_ROOT / script_relpath).read_text(encoding="utf-8")
                ast.parse(source)


if __name__ == "__main__":
    unittest.main()
