from __future__ import annotations

import ast
from pathlib import Path
import tempfile
import unittest

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]


def load_function(script_relpath: str, function_name: str, namespace: dict | None = None):
    script_path = REPO_ROOT / script_relpath
    source = script_path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(script_path))

    fn_nodes = [
        node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == function_name
    ]
    if not fn_nodes:
        raise AssertionError(f"Function {function_name!r} not found in {script_relpath}")

    module = ast.Module(body=fn_nodes, type_ignores=[])
    compiled = compile(module, filename=str(script_path), mode="exec")
    env: dict = {}
    if namespace:
        env.update(namespace)
    exec(compiled, env)
    return env[function_name]


class TestScriptHelpers(unittest.TestCase):
    def test_gaussian_has_expected_peak_and_symmetry(self):
        gaussian = load_function("scripts/plot_dang.py", "gaussian", namespace={"np": np})

        amp = 7.5
        mu = -3.0
        std = 0.2
        center = gaussian(mu, amp, mu, std)
        left = gaussian(mu - 0.1, amp, mu, std)
        right = gaussian(mu + 0.1, amp, mu, std)

        self.assertTrue(np.isclose(center, amp))
        self.assertTrue(np.isclose(left, right))
        self.assertLess(left, center)

    def test_getminmax_returns_95_percent_interval(self):
        getminmax = load_function("scripts/plot_dang.py", "getminmax", namespace={"np": np})

        values = np.arange(0, 100)
        lo, hi = getminmax(values)

        self.assertTrue(np.isclose(lo, np.percentile(values, 2.5)))
        self.assertTrue(np.isclose(hi, np.percentile(values, 97.5)))

    def test_nearest_returns_closest_index(self):
        nearest = load_function("scripts/dust_amp_plot.py", "nearest", namespace={"np": np})
        values = np.array([1.0, 5.0, 8.0, 10.0])

        self.assertEqual(nearest(values, 7.2), 2)
        self.assertEqual(nearest(values, 4.9), 1)

    def test_make_mean_maps_read_params_parses_minimal_paramfile(self):
        read_params = load_function("scripts/make_mean_maps.py", "read_params")

        with tempfile.TemporaryDirectory() as tmpdir:
            param_file = Path(tmpdir) / "param_test.txt"
            param_file.write_text(
                "\n".join(
                    [
                        "NUMBAND=2",
                        "NUMGIBBS=00123",
                        "MASKFILE='mask.fits'",
                        "BAND_LABEL001='bp_030'",
                        "BAND_FREQ001=28.4",
                        "BAND_LABEL002='bp_044'",
                        "BAND_FREQ002=44.1",
                    ]
                ),
                encoding="utf-8",
            )

            labels, freq, numgibbs, maskfile = read_params(str(param_file))

        self.assertEqual(labels, ["'bp_030'", "'bp_044'"])
        self.assertEqual(freq, [28.4, 44.1])
        self.assertEqual(numgibbs, 123)
        self.assertEqual(maskfile, "'mask.fits'")


if __name__ == "__main__":
    unittest.main()
