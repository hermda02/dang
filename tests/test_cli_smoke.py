from __future__ import annotations

from pathlib import Path
import subprocess
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts/metrop_test.py"


class TestCliSmoke(unittest.TestCase):
    def test_metrop_help_prints_usage(self):
        result = subprocess.run(
            ["python3", str(SCRIPT), "--help"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=False,
        )

        output = (result.stdout or "") + (result.stderr or "")
        self.assertIn("Usage:", output)

    def test_metrop_no_args_prints_usage(self):
        result = subprocess.run(
            ["python3", str(SCRIPT)],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=False,
        )

        output = (result.stdout or "") + (result.stderr or "")
        self.assertIn("Usage:", output)


if __name__ == "__main__":
    unittest.main()
