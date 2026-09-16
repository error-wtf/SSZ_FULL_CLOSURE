import subprocess
import sys
from pathlib import Path


def test_constructive_audit_passes():
    root = Path(__file__).resolve().parents[2]
    p = subprocess.run(
        [
            sys.executable,
            str(root / "ssz_p5_full_closure_auditor.py"),
            "--data-dir",
            str(root),
            "--quick",
        ],
        capture_output=True,
        text=True,
    )
    assert p.returncode == 0, p.stdout + p.stderr
    assert "FULL_CONSTRUCTIVE_CLOSURE" in p.stdout and "PASS" in p.stdout
