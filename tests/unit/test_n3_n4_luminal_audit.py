"""Regression: the real frozen N2 artifact cannot open a scientific gate."""
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PROFILE = ROOT / 'data/generated/spectral/N2_LUMINAL_BACKGROUND_PROFILE_V2.csv'


def test_real_runner_blocks_invalid_background_and_preserves_all_rows(tmp_path):
    before = hashlib.sha256(PROFILE.read_bytes()).hexdigest()
    out = tmp_path / 'decision.json'
    result = subprocess.run(
        [sys.executable, str(ROOT / 'tools/run_n3_n4_luminal_decision.py'),
         '--output', str(out)], cwd=ROOT, env={**os.environ, 'PYTHONPATH': str(ROOT / 'src')},
        capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stderr
    report = json.loads(out.read_text())
    assert report['verdict'] == 'BLOCKED'
    assert report['physical_qnm_claim_allowed'] is False
    assert report['profile_rows'] == report['diagnostic_rows'] == 200
    assert report['n3']['status'] == 'NOT_EVALUABLE'
    assert report['n4']['status'] == 'NOT_EVALUABLE'
    assert report['n4']['executed_L'] == []
    assert report['n4']['requested_L'] == list(range(6, 1001))
    assert not out.with_suffix('.stream41.csv').exists()
    assert hashlib.sha256(PROFILE.read_bytes()).hexdigest() == before


def load_runner():
    import importlib.util
    path = ROOT / 'tools/run_n3_n4_luminal_decision.py'
    spec = importlib.util.spec_from_file_location('n3_audit', path)
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    return runner


def test_kinetic_precheck_rejects_partial_domain_even_if_positive():
    runner = load_runner()
    assert runner.kinetic_precheck([1.0] * 54, expected_rows=200,
                                    certified_error_bound=0.0)['status'] == 'NOT_EVALUABLE'
