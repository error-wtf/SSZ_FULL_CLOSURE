"""P2 positive controls (injection/recovery) + P3 pipeline freeze."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.observations.injection_recovery import (  # noqa: E402
    injection_recovery_trial,
    inject_signal,
    positive_control_battery,
    recover_frequency_and_width,
)
from ssz_p5.observations.pipeline_freeze import (  # noqa: E402
    code_file_hashes,
    freeze_pipeline_record,
    verify_pipeline_record,
)


def test_positive_control_battery_all_recovered():
    r = positive_control_battery()
    assert r["all_recovered"], (
        f"only {r['n_recovered']}/{r['n_trials']} injections recovered")
    assert r["n_trials"] == 8


def test_all_control_frequencies_below_nyquist():
    """The control set must stay below Nyquist (dt=1/512 -> 256 Hz);
    an aliased injection would test the alias, not the finder."""
    rng = np.random.default_rng(0)
    dt = 1.0 / 512.0
    for f in (20.0, 55.0, 133.0, 200.0):
        assert f < 256.0
        tr = injection_recovery_trial(rng, dt, 32.0, f, 1.0, 4.0)
        assert tr["recovered"]


def test_injection_changes_the_series():
    rng = np.random.default_rng(1)
    counts = rng.poisson(50.0, 4096).astype(float)
    injected = inject_signal(counts, 1 / 512.0, 100.0, 1.0, 5.0)
    assert not np.array_equal(counts, injected)


def test_recovery_locator_is_deterministic():
    rng = np.random.default_rng(2)
    counts = rng.poisson(50.0, 4096).astype(float)
    injected = inject_signal(counts, 1 / 512.0, 100.0, 1.0, 5.0)
    a = recover_frequency_and_width(injected, 1 / 512.0, 100.0)
    b = recover_frequency_and_width(injected, 1 / 512.0, 100.0)
    assert a == b


# ---------------- pipeline freeze ----------------

def _record():
    return freeze_pipeline_record(
        analysis_config={"psd_normalization": "leahy",
                         "significance_threshold_sigma": 4.0},
        code_files={"tools/bridge_stage_a_observed_modes.py": "a" * 64},
        git_commit="deadbeefcafe",
        gti_policy={"min_continuous_s": 8.0, "drop_short_gaps": True},
        negative_controls={"spurious_55hz": "REJECTED"},
        positive_controls={"battery": "8/8"},
    )


def test_freeze_verify_round_trip():
    rec, digest = _record()
    assert rec["sha256"] == digest
    assert verify_pipeline_record(rec)


def test_tampered_record_fails_verification():
    rec, _ = _record()
    tampered = dict(rec, gti_policy={"min_continuous_s": 1.0,
                                     "drop_short_gaps": True})
    assert not verify_pipeline_record(tampered)


def test_wrong_version_fails_verification():
    rec, _ = _record()
    bad = dict(rec, record_version="FROZEN_PIPELINE_RECORD_V2")
    assert not verify_pipeline_record(bad)


def test_code_file_hashes_reports_absent(tmp_path):
    h = code_file_hashes(tmp_path, ["nope.py"])
    assert h["nope.py"] == "ABSENT"
