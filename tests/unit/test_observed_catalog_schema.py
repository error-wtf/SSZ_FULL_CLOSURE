"""OBSERVED catalog builder must produce schema-conformant documents;
freeze/verify round-trip must be consistent (contract rule 5)."""
import json
import sys
from pathlib import Path

import jsonschema
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.observations.blind_catalog import (  # noqa: E402
    build_observed_catalog,
    freeze_observed_catalog,
    verify_frozen_catalog,
)

SCHEMA = json.loads(
    (ROOT / "schemas/OBSERVED_MODE_CATALOG.schema.json").read_text())

CONFIG = {
    "psd_normalization": "leahy",
    "freq_range_hz": [0.1, 2048.0],
    "energy_bands_kev": [0.5, 2.0, 10.0],
    "significance_threshold_sigma": 4.0,
    "trial_correction": "bonferroni",
    "frozen_before_first_psd": True,
}


def _doc():
    return build_observed_catalog(
        target="MAXI J1820+070",
        obsids=["5200120403"],
        total_exposure_s=1896.0,
        analysis_config=CONFIG,
        candidates=[{
            "f_hz": 55.6, "sigma_f_hz": 0.5, "amplitude_rms": 0.01,
            "covariance": [[1e-6, 0.0], [0.0, 1e-4]],
            "significance_sigma": 2.1, "energy_dependence": "none",
            "obs_occurrence": 1,
        }],
        instrumental_controls=[{
            "f_hz": 55.6, "predicted_by": "gti_gap_comb",
            "verdict": "RECOVERED_AND_REJECTED",
            "rejection_reason": "on gap comb and unstable in continuous segments",
        }],
    )


def test_unfrozen_document_rejected_by_schema():
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(_doc(), SCHEMA)


def test_frozen_document_validates():
    frozen, digest = freeze_observed_catalog(_doc())
    assert frozen["frozen"] is True
    assert frozen["sha256"] == digest
    jsonschema.validate(frozen, SCHEMA)


def test_freeze_verify_round_trip():
    frozen, _ = freeze_observed_catalog(_doc())
    assert verify_frozen_catalog(frozen)
    tampered = dict(frozen)
    tampered["modes"] = [dict(frozen["modes"][0], f_hz=99.0)]
    assert not verify_frozen_catalog(tampered)


def test_double_freeze_refused():
    frozen, _ = freeze_observed_catalog(_doc())
    with pytest.raises(ValueError):
        freeze_observed_catalog(frozen)


def test_artifact_mode_rejected_by_schema():
    bad = freeze_observed_catalog(_doc())[0]
    bad["modes"][0]["energy_dependence"] = "purple"
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad, SCHEMA)
