"""OBSERVED/PREDICTED catalog schemas must accept conforming documents and
reject malformed ones (contract rule 5/6 enforceable validation)."""
import json
import sys
from pathlib import Path

import jsonschema
import pytest

ROOT = Path(__file__).resolve().parents[2]
OBS_SCHEMA = json.loads((ROOT / "schemas/OBSERVED_MODE_CATALOG.schema.json").read_text())
PRED_SCHEMA = json.loads((ROOT / "schemas/PREDICTED_MODE_CATALOG.schema.json").read_text())


def test_smoke_stage_a_observed_validates():
    p = ROOT / "data/observed/SMOKE_5200120403_STAGE_A.json"
    if not p.exists():
        pytest.skip("smoke artifact not present")
    doc = json.loads(p.read_text())
    cat = doc.get("catalog")
    if not isinstance(cat, dict):
        pytest.skip("smoke catalog is a placeholder string, not a document")
    jsonschema.validate(cat, OBS_SCHEMA)


def test_predicted_catalog_builder_freeze_validates():
    sys.path.insert(0, str(ROOT / "src"))
    from ssz_p5.spectroscopy.physical_units import FrozenMassPrior
    from ssz_p5.spectroscopy.theory_catalog import (
        build_predicted_catalog,
        freeze_predicted_catalog,
        verify_frozen_catalog,
    )

    prior = FrozenMassPrior("unit", 10.0, 1.0, "test")
    rows = [{"omega_bar_re": 0.1, "observable": "q1", "u": 0.65, "L": 6,
             "omega2": 1.2, "cluster_residue_weight": 0.9, "ipr": 0.02,
             "localization_class": "EXTENDED"}]
    doc = build_predicted_catalog(
        rows, prior,
        member_hash="0" * 64,
        healthy_window={"u_min": 0.62, "u_max": 0.70,
                        "L_values": [6, 12, 20], "controls_all_pass": True},
        observables=["q1"])
    # unfrozen document must NOT validate (frozen must be const true):
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(doc, PRED_SCHEMA)
    frozen, digest = freeze_predicted_catalog(doc)
    assert frozen["frozen"] is True
    assert frozen["sha256"] == digest
    assert verify_frozen_catalog(frozen)
    jsonschema.validate(frozen, PRED_SCHEMA)


def test_predicted_schema_rejects_unfrozen():
    bad = {"catalog": "PREDICTED_MODE_CATALOG", "version": 1, "frozen": False,
           "sha256": "", "member_hash": "x", "modes": [],
           "healthy_window": {"u_min": 0.0, "u_max": 1.0, "L_values": [],
                              "controls_all_pass": True},
           "blindness_assert": {"did_not_read_observed_catalog": True,
                                "enforced_by": "t"}}
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad, PRED_SCHEMA)
