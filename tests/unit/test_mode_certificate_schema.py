"""MODE_CERTIFICATE_V1 output must validate against its JSON schema."""
import json
import sys
from pathlib import Path

import jsonschema
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from ssz_p5.qnm.mode_certificate import certify_mode  # noqa: E402

SCHEMA = json.loads(
    (Path(__file__).resolve().parents[2] / "schemas/MODE_CERTIFICATE.schema.json")
    .read_text()
)


def _cert(**overrides):
    n = 10
    base = dict(
        L=6,
        mode_index=0,
        kinetic_norms=[1.0],
        tracking_overlaps=[0.99],
        relative_gaps=[0.1],
        omega_values=[1.0, 1.0],
        observable_shares=[0.3],
        d2_value=0.9,
        scale_control_ok=True,
        basis_control_ok=True,
    )
    base.update(overrides)
    return certify_mode(**base)


def test_certified_certificate_validates():
    jsonschema.validate(_cert().to_json(), SCHEMA)


def test_rejected_certificate_validates():
    c = _cert(kinetic_norms=[-1.0])
    jsonschema.validate(c.to_json(), SCHEMA)
    assert c.verdict == "REJECTED_V1"


def test_not_evaluable_certificate_validates():
    c = _cert(d2_value=float("nan"), skip_axes=("A6",))
    jsonschema.validate(c.to_json(), SCHEMA)
    assert c.verdict == "NOT_EVALUABLE_V1"


def test_invalid_version_rejected_by_schema():
    bad = _cert().to_json()
    bad["certificate_version"] = "MODE_CERTIFICATE_V2"
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad, SCHEMA)


def test_bad_verdict_enum_rejected():
    bad = _cert().to_json()
    bad["verdict"] = "SUPER_CERTIFIED"
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad, SCHEMA)
