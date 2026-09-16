import json

import numpy as np
import pytest

from ssz_p5.coefficients.identities import a5_holonomic, require_v12
from ssz_p5.config import SLOT_NAMES, repo_root
from ssz_p5.jets.jet9d8 import derivative
from ssz_p5.provenance.blacklist import require_production_input
from ssz_p5.provenance.manifest import verify_sha256s
from ssz_p5.qnm.gate import require_direct_krgm_certificate


@pytest.mark.parametrize(
    "name",
    [
        "data/diagnostic/ssz_p5_SELECTED_41STREAM_V2_2026-09-16.csv",
        "archive/full_working_snapshot/ssz_p5_SELECTED_41STREAM_V2_2026-09-16.csv",
    ],
)
def test_blacklisted_actual_archive_rejected(name):
    with pytest.raises(ValueError):
        require_production_input(repo_root() / name)


def test_forged_boolean_certificate_rejected(tmp_path):
    p = tmp_path / "cert.json"
    p.write_text(json.dumps({"pass": True}))
    with pytest.raises(RuntimeError):
        require_direct_krgm_certificate(p)


def test_corrupt_hash_rejected(tmp_path):
    (tmp_path / "input").write_text("changed")
    (tmp_path / "SHA256SUMS").write_text("0" * 64 + "  input\n")
    assert verify_sha256s(tmp_path) == ["hash:input"]


def test_v12_sign():
    h = np.ones(9)
    v6 = np.ones(9)
    require_v12({"v6": v6, "v12": -v6 / 2}, h)
    with pytest.raises(ValueError):
        require_v12({"v6": v6, "v12": v6 / 2}, h)


def test_a5_requires_second_derivative():
    r = np.linspace(0.2, 2, 41)
    z = np.zeros_like(r)
    got = a5_holonomic(r, r**3, z, z, z, z)
    np.testing.assert_allclose(got, -6 * r, atol=1e-9)
    assert not np.allclose(got, z)


@pytest.mark.parametrize(
    "x,y",
    [(np.ones(9), np.ones(9)), (np.arange(9.0), np.full(9, np.nan)), (np.arange(3.0), np.ones(3))],
)
def test_jet_refuses_invalid_input(x, y):
    with pytest.raises(ValueError):
        derivative(x, y)


def test_schema_count():
    assert len(SLOT_NAMES) == len(set(SLOT_NAMES)) == 41
