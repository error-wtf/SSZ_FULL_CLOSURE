import json
import shutil

import pytest

from ssz_p5.config import DEFAULT_L, repo_root
from ssz_p5.provenance.manifest import sha256
from ssz_p5.qnm.gate import REQUIRED_GATES, require_direct_krgm_certificate


def synthetic_certificate(tmp_path):
    # Synthetic fixture tests the guard, never a scientific certificate.
    (tmp_path / "schemas").mkdir()
    shutil.copyfile(
        repo_root() / "schemas/direct_global_krgm_certificate.schema.json",
        tmp_path / "schemas/direct_global_krgm_certificate.schema.json",
    )
    for name in ("PRODUCTION_BLACKLIST.json", "SSZ_P5_REGIONAL_PRODUCTION_MEMBER_2026-09-17.json"):
        shutil.copyfile(repo_root() / name, tmp_path / name)
    artifact = tmp_path / "synthetic.csv"
    artifact.write_text("test fixture only\n")
    data = dict(
        pass_=True,
        release="2026-09-16",
        action_sha256=sha256(tmp_path / "SSZ_P5_REGIONAL_PRODUCTION_MEMBER_2026-09-17.json"),
        coefficient_stream_sha256=sha256(artifact),
        L_values=list(DEFAULT_L),
        gates=[
            dict(name=n, status="PASS", criterion="fixture", evidence="fixture")
            for n in sorted(REQUIRED_GATES)
        ],
        artifacts=[dict(path="synthetic.csv", sha256=sha256(artifact))],
    )
    data["pass"] = data.pop("pass_")
    path = tmp_path / "certificate.json"
    path.write_text(json.dumps(data))
    return path, data


def test_schema_and_hash_guard_accepts_complete_synthetic_fixture(tmp_path):
    path, _ = synthetic_certificate(tmp_path)
    assert require_direct_krgm_certificate(path, tmp_path)["pass"] is True


@pytest.mark.parametrize("failure", ["artifact_hash", "action_hash", "gate", "multipole", "escape"])
def test_schema_valid_forgery_is_rejected(tmp_path, failure):
    path, data = synthetic_certificate(tmp_path)
    if failure == "artifact_hash":
        (tmp_path / "synthetic.csv").write_text("tampered")
    elif failure == "action_hash":
        data["action_sha256"] = "0" * 64
    elif failure == "gate":
        data["gates"].pop()
    elif failure == "multipole":
        data["L_values"].remove(1000)
    elif failure == "escape":
        data["artifacts"][0]["path"] = "../outside.csv"
    path.write_text(json.dumps(data))
    with pytest.raises(RuntimeError):
        require_direct_krgm_certificate(path, tmp_path)
