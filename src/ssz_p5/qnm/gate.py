"""A claimed PASS cannot enable QNM without schema, gate and file provenance checks."""

import json
from pathlib import Path

import jsonschema

from ..config import DEFAULT_L, repo_root
from ..provenance.blacklist import require_production_input
from ..provenance.manifest import sha256

REQUIRED_GATES = {
    "input_provenance",
    "background",
    "slots41",
    "a5",
    "v12",
    "auxiliary",
    "epsilon_y",
    "constraints",
    "operator_structure",
    "kinetic",
    "radial",
    "angular",
    "pure_limits",
    "center",
    "high_L_mass",
    "interfaces",
    "jet_convergence",
}


def require_direct_krgm_certificate(path: Path, root: Path | None = None) -> dict:
    root = (root or repo_root()).resolve()
    try:
        data = json.loads(path.read_text())
        schema = json.loads(
            (root / "schemas/direct_global_krgm_certificate.schema.json").read_text()
        )
        jsonschema.validate(data, schema)
        action = root / "SSZ_P5_HSVT_ACTION_MEMBER_2026-09-16.json"
        if data["action_sha256"] != sha256(action):
            raise ValueError("action hash mismatch")
        if not set(DEFAULT_L).issubset(data["L_values"]):
            raise ValueError("missing required multipoles")
        names = {g["name"] for g in data["gates"]}
        if not REQUIRED_GATES.issubset(names):
            raise ValueError("missing required direct-export gates")
        hashes = set()
        for rec in data["artifacts"]:
            file = require_production_input(root / rec["path"], root)
            actual = sha256(file)
            if actual != rec["sha256"]:
                raise ValueError("artifact hash mismatch")
            hashes.add(actual)
        if data["coefficient_stream_sha256"] not in hashes:
            raise ValueError("coefficient stream not bound to an artifact")
        return data
    except (OSError, ValueError, KeyError, TypeError, jsonschema.ValidationError) as exc:
        raise RuntimeError(
            f"QNM disabled: invalid or absent direct-global-KRGM certificate: {exc}"
        ) from exc
