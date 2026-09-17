"""Fail-closed validation of direct products and same-operator spectral evidence.

Archive witnesses and a syntactically valid certificate alone cannot pass this gate.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

from ..coefficients.schema41 import validate_41_schema
from ..config import DEFAULT_L, SLOT_NAMES
from ..geometry.p5 import from_frame
from ..policy import numerical_policy
from ..production.member import validate_stream_regions
from ..provenance.blacklist import require_production_input
from ..provenance.manifest import sha256
from ..qnm.gate import require_direct_krgm_certificate
from ..reducer.canonical import reduce_profile
from ..stability.finite_l import stability_diagnostics
from ..types import Coefficients41, ConstraintPivots, ReducedOperator

DIRECT = "data/certificates/DIRECT_GLOBAL_KRGM_CERTIFICATE.json"
SPECTRAL = "data/certificates/COUPLED_SPECTRAL_CERTIFICATE.json"
ABSOLUTE = "ABSOLUTE_FULL_CLOSURE.json"
CONVERGENCE_AXES = {
    "radial_resolution",
    "outer_domain",
    "asymptotic_order",
    "root_tolerance",
    "basis_normalization",
    "branch_continuation",
}
ABSOLUTE_GATES = {
    "p5_geometry",
    "action_member",
    "global_direct_41of41",
    "constraint_reduction",
    "direct_global_krgsm",
    "finite_l_stability",
    "spectral_solver",
    "qnm_convergence",
    "absolute_full_closure",
}


def bound_file(root, record):
    path = require_production_input(root / record["path"], root)
    if sha256(path) != record["sha256"]:
        raise ValueError(f"artifact hash mismatch: {record['path']}")
    return path


def verify_direct_products(root: Path):
    root = root.resolve()
    path = root / DIRECT
    cert = require_direct_krgm_certificate(path, root)
    records = cert["artifacts"]
    sources = [r for r in records if r["sha256"] == cert["coefficient_stream_sha256"]]
    if len(sources) != 1:
        raise ValueError("one unambiguous direct coefficient stream is required")
    frame = pd.read_csv(bound_file(root, sources[0]))
    provenance = {"region", "patch_id", "source_action_id", "source_sha256", "generator_version"}
    if not provenance.issubset(frame.columns) or frame[list(provenance)].isna().any().any():
        raise ValueError("direct stream lacks row-level action provenance")
    bg = from_frame(frame.rename(columns={"r": "x"}) if "x" not in frame else frame)
    c = Coefficients41(bg, {s: frame[s].to_numpy(float) for s in SLOT_NAMES}, sha256(path))
    if any(g.status != "PASS" for g in validate_41_schema(c)):
        raise ValueError("invalid direct 41-slot stream")
    validate_stream_regions(frame)
    operators = {}
    for record in records:
        if not record["path"].endswith(".npz"):
            continue
        with np.load(bound_file(root, record), allow_pickle=False) as data:
            L = int(data["L"].item())
            if L in operators:
                raise ValueError(f"duplicate direct operator L={L}")
            np.testing.assert_array_equal(data["r"], bg.r)
            pivots = ConstraintPivots(*(data[k] for k in ("Dh1", "DeltaV", "pivotA0")))
            op = ReducedOperator(
                L,
                bg.r,
                *(data[k] for k in ("K", "R", "G", "S", "M")),
                pivots,
                cert["coefficient_stream_sha256"],
            )
            replay = reduce_profile(c, L)
            tolerance = numerical_policy()["matrix_symmetry_scaled"]
            for key in ("K", "R", "G", "S", "M"):
                expected = getattr(replay, key)
                actual = getattr(op, key)
                error = np.abs(actual - expected) / np.maximum(1.0, np.abs(expected))
                if np.max(error) > tolerance:
                    raise ValueError(
                        f"exported {key} does not reproduce from the 41 stream at L={L}"
                    )
            result = stability_diagnostics(op)
            if not result["pass"]:
                raise ValueError(f"finite-L stability failed: {result}")
            operators[L] = record["sha256"]
    if not set(DEFAULT_L).issubset(operators):
        raise ValueError("missing required direct finite-L matrix products")
    return cert, operators


def verify_spectral_evidence(root, direct_hash, operators):
    path = root / SPECTRAL
    cert = json.loads(path.read_text())
    if cert.get("pass") is not True or cert.get("direct_certificate_sha256") != direct_hash:
        raise ValueError("spectral result is not bound to the passing direct operator")
    if cert.get("operator_sha256_by_L") != {str(k): v for k, v in operators.items()}:
        raise ValueError("QNM must use the same direct matrices at every L")
    if cert.get("method") not in {"compactified_jost", "ecs", "continued_fraction"}:
        raise ValueError("unaccepted resonance method")
    table = pd.read_csv(bound_file(root, cert["convergence_artifact"]))
    required = {
        "branch",
        "L",
        "axis",
        "setting",
        "omega_re",
        "omega_im",
        "residual",
        "difference",
        "tolerance",
    }
    if not required.issubset(table.columns) or table.empty or table.isna().any().any():
        raise ValueError("missing convergence evidence")
    numeric = table[["omega_re", "omega_im", "residual", "difference", "tolerance"]].to_numpy(float)
    if not np.isfinite(numeric).all() or (table.tolerance <= 0).any():
        raise ValueError("invalid convergence numerics")
    if (table.residual < 0).any() or (table.difference < 0).any():
        raise ValueError("negative error measure")
    if (table.difference > table.tolerance).any() or (table.residual > table.tolerance).any():
        raise ValueError("spectral convergence failed")
    axes = CONVERGENCE_AXES | ({"ecs_angle"} if cert["method"] == "ecs" else set())
    for _, group in table.groupby(["L", "branch"]):
        if not axes.issubset(set(group.axis)):
            raise ValueError("missing branch convergence axis")
        if any(group.loc[group.axis == a, "setting"].nunique() < 2 for a in axes):
            raise ValueError("one numerical setting cannot establish convergence")
    scan = cert["unstable_mode_search"]
    for key in ("real_interval", "imaginary_interval"):
        bounds = np.asarray(scan[key], float)
        if bounds.shape != (2,) or not np.isfinite(bounds).all() or bounds[1] <= bounds[0]:
            raise ValueError("invalid spectral search domain")
    if scan["imaginary_interval"][1] <= 0:
        raise ValueError("no search for Im(omega)>0")
    resolution = np.asarray(scan["resolution"], float)
    if resolution.shape != (2,) or not np.isfinite(resolution).all() or (resolution < 2).any():
        raise ValueError("missing spectral search resolution")
    bound_file(root, scan["artifact"])
    return cert


def verify_absolute_closure(root: Path):
    root = root.resolve()
    _, operators = verify_direct_products(root)
    direct_hash = sha256(root / DIRECT)
    verify_spectral_evidence(root, direct_hash, operators)
    cert = json.loads((root / ABSOLUTE).read_text())
    if any(cert.get(k) != "PASS" for k in ABSOLUTE_GATES):
        raise ValueError("a mandatory absolute-closure gate is absent or not PASS")
    if cert.get("direct_certificate_sha256") != direct_hash:
        raise ValueError("absolute certificate has stale direct operator")
    if cert.get("spectral_certificate_sha256") != sha256(root / SPECTRAL):
        raise ValueError("absolute certificate has stale spectral evidence")
    if not (root / "ABSOLUTE_FULL_CLOSURE.md").is_file():
        raise ValueError("missing readable absolute certificate")
    return cert
