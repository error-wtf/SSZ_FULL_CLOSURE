#!/usr/bin/env python3
"""Rebuild the current electric member and audit all available spectral layers."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from ssz_p5.numerics import module
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central
from ssz_p5.qnm.principal_tracking import analyze_principal_tracking

ROOT = Path(__file__).resolve().parents[1]
MEMBER_CSV = ROOT / "data/generated/phase2_q2/ELECTRIC_PRODUCTION_MEMBER_CURRENT.csv"
MEMBER_MANIFEST = ROOT / "data/generated/phase2_q2/ELECTRIC_PRODUCTION_MEMBER_CURRENT.json"
DEFAULT_OUT = ROOT / "data/generated/spectral_selection/FULL_AVAILABLE_SPECTRAL_AUDIT.json"
DEFAULT_L = (6, 12, 20, 42, 110, 420, 1000)
J_P5 = 0.006114458428


def frozen_member_hash() -> str:
    return hashlib.sha256(MEMBER_CSV.read_bytes()).hexdigest()


def builder_replay(build, atol: float = 1e-10):
    """Compare the numerical builder to the hash-pinned stream.

    The release policy itself uses an absolute numerical tolerance for the
    regenerated floating-point stream; the immutable CSV bytes carry the
    member hash.  Requiring a newly integrated ODE trajectory to reproduce the
    CSV bit-for-bit would be stronger than the repository's own G05/G20 gate
    and is not scientifically justified.
    """
    frozen = pd.read_csv(MEMBER_CSV, index_col=0).drop(columns=["A0"])
    fresh = (
        build.action[frozen.columns]
        .sort_values("u")
        .reset_index(drop=True)
    )
    if len(fresh) != len(frozen):
        raise RuntimeError("current member builder row count mismatch")

    per_column = {}
    worst = (0.0, None)
    for col in frozen.columns:
        a = fresh[col].to_numpy(float)
        b = frozen[col].to_numpy(float)
        delta = float(np.max(np.abs(a - b)))
        per_column[col] = delta
        if delta > worst[0]:
            worst = (delta, col)
    return {
        "atol": atol,
        "max_abs_delta": worst[0],
        "worst_column": worst[1],
        "pass": bool(worst[0] <= atol),
        "per_column_max_abs_delta": per_column,
    }

def maxwell_threshold():
    ell = 1
    while math.sqrt(ell * (ell + 1)) * J_P5 < math.pi / 2:
        ell += 1
    return {
        "J_P5": J_P5,
        "threshold_ell": ell,
        "I_at_threshold": math.sqrt(ell * (ell + 1)) * J_P5,
        "I_previous": math.sqrt((ell - 1) * ell) * J_P5,
        "pi_over_2": math.pi / 2,
    }


def audit_current_member():
    manifest = json.loads(MEMBER_MANIFEST.read_text())
    frozen_hash = frozen_member_hash()
    if frozen_hash != manifest["member_hash"]:
        raise RuntimeError(
            f"frozen member hash mismatch: {frozen_hash} != {manifest['member_hash']}"
        )

    build = build_onshell_central(ROOT)
    replay = builder_replay(build)
    if not replay["pass"]:
        raise RuntimeError(
            "current member builder exceeds release replay tolerance: "
            f"{replay['max_abs_delta']} > {replay['atol']} "
            f"in {replay['worst_column']}"
        )
    rows = len(build.action)
    if rows != manifest["rows"]:
        raise RuntimeError("current member row count mismatch")

    d = build.direct41.copy().sort_values("x").reset_index(drop=True)
    u = d.u.to_numpy(float)
    prod = (u > 0.62) & (u < 0.70)
    if int(prod.sum()) != int(manifest["domain"]["production_rows"]):
        raise RuntimeError("production-window row count mismatch")

    reducer = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    per_l = {}
    for L in DEFAULT_L:
        op = reducer.canonical_audit(d, L)
        K = np.asarray(op["K"], float)[prod]
        G = np.asarray(op["G"], float)[prod]
        scan = analyze_principal_tracking(K, G)
        eig_k = np.linalg.eigvalsh((K + np.swapaxes(K, 1, 2)) / 2)
        per_l[str(L)] = {
            "min_eig_K": float(np.min(eig_k)),
            "min_cr2": float(np.min(scan["values"])),
            "min_adjacent_overlap": scan["min_adjacent_overlap"],
            "min_adjacent_overlap_per_branch": scan["min_adjacent_overlap_per_branch"],
            "naive_sort_relabel_steps": scan["naive_sort_relabel_steps"],
            "dominant_component_changes": scan["dominant_component_changes"],
            "endpoint_abs_overlap_matrix": scan["endpoint_abs_overlap_matrix"].tolist(),
            "principal_tracking_status": (
                "SMOOTH"
                if scan["min_adjacent_overlap"] > 0.99
                else "DISCONTINUITY_DETECTED"
            ),
        }

    return {
        "member_hash_expected": manifest["member_hash"],
        "frozen_csv_sha256": frozen_hash,
        "frozen_hash_verification": "PASS",
        "builder_numerical_replay": replay,
        "member_rows": rows,
        "production_window": manifest["domain"]["production_window"],
        "production_rows": int(prod.sum()),
        "per_L": per_l,
        "maxwell_testfield": maxwell_threshold(),
    }


def global_spectral_status():
    direct_paths = [
        ROOT / "data/certificates/DIRECT_GLOBAL_KRGM_CERTIFICATE.json",
        ROOT / "data/certificates/SSZ_P5_DIRECT_GLOBAL_KRGM_CERTIFICATE.json",
    ]
    present = [str(p.relative_to(ROOT)) for p in direct_paths if p.is_file()]
    return {
        "direct_global_certificate_present": present,
        "status": (
            "READY_FOR_GLOBAL_SPECTRAL_SOLVE"
            if present
            else "SPECTRAL_SELECTION_NOT_YET_EVALUABLE"
        ),
        "blocker": (
            None
            if present
            else "MISSING_CURRENT_MEMBER_CENTER_TO_INFINITY_DIRECT_KRGM"
        ),
        "note": (
            "The frozen current electric member covers u in [0.61,0.71]. "
            "No current-member center-to-infinity direct K/R/G/S/M certificate "
            "is exposed by main. Local principal analysis is therefore not a "
            "substitute for global modal residues."
        ),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, default=DEFAULT_OUT)
    ns = ap.parse_args()

    result = {
        "schema_version": "1.0",
        "current_member": audit_current_member(),
        "global_spectral": global_spectral_status(),
    }
    ns.output.parent.mkdir(parents=True, exist_ok=True)
    ns.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
