#!/usr/bin/env python3
"""Fresh full-available-domain kinetic/radial audit of the current electric member.

This does NOT claim center-to-infinity/global closure.  It deliberately tests the
entire native domain exposed by build_onshell_central(), in addition to the
registered production window, so a local G40 PASS cannot be silently promoted
into a global spectral prerequisite.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.config import DEFAULT_L  # noqa: E402
from ssz_p5.numerics import module  # noqa: E402
from ssz_p5.production.electric_hybrid_onshell_central import (  # noqa: E402
    build_onshell_central,
)

MEMBER_CSV = ROOT / "data/generated/phase2_q2/ELECTRIC_PRODUCTION_MEMBER_CURRENT.csv"
MEMBER_JSON = ROOT / "data/generated/phase2_q2/ELECTRIC_PRODUCTION_MEMBER_CURRENT.json"
OUT = ROOT / "data/generated/spectral/CURRENT_MEMBER_FULL_AVAILABLE_DOMAIN_AUDIT.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def analyze_region(K: np.ndarray, G: np.ndarray, u: np.ndarray, mask: np.ndarray) -> dict:
    inds = np.flatnonzero(mask)
    if not len(inds):
        raise ValueError("empty audit region")

    Ks = (K + K.swapaxes(1, 2)) / 2.0
    Gs = (G + G.swapaxes(1, 2)) / 2.0
    ke = np.linalg.eigvalsh(Ks)
    row_min = ke[:, 0]

    local = row_min[inds]
    q = int(np.argmin(local))
    i = int(inds[q])

    cr2 = np.full(len(inds), np.nan, float)
    for j, idx in enumerate(inds):
        w, U = np.linalg.eigh(Ks[idx])
        if np.min(w) <= 0:
            continue
        invsqrt = U @ np.diag(1.0 / np.sqrt(w)) @ U.T
        C = invsqrt @ Gs[idx] @ invsqrt
        cr2[j] = np.linalg.eigvalsh((C + C.T) / 2.0)[0]

    positive = np.isfinite(cr2)
    return {
        "rows": int(len(inds)),
        "u_min": float(np.min(u[inds])),
        "u_max": float(np.max(u[inds])),
        "min_K": float(row_min[i]),
        "u_at_min_K": float(u[i]),
        "negative_K_rows": int(np.sum(local <= 0)),
        "first_negative_K_u": (
            float(u[inds[np.flatnonzero(local <= 0)[0]]])
            if np.any(local <= 0)
            else None
        ),
        "min_cr2_where_K_positive": (
            float(np.min(cr2[positive])) if np.any(positive) else None
        ),
        "negative_cr2_rows_where_defined": (
            int(np.sum(cr2[positive] <= 0)) if np.any(positive) else None
        ),
        "K_positive_everywhere": bool(np.all(local > 0)),
        "radial_positive_where_defined": bool(
            np.all(cr2[positive] > 0) if np.any(positive) else False
        ),
    }


def main() -> int:
    meta = json.loads(MEMBER_JSON.read_text())
    expected_hash = meta["member_hash"]
    actual_hash = sha256(MEMBER_CSV)
    if actual_hash != expected_hash:
        raise RuntimeError(
            f"member hash mismatch: expected {expected_hash}, got {actual_hash}"
        )

    build = build_onshell_central(ROOT)
    d = build.direct41.sort_values("x").reset_index(drop=True)
    u = d.u.to_numpy(float)

    production = (u > 0.62) & (u < 0.70)
    full = np.ones(len(d), dtype=bool)
    outer_tail = u >= 0.70
    inner_tail = u <= 0.62

    reducer = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    per_l = {}
    for L in DEFAULT_L:
        a = reducer.canonical_audit(d, int(L))
        K = np.asarray(a["K"], float)
        G = np.asarray(a["G"], float)
        per_l[str(L)] = {
            "production_window": analyze_region(K, G, u, production),
            "full_available_domain": analyze_region(K, G, u, full),
            "outer_tail_u_ge_0p70": (
                analyze_region(K, G, u, outer_tail) if np.any(outer_tail) else None
            ),
            "inner_tail_u_le_0p62": (
                analyze_region(K, G, u, inner_tail) if np.any(inner_tail) else None
            ),
            "max_K_asymmetry": float(np.max(np.abs(K - K.swapaxes(1, 2)))),
            "max_G_asymmetry": float(np.max(np.abs(G - G.swapaxes(1, 2)))),
        }

    production_pass = all(
        x["production_window"]["K_positive_everywhere"]
        and x["production_window"]["radial_positive_where_defined"]
        for x in per_l.values()
    )
    full_k_pass = all(
        x["full_available_domain"]["K_positive_everywhere"] for x in per_l.values()
    )
    full_radial_pass = all(
        x["full_available_domain"]["radial_positive_where_defined"]
        for x in per_l.values()
    )

    if not full_k_pass:
        verdict = "CURRENT_MEMBER_FULL_AVAILABLE_DOMAIN_K_FAIL"
        blocker = "KINETIC_GHOST_OUTSIDE_REGISTERED_PRODUCTION_WINDOW"
        downstream = "BLOCKED_BEFORE_GLOBAL_SPECTRAL_EXPORT"
    elif not full_radial_pass:
        verdict = "CURRENT_MEMBER_FULL_AVAILABLE_DOMAIN_RADIAL_FAIL"
        blocker = "RADIAL_GRADIENT_INSTABILITY_OUTSIDE_REGISTERED_PRODUCTION_WINDOW"
        downstream = "BLOCKED_BEFORE_GLOBAL_SPECTRAL_EXPORT"
    else:
        verdict = "CURRENT_MEMBER_FULL_AVAILABLE_DOMAIN_PASS"
        blocker = "MISSING_CENTER_TO_INFINITY_DIRECT_GLOBAL_KRGM"
        downstream = "LOCAL_NATIVE_DOMAIN_CLEAR_GLOBAL_EXPORT_STILL_REQUIRED"

    payload = {
        "member": meta.get("name"),
        "member_hash": expected_hash,
        "member_csv_sha256": actual_hash,
        "scope": (
            "Fresh canonical reducer audit on the entire available current-member "
            "domain only; this is not a center-to-infinity certificate."
        ),
        "available_domain": [float(np.min(u)), float(np.max(u))],
        "registered_production_window": [0.62, 0.70],
        "required_L": list(DEFAULT_L),
        "per_L": per_l,
        "registered_production_window_pass": production_pass,
        "full_available_domain_K_pass": full_k_pass,
        "full_available_domain_radial_pass": full_radial_pass,
        "verdict": verdict,
        "blocker": blocker,
        "downstream_spectral_semantics": downstream,
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    print(json.dumps(payload, indent=2, allow_nan=False))

    # CI success means the audit executed and reproduced the registered local
    # production-window PASS. A physical failure outside that window is a
    # scientific result, not a CI infrastructure failure.
    return 0 if production_pass else 5


if __name__ == "__main__":
    raise SystemExit(main())
