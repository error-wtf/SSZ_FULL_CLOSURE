#!/usr/bin/env python3
"""Current-member spectral preflight.

Rebuild the frozen electric production member from source, replay the certified
production-window principal gates, then extend the same reducer diagnostically
over the member's full stored u-domain. This does NOT manufacture a global
center-to-infinity certificate or a QNM spectrum.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ssz_p5.config import DEFAULT_L
from ssz_p5.numerics import module
from ssz_p5.production.electric_hybrid_onshell_central import (
    build_onshell_central,
    principal_audit,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/generated/spectral/CURRENT_MEMBER_SPECTRAL_PREFLIGHT_2026-09-29.json"


def generalized_min(Ks, Gs):
    vals = np.full(len(Ks), np.nan)
    for i in range(len(Ks)):
        w, U = np.linalg.eigh(Ks[i])
        if np.min(w) <= 0:
            continue
        invsqrt = U @ np.diag(1.0 / np.sqrt(w)) @ U.T
        C = invsqrt @ Gs[i] @ invsqrt
        vals[i] = np.linalg.eigvalsh((C + C.T) / 2.0)[0]
    return vals


def main() -> int:
    build, production_rows = principal_audit(ROOT)
    d = build.direct41.sort_values("x").reset_index(drop=True)
    red = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")

    u = d.u.to_numpy(float)
    masks = {
        "production_window": (u > 0.62) & (u < 0.70),
        "stored_full_domain": np.ones(len(d), dtype=bool),
        "inner_ring_band": (u >= 0.703) & (u <= 0.7095),
    }

    per_L = {}
    for L in DEFAULT_L:
        a = red.canonical_audit(d, int(L))
        K = np.asarray(a["K"], float)
        G = np.asarray(a["G"], float)
        Ks = (K + K.swapaxes(1, 2)) / 2.0
        Gs = (G + G.swapaxes(1, 2)) / 2.0
        ke = np.linalg.eigvalsh(Ks)
        cr2 = generalized_min(Ks, Gs)

        scopes = {}
        for name, mask in masks.items():
            inds = np.flatnonzero(mask)
            rowmin = ke[inds, 0]
            j = int(inds[np.argmin(rowmin)])
            finite_cr = cr2[inds][np.isfinite(cr2[inds])]
            scopes[name] = {
                "rows": int(len(inds)),
                "min_eig_K": float(ke[j, 0]),
                "u_at_min_K": float(u[j]),
                "mode_index_at_min_K": int(np.argmin(ke[j])),
                "negative_K_rows": int(np.sum(rowmin <= 0)),
                "min_cr2_where_K_positive": (
                    float(np.min(finite_cr)) if finite_cr.size else None
                ),
                "nonpositive_cr2_rows_where_defined": int(
                    np.sum(finite_cr <= 0)
                ) if finite_cr.size else None,
            }
        per_L[str(int(L))] = scopes

    payload = {
        "member": "ELECTRIC_PRODUCTION_MEMBER_CURRENT",
        "scope_note": (
            "production_window is certified usage; stored_full_domain and "
            "inner_ring_band are diagnostic extensions of the same rebuilt "
            "member/reducer and are NOT a center-to-infinity direct-global certificate"
        ),
        "u_domain": [float(u.min()), float(u.max())],
        "production_principal_audit": production_rows,
        "per_L": per_L,
        "production_all_pass": bool(all(r["pass"] for r in production_rows)),
        "stored_full_domain_all_K_positive": bool(
            all(v["stored_full_domain"]["negative_K_rows"] == 0 for v in per_L.values())
        ),
        "inner_ring_band_all_K_positive": bool(
            all(v["inner_ring_band"]["negative_K_rows"] == 0 for v in per_L.values())
        ),
        "direct_global_krgm_certificate_present": (
            ROOT / "data/certificates/DIRECT_GLOBAL_KRGM_CERTIFICATE.json"
        ).is_file(),
        "spectral_verdict": "PREFLIGHT_ONLY_NOT_GLOBAL_QNM",
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    print(json.dumps(payload, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
