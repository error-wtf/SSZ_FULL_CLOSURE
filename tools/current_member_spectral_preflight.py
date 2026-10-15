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

        # Locate the first loss of a healthy radial principal sector when
        # moving inward (increasing u) from the certified production window.
        order = np.argsort(u)
        us = u[order]
        k0 = ke[order, 0]
        c0 = cr2[order]
        k_cross = None
        for q in range(1, len(us)):
            if k0[q - 1] > 0 and k0[q] <= 0:
                # Linear zero interpolation is diagnostic only.
                frac = k0[q - 1] / (k0[q - 1] - k0[q])
                k_cross = float(us[q - 1] + frac * (us[q] - us[q - 1]))
                break
        cr_cross = None
        last = None
        for q in range(len(us)):
            if not np.isfinite(c0[q]):
                continue
            if last is not None and c0[last] > 0 and c0[q] <= 0:
                frac = c0[last] / (c0[last] - c0[q])
                cr_cross = float(us[last] + frac * (us[q] - us[last]))
                break
            last = q
        ring_samples = {}
        for label, ur in {
            "outer_ring": 2.0 / 3.0,
            "inner_stable_ring": 0.7061346,
        }.items():
            nearest = int(np.argmin(np.abs(u - ur)))
            # Linear interpolation of the symmetric K matrix on the dense u-grid.
            su = np.argsort(u)
            Kui = np.empty((Ks.shape[1], Ks.shape[2]))
            Gui = np.empty((Gs.shape[1], Gs.shape[2]))
            for aa in range(Ks.shape[1]):
                for bb in range(Ks.shape[2]):
                    Kui[aa, bb] = np.interp(ur, u[su], Ks[su, aa, bb])
                    Gui[aa, bb] = np.interp(ur, u[su], Gs[su, aa, bb])
            kwe = np.linalg.eigvalsh((Kui + Kui.T) / 2.0)
            if np.min(kwe) > 0:
                w, U = np.linalg.eigh((Kui + Kui.T) / 2.0)
                invsqrt = U @ np.diag(1.0 / np.sqrt(w)) @ U.T
                C = invsqrt @ ((Gui + Gui.T) / 2.0) @ invsqrt
                cr_ring = float(np.linalg.eigvalsh((C + C.T) / 2.0)[0])
            else:
                cr_ring = None
            ring_samples[label] = {
                "u_target": float(ur),
                "nearest_grid_u": float(u[nearest]),
                "nearest_grid_delta_u": float(abs(u[nearest] - ur)),
                "min_eig_K_interpolated": float(kwe[0]),
                "mode_index_interpolated": int(np.argmin(kwe)),
                "min_cr2_if_K_positive": cr_ring,
            }

        per_L[str(int(L))] = {
            **scopes,
            "first_inward_K_zero_u_diagnostic": k_cross,
            "first_inward_cr2_zero_u_diagnostic": cr_cross,
            "ring_samples": ring_samples,
        }

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
