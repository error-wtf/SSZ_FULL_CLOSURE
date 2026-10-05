#!/usr/bin/env python3
"""Track eigenvalue branches ACROSS L and certify L-invariant families.

The local principal solver G(r,L) v = lambda(r,L) K(r,L) v is run at
several L for the same radial window.  Branches whose relative frequency
spread over L is below the declared threshold are recorded as
L-INVARIANT CANDIDATE FAMILY (a structural-selection candidate — NOT a
QNM claim; the global operator later reproduces or falsifies them).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ssz_p5.config import DEFAULT_L  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def sym(a: np.ndarray) -> np.ndarray:
    """Symmetrize a (rows, m, m) stack or a single (m, m) matrix."""
    if a.ndim == 3:
        return 0.5 * (a + np.swapaxes(a, 1, 2))
    return 0.5 * (a + a.T)


def solve_modes(K: np.ndarray, G: np.ndarray):
    import scipy.linalg as sla
    vals, vecs = sla.eigh(sym(G), sym(K))
    return vals, vecs, sym(K)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--u-min", type=float, default=0.62)
    ap.add_argument("--u-max", type=float, default=0.70)
    ap.add_argument("--n-u", type=int, default=4)
    ap.add_argument("--spread-threshold", type=float, default=1e-9,
                    help="max relative omega spread over L for invariance")
    ap.add_argument("--output", type=Path,
                    default=ROOT / "data/generated/spectral/MODE_FAMILY_INVARIANCE.json")
    args = ap.parse_args()

    from ssz_p5.numerics import module  # noqa: E402
    from ssz_p5.production.electric_hybrid_onshell_central import (  # noqa: E402
        build_onshell_central,
    )

    build = build_onshell_central(ROOT)
    stream = build.direct41.sort_values("x").reset_index(drop=True)
    reducer = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    uall = stream.u.to_numpy(float)
    ids = np.flatnonzero((uall > args.u_min) & (uall < args.u_max))
    sel = ids[np.unique(np.linspace(0, len(ids) - 1,
                                    args.n_u).round().astype(int))]
    us = uall[sel]          # ACTUAL grid points, not a synthetic linspace
    idx = [int(np.argmin(np.abs(uall - u))) for u in us]

    per_L = {}
    for L in DEFAULT_L:
        audit = reducer.canonical_audit(stream, int(L))
        Kall = sym(np.asarray(audit["K"], float))
        Gall = sym(np.asarray(audit["G"], float))
        omegas_at_us = []
        for i in idx:
            vals, _, _ = solve_modes(Kall[i], Gall[i])
            omegas_at_us.append(np.sqrt(np.abs(vals)))
        u_used = [float(uall[i]) for i in idx]
        per_L[int(L)] = {"u_used": u_used,
                         "omega_rows": omegas_at_us}

    # match branches across L by frequency proximity at each u
    ref_L = int(DEFAULT_L[0])
    ref = np.asarray(per_L[ref_L]["omega_rows"])
    families = []
    # Branch matching across L by frequency proximity at each u (eigh
    # ordering can permute between L; nearest-frequency is the tracker).
    seed = np.asarray(per_L[int(DEFAULT_L[0])]["omega_rows"])
    for branch in range(seed.shape[1]):
        omega_means, spread_rows = [], []
        for ui in range(seed.shape[0]):
            seed_om = float(seed[ui][branch])
            col = []
            for L in DEFAULT_L:
                om = np.asarray(per_L[int(L)]["omega_rows"][ui], float)
                j = int(np.argmin(np.abs(om - seed_om)))
                col.append(float(om[j]))
            col = np.asarray(col, float)
            spread = float((col.max() - col.min()) / max(abs(col.mean()), 1e-300))
            spread_rows.append(spread)
            omega_means.append(float(col.mean()))
        # per-u families: a branch is L-invariant AT a radius if that
        # radius's spread passes.  Branch identity is tracked by the
        # frequency VALUE (fingerprint), not by eigh ordering — ordering
        # permutes across L and u.  Each surviving (branch, u) is one
        # L-invariant family member.
        for ui, spread in enumerate(spread_rows):
            if spread <= args.spread_threshold:
                families.append({
                    "branch": int(branch),
                    "u": float(us[ui]),
                    "omega_L_mean": omega_means[ui],
                    "rel_spread": spread,
                })

    result = {
        "audit": "MODE_FAMILY_INVARIANCE_V1",
        "claim_level": ("L-invariant local principal family candidate — "
                        "NOT a QNM claim; global operator decides later"),
        "spread_threshold": args.spread_threshold,
        "L_values": [int(L) for L in DEFAULT_L],
        "n_families_invariant": len(families),
        "families": families,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=1) + "\n")
    print(f"written: {args.output}")
    ref0 = np.asarray(per_L[int(DEFAULT_L[0])]["omega_rows"])
    all_spreads = []
    for branch in range(ref0.shape[1]):
        col = np.asarray([per_L[int(L)]["omega_rows"][0][branch]
                          for L in DEFAULT_L], float)
        all_spreads.append(float((col.max() - col.min())
                                 / max(abs(col.mean()), 1e-300)))
    if all_spreads:
        best = int(np.argmin(all_spreads))
        print(f"best branch: {best}, spread {all_spreads[best]:.3e}")
        print(f"spread distribution: min={min(all_spreads):.2e} "
              f"median={float(np.median(all_spreads)):.2e} "
              f"max={max(all_spreads):.2e}")
    print(f"L-invariant (branch, u) family members: {len(families)}")
    for f in families:
        print(f"  branch {f['branch']:>3} @ u={f['u']:.6f}: "
              f"omega={f['omega_L_mean']:.9f} spread={f['rel_spread']:.2e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
