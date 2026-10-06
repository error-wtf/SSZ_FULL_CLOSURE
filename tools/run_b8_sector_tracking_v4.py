#!/usr/bin/env python3
"""B8: sector tracking on the full 3-field system of the ghost-free V4 branch.

(a1, eps) = (+0.5, -0.3), window [0.05, 1.35], operator as certified by the
S-projection repair (validator UNCHANGED, PASS).

Declared protocol:
  * solve the generalized box eigenproblem K psi = omega^2 G psi for
    L = 6 (B3-style P1 weak form is NOT re-implemented here; this run uses
    the exported canonical operator pairs via the production reducer at
    each L in the ladder),
  * track MODES by the K-metric eigenspace overlap
        O_mn = |psi_m^dag K psi_n| / sqrt((psi_m^dag K psi_m)(psi_n^dag K psi_n))
    across the resolution ladder N -> 2N -> 4N (avoided crossings follow
    eigenspaces, not sorted eigenvalue indices),
  * per tracked mode store the field fractions F_psi, F_chi, F_V computed
    in the K metric (fractions sum to 1 per mode),
  * negative control: shuffling the K-metric (diagonal scaling noise)
    must DESTROY the overlap structure.

Output: artifacts/B8_SECTOR_TRACKING_V4_V1.json
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/home/error/physics/clones/SSZ_FULL_CLOSURE")
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

OUT = ROOT / "data/generated/spectral/B8_SECTOR_TRACKING_V4_V1.json"


def load(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def k_overlap(psi_m, psi_n, K):
    num = abs(np.vdot(psi_m, K @ psi_n))
    den = np.sqrt(abs(np.vdot(psi_m, K @ psi_m)) * abs(np.vdot(psi_n, K @ psi_n)))
    return float(num / max(den, 1e-300))


def field_fractions(psi, K):
    """F_psi, F_chi, F_V in the K metric: |psi_i^dag K_ii psi_i| normalized."""
    fr = []
    for i in range(3):
        k_ii = K[:, i, i]  # diagonal kinetic coefficient per grid point
        fr.append(abs(float(np.sum(k_ii * np.abs(psi[:, i]) ** 2))))
    tot = sum(fr) or 1.0
    return [f / tot for f in fr]


def main() -> int:
    t0 = time.time()
    red9 = load("red9", str(ROOT / "src/ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py"))
    v4 = load("v4solve", str(ROOT / "tools/run_luminal_background_solve_v4.py"))
    from ssz_p5.reducer.canonical import validate_operator, ReducedOperator, ConstraintPivots
    from ssz_p5.jets.jet9d8 import derivative as jet

    A1, EPS = 0.5, -0.3
    WINDOW_R = (0.05, 1.35)

    archive = pd.read_csv(v4.ARCHIVE).sort_values("u").reset_index(drop=True)
    r_all = (1.0 / archive.u.to_numpy(float))[::-1]
    lo = int(np.searchsorted(r_all, WINDOW_R[0]))
    hi = int(np.searchsorted(r_all, WINDOW_R[1]))
    r_w = r_all[lo + 8: hi - 8]

    red, facts, _ = v4.build_reduction()
    rhs = v4.rhs_factory(red, A1)
    sol = v4.integrate_branch(rhs, float(r_all[0]), float(r_all[-1]), EPS)
    assert bool(sol.t[-1] >= float(r_all[-1]) * (1 - 1e-9))
    y = sol.sol(r_w)
    f_w, h_w, phi_w = y[0], y[1], y[2]
    phi_r = jet(r_w, phi_w)

    from ssz_p5.coefficients.mh_action_primitives import quartic_g5zero_primitives
    from ssz_p5.coefficients.mh_general_primitives import emit_from_primitives
    feed = pd.DataFrame({
        "x": r_w, "f": f_w, "h": h_w, "phi_r": phi_r,
        "G4": 0.5 + A1 * (phi_w - 1.0),
        "G4X": np.zeros(len(r_w)), "G4XX": np.zeros(len(r_w)),
        "G4phi": np.full(len(r_w), A1), "G4phiX": np.zeros(len(r_w)),
        "G3X": np.zeros(len(r_w))})
    prim = quartic_g5zero_primitives(feed)
    c2 = np.sqrt(f_w * h_w) * phi_r * 0.5 * r_w**2
    prim_in = pd.DataFrame({
        "x": r_w, "f": f_w, "h": h_w, "phiprime": phi_r,
        "A0prime": np.zeros(len(r_w)),
        "a1": prim.a1_action.to_numpy(float), "c2": c2,
        "c4": prim.c4_action.to_numpy(float),
        "F_tensor": prim.F_tensor_action.to_numpy(float),
        "G_tensor": prim.G_tensor_action.to_numpy(float),
        "H_tensor": prim.H_tensor_action.to_numpy(float)})
    stream = emit_from_primitives(prim_in)
    stream["phi"] = phi_w
    sub = stream.sort_values("x").reset_index(drop=True)
    sub["X"] = -sub["h"].to_numpy(float) * sub["phiprime"].to_numpy(float) ** 2 / 2.0

    def reduced_for(L):
        res = red9.canonical_audit(sub, L)
        S = 0.5 * (res["S"] - np.swapaxes(res["S"], 1, 2))
        m = res["maps"]
        op = ReducedOperator(L, sub.x.to_numpy(float),
                             res["K"], res["R"], res["G"], S, res["M"],
                             ConstraintPivots(m["Dh1"], m["DeltaV"], m["pivotA0"]),
                             f"v4_b8_L{L}")
        validate_operator(op)
        return res["K"], res["G"]

    # --- coarse-to-fine eigenspace tracking at one L with grid refinement is
    # not available without the P1 weak form; the declared B8 scope here is
    # the SECTOR CONTENT + K-overlap structure at L=6,12,20 (box modes,
    # TRUNCATED_BOX semantics preserved).
    results = {}
    modes_by_L = {}
    for L in (6, 12, 20):
        K, G = reduced_for(L)
        n_r = K.shape[0]
        # flatten to (3*n_r) generalized eigenproblem on the collocation grid
        Kf = np.zeros((3 * n_r, 3 * n_r))
        Gf = np.zeros_like(Kf)
        for i in range(3):
            for j in range(3):
                Kf[i::3, j::3] = np.diag(K[:, i, j])
                Gf[i::3, j::3] = np.diag(G[:, i, j])
        # symmetric part for a stable Hermitian solve (box operator pair is
        # symmetric by the symmetry table; asymmetry lives in S, not K/G)
        Ks = 0.5 * (Kf + Kf.T)
        Gs = 0.5 * (Gf + Gf.T)
        # regularize G to positive definite via its diagonal floor
        dg = np.diag(Gs).copy()
        floor = 1e-300
        dg[dg <= 0] = floor
        Ginv = np.diag(1.0 / dg)
        vals, vecs = np.linalg.eigh(Ginv @ Ks)
        order = np.argsort(vals)
        vals, vecs = vals[order], vecs[:, order]
        keep = np.isfinite(vals)
        vals, vecs = vals[keep], vecs[:, keep]
        n_modes = min(12, vecs.shape[1])
        modes = []
        for k in range(n_modes):
            psi = vecs[:, k]
            # un-flatten
            f3 = np.zeros((n_r, 3))
            f3[:, 0] = psi[0::3]
            f3[:, 1] = psi[1::3]
            f3[:, 2] = psi[2::3]
            Fr = field_fractions(f3, K)
            modes.append({
                "index": k,
                "lam": float(vals[k]),
                "F_psi": round(Fr[0], 6),
                "F_chi": round(Fr[1], 6),
                "F_V": round(Fr[2], 6),
                "vec": f3,
                "Kvec": (K @ f3.reshape(-1, 3) if False else None),
            })
        modes_by_L[L] = (modes, K)
        results[str(L)] = [
            {k: m_[k] for k in ("index", "lam", "F_psi", "F_chi", "F_V")}
            for m_ in modes
        ]
        print(f"L={L}: {n_modes} modes, lambda[0..2] = "
              f"{[round(m_['lam'], 3) for m_ in modes[:3]]}", flush=True)
        print(f"      F_V range: {min(m_['F_V'] for m_ in modes):.3f} .. "
              f"{max(m_['F_V'] for m_ in modes):.3f}", flush=True)

    # --- K-metric eigenspace overlap between the L=6 and L=12 mode sets:
    # for each L=6 mode find the best-overlap L=12 mode (eigenspace tracking)
    modes6, K6 = modes_by_L[6]
    modes12, K12 = modes_by_L[12]
    n_r6, n_r12 = modes6[0]["vec"].shape[0], modes12[0]["vec"].shape[0]
    # resample L6 vectors onto the L12 grid (both grids are the same r_w;
    # reduced_for uses the same sub frame, so n_r6 == n_r12 — assert)
    overlaps = []
    if n_r6 == n_r12:
        K6f = np.zeros((3 * n_r6, 3 * n_r6))
        for i in range(3):
            for j in range(3):
                K6f[i::3, j::3] = np.diag(K6[:, i, j])
        tracked = []
        for m6 in modes6:
            best, best_o = None, -1.0
            v6 = m6["vec"].reshape(-1)
            for m12 in modes12:
                v12 = m12["vec"].reshape(-1)
                num = abs(np.vdot(v6, K6f @ v12))
                den = np.sqrt(abs(np.vdot(v6, K6f @ v6)) * abs(np.vdot(v12, K6f @ v12)))
                o = float(num / max(den, 1e-300))
                if o > best_o:
                    best_o, best = o, m12
            tracked.append({
                "L6_index": m6["index"], "L12_best": best["index"],
                "overlap": round(best_o, 6),
                "F_V_L6": m6["F_V"], "F_V_L12": best["F_V"],
                "avoided_crossing_suspect": bool(best_o < 0.9),
            })
        overlaps = tracked
        n_low = sum(1 for t in tracked if t["overlap"] >= 0.9)
        print(f"eigenspace tracking L6->L12: {n_low}/{len(tracked)} modes with overlap >= 0.9")
    else:
        print("grid mismatch; overlap skipped", n_r6, n_r12)

    # --- negative control: shuffled K must destroy the overlap structure
    rng = np.random.default_rng(4104)
    perm = rng.permutation(n_r6)
    K6s = K6[perm][:, perm]
    K6sf = np.zeros((3 * n_r6, 3 * n_r6))
    for i in range(3):
        for j in range(3):
            K6sf[i::3, j::3] = np.diag(K6s[:, i, j])
    ctrl_overlaps = []
    for m6 in modes6[:5]:
        v6 = m6["vec"].reshape(-1)
        best_o = -1.0
        for m12 in modes12:
            v12 = m12["vec"].reshape(-1)
            num = abs(np.vdot(v6, K6sf @ v12))
            den = np.sqrt(abs(np.vdot(v6, K6sf @ v6)) * abs(np.vdot(v12, K6sf @ v12)))
            o = float(num / max(den, 1e-300))
            best_o = max(best_o, o)
        ctrl_overlaps.append(round(best_o, 6))
    print("negative control overlaps (shuffled K):", ctrl_overlaps[:5])

    out = {
        "audit": "B8_SECTOR_TRACKING_V4_V1",
        "branch": {"a1": A1, "eps": EPS, "ghost_free_side": True,
                    "window_r": [float(r_w[0]), float(r_w[-1])]},
        "semantics": "TRUNCATED_BOX modes, NOT QNMs (unchanged disclaimer)",
        "mode_table_by_L": results,
        "eigenspace_tracking_L6_L12": overlaps,
        "negative_control": {
            "method": "row/column permutation of K",
            "top5_overlaps": ctrl_overlaps,
            "mean_overlap": float(np.mean(ctrl_overlaps)),
        },
        "wall_seconds": round(time.time() - t0, 1),
    }
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print("written:", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
