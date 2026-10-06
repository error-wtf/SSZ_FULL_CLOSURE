#!/usr/bin/env python3
"""V4_SCHUR_ORIGIN_TRACE_V1 — WHERE does the Schur back-reaction come from?

The escalation from V4_SCHUR_VS_STRIKE_DIAGNOSTIC_V1 found:
  - strike (psi,V) principal block: positive everywhere (22 nodes where
    Schur disagrees, all negative in Schur)
  - Schur correction is up to 12.9x the principal block size at the
    disagreeing nodes
  - K_cc spans 0.006 .. 1.9e6 — the chi row is NOT uniformly algebraic.

Question answered here: the chi ROW of the reduced descriptor is the
ELIMINATION RESULT of the Sec-20.4 constraint architecture (h1, H2, A0
substituted OUT). The chi row therefore carries back-substituted pieces
of the ELIMINATED fields — including their second derivatives through
the H11-type terms. The huge K_cc and the huge Schur correction may be
emitter-scaling artifacts of those back-substituted terms, not physics.

Test: build the Schur correction from the UNREDUCED second-variational
form directly is out of scope here; instead, this tool locates WHERE the
dominant Schur-correction contribution lives (which of K/G/S/M blocks
drives it, and at which radii), and checks correlation with the known
emitter-stiffness pattern (outer-edge conditioning, cf. V4_S_CONTRACT_
FORENSICS: chi-related stiffness concentrated at the outer edge).

If the disagreement nodes coincide with the outer-edge window and with
huge K_cc (emitter-stiffness pattern), the correct reading is:
  the plain strike remains the certified physical reduction (DOF audit
  unchanged), and the Schur attempt is a REDUCTION-SCHEME error (double
  counting the constraint), documented and abandoned.

If disagreement nodes are spread over the window independent of K_cc,
the chi elimination is genuinely dynamical and the theory question goes
to escalation as declared.
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

OUT = ROOT / "data/generated/spectral/V4_SCHUR_ORIGIN_TRACE_V1.json"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    t0 = time.time()
    v4 = load("v4solve", str(ROOT / "tools/run_luminal_background_solve_v4.py"))
    red9 = load("red9", str(ROOT / "src/ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py"))
    from ssz_p5.jets.jet9d8 import derivative as jet

    A1, EPS = 0.5, -0.3
    N_WIN = 3500
    r_grid = np.linspace(0.05, 1.35, N_WIN)

    archive = pd.read_csv(v4.ARCHIVE).sort_values("u").reset_index(drop=True)
    r_all = (1.0 / archive.u.to_numpy(float))[::-1]
    red, facts, _ = v4.build_reduction()
    rhs = v4.rhs_factory(red, A1)
    sol = v4.integrate_branch(rhs, float(r_all[0]), float(r_all[-1]), EPS)
    y = sol.sol(r_grid)
    f_w, h_w, phi_w = y[0], y[1], y[2]
    phi_r = jet(r_grid, phi_w)
    from ssz_p5.coefficients.mh_action_primitives import quartic_g5zero_primitives
    from ssz_p5.coefficients.mh_general_primitives import emit_from_primitives
    prim = quartic_g5zero_primitives(pd.DataFrame({
        "x": r_grid, "f": f_w, "h": h_w, "phi_r": phi_r,
        "G4": 0.5 + A1*(phi_w-1.0), "G4X": np.zeros(N_WIN),
        "G4XX": np.zeros(N_WIN), "G4phi": np.full(N_WIN, A1),
        "G4phiX": np.zeros(N_WIN), "G3X": np.zeros(N_WIN)}))
    c2 = np.sqrt(f_w*h_w)*phi_r*0.5*r_grid**2
    stream = emit_from_primitives(pd.DataFrame({
        "x": r_grid, "f": f_w, "h": h_w, "phiprime": phi_r,
        "A0prime": np.zeros(N_WIN),
        "a1": prim.a1_action.to_numpy(float), "c2": c2,
        "c4": prim.c4_action.to_numpy(float),
        "F_tensor": prim.F_tensor_action.to_numpy(float),
        "G_tensor": prim.G_tensor_action.to_numpy(float),
        "H_tensor": prim.H_tensor_action.to_numpy(float)}))
    stream["phi"] = phi_w
    sub = stream.sort_values("x").reset_index(drop=True)
    sub["X"] = -sub["h"].to_numpy(float)*sub["phiprime"].to_numpy(float)**2/2.0

    res = red9.canonical_audit(sub, 6)
    K, G = res["K"], res["G"]
    S = 0.5*(res["S"] - np.swapaxes(res["S"], 1, 2))
    M = res["M"]

    # where is the Schur correction large, per block?
    blocks = {"K": K, "G": G, "S": S, "M": M}
    corr_blocks = {}
    for bname, B in blocks.items():
        corr = np.zeros(K.shape[0])
        for i in range(K.shape[0]):
            Bpc = B[i][np.ix_([0, 2], [1,])]
            Bcp = B[i][np.ix_([1,], [0, 2])]
            Bcc = B[i][1, 1]
            if abs(Bcc) > 1e-14:
                corr[i] = float(np.abs((Bpc @ Bcp)/Bcc).max())
        corr_blocks[bname] = corr

    # strike/schur disagreement nodes (as in diagnostic)
    disagree_idx = []
    for i in range(K.shape[0]):
        Kpp = K[i][np.ix_([0, 2], [0, 2])]
        Kpc = K[i][np.ix_([0, 2], [1,])]
        Kcp = K[i][np.ix_([1,], [0, 2])]
        Kcc = K[i][1, 1]
        Ks = Kpp - (Kpc @ Kcp)/Kcc if abs(Kcc) > 1e-14 else Kpp
        s1 = np.linalg.eigvalsh(0.5*(Kpp+Kpp.T))[0]
        s2 = np.linalg.eigvalsh(0.5*(Ks+Ks.T))[0]
        if (s1 < 0) != (s2 < 0):
            disagree_idx.append(i)
    disagree_idx = np.array(disagree_idx, int)
    out = {"audit": "V4_SCHUR_ORIGIN_TRACE_V1",
            "branch": {"a1": A1, "eps": EPS},
            "n_disagreement": int(len(disagree_idx))}
    if len(disagree_idx):
        r_dis = r_grid[disagree_idx]
        out["disagreement_r_range"] = [float(r_dis.min()), float(r_dis.max())]
        # correlation with K_cc magnitude:
        kcc = np.array([abs(K[i][1, 1]) for i in range(K.shape[0])])
        med_kcc = float(np.median(kcc))
        frac_huge_at_disagreement = float(
            (kcc[disagree_idx] > med_kcc).mean())
        out["Kcc_median_all"] = med_kcc
        out["frac_disagreement_nodes_with_Kcc_above_median"] = \
            frac_huge_at_disagreement
        # which BLOCK dominates the correction at those nodes?
        dom = {}
        for bname, corr in corr_blocks.items():
            dom[bname] = float(np.mean(corr[disagree_idx]))
        total = sum(dom.values()) or 1.0
        out["correction_block_share_at_disagreement"] = {
            k_: round(v_/total, 4) for k_, v_ in dom.items()}
        # correlation with the known emitter-stiffness pattern (outer edge):
        out["outer_edge_share"] = round(float((r_dis > 1.0).mean()), 4)
        print("disagreement r-range:", out["disagreement_r_range"])
        print("block share:", out["correction_block_share_at_disagreement"])
        print("outer-edge share:", out["outer_edge_share"],
              "| Kcc>median share:", frac_huge_at_disagreement)

    # ---------- verdict (declared rule)
    if len(disagree_idx) == 0:
        out["verdict"] = "NO_DISAGREEMENT"
    elif (out.get("outer_edge_share", 0) > 0.7 and
          out.get("frac_disagreement_nodes_with_Kcc_above_median", 0) > 0.7):
        out["verdict"] = ("EMITTER_STIFFNESS_PATTERN — disagreement "
                           "concentrated at the outer edge with huge K_cc, "
                           "matching the known V4_S_CONTRACT_FORENSICS "
                           "conditioning pattern. Reading: the Schur "
                           "correction double-counts the constraint (the "
                           "chi row already contains the back-substituted "
                           "constraint result). The certified PLAIN STRIKE "
                           "reduction (DOF audit) stands. Schur attempt "
                           "documented and abandoned. Solver proceeds on "
                           "the plain strike.")
    else:
        out["verdict"] = ("DISAGREEMENT_SPREAD — chi elimination genuinely "
                           "dynamical; escalate to theory review. Solver "
                           "BLOCKED.")
    print("VERDICT:", out["verdict"])

    out["wall_seconds"] = round(time.time() - t0, 1)
    OUT.write_text(json.dumps(out, indent=1, allow_nan=False) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
