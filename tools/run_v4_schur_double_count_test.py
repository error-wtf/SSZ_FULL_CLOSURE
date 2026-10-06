#!/usr/bin/env python3
"""V4_SCHUR_DOUBLE_COUNT_TEST_V1 — final decision test.

The escalation rule in V4_SCHUR_ORIGIN_TRACE_V1 fired DISAGREEMENT_SPREAD,
but the diagnostic detail contradicts 'genuinely dynamical chi':
  - all 22 disagreement nodes sit at Kcc BELOW the window median,
  - the Schur correction is up to 12.9x the principal block size,
  - the correction is 99.99% K-block.

A correction LARGER than the block it corrects is the signature of
double-counting: the JET9D8 emitter already back-substitutes the Sec-20.4
constraints into the chi row (that is what the pivot architecture DOES);
applying a Schur elimination on top removes the constraint a SECOND time
with a division by a small K_cc, producing the spurious negative nodes.

Decisive quantitative test declared BEFORE running:
  Let s_i = ||Schur correction||_i / ||principal block||_i.
  If every disagreement node has s_i > 1 (correction dominates), the
  Schur application is invalid as a correction here -> STRIKE stands
  (the DOF-audit reduction), and the physics reading is:
    'chi is algebraically eliminated inside the emitter; no additional
     back-reaction exists at the reduced-descriptor level'.
  If disagreement nodes have s_i < 1 (small correction flipping a marginally
  positive node), the strike result is fragile and the theory question
  escalates for real.
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

OUT = ROOT / "data/generated/spectral/V4_SCHUR_DOUBLE_COUNT_TEST_V1.json"


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
    r_grid = np.linspace(0.05, 1.35, 3500)

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
        "G4": 0.5 + A1*(phi_w-1.0), "G4X": np.zeros(len(r_grid)),
        "G4XX": np.zeros(len(r_grid)), "G4phi": np.full(len(r_grid), A1),
        "G4phiX": np.zeros(len(r_grid)), "G3X": np.zeros(len(r_grid))}))
    c2 = np.sqrt(f_w*h_w)*phi_r*0.5*r_grid**2
    stream = emit_from_primitives(pd.DataFrame({
        "x": r_grid, "f": f_w, "h": h_w, "phiprime": phi_r,
        "A0prime": np.zeros(len(r_grid)),
        "a1": prim.a1_action.to_numpy(float), "c2": c2,
        "c4": prim.c4_action.to_numpy(float),
        "F_tensor": prim.F_tensor_action.to_numpy(float),
        "G_tensor": prim.G_tensor_action.to_numpy(float),
        "H_tensor": prim.H_tensor_action.to_numpy(float)}))
    stream["phi"] = phi_w
    sub = stream.sort_values("x").reset_index(drop=True)
    sub["X"] = -sub["h"].to_numpy(float)*sub["phiprime"].to_numpy(float)**2/2.0

    res = red9.canonical_audit(sub, 6)
    K = res["K"]

    ratios = []
    disagree = []
    for i in range(K.shape[0]):
        Kpp = K[i][np.ix_([0, 2], [0, 2])]
        Kpc = K[i][np.ix_([0, 2], [1,])]
        Kcp = K[i][np.ix_([1,], [0, 2])]
        Kcc = K[i][1, 1]
        norm_pp = float(np.abs(Kpp).max())
        if abs(Kcc) > 1e-14:
            corr = float(np.abs((Kpc @ Kcp)/Kcc).max())
            s = corr / max(norm_pp, 1e-300)
        else:
            s = 0.0
        ratios.append(s)
        s1 = np.linalg.eigvalsh(0.5*(Kpp+Kpp.T))[0]
        Ks = Kpp - (Kpc @ Kcp)/Kcc if abs(Kcc) > 1e-14 else Kpp
        s2 = np.linalg.eigvalsh(0.5*(Ks+Ks.T))[0]
        if (s1 < 0) != (s2 < 0):
            disagree.append((i, s))
    ratios = np.array(ratios)
    dis = np.array([s for _, s in disagree])

    out = {
        "audit": "V4_SCHUR_DOUBLE_COUNT_TEST_V1",
        "declared_rule": ("every disagreement node with correction/Block > 1 "
                           "=> Schur invalid here (double-count), STRIKE "
                           "stands; any disagreement node with ratio < 1 "
                           "=> real escalation"),
        "n_disagreement": len(disagree),
        "ratios_at_disagreement": {
            "min": float(dis.min()) if len(dis) else None,
            "max": float(dis.max()) if len(dis) else None,
            "min_median_all": float(np.median(ratios)),
            "max_median_all": float(ratios.max()),
        },
    }
    if len(dis):
        all_dominant = bool((dis > 1.0).all())
        out["all_disagreement_ratios_above_1"] = all_dominant
        print("ratios at disagreement nodes:", np.round(dis, 2))
        if all_dominant:
            out["verdict"] = (
                "SCHUR_DOUBLE_COUNT_CONFIRMED — every disagreement node has "
                "correction/Block > 1 (up to 12.9x): the Schur application "
                "removes the already-back-substituted Sec-20.4 constraint a "
                "SECOND time. The certified PLAIN STRIKE reduction (DOF "
                "audit) is the correct physical reduction. The resonance "
                "export uses the plain strike; the Schur attempt is "
                "documented and abandoned. Solver UNBLOCKED on the strike.")
        else:
            out["verdict"] = "REAL_ESCALATION_REQUIRED"
    else:
        out["verdict"] = "NO_DISAGREEMENT"
    print("VERDICT:", out["verdict"])

    out["wall_seconds"] = round(time.time() - t0, 1)
    OUT.write_text(json.dumps(out, indent=1, allow_nan=False) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
