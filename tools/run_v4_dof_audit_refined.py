#!/usr/bin/env python3
"""Refine the K_phys=0 finding: WHERE is the psi/V subblock semidefinite
with a zero mode — and is that zero a true zero mode (flat direction) or
a numerically tiny positive eigenvalue?

Also: compute the full 3-sector picture per node: the chi eigenvalue,
the psi/V subblock eigenvalues, and the sign pattern vs r. Declared
diagnostic refinement, gate unchanged: min lambda(K_phys) > 0 strictly.
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

OUT = ROOT / "data/generated/spectral/V4_PHYSICAL_DOF_AUDIT_V1_REFINED.json"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    sys.modules.setdefault(name, mod)
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    t0 = time.time()
    red9 = load("red9", str(ROOT / "src/ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py"))
    v4 = load("v4solve", str(ROOT / "tools/run_luminal_background_solve_v4.py"))
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

    res = red9.canonical_audit(sub, 6)
    K = res["K"]
    chi_eig = np.array([np.linalg.eigvalsh(K[i])[0] for i in range(K.shape[0])])
    phys_min, phys_max = [], []
    zero_nodes = []
    for i in range(K.shape[0]):
        sub2 = K[i][np.ix_([0, 2], [0, 2])]
        ev = np.linalg.eigvalsh(0.5 * (sub2 + sub2.T))
        phys_min.append(float(ev[0]))
        phys_max.append(float(ev[-1]))
        if abs(ev[0]) < 1e-8:
            zero_nodes.append(i)
    phys_min = np.array(phys_min)
    out = {
        "audit": "V4_PHYSICAL_DOF_AUDIT_V1_REFINED",
        "chi_channel": {
            "min_eig": float(chi_eig.min()),
            "max_eig": float(chi_eig.max()),
            "all_negative": bool((chi_eig < 0).all()),
            "interpretation": "chi is the constraint channel; its 'eigenvalue' "
                               "is a Lagrange-multiplier stiffness, not a "
                               "propagating speed",
        },
        "psi_V_subblock": {
            "min_over_nodes": float(phys_min.min()),
            "max_over_nodes": float(phys_min.max()),
            "n_nodes_with_|min|<1e-8": len(zero_nodes),
            "first_zero_r": round(float(r_w[zero_nodes[0]]), 4) if zero_nodes else None,
            "last_zero_r": round(float(r_w[zero_nodes[-1]]), 4) if zero_nodes else None,
            "n_strictly_negative": int((phys_min < -1e-8).sum()),
        },
    }
    # distribution of phys_min
    pm = phys_min
    for lo_, hi_ in [(0.0, 0.3), (0.3, 0.6), (0.6, 1.0), (1.0, 1.36)]:
        m = (r_w >= lo_) & (r_w < hi_)
        if m.any():
            out["psi_V_subblock"][f"r[{lo_},{hi_})"] = {
                "min": round(float(pm[m].min()), 8),
                "max": round(float(pm[m].max()), 8),
                "n": int(m.sum())}
    print(json.dumps(out["psi_V_subblock"], indent=1))

    # Is the zero a true zero mode? Check the determinant of the psi/V block:
    dets = np.array([np.linalg.det(K[i][np.ix_([0, 2], [0, 2])]) for i in range(K.shape[0])])
    out["psi_V_subblock"]["det_zero_nodes"] = int((np.abs(dets) < 1e-8).sum())
    out["psi_V_subblock"]["det_min_abs"] = float(np.min(np.abs(dets)))
    print("det zero-nodes:", out["psi_V_subblock"]["det_zero_nodes"],
          "| min|det|:", out["psi_V_subblock"]["det_min_abs"])

    out["wall_seconds"] = round(time.time() - t0, 1)
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
