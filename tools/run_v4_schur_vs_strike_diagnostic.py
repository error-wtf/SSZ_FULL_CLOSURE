#!/usr/bin/env python3
"""V4_SCHUR_VS_STRIKE_DIAGNOSTIC_V1 — is the Schur back-reaction of the
chi constraint physically required, or does it corrupt the (psi,V) blocks?

Background: V4_PHYSICAL_DOF_AUDIT_V1 certified K_phys > 0 using the plain
principal (psi,V) subblock. The new Schur-reduced export (with chi
back-reaction -K_pc K_cc^-1 K_cp) turns K_phys NEGATIVE on parts of the
window. One of the two reductions is wrong for this structure.

This diagnostic answers, per node in the window:
  1. Is the chi equation algebraic (K_cc ~ 0) or dynamical (K_cc large)?
  2. How large is the Schur correction relative to the principal block?
  3. Does the chi row actually contain chi second-derivatives (i.e. is
     chi a real DOF that the DOF audit misclassified), or is the large
     K_cc an artifact of the emitter's coefficient scaling?

Decision output:
  SCHUR_REQUIRED        -> DOF audit needs revision, solver uses Schur
  ALGEBRAIC_STRIKE_OK   -> naive strike is exact, export uses it
  AMBIGUOUS             -> stop, escalate to theory review
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

OUT = ROOT / "data/generated/spectral/V4_SCHUR_VS_STRIKE_DIAGNOSTIC_V1.json"


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

    # ---------- per-node diagnostics in the window
    kcc_over_r = []
    schur_rel = []
    strike_min = []
    schur_min = []
    for i in range(K.shape[0]):
        Kcc = K[i][1, 1]
        kcc_over_r.append(float(abs(Kcc)))
        # Schur correction on K:
        Kpc = K[i][np.ix_([0, 2], [1,])]
        Kcp = K[i][np.ix_([1,], [0, 2])]
        Kpp = K[i][np.ix_([0, 2], [0, 2])]
        if abs(Kcc) > 1e-14:
            corr = (Kpc @ Kcp) / Kcc
            Ks = Kpp - corr
            schur_rel.append(float(np.abs(corr).max() /
                                    max(np.abs(Kpp).max(), 1e-300)))
        else:
            Ks = Kpp
            schur_rel.append(0.0)
        strike_min.append(float(np.linalg.eigvalsh(0.5*(Kpp+Kpp.T))[0]))
        schur_min.append(float(np.linalg.eigvalsh(0.5*(Ks+Ks.T))[0]))
    kcc_over_r = np.array(kcc_over_r)
    schur_rel = np.array(schur_rel)
    strike_min = np.array(strike_min)
    schur_min = np.array(schur_min)

    out = {
        "audit": "V4_SCHUR_VS_STRIKE_DIAGNOSTIC_V1",
        "branch": {"a1": A1, "eps": EPS},
        "window": [0.05, 1.35], "n_nodes": int(K.shape[0]),
        "Kcc_absolute": {
            "min": float(kcc_over_r.min()),
            "max": float(kcc_over_r.max()),
            "median": float(np.median(kcc_over_r)),
            "interpretation": ("K_cc spans orders of magnitude; if it is "
                                " LARGE the chi row carries chi-second-"
                                "derivative terms and chi is NOT purely "
                                "algebraic in this reduced form"),
        },
        "schur_correction_relative_max": float(schur_rel.max()),
        "schur_correction_relative_median": float(np.median(schur_rel)),
        "strike_min_eigenvalue": {
            "min": float(strike_min.min()),
            "n_negative": int((strike_min < 0).sum()),
        },
        "schur_min_eigenvalue": {
            "min": float(schur_min.min()),
            "n_negative": int((schur_min < 0).sum()),
        },
    }
    print(json.dumps({k: out[k] for k in
                       ("Kcc_absolute", "schur_correction_relative_max",
                        "strike_min_eigenvalue", "schur_min_eigenvalue")},
                       indent=1))

    # ---------- correlation: where strike and schur disagree, how big is Kcc?
    disagree = (strike_min < 0) != (schur_min < 0)
    out["disagreement_nodes"] = int(disagree.sum())
    if disagree.any():
        out["disagreement_Kcc_range"] = [
            float(kcc_over_r[disagree].min()),
            float(kcc_over_r[disagree].max())]

    # ---------- verdict
    # If Kcc is systematically huge (>> principal block), the chi row is
    # dynamically stiff, the Schur back-reaction is REQUIRED, and the
    # negative Schur eigenvalues mean the (psi,V) subspace alone is not
    # a closed positive system -> escalate.
    # If Kcc is tiny, strike is exact and the export should use strike.
    ratio_huge = kcc_over_r.max() / max(np.abs(K[[0, 2]][0][0]).max(), 1e-300)
    out["Kcc_to_principal_scale_ratio"] = float(ratio_huge)
    if out["disagreement_nodes"] == 0:
        out["verdict"] = "STRIKE_AND_SCHUR_AGREE"
    elif out["schur_min_eigenvalue"]["n_negative"] > 0 and \
            out["strike_min_eigenvalue"]["n_negative"] == 0:
        out["verdict"] = "SCHUR_CORRUPTS_STRIKE_HEALTHY — AMBIGUOUS, escalate"
        out["escalation_note"] = (
            "The DOF audit certified the principal (psi,V) strike as "
            "positive; the Schur back-reaction turns it negative. Either "
            "the chi elimination is algebraic (strike exact, Schur "
            "double-counts the constraint) or the chi channel is dynamical "
            "(the DOF audit's sector attribution needs revision). This is "
            "a THEORY question, not a numerics question. Solver BLOCKED "
            "until resolved.")
    print("VERDICT:", out["verdict"])

    out["wall_seconds"] = round(time.time() - t0, 1)
    OUT.write_text(json.dumps(out, indent=1, allow_nan=False) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
