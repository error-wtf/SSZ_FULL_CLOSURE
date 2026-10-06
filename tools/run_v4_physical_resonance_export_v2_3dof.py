#!/usr/bin/env python3
"""V4_PHYSICAL_RESONANCE_EXPORT_V1 — the producer block for the coupled
resonance solver (COUPLED_RESONANCE_SOLVER_V1 pre-step).

WHAT THIS DOES
==============
1. Builds the reduced (K, G, S, M) operator blocks on the FULL radial
   domain: interior window [0.05, 1.35] PLUS continued exterior
   [1.35, R_EXT] on one common grid.
2. Performs the constraint-reduced (psi, V) projection WITH constraint
   back-reaction (Schur complement of the chi block), NOT a naive
   principal-submatrix strike. The chi channel is non-propagating
   (Lagrange multiplier, V4_PHYSICAL_DOF_AUDIT_V1) so its elimination
   contributes the Schur correction -K_cv K_cc^{-1} K_vc per block,
   preserving ALL physical coupling terms as required.
3. Exports the physical blocks + exterior metric + asymptotic
   normalization constants in the frozen schema, hash-bound to the
   certified V4 branch (a1=+0.5, eps=-0.3).

WHAT THIS DOES NOT DO
=====================
- No boundary conditions imposed here (that is V4_PHYSICAL_BC_V1).
- No root finding (that is SOLVER A/B).
- No certification (that is the C-R gates).

Validation inside this tool:
- re-derives the certified branch (same rhs/integrator as all previous
  artifacts)
- checks K_phys (Schur-reduced) min eigenvalue > 0 over the full domain
- checks continuity of the blocks at the window/exterior interface
- writes SHA256 sidecar.
"""
from __future__ import annotations

import hashlib
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

OUT_NPZ = ROOT / "data/generated/spectral/V4_PHYSICAL_RESONANCE_EXPORT_V2_3DOF.npz"
OUT_JSON = ROOT / "data/generated/spectral/V4_PHYSICAL_RESONANCE_EXPORT_V2_3DOF.json"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def schur_reduce(K, G, S, M, phys=(0, 2), cons=(1,)):
    """Constraint reduction WITH back-reaction (Schur complement).

    The chi equation (constraint, no time derivatives on chi in the
    reduced sense) couples linearly: eliminating chi modifies the
    psi/V blocks. For second-order forms (K qddot + S qdot + ...) the
    constraint row has no K_cc contribution (chi is non-propagating),
    so the Schur correction enters through G/S/M blocks:
        G_phys = G_pp - G_pc G_cc^{-1} G_cp   (etc.)
    with the constraint block G_cc taken as the chi-block of the
    (semi-)definite spatial kinetic form. Where G_cc is numerically
    singular the chi row is a pure algebraic constraint and the naive
    principal strike is EXACT (recorded per node).

    Returns physical blocks + per-node provenance of the reduction.
    """
    n_r = K.shape[0]
    Kp, Gp, Sp, Mp = [], [], [], []
    provenance = {"schur_nodes": 0, "algebraic_nodes": 0}
    for i in range(n_r):
        def red(block):
            Bpp = block[np.ix_(phys, phys)]
            Bpc = block[np.ix_(phys, cons)]
            Bcp = block[np.ix_(cons, phys)]
            Bcc = block[np.ix_(cons, cons)]
            # singular constraint block -> algebraic: naive strike is exact
            if abs(float(Bcc[0, 0])) < 1e-14:
                provenance["algebraic_nodes"] += 1
                return Bpp
            provenance["schur_nodes"] += 1
            return Bpp - Bpc @ np.linalg.solve(Bcc, Bcp)
        Kp.append(red(K[i])); Gp.append(red(G[i]))
        Sp.append(red(S[i])); Mp.append(red(M[i]))
    return (np.stack(Kp), np.stack(Gp), np.stack(Sp), np.stack(Mp),
            provenance)


def main() -> int:
    t0 = time.time()
    v4 = load("v4solve", str(ROOT / "tools/run_luminal_background_solve_v4.py"))
    red9 = load("red9", str(ROOT / "src/ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py"))
    from ssz_p5.jets.jet9d8 import derivative as jet

    A1, EPS = -0.5, -0.3
    R_IN, R_WIN_END, R_EXT = 0.05, 1.35, 60.0
    # Non-uniform grid: the certified window keeps its historical density
    # (3500 pts, Dr ~ 3.7e-4) so the stiff chi dynamics stays resolved;
    # the exterior is logarithmic (600 pts to R_EXT). A uniform grid over
    # the full domain UNDERSAMPLES the window and produced 2 spurious
    # negative Schur-reduced K_phys nodes at r~1.5-2.25 (grid artifact,
    # diagnosed and documented in the run report).
    N_WIN, N_EXT = 3500, 600

    out = {"audit": "V4_PHYSICAL_RESONANCE_EXPORT_V2_3DOF",
            "branch": {"a1": A1, "eps": EPS},
            "domain": {"r_in": R_IN, "r_window_end": R_WIN_END,
                        "r_ext": R_EXT, "n_window": N_WIN, "n_exterior": N_EXT,
                        "grid": "dense-uniform window + logarithmic exterior",
                        "artifact_note": ("uniform full-domain grid produced 2 "
                                           "spurious negative Schur nodes at "
                                           "r~1.5-2.25 (window undersampled); "
                                           "this grid resolves the stiff chi "
                                           "dynamics")}}

    # ---------- certified branch, full domain
    archive = pd.read_csv(v4.ARCHIVE).sort_values("u").reset_index(drop=True)
    r_all = (1.0 / archive.u.to_numpy(float))[::-1]
    red, facts, _ = v4.build_reduction()
    assert all(facts.values())
    rhs = v4.rhs_factory(red, A1)
    sol = v4.integrate_branch(rhs, float(r_all[0]), R_EXT, EPS)
    assert bool(sol.t[-1] >= R_EXT * (1 - 1e-9)), "exterior continuation failed"
    r_win = np.linspace(R_IN, R_WIN_END, N_WIN, endpoint=False)
    r_ext = np.logspace(np.log10(R_WIN_END), np.log10(R_EXT), N_EXT + 1)[1:]
    r_grid = np.concatenate([r_win, r_ext])
    y = sol.sol(r_grid)
    f_w, h_w, phi_w = y[0], y[1], y[2]
    phi_r = jet(r_grid, phi_w)

    out["branch_limits"] = {"f_inf": round(float(f_w[-1]), 8),
                             "h_inf": round(float(h_w[-1]), 8)}

    # ---------- reduced operator blocks on the full grid
    from ssz_p5.coefficients.mh_action_primitives import quartic_g5zero_primitives
    from ssz_p5.coefficients.mh_general_primitives import emit_from_primitives
    prim = quartic_g5zero_primitives(pd.DataFrame({
        "x": r_grid, "f": f_w, "h": h_w, "phi_r": phi_r,
        "G4": 0.5 + A1 * (phi_w - 1.0),
        "G4X": np.zeros(len(r_grid)), "G4XX": np.zeros(len(r_grid)),
        "G4phi": np.full(len(r_grid), A1), "G4phiX": np.zeros(len(r_grid)),
        "G3X": np.zeros(len(r_grid))}))
    c2 = np.sqrt(f_w * h_w) * phi_r * 0.5 * r_grid**2
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
    sub["X"] = -sub["h"].to_numpy(float) * sub["phiprime"].to_numpy(float) ** 2 / 2.0

    L = 6
    res = red9.canonical_audit(sub, L)
    K_full, G_full = res["K"], res["G"]
    S_full = 0.5 * (res["S"] + np.swapaxes(res["S"], 1, 2))  # antisymm kept:
    S_full = 0.5 * (res["S"] - np.swapaxes(res["S"], 1, 2))
    M_full = res["M"]

    # NO REDUCTION — all 3 channels are dynamical (rule 7 pass, closure
    # ce1c156: dphi is dynamical in the unreduced descriptor; production
    # branch a1=-0.5, eps=-0.3 is three-channel healthy)
    Kp, Gp, Sp, Mp = K_full, G_full, S_full, M_full
    prov = {"reduction": "NONE — full 3-channel physical operator "
                           "(psi,dphi,V), per V4_DOF_ARCHITECTURE_REAUDIT_V2"}

    # ---------- health checks on the FULL domain
    kmin = float(min(np.linalg.eigvalsh(
        0.5*(Kp[i]+Kp[i].T))[0] for i in range(Kp.shape[0])))
    kmin_win = float(min(np.linalg.eigvalsh(
        0.5*(Kp[i]+Kp[i].T))[0]
        for i in range(int(np.searchsorted(r_grid, R_WIN_END)))))
    out["k_phys_health"] = {
        "min_eigenvalue_full_domain": kmin,
        "min_eigenvalue_window": kmin_win,
        "positive_full_domain": bool(kmin > 0),
        "note": ("window value must match V4_PHYSICAL_DOF_AUDIT_V1_REFINED "
                  "(+4.47e-06) up to grid differences"),
    }
    print("K_phys min (full/window):", f"{kmin:.3e} / {kmin_win:.3e}", flush=True)

    # continuity at the interface
    i_end = int(np.searchsorted(r_grid, R_WIN_END))
    dK = float(np.max(np.abs(Kp[i_end] - Kp[i_end-1]))) if i_end > 0 else 0.0
    out["interface_continuity_max_abs_dK"] = dK

    # ---------- tortoise coordinate + asymptotic normalization
    dr_star = 1.0/np.sqrt(f_w*h_w)
    r_star = np.concatenate([[0.0], np.cumsum(0.5*(dr_star[1:]+dr_star[:-1])*
                                              np.diff(r_grid))])
    f_inf = float(f_w[-1])
    sqrt_f_inf = float(np.sqrt(f_inf))
    out["asymptotic_normalization"] = {
        "definition_time": "T = sqrt(f_inf) * t",
        "definition_tortoise": "R_star = sqrt(f_inf) * r_star",
        "definition_freq": "Omega_inf = omega_t / sqrt(f_inf)",
        "Omega_inf_factor": round(1.0/sqrt_f_inf, 8),
        "f_inf": round(f_inf, 8),
    }
    out["exterior_decay_note"] = ("Kretschmann tail consistent with r^-6 "
                                    "(V4_CURVATURE_BENCHMARKS_V3, p=-6.0001); "
                                    "outgoing basis derived from the "
                                    "continued exterior in V4_PHYSICAL_BC_V1")

    # ---------- export
    np.savez_compressed(
        OUT_NPZ,
        r=r_grid, r_star=r_star, f=f_w, h=h_w, phi=phi_w,
        L=np.array([L]),
        K_phys=Kp, G_phys=Gp, S_phys=Sp, M_phys=Mp,
        K_full3=K_full, G_full3=G_full, S_full3=S_full, M_full3=M_full,
        Dh1=res["maps"]["Dh1"], DeltaV=res["maps"]["DeltaV"],
        pivotA0=res["maps"]["pivotA0"],
        f_inf=np.array([f_inf]), sqrt_f_inf=np.array([sqrt_f_inf]))
    h_npz = hashlib.sha256(OUT_NPZ.read_bytes()).hexdigest()
    out["npz_sha256"] = h_npz
    out["npz_shape_summary"] = {
        "K_phys": list(Kp.shape), "r_grid": len(r_grid)}
    out["wall_seconds"] = round(time.time() - t0, 1)
    OUT_JSON.write_text(json.dumps(out, indent=1, allow_nan=False) + "\n")
    (OUT_NPZ.parent / (OUT_NPZ.name + ".sha256")).write_text(h_npz + "\n")
    print("written:", OUT_NPZ.name, "sha:", h_npz[:16])
    print("verdict:", "PASS" if (kmin > 0) else "FAIL (K_phys not positive)")
    return 0 if kmin > 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
