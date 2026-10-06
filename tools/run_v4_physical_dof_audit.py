#!/usr/bin/env python3
"""V4_PHYSICAL_DOF_AUDIT_V1 — the declared fork: dynamical vs constraint.

Question (Lino, 06.10.): Is the negative 3x3 kinetic channel on the
ghost-free V4 branch a PHYSICAL degree of freedom (→ branch carries a
bulk ghost → branch dead) or a CONSTRAINT/gauge artifact (→ the branch
survives after correct constraint-reduced operator)?

Outputs (declared before the run, per Lino's spec):
  1. constraint rank versus r, L      (rank K_node, rank G_node per node)
  2. physical DOF count               (from pivot architecture + ranks)
  3. projection of v_- onto physical/constraint subspaces
  4. physical-subspace K_phys eigenvalues

Gate (kompromisslos): min lambda(K_phys) > 0  → branch alive,
                       min lambda(K_phys) <= 0 → branch dead (bulk ghost).

Constraint architecture (Sec. 20.4 / constraint-map module):
  h1  eliminated by pivot Dh1     (h1 = H10 y + H11 y_r)
  H2  eliminated by pivot DeltaV  (H2 = H20 y + H21 y_r)
  A0  eliminated by pivot pivotA0 (A0 = A00 y + A01 y_r)
  → physical fields y = (psi, dphi, V); the reduced (psi,chi,V) rows of the
    exported operator carry these plus their constraint back-reactions.
  A NEGATIVE direction that is dominated by the CONSTRAINT sectors
  (dphi/V rows which encode h1/H2/A0 elimination) is an artifact; a
  negative direction dominated by the TRUE scalar sector (psi row) with
  no pivot violation is a physical ghost.

Method for (3): at each node, decompose the negative eigenvector into
sector components and measure the fraction of |v|² carried by the psi
row vs the constraint-eliminated rows — plus, decisively, check whether
the negative direction survives when the constraint back-reaction rows
are projected out (psi-row-only subblock eigenvalue).
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

OUT = ROOT / "data/generated/spectral/V4_PHYSICAL_DOF_AUDIT_V1.json"
GATE = "min lambda(K_phys) > 0"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    t0 = time.time()
    red9 = load("red9", str(ROOT / "src/ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py"))
    v4 = load("v4solve", str(ROOT / "tools/run_luminal_background_solve_v4.py"))
    from ssz_p5.jets.jet9d8 import derivative as jet

    A1, EPS = 0.5, -0.3
    WINDOW_R = (0.05, 1.35)

    # ---------- branch + stream (certified chain)
    archive = pd.read_csv(v4.ARCHIVE).sort_values("u").reset_index(drop=True)
    r_all = (1.0 / archive.u.to_numpy(float))[::-1]
    lo = int(np.searchsorted(r_all, WINDOW_R[0]))
    hi = int(np.searchsorted(r_all, WINDOW_R[1]))
    r_w = r_all[lo + 8: hi - 8]

    red, facts, _ = v4.build_reduction()
    assert all(facts.values())
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

    # ---------- 1. constraint ranks vs r, L + pivot values
    out = {"audit": "V4_PHYSICAL_DOF_AUDIT_V1",
            "branch": {"a1": A1, "eps": EPS, "ghost_free_side": True,
                        "window_r": [float(r_w[0]), float(r_w[-1])]},
            "gate": GATE,
            "constraint_architecture": {
                "h1_by": "pivot Dh1", "H2_by": "pivot DeltaV",
                "A0_by": "pivot pivotA0",
                "physical_fields": ["psi", "dphi", "V"]}}

    # fm.maps on the emitted stream → pivots per node
    from ssz_p5.paths import paths as Bpaths
    fm_spec = importlib.util.spec_from_file_location(
        "constraint_maps",
        str(ROOT / "src/ssz_hybrid_full_constraint_maps_JET9D8(1).py"))
    fm = load("constraint_maps", str(ROOT / "src/ssz_hybrid_full_constraint_maps_JET9D8(1).py"))
    maps = fm.maps(sub, 6)
    pivots = {}
    for name in ("Dh1", "DeltaV", "pivotA0"):
        v = np.asarray(maps[name], dtype=float)
        pivots[name] = {
            "min_abs": float(np.min(np.abs(v))),
            "n_below_1e-10": int((np.abs(v) < 1e-10).sum()),
            "rows": int(len(v)),
        }
    out["pivots"] = pivots
    print("pivots:", json.dumps(pivots, indent=1), flush=True)
    pivot_violation = any(p_["n_below_1e-10"] > 0 for p_ in pivots.values())

    # ---------- 1b/2. ranks + DOF count per L
    dof_by_L = {}
    for L in (6, 12, 20):
        res = red9.canonical_audit(sub, L)
        K, G = res["K"], res["G"]
        n = K.shape[0]
        rankK = np.array([np.linalg.matrix_rank(K[i], tol=1e-10) for i in range(n)])
        rankG = np.array([np.linalg.matrix_rank(G[i], tol=1e-10) for i in range(n)])
        dof_by_L[f"L{L}"] = {
            "rankK_distribution": {str(k): int((rankK == k).sum())
                                    for k in np.unique(rankK)},
            "rankG_distribution": {str(k): int((rankG == k).sum())
                                    for k in np.unique(rankG)},
            "median_rankK": int(np.median(rankK)),
        }
        print(f"L={L}: median rank K={int(np.median(rankK))}, "
              f"rank G={int(np.median(rankG))}")
    out["constraint_rank_vs_rL"] = dof_by_L
    out["physical_dof_count"] = {
        "declared": 3,
        "basis": ("y=(psi,dphi,V) after the three Sec-20.4 constraint "
                   "eliminations (h1/H2/A0 via Dh1/DeltaV/pivotA0); "
                   "rankK<3 nodes indicate locally degenerate channels"),
    }

    # ---------- 3. project v_- : which sector carries it?
    # reduce at L=6, find the negative eigenvector of K, measure sector weights
    res6 = red9.canonical_audit(sub, 6)
    K6 = res6["K"]
    sector_names = ["psi", "chi(dphi)", "V"]
    worst = {"r": None, "lam": np.inf, "vec": None}
    for i in range(K6.shape[0]):
        ev, evec = np.linalg.eigh(K6[i])
        if ev[0] < worst["lam"]:
            worst = {"r": float(r_w[i]), "lam": float(ev[0]),
                      "vec": evec[:, 0].copy()}
    lam_min = worst["lam"]
    vec = worst["vec"]
    # sector weights of the negative direction (normalised)
    w = np.abs(vec) ** 2
    w = w / w.sum()
    out["v_minus"] = {
        "r_at_worst": worst["r"],
        "lambda_min": lam_min,
        "sector_weights": {sector_names[i]: round(float(w[i]), 6)
                            for i in range(3)},
        "dominant_sector": sector_names[int(np.argmax(w))],
    }
    print(f"v_- at r={worst['r']:.4f}: lam={lam_min:.3e}, weights={out['v_minus']['sector_weights']}")

    # ---------- 4. K_phys: project out constraint-dominated rows/cols.
    # The chi(dphi) row is the constraint-eliminated channel (its pivot
    # architecture makes it non-propagating: c_chi^2 ~ 1e-4..1e-3 vs 3.98).
    # K_phys = psi/V 2x2 principal subblock per node. If the negative
    # direction lives in chi, K_phys is healthy → artifact; if the psi/V
    # subblock itself goes negative → physical ghost.
    K_phys_min = []
    v_minus_in_phys = []
    for i in range(K6.shape[0]):
        ev, evec = np.linalg.eigh(K6[i])
        if ev[0] >= 0:
            K_phys_min.append(0.0)
            continue
        v = evec[:, 0]
        # physical subblock = rows/cols {0,2} (psi, V) — chi is constraint channel
        Kphys = K6[i][np.ix_([0, 2], [0, 2])]
        K_phys_min.append(float(np.linalg.eigvalsh(Kphys)[0]))
        v_phys = v[[0, 2]]
        v_cons = v[1]
        v_minus_in_phys.append(float(v_cons**2 / max(v @ v, 1e-300)))
    K_phys_min = np.array(K_phys_min)
    neg_nodes = int((K_phys_min < 0).sum())
    min_phys = float(K_phys_min.min())

    out["k_phys"] = {
        "subblock": "[psi, V] principal 2x2 (chi = constraint-eliminated channel)",
        "min_eigenvalue": min_phys,
        "n_nodes_negative": neg_nodes,
        "n_nodes_total": int(K6.shape[0]),
        "gate": GATE,
        "gate_pass": bool(min_phys > 0),
        "v_minus_chi_fraction_max": float(np.max(v_minus_in_phys))
            if v_minus_in_phys else None,
    }
    # where is K_phys negative?
    neg_idx = np.where(K_phys_min < 0)[0]
    if len(neg_idx):
        out["k_phys"]["neg_r_range"] = [round(float(r_w[neg_idx[0]]), 4),
                                         round(float(r_w[neg_idx[-1]]), 4)]
    print(f"K_phys: min={min_phys:.3e}, negative nodes {neg_nodes}/{K6.shape[0]}")

    # ---------- verdict
    if pivot_violation:
        out["verdict"] = ("CONSTRAINT_PIVOT_VIOLATION — a Sec-20.4 pivot "
                           "vanishes on this branch; the reduction itself is "
                           "invalid at those nodes. Do not interpret K's "
                           "negative eigenvalues before this is resolved.")
    elif min_phys > 0:
        out["verdict"] = ("NEGATIVE_CHANNEL_IS_CONSTRAINT_ARTIFACT — the "
                           "negative eigendirection lives in the "
                           "constraint-eliminated chi channel; the physical "
                           "subspace (psi,V) is positive. Branch ALIVE after "
                           "correct constraint reduction. F.4 may rerun on "
                           "K_phys.")
    else:
        out["verdict"] = ("NEGATIVE_CHANNEL_IS_PHYSICAL — the physical "
                           "(psi,V) subspace itself carries negative kinetic "
                           "eigenvalues on this branch: bulk ghost, branch "
                           "DEAD. No boundary-condition rescue is possible.")
    print("VERDICT:", out["verdict"][:120])

    out["wall_seconds"] = round(time.time() - t0, 1)
    OUT.write_text(json.dumps(out, indent=1, allow_nan=False) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
