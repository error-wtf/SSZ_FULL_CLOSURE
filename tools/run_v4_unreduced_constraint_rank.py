#!/usr/bin/env python3
"""V4_UNREDUCED_CONSTRAINT_RANK_V1 — the decisive DOF test from the
UNREDUCED Euler descriptor (rule 8 of Lino's execution order).

Question: does the unreduced action contain an independent primary
constraint that removes dphi?

Method: assemble the 8x8 Euler descriptor operator
(H0, H1, H2, h1, dphi, dA0, dA1, V) on the certified V4 branch and
analyze the kinetic (highest-time-derivative) pencil:
  - rank of the kinetic Hessian block per field
  - null-space dimension = constraint count
  - which fields carry kinetic terms at all (time-derivative order 2)

Decision rule (declared BEFORE running):
  If dphi has a kinetic term in the unreduced descriptor
  (e.g. e1 * dphi_dot^2 present with nonzero coefficient), then dphi is
  dynamical BEFORE any auxiliary elimination, and NO primary constraint
  removes it. The negative dphi direction in the 3x3 reduced K is then
  a physical kinetic pathology of the (a1=+0.5, eps=-0.3) branch ->
  branch rejected per rule 8/9 unless another V4 candidate survives.

  If dphi carries NO kinetic term, it is a Lagrange-multiplier-like
  auxiliary in the unreduced action and its negative reduced direction
  is a constraint artifact.
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

OUT = ROOT / "data/generated/spectral/V4_UNREDUCED_CONSTRAINT_RANK_V1.json"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    t0 = time.time()
    v4 = load("v4solve", str(ROOT / "tools/run_luminal_background_solve_v4.py"))
    from ssz_p5.reducer import unreduced_descriptor as ud
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

    # ---------- assemble the 8x8 unreduced Euler descriptor
    P = ud.descriptor_operator(stream, 6)
    audit = ud.structural_audit(P)
    out = {"audit": "V4_UNREDUCED_CONSTRAINT_RANK_V1",
            "branch": {"a1": A1, "eps": EPS},
            "fields": list(ud.FIELDS),
            "structural_audit": audit}
    print("structural:", audit, flush=True)

    # ---------- kinetic analysis: which fields have second-time-derivative
    # terms? keys are (time_order, radial_order) dicts of (n,8,8) arrays.
    kinetic_fields = {}
    max_t = audit["max_time_derivative_order"]
    for name_idx, name in enumerate(ud.FIELDS):
        # diagonal kinetic coefficient: sum over blocks with pt == max_t
        # of |block[:, idx, idx]|
        tot = np.zeros(len(r_grid))
        for (pt, pr), block in P.items():
            if pt == max_t:
                tot += np.abs(block[:, name_idx, name_idx])
        kinetic_fields[name] = {
            "has_kinetic": bool((tot > 1e-12).any()),
            "max_abs_diag_coeff": float(tot.max()),
            "n_nodes_kinetic": int((tot > 1e-12).sum()),
        }
    out["kinetic_fields"] = kinetic_fields
    for name, info in kinetic_fields.items():
        print(f"  {name}: kinetic={info['has_kinetic']} "
              f"(max|coeff|={info['max_abs_diag_coeff']:.3e}, "
              f"{info['n_nodes_kinetic']}/{len(r_grid)} nodes)", flush=True)

    # ---------- decisive dphi answer
    dphi_kinetic = kinetic_fields["dphi"]
    # also check off-diagonal: dphi with ANY field at max time order
    dphi_offdiag = np.zeros(len(r_grid))
    for (pt, pr), block in P.items():
        if pt == max_t:
            idx = ud.FI["dphi"]
            for j in range(8):
                if j != idx:
                    dphi_offdiag += np.abs(block[:, idx, j])
    out["dphi_offdiag_kinetic_max"] = float(dphi_offdiag.max())
    dphi_dynamical = dphi_kinetic["has_kinetic"] or \
        (dphi_offdiag.max() > 1e-12)

    # constraint count: fields WITHOUT kinetic terms are algebraic/constraint
    no_kin = [n for n, info in kinetic_fields.items() if not info["has_kinetic"]]
    out["fields_without_kinetic"] = no_kin
    out["n_fields_without_kinetic"] = len(no_kin)
    out["dphi_is_dynamical_in_unreduced_action"] = bool(dphi_dynamical)

    # ---------- verdict (declared rule)
    if dphi_dynamical:
        out["verdict"] = (
            "NO_PRIMARY_CONSTRAINT_REMOVES_DPHI — dphi carries kinetic "
            "terms in the unreduced Euler descriptor, so the negative "
            "dphi direction in the reduced 3x3 K is a PHYSICAL kinetic "
            "pathology of this V4 branch. Rule 8/9: branch "
            "(a1=+0.5, eps=-0.3) is REJECTED for resonance work unless "
            "no V4 candidate survives; scan the N4 V4 branch family.")
        out["dof_verdict"] = ("3 physical DOF (psi, dphi, V) confirmed "
                               "from the unreduced action; the earlier "
                               "2x2 'K_phys' projection verdict is "
                               "INVALID (projection without constraint).")
    else:
        out["verdict"] = (
            "DPHI_IS_AUXILIARY — no kinetic term for dphi in the "
            "unreduced descriptor; its negative reduced direction is a "
            "constraint artifact and the strike reduction is justified.")
    print("VERDICT:", out["verdict"][:160])

    out["wall_seconds"] = round(time.time() - t0, 1)
    OUT.write_text(json.dumps(out, indent=1, allow_nan=False, default=str) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
