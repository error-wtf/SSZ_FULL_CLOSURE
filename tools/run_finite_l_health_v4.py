#!/usr/bin/env python3
"""F.4: finite-L health rerun on a V4 (varying-G4) luminal branch.

Consumes a certified V4 branch (default: the ghost-free side a1*eps < 0,
K >= 0 everywhere per the V4 N4 run), emits the 41-slot stream through the
VERIFIED chain, reduces per L in DEFAULT_L, and runs the fail-closed
finite-L stability diagnostics (SSZ-B6 pre-stage).

Honest scope: this is NOT a QNM claim and NOT a bridge certificate.
QNM_SPECTROSCOPY_PATH_STATUS stays qnm_claims_allowed=false until the
physical-BC layer (B6) is healthy.  This tool only decides the H1
'finite-L health pass' clause of data/diagnostic/N1_N4_DECISION_RULES.json
on the declared V4 branch.

Usage:
    python tools/run_finite_l_health_v4.py \
        --a1 0.5 --eps -0.5 [--window-r 0.05 1.35] [--out NAME]
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.coefficients.mh_action_primitives import (  # noqa: E402
    quartic_g5zero_primitives,
)
from ssz_p5.coefficients.mh_general_primitives import (  # noqa: E402
    emit_from_primitives,
)
from ssz_p5.config import DEFAULT_L, SLOT_NAMES  # noqa: E402
from ssz_p5.geometry.p5 import from_frame  # noqa: E402
from ssz_p5.provenance.manifest import sha256  # noqa: E402
from ssz_p5.reducer.canonical import reduce_profile  # noqa: E402
from ssz_p5.stability.finite_l import stability_diagnostics  # noqa: E402
from ssz_p5.types import Coefficients41  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "v4solver", ROOT / "tools" / "run_luminal_background_solve_v4.py"
)
_v4 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_v4)

ARCHIVE = ROOT / (
    "data/production/ssz_p5_F1b_G4XX_transverse_core_candidate_2026-09-14.csv"
)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--a1", type=float, required=True)
    ap.add_argument("--eps", type=float, required=True)
    ap.add_argument("--window-r", type=float, nargs=2, default=(0.05, 1.35))
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    if args.a1 * args.eps >= 0:
        print(
            "REFUSED: this tool decides the H1 health clause and therefore "
            "requires the ghost-free side a1*eps < 0 (K >= 0 measured in the "
            "V4 N4 run); the ghost side goes through H3, not F.4"
        )
        return 2

    archive = pd.read_csv(ARCHIVE).sort_values("u").reset_index(drop=True)
    r_all = (1.0 / archive.u.to_numpy(float))[::-1]
    r_min, r_core = float(r_all[0]), float(r_all[-1])

    red, facts, _ = _v4.build_reduction()
    sol = _v4.integrate_branch(
        _v4.rhs_factory(red, args.a1), r_min, r_core, args.eps
    )
    reached = bool(sol.t[-1] >= r_core * (1 - 1e-9))
    if not reached:
        print(f"REFUSED: branch did not reach the core edge (t={sol.t[-1]})")
        return 2
    y = sol.sol(r_all)
    f_p, h_p, phi_p = y[0], y[1], y[2]

    # interior window (declared margin, jet9d8 edge stencils)
    lo = int(np.searchsorted(r_all, args.window_r[0]))
    hi = int(np.searchsorted(r_all, args.window_r[1]))
    M = 8
    win = slice(lo + M, hi - M)
    r_w = r_all[win]
    f_w, h_w, phi_w = f_p[win], h_p[win], phi_p[win]
    # u_w = 1.0 / r_w  # (emitter supplies u internally)

    from ssz_p5.jets.jet9d8 import derivative as jet

    phi_r = jet(r_w, phi_w)
    feed = pd.DataFrame({
        "x": r_w, "f": f_w, "h": h_w, "phi_r": phi_r,
        "G4": 0.5 + args.a1 * (phi_w - 1.0),
        "G4X": np.zeros(len(r_w)), "G4XX": np.zeros(len(r_w)),
        "G4phi": np.full(len(r_w), args.a1),
        "G4phiX": np.zeros(len(r_w)), "G3X": np.zeros(len(r_w)),
    })
    prim = quartic_g5zero_primitives(feed)
    c2 = np.sqrt(f_w * h_w) * phi_r * 0.5 * r_w**2
    prim_in = pd.DataFrame({
        "x": r_w, "f": f_w, "h": h_w, "phiprime": phi_r,
        "A0prime": np.zeros(len(r_w)),
        "a1": prim.a1_action.to_numpy(float),
        "c2": c2,
        "c4": prim.c4_action.to_numpy(float),
        "F_tensor": prim.F_tensor_action.to_numpy(float),
        "G_tensor": prim.G_tensor_action.to_numpy(float),
        "H_tensor": prim.H_tensor_action.to_numpy(float),
    })
    stream = emit_from_primitives(prim_in)
    # the emitter's internal phi-column integration is NaN-prone; overwrite
    # with the solver's own phi(r) (exact state) and re-sort
    stream["phi"] = phi_w
    sub = stream.sort_values("x").reset_index(drop=True)
    # X (static radial convention) is reconstructed from the branch identity
    # h*phi_r^2 + 2X = 0 (geometry/p5.py) if the emitter does not provide it
    if "X" not in sub.columns:
        sub["X"] = -sub["h"].to_numpy(float) * (
            sub["phiprime"].to_numpy(float) ** 2
        ) / 2.0
    for col in ("x", "u", "phi", "f", "h", "phiprime", "A0prime", "X"):
        if col not in sub.columns:
            raise ValueError(f"emitter stream lacks {col}")

    c = Coefficients41(
        from_frame(sub),
        {s: sub[s].to_numpy(float) for s in SLOT_NAMES},
        f"v4_branch_a1_{args.a1}_eps_{args.eps}",
    )
    health = {}
    all_pass = True
    for L in DEFAULT_L:
        op = reduce_profile(c, int(L))
        diag = stability_diagnostics(op)
        health[str(int(L))] = diag
        all_pass = all_pass and bool(diag["pass"])
        print(
            f"L={int(L)}: min_eig_K={diag['min_eig_K']:.3e} "
            f"min_radial={diag['min_radial']:.3e} pass={diag['pass']}"
        )

    result = {
        "audit": "F4_FINITE_L_HEALTH_V4",
        "branch": {"a1": args.a1, "eps": args.eps,
                   "ghost_free_side": args.a1 * args.eps < 0},
        "window_r": [float(r_w[0]), float(r_w[-1])],
        "window_rows": int(len(r_w)),
        "K_on_branch": "K >= 0 everywhere (V4 N4 measurement, ghost-free side)",
        "finite_l_health": health,
        "all_pass": bool(all_pass),
        "h1_health_clause": (
            "PASS" if all_pass else "FAIL"
        ),
        "claim_scope": (
            "decides the finite-L health clause of N1_N4_DECISION_RULES "
            "H1 only; NO QNM claim, qnm_claims_allowed stays false"
        ),
    }
    out_name = args.out or (
        f"F4_FINITE_L_HEALTH_V4_a1_{args.a1}_eps_{args.eps}.json"
        .replace("-", "m").replace(".", "p").replace("pp", "p")
    )
    out_path = ROOT / "data/generated/spectral" / out_name
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=1, allow_nan=False) + "\n")
    (out_path.parent / (out_path.name + ".sha256")).write_text(
        sha256(out_path) + "\n"
    )
    print("all_pass:", all_pass, "->", out_path.relative_to(ROOT))
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
