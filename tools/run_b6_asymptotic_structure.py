#!/usr/bin/env python3
"""B6 step 1: asymptotic analysis of the V4 operator (declared diagnostic).

Before choosing physical BCs we need the asymptotic structure of the
coupled radial operator at BOTH ends of the ghost-free V4 branch
(a1=+0.5, eps=-0.3):
  * outer edge (r -> 1.41): where the negative kinetic channel grows —
    what is the local wave equation?  Plane-wave ansatz gives the local
    dispersion and the local wave speeds per sector.
  * inner edge (r -> ~0.01): same.
From the local speeds c_n^2(r) = eig(K/G) we can decide which boundary is
wavelike (absorbing/outgoing possible) and which is not.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/home/error/physics/clones/SSZ_FULL_CLOSURE")
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

OUT = ROOT / "data/generated/spectral/B6_ASYMPTOTIC_STRUCTURE_V1.json"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    v4 = load("v4solve", str(ROOT / "tools/run_luminal_background_solve_v4.py"))
    from ssz_p5.jets.jet9d8 import derivative as jet

    archive = pd.read_csv(v4.ARCHIVE).sort_values("u").reset_index(drop=True)
    r_all = (1.0 / archive.u.to_numpy(float))[::-1]
    # full domain profile (not window) for the global picture
    red, facts, _ = v4.build_reduction()
    rhs = v4.rhs_factory(red, 0.5)
    sol = v4.integrate_branch(rhs, float(r_all[0]), float(r_all[-1]), -0.3)
    assert bool(sol.t[-1] >= float(r_all[-1]) * (1 - 1e-9))
    y = sol.sol(r_all)
    f, h, phi = y[0], y[1], y[2]

    out = {"audit": "B6_ASYMPTOTIC_STRUCTURE_V1",
            "branch": {"a1": 0.5, "eps": -0.3},
            "domain": [float(r_all[0]), float(r_all[-1])], "rows": int(len(r_all))}

    # local metric functions at both ends and midpoint
    regions = {"inner_core": (r_all[0], 0.05),
                "interior": (0.05, 0.6),
                "outer": (0.6, 1.2),
                "outer_edge": (1.2, float(r_all[-1]))}
    stats = {}
    for name, (lo, hi) in regions.items():
        m = (r_all >= lo) & (r_all < hi)
        stats[name] = {
            "f": [round(float(f[m].min()), 5), round(float(f[m].max()), 5)],
            "h": [round(float(h[m].min()), 5), round(float(h[m].max()), 5)],
            "phi": [round(float(phi[m].min()), 5), round(float(phi[m].max()), 5)],
            "rows": int(m.sum()),
        }
    out["profile_regions"] = stats
    print(json.dumps(stats, indent=1))

    # Tortoise coordinate behaviour: dr*/dr = 1/sqrt(f h) (static spherical)
    sqfh = np.sqrt(f * h)
    inv_sqfh = 1.0 / sqfh
    out["tortoise"] = {
        "drstar_dr_inner": round(float(inv_sqfh[0]), 4),
        "drstar_dr_outer": round(float(inv_sqfh[-1]), 4),
        "monotone": bool(np.all(np.diff(inv_sqfh) >= -1e-12)),
    }
    # local wave speed per sector from G and K diagonals on the exported window
    arrays = dict(np.load(ROOT / "data/generated/spectral/V4_GHOSTFREE_BRANCH_EXPORT_V1.npz"))
    for L in (6, 12, 20):
        r_op = arrays[f"{L}_r"]
        G, K = arrays[f"{L}_G"], arrays[f"{L}_K"]
        # local speeds: eig of diag-blocks K_ii / G_ii per node
        spd = np.zeros((len(r_op), 3))
        for i in range(3):
            kii = K[:, i, i]
            gii = G[:, i, i]
            with np.errstate(divide="ignore", invalid="ignore"):
                s = kii / gii
            s = np.where(np.isfinite(s), s, np.nan)
            spd[:, i] = s
        spd = np.abs(spd)
        seg = {"inner": slice(0, len(r_op)//3),
                "middle": slice(len(r_op)//3, 2*len(r_op)//3),
                "outer": slice(2*len(r_op)//3, len(r_op))}
        out[f"local_speed2_L{L}"] = {
            seg_name: {f"sector_{i}": (round(float(np.nanmedian(spd[sl, i])), 4)
                                        if np.isfinite(spd[sl, i]).any() else None)
                        for i in range(3)}
            for seg_name, sl in seg.items()
        }
    print(json.dumps(out["local_speed2_L6"], indent=1)[:600])

    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print("written:", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
