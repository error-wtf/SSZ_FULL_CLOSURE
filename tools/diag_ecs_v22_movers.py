#!/usr/bin/env python3
"""ECS_V2_2 DEV DIAG 2: find the V-DEPENDENT eigenvalue (physical root)."""
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, "/home/error/physics/clones/SSZ_FULL_CLOSURE/tools")
from run_ecs_square_well_v22 import build_grid, build_h, load_ref

REF = load_ref()
LAM_REF = REF * REF


def spectrum(v0, h, k, theta=45.0, tail=12.0, sigma=None):
    z = build_grid(h, theta, tail)
    H, _ = build_h(z, v0, 1.0)
    lam, _ = spla.eigs(H, k=k, sigma=LAM_REF if sigma is None else sigma,
                       which="LM", return_eigenvectors=True)
    return np.array(lam)


h = 0.005
for k in (60,):
    lam_w = spectrum(-2.5, h, k)
    lam_b = spectrum(+2.5, h, k)
    print(f"k={k}: |well|={len(lam_w)} |barrier|={len(lam_b)}")
    # movers: for each well eigenvalue, distance to nearest barrier one
    movers = []
    for l in lam_w:
        d = np.min(np.abs(lam_b - l))
        movers.append((d, l))
    movers.sort(key=lambda t: -t[0])
    print("top movers (potential-sensitive):")
    for d, l in movers[:10]:
        w = np.sqrt(complex(l))
        if w.imag > 0:
            w = -w
        print(f"  lam={l.real:9.4f}{l.imag:+9.4f}j  omega={w.real:8.4f}"
              f"{w.imag:+8.4f}j  mover_dist={d:.4f}  "
              f"dist_to_REF={abs(w-REF)/abs(REF):.4e}")
