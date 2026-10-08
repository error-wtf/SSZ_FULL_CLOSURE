#!/usr/bin/env python3
"""ECS_V2_2 DEV DIAG 4: tracked-root discriminator (well vs barrier).

The rotated-continuum band exists for EVERY V0, so nearest-neighbour
mover distances are meaningless. The physical root is the eigenvalue
whose h->0 limit is REF for the well and is NOT for the barrier.
"""
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, "/home/error/physics/clones/SSZ_FULL_CLOSURE/tools")
from run_ecs_square_well_v22 import load_ref

REF = load_ref()
LAM_REF = REF * REF


def grid(h, theta, tail, r_match, a=1.0):
    m = int(round(1.0 / h))
    n_in = int(round(r_match / h))
    r1 = np.arange(0, n_in + 1) * h
    s = np.arange(1, int(round(tail / h)) + 1) * h
    th = np.radians(theta)
    return np.concatenate([r1.astype(complex), r_match + np.exp(1j*th)*s])


def ham(z, v0, a=1.0):
    n = len(z)
    rows, cols, vals = [], [], []
    for i in range(1, n - 1):
        hl = z[i] - z[i-1]; hr = z[i+1] - z[i]
        c = 2.0 / (hl * hr * (hl + hr))
        rows += [i, i, i]; cols += [i-1, i, i+1]
        vals += [c*hr, -c*(hl+hr), c*hl]
    h2 = z[2] - z[0]
    rows += [0,0,0]; cols += [0,1,2]
    vals += [-3.0/h2, 4.0/h2, -1.0/h2]
    rows.append(n-1); cols.append(n-1); vals.append(1.0+0j)
    D2 = sp.csc_matrix((vals, (rows, cols)), shape=(n, n))
    V = np.where(np.real(z) <= a, v0, 0.0)
    V[-1] = 0.0
    return ((-D2) + sp.diags(V.astype(complex))).tocsc()


def tracked(v0, h, r_match, theta=45.0, tail=12.0, k=40):
    z = grid(h, theta, tail, r_match)
    H = ham(z, v0)
    lam, _ = spla.eigs(H, k=k, sigma=LAM_REF, which="LM",
                       return_eigenvectors=True)
    best = min(lam, key=lambda l: abs(complex(l) - LAM_REF))
    w = np.sqrt(complex(best))
    if w.imag > 0:
        w = -w
    return w


for rm in (2.0, 6.0):
    print(f"===== r_match={rm}")
    for v0, tag in ((-2.5, "WELL"), (+2.5, "BARR")):
        for h in (0.005, 0.0025):
            w = tracked(v0, h, rm)
            print(f"  {tag} h={h:<7} omega={w.real:.6f}{w.imag:+.6f}j  "
                  f"dist_REF={abs(w-REF)/abs(REF):.4e}")
