#!/usr/bin/env python3
"""ECS_V2_2 DEV DIAG 5: DENSE spectrum — no ARPACK selection bias.
Does the operator actually HAVE the resonance eigenvalue?"""
import sys
import numpy as np
import scipy.sparse as sp
from scipy.linalg import eig as dense_eig

sys.path.insert(0, "/home/error/physics/clones/SSZ_FULL_CLOSURE/tools")
from run_ecs_square_well_v22 import load_ref

REF = load_ref()
LAM_REF = REF * REF
print("lambda_ref =", LAM_REF)


def grid(h, theta, tail, r_match, a=1.0):
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


h = 0.05
for v0, tag in ((-2.5, "WELL"), (+2.5, "BARR")):
    z = grid(h, 45.0, 12.0, 6.0)
    H = ham(z, v0).toarray()
    lam = dense_eig(H, right=False)
    if isinstance(lam, tuple):
        lam = lam[0]
    d = np.abs(lam - LAM_REF)
    order = np.argsort(d)
    print(f"--- {tag} (n={len(lam)}): 5 nearest to lambda_ref")
    for j in order[:5]:
        l = lam[j]
        print(f"    lam={l.real:10.5f}{l.imag:+10.5f}j  dist={d[j]:.4e}")
