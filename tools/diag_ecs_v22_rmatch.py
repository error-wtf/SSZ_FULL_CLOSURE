#!/usr/bin/env python3
"""ECS_V2_2 DEV DIAG 3: r_match scan — rotation must start right behind
the potential (classic ECS prescription), not at r_match=6."""
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, "/home/error/physics/clones/SSZ_FULL_CLOSURE/tools")
from run_ecs_square_well_v22 import load_ref

REF = load_ref()
LAM_REF = REF * REF


def build_grid2(h, theta, tail, r_match, a=1.0):
    m = int(round(1.0 / h))
    n_in = int(round(r_match / h))
    assert abs(n_in * h - r_match) < 1e-12
    assert abs(a * m % 1) < 1e-9
    r1 = np.arange(0, n_in + 1) * h
    s = np.arange(1, int(round(tail / h)) + 1) * h
    th = np.radians(theta)
    z2 = r_match + np.exp(1j * th) * s
    return np.concatenate([r1.astype(complex), z2])


def build_h3(z, v0, a=1.0):
    n = len(z)
    rows, cols, vals = [], [], []
    for i in range(1, n - 1):
        hl = z[i] - z[i - 1]
        hr = z[i + 1] - z[i]
        c = 2.0 / (hl * hr * (hl + hr))
        rows += [i, i, i]; cols += [i - 1, i, i + 1]
        vals += [c * hr, -c * (hl + hr), c * hl]
    h2 = z[2] - z[0]
    rows += [0, 0, 0]; cols += [0, 1, 2]
    vals += [-3.0 / h2, 4.0 / h2, -1.0 / h2]
    rows.append(n - 1); cols.append(n - 1); vals.append(1.0 + 0j)
    D2 = sp.csc_matrix((vals, (rows, cols)), shape=(n, n))
    V = np.where(np.real(z) <= a, v0, 0.0)
    V[-1] = 0.0
    return ((-D2) + sp.diags(V.astype(complex))).tocsc()


def movers_for(r_match, h=0.005, theta=45.0, tail=12.0, k=40):
    out = []
    res = {}
    for tag, v0 in (("well", -2.5), ("barr", +2.5)):
        z = build_grid2(h, theta, tail, r_match)
        H = build_h3(z, v0)
        lam, _ = spla.eigs(H, k=k, sigma=LAM_REF, which="LM",
                           return_eigenvectors=True)
        res[tag] = lam
    md = []
    for l in res["well"]:
        d = np.min(np.abs(res["barr"] - l))
        md.append((d, l))
    md.sort(key=lambda t: -t[0])
    return md


for rm in (6.0, 3.0, 2.0, 1.5, 1.2, 1.1, 1.05):
    md = movers_for(rm)
    print(f"--- r_match={rm}: top-3 movers")
    for d, l in md[:3]:
        w = np.sqrt(complex(l))
        if w.imag > 0:
            w = -w
        print(f"    lam={l.real:9.4f}{l.imag:+9.4f}j  omega={w.real:8.4f}"
              f"{w.imag:+8.4f}j  mover={d:.4f}  dist_REF={abs(w-REF)/abs(REF):.4e}")
