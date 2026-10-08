#!/usr/bin/env python3
"""ECS_V2_2 DEV DIAG 6: dense well-vs-barrier MOVERS — where is the
interior-localized (potential-sensitive) resonance eigenvalue?"""
import sys
import numpy as np
import scipy.sparse as sp
from scipy.linalg import eig as dense_eig

sys.path.insert(0, "/home/error/physics/clones/SSZ_FULL_CLOSURE/tools")
from run_ecs_square_well_v22 import load_ref

REF = load_ref()
LAM_REF = REF * REF


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
z = grid(h, 45.0, 12.0, 6.0)
lams = {}
for v0, tag in ((-2.5, "well"), (+2.5, "barr")):
    H = ham(z, v0).toarray()
    lam = dense_eig(H, right=False)
    if isinstance(lam, tuple):
        lam = lam[0]
    lams[tag] = lam
    print(tag, "n =", len(lam))

# potential-sensitive pairs: greedy nearest matching
bw = list(lams["well"])
bb = list(lams["barr"])
pairs = []
for l in bw:
    j = int(np.argmin([abs(l - m) for m in bb]))
    pairs.append((abs(l - bb[j]), l, bb[j]))
    bb.pop(j)
pairs.sort(key=lambda t: -t[0])
print("top-10 movers (well vs barrier):")
for d, lw, lb in pairs[:10]:
    w = np.sqrt(complex(lw))
    if w.imag > 0:
        w = -w
    print(f"  dλ={d:8.4f}  well λ={lw.real:9.4f}{lw.imag:+9.4f}j  "
          f"ω={w.real:8.4f}{w.imag:+8.4f}j  dist_REF={abs(w-REF)/abs(REF):.4e}")

print()
print("distance of lambda_ref to nearest WELL eigenvalue:",
      float(np.min(np.abs(lams["well"] - LAM_REF))))
