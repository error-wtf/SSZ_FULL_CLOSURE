#!/usr/bin/env python3
"""ECS_V2_2 DEV DIAGNOSTIC (pre-certification, not an artifact run).

Questions:
  Q1. Which eigenvalue is the TRUE square-well resonance? Test: eigenvalue
      RESPONSE to V0 (well -2.5 vs barrier +2.5 vs empty 0). The physical
      resonance moves; grid artifacts do not.
  Q2. Convergence order: is the O(h^0.6) slope a V-jump artifact?
      Compare node-sampled V vs cell-averaged V at the discontinuity.
  Q3. Window size: does k=14 hide the true root behind continuum states?
"""
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, "/home/error/physics/clones/SSZ_FULL_CLOSURE/tools")
from run_ecs_square_well_v22 import build_grid, load_ref

REF = load_ref()
LAM_REF = REF * REF


def build_h2(z, v0, a, cell_avg=False):
    n = len(z)
    rows, cols, vals = [], [], []
    for i in range(1, n - 1):
        hl = z[i] - z[i - 1]
        hr = z[i + 1] - z[i]
        c = 2.0 / (hl * hr * (hl + hr))
        rows += [i, i, i]
        cols += [i - 1, i, i + 1]
        vals += [c * hr, -c * (hl + hr), c * hl]
    h2 = z[2] - z[0]
    rows += [0, 0, 0]
    cols += [0, 1, 2]
    vals += [-3.0 / h2, 4.0 / h2, -1.0 / h2]
    rows.append(n - 1); cols.append(n - 1); vals.append(1.0 + 0j)
    D2 = sp.csc_matrix((vals, (rows, cols)), shape=(n, n))
    if cell_avg:
        V = np.zeros(n)
        left = np.real(z) - (z[1] - z[0]) / 2 > a   # cell mid in exterior?
        right = np.real(z) + (z[1] - z[0]) / 2 <= a
        V[right & ~left] = v0          # full cell inside well
        V[~right & ~left] = v0 / 2.0   # jump node: half cell
        V[0] = v0
    else:
        V = np.where(np.real(z) <= a, v0, 0.0)
    V[-1] = 0.0
    H = (-D2) + sp.diags(V.astype(complex))
    return H.tocsc()


def nearest_roots(v0, h, cell_avg, k=40):
    z = build_grid(h, 45.0, 12.0)
    H = build_h2(z, v0, 1.0, cell_avg)
    lam, _ = spla.eigs(H, k=k, sigma=LAM_REF, which="LM",
                       return_eigenvectors=True)
    ws = []
    for l in lam:
        w = np.sqrt(complex(l))
        if w.imag > 0:
            w = -w
        if w.real > 0 and w.imag < 0:
            ws.append(w)
    return sorted(ws, key=lambda w: abs(w - REF))


print("=== Q1: V-response at h=0.005 (cell_avg=False) ===")
sets = {}
for tag, v0 in (("well-2.5", -2.5), ("empty 0.0", 0.0), ("barrier+2.5", 2.5)):
    ws = nearest_roots(v0, 0.005, False, k=40)
    sets[tag] = ws
    print(f"--- {tag}: 5 nearest (Re>0, Im<0)")
    for w in ws[:5]:
        print(f"    {w.real:.6f} {w.imag:+.6f}j  dist {abs(w-REF)/abs(REF):.4e}")

print()
print("=== Q2: convergence, cell_avg True vs False (V0=-2.5) ===")
for cell in (False, True):
    print(f"--- cell_avg={cell}")
    for h in (0.02, 0.01, 0.005, 0.0025):
        ws = nearest_roots(-2.5, h, cell, k=40)
        w = ws[0] if ws else None
        if w is not None:
            print(f"    h={h:<7} omega={w.real:.6f} {w.imag:+.6f}j  "
                  f"rel_err={abs(w-REF)/abs(REF):.4e}")
        else:
            print(f"    h={h:<7} NO CANDIDATE")
