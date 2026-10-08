#!/usr/bin/env python3
"""REPRODUCER for the ARPACK phantom-pair finding (F1 of ECS_V2_2).

scipy.sparse.linalg.eigs(A, M=L, sigma=lam_target) on a COMPLEX
generalized pencil can return pairs (lam, v) that do NOT satisfy
(A - lam*L) v = 0 — they cluster around sigma by construction of the
shift-invert transform but are not spectrum members.

This script demonstrates it on the ECS square-well FEM pencil and is
the permanent reproducer/detector: ALWAYS validate pairs before use.
"""
import sys

import numpy as np
import scipy.sparse.linalg as spla
from scipy.linalg import eig as dense_eig

sys.path.insert(0, "/home/error/physics/clones/SSZ_FULL_CLOSURE/tools")
from run_ecs_square_well_v22 import fem_contour, load_ref

REF = load_ref()
LAM = REF * REF

A, L, _ = fem_contour(45.0, 0.05, -2.5, 12.0)
Ad, Ld = A.toarray(), L.toarray()

lam_d, vec_d = dense_eig(Ad, Ld, right=True)
dres = min(np.linalg.norm(Ad@vec_d[:, j] - lam_d[j]*(Ld@vec_d[:, j]))
           / np.linalg.norm(Ld@vec_d[:, j]) for j in range(lam_d.size))
d_near = float(np.min(np.abs(lam_d - LAM)))
print(f"dense: {lam_d.size} pairs, max resid ~ {dres:.2e}, "
      f"nearest to sigma: {d_near:.3e}")

lam_a, vec_a = spla.eigs(A, M=L, k=40, sigma=LAM, which="LM",
                         return_eigenvectors=True)
bad = 0
for j in range(lam_a.size):
    lam_j = complex(lam_a[j])
    v = vec_a[:, j]
    r = np.linalg.norm(A@v - lam_j*(L@v)) / np.linalg.norm(L@v)
    in_dense = np.min(np.abs(lam_d - lam_j))
    if r > 1e-6 or in_dense > 1e-6:
        bad += 1
print(f"ARPACK: {lam_a.size} pairs returned, {bad} PHANTOM "
      f"(resid>1e-6 or absent from dense spectrum)")
assert bad == lam_a.size, "expected all-phantom on this pencil"
print("REPRODUCED: do not trust scipy generalized shift-invert pairs "
      "without residual validation.")
