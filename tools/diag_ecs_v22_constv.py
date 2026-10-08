#!/usr/bin/env python3
"""ECS_V2_2 DEV DIAG 9: constant-complex-V known-answer test on the
contour. Exact eigenvalues: lam_n = wp^2 + (n pi / L_c)^2 with complex
contour length L_c = r_match + e^{i theta} * tail. FD must reproduce
them to O(h^2) or the machinery is broken (not the concept)."""
import sys
import numpy as np
import scipy.sparse as sp
from scipy.linalg import eig as dense_eig

r_match, tail, theta = 6.0, 12.0, 45.0
wp = complex(1.6 - 0.9j)
h = 0.01

n_in = int(round(r_match / h))
r1 = np.arange(0, n_in + 1) * h
s = np.arange(1, int(round(tail / h)) + 1) * h
z = np.concatenate([r1.astype(complex), r_match + np.exp(1j*np.radians(theta))*s])
n = len(z)

rows, cols, vals = [], [], []
for i in range(1, n-1):
    hl = z[i]-z[i-1]; hr = z[i+1]-z[i]
    c = 2.0/(hl*hr*(hl+hr))
    rows += [i,i,i]; cols += [i-1,i,i+1]
    vals += [c*hr, -c*(hl+hr), c*hl]
h2 = z[2]-z[0]
rows += [0]; cols += [0]; vals += [1.0+0j]          # Dirichlet center
rows += [n-1]; cols += [n-1]; vals += [1.0+0j]      # Dirichlet end
D2 = sp.csc_matrix((vals, (rows, cols)), shape=(n, n))
H = ((-D2) + sp.diags(np.full(n, wp))).tocsc()

lam = dense_eig(H.toarray(), right=False)
if isinstance(lam, tuple):
    lam = lam[0]

Lc = r_match + np.exp(1j*np.radians(theta))*tail
print("L_c =", Lc, " |L_c| =", abs(Lc))
err = []
for nn in range(1, 8):
    exact = wp**2 + (nn*np.pi/Lc)**2
    d = np.abs(lam - exact)
    j = int(np.argmin(d))
    print(f"n={nn}: exact={exact.real:9.5f}{exact.imag:+9.5f}j  "
          f"fd={lam[j].real:9.5f}{lam[j].imag:+9.5f}j  |err|={d[j]:.3e}")
    err.append(d[j])
print("max err:", max(err))
