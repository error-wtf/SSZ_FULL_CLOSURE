#!/usr/bin/env python3
"""ECS_V2_2 DEV DIAG 8: exact-eigenvector residual on the rotated grid.
H psi_0 with psi_0(z) = exp(-z^2/4) must give 0.5*psi_0. Find the bug."""
import sys
import numpy as np
import scipy.sparse as sp


def ham_general(z, vfun):
    n = len(z)
    rows, cols, vals = [], [], []
    for i in range(1, n - 1):
        hl = z[i] - z[i-1]; hr = z[i+1] - z[i]
        c = 2.0 / (hl * hr * (hl + hr))
        rows += [i, i, i]; cols += [i-1, i, i+1]
        vals += [c*hr, -c*(hl+hr), c*hl]
    rows += [0]; cols += [0]; vals += [1.0+0j]
    rows += [n-1]; cols += [n-1]; vals += [1.0+0j]
    D2 = sp.csc_matrix((vals, (rows, cols)), shape=(n, n))
    V = np.asarray(vfun(z), complex)
    return ((-D2) + sp.diags(V)).tocsc()


x = np.linspace(-15.0, 15.0, 301)
th = np.radians(20.0)
z = x * np.exp(1j*th)
psi0 = np.exp(-z**2/4.0)
H = ham_general(z, lambda zz: zz**2/4.0)
r = H @ psi0 - 0.5*psi0
# exclude the two Dirichlet boundary rows from the residual
print("interior max |r|:", np.max(np.abs(r[1:-1])))
print("interior max |r| / 0.5:",
      np.max(np.abs(r[1:-1]))/0.5)
# inspect single-row stencil response:
i = 150  # x=0
print("row150 r:", r[i])
# analytic 2nd derivative of exp(-z^2/4): (z^2/4 - 1/2) e^{-z^2/4}
d2_exact = (z**2/4.0 - 0.5)*psi0
# FD d2 from the stencil coefficients at row i:
hl = z[i]-z[i-1]; hr = z[i+1]-z[i]
c = 2.0/(hl*hr*(hl+hr))
fd2 = c*hr*psi0[i-1] - c*(hl+hr)*psi0[i] + c*hl*psi0[i+1]
print("FD d2  :", fd2)
print("exa d2 :", d2_exact[i])
print("ratio  :", fd2/d2_exact[i])
