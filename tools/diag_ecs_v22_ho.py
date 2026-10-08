#!/usr/bin/env python3
"""ECS_V2_2 DEV DIAG 7: end-to-end machinery validation.

(A) H(0) on a real grid vs harmonic oscillator -psi'' + x^2/4 on (-L, L):
    known eigenvalues 0.5, 1.5, 2.5, ...
(B) same H on a COMPLEX-ROTATED grid about x=0 (theta): known spectrum
    (n+1/2) e^{-2i theta}.  This validates the complex FD grid, the
    one-sided center row and the Dirichlet truncation end-to-end.
"""
import sys
import numpy as np
import scipy.sparse as sp
from scipy.linalg import eig as dense_eig


def sym_grid(L, n):
    return np.linspace(-L, L, n).astype(complex)


def ham_general(z, vfun, bc0="dirichlet", bc1="dirichlet"):
    n = len(z)
    rows, cols, vals = [], [], []
    for i in range(1, n - 1):
        hl = z[i] - z[i-1]; hr = z[i+1] - z[i]
        c = 2.0 / (hl * hr * (hl + hr))
        rows += [i, i, i]; cols += [i-1, i, i+1]
        vals += [c*hr, -c*(hl+hr), c*hl]
    h2 = z[2] - z[0]
    # left BC
    if bc0 == "neumann":
        rows += [0,0,0]; cols += [0,1,2]
        vals += [-3.0/h2, 4.0/h2, -1.0/h2]
    else:
        rows += [0]; cols += [0]; vals += [1.0+0j]
    if bc1 == "neumann":
        hm = z[-1] - z[-3]
        rows += [n-1,n-1,n-1]; cols += [n-1,n-2,n-3]
        vals += [3.0/hm, -4.0/hm, 1.0/hm]
    else:
        rows += [n-1]; cols += [n-1]; vals += [1.0+0j]
    D2 = sp.csc_matrix((vals, (rows, cols)), shape=(n, n))
    V = np.asarray(vfun(z), complex)
    return ((-D2) + sp.diags(V)).tocsc()


print("=== A: theta=0 harmonic oscillator, Dirichlet, L=15, n=301 ===")
z = sym_grid(15.0, 301)
H = ham_general(z, lambda x: x**2/4.0)
lam = dense_eig(H.toarray(), right=False)
if isinstance(lam, tuple):
    lam = lam[0]
lam = np.sort(lam.real)
print("first 5:", lam[:5])

print()
print("=== B: complex-rotated grid (rotation about 0), theta=20deg ===")
th = np.radians(20.0)
zrot = z * np.exp(1j*th)     # z = e^{i theta} x  (full complex scaling)
H = ham_general(zrot, lambda x: x**2/4.0)
lam = dense_eig(H.toarray(), right=False)
if isinstance(lam, tuple):
    lam = lam[0]
expect = (np.arange(6) + 0.5) * np.exp(-2j*th)
order = np.argsort(np.abs(lam))
print("5 smallest |lambda|:")
for j in order[:5]:
    print(f"   lam={lam[j].real:9.5f}{lam[j].imag:+9.5f}j")
print("expected:", np.round(expect, 5))
