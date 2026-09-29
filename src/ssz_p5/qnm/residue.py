"""Fail-closed observable-residue definitions for open/QNM systems.

A local positive K-density is useful for self-adjoint bound-state problems but
is not, by itself, a physical QNM residue.  Open systems require a retarded
Green function or an equivalent left/right mode normalization.
"""
from __future__ import annotations

import numpy as np


def pole_residue(
    right: np.ndarray,
    left: np.ndarray,
    dL_domega: np.ndarray,
    observable: np.ndarray | None = None,
) -> complex:
    """Biorthogonal residue factor at an isolated simple pole.

    denominator = left^H (dL/domega) right.
    If observable is supplied, the numerator is the observable-projected local
    right/left product.  This helper intentionally does not infer an observable.
    """
    r = np.asarray(right, complex)
    l = np.asarray(left, complex)
    d = np.asarray(dL_domega, complex)
    denom = np.vdot(l, d @ r)
    if abs(denom) == 0:
        raise ValueError("singular QNM pole normalization")
    if observable is None:
        numer = np.vdot(l, r)
    else:
        o = np.asarray(observable, complex)
        numer = np.vdot(l, o @ r)
    return numer / denom


def spectral_density_from_green(green_diag: np.ndarray, observable: np.ndarray | None = None):
    """rho = -(1/pi) Im Tr[O G^R O^H] for diagonal-position Green blocks."""
    g = np.asarray(green_diag, complex)
    if g.ndim < 2 or g.shape[-1] != g.shape[-2]:
        raise ValueError("green_diag must end in square field blocks")
    if observable is None:
        projected = g
    else:
        o = np.asarray(observable, complex)
        projected = o @ g @ o.conj().T
    return -np.imag(np.trace(projected, axis1=-2, axis2=-1)) / np.pi
