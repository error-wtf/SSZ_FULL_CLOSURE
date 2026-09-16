from __future__ import annotations

import math

import numpy as np
from numpy.typing import ArrayLike


def derivative(
    x: ArrayLike, y: ArrayLike, order: int = 1, *, window: int = 9, degree: int = 8
) -> np.ndarray:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if x.ndim != 1 or y.ndim != 1 or len(x) != len(y):
        raise ValueError("x and y must be 1-D arrays of equal length")
    if order < 0 or order > degree:
        raise ValueError("order must be between zero and degree")
    if window < degree + 1 or len(x) < window:
        raise ValueError("insufficient stencil")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("nonfinite jet input")
    dx = np.diff(x)
    if not (np.all(dx > 0) or np.all(dx < 0)):
        raise ValueError("grid must be strictly monotone")
    n = len(x)
    out = np.empty(n, float)
    half = window // 2
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, lo + window)
        lo = max(0, hi - window)
        xx = x[lo:hi]
        yy = y[lo:hi]
        x0 = x[i]
        scale = float(np.max(np.abs(xx - x0))) or 1.0
        z = (xx - x0) / scale
        deg = min(degree, len(xx) - 1)
        coeff = np.polynomial.polynomial.polyfit(z, yy, deg)
        out[i] = 0.0 if order > deg else math.factorial(order) * coeff[order] / scale**order
    return out


_JET_CACHE = {}


def _jet_weights(r, order, window=9, degree=8):
    r = np.asarray(r, float)
    key = (r.shape, r.tobytes(), int(order), int(window), int(degree))
    if key in _JET_CACHE:
        return _JET_CACHE[key]
    n = len(r)
    half = window // 2
    inds = np.empty((n, window), dtype=int)
    weights = np.zeros((n, window), float)
    import math

    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, lo + window)
        lo = max(0, hi - window)
        idx = np.arange(lo, hi, dtype=int)
        xx = r[idx]
        x0 = r[i]
        scale = float(np.max(np.abs(xx - x0))) or 1.0
        z = (xx - x0) / scale
        deg = min(degree, len(idx) - 1)
        V = np.polynomial.polynomial.polyvander(z, deg)
        # coef = solve(V,y), hence derivative weights are row(order) of inv(V)
        if order <= deg:
            rhs = np.zeros(deg + 1)
            rhs[order] = math.factorial(order) / scale**order
            w = np.linalg.solve(V.T, rhs)
        else:
            w = np.zeros(len(idx))
        inds[i, : len(idx)] = idx
        if len(idx) < window:
            inds[i, len(idx) :] = idx[-1]
        weights[i, : len(idx)] = w
    _JET_CACHE[key] = (inds, weights)
    return inds, weights


def profile_derivative(r, a, order, window=9, degree=8):
    a = np.asarray(a, float)
    if order == 0:
        return a
    inds, w = _jet_weights(r, order, window, degree)
    vals = a[inds]  # n,window,...
    return np.einsum("nw,nw...->n...", w, vals, optimize=True)
