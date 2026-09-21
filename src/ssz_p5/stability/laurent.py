"""Small truncated Laurent-series algebra for the SSZ high-L reducers.

The asymptotic variable is eps = 1/L.  Coefficients may be numpy arrays; eps is
independent of radius, so radial differentiation acts coefficient-wise.

This module deliberately does *not* infer physics from large numerical L.  It is
used to propagate formal Laurent coefficients through the same algebraic maps
as the finite-L common-action reducer.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, Iterable

import numpy as np


def _mul_arrays(a, b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    if a.ndim==1 and b.ndim>1 and a.shape[0]==b.shape[0]:
        a=a.reshape((a.shape[0],)+(1,)*(b.ndim-1))
    if b.ndim==1 and a.ndim>1 and b.shape[0]==a.shape[0]:
        b=b.reshape((b.shape[0],)+(1,)*(a.ndim-1))
    return a*b


@dataclass(frozen=True)
class LaurentSeries:
    coeffs: Dict[int, np.ndarray]
    pmin: int = -6
    pmax: int = 14

    def __post_init__(self):
        clean: Dict[int, np.ndarray] = {}
        for p, a in self.coeffs.items():
            p = int(p)
            if self.pmin <= p <= self.pmax:
                clean[p] = np.asarray(a, dtype=float)
        object.__setattr__(self, "coeffs", clean)

    @classmethod
    def const(cls, a, *, pmin=-6, pmax=14) -> "LaurentSeries":
        return cls({0: np.asarray(a, dtype=float)}, pmin=pmin, pmax=pmax)

    @classmethod
    def monomial(cls, p: int, a, *, pmin=-6, pmax=14) -> "LaurentSeries":
        return cls({int(p): np.asarray(a, dtype=float)}, pmin=pmin, pmax=pmax)

    def copy(self) -> "LaurentSeries":
        return LaurentSeries({p: np.array(a, copy=True) for p, a in self.coeffs.items()}, self.pmin, self.pmax)

    def powers(self) -> list[int]:
        return sorted(self.coeffs)

    def coeff(self, p: int, like=None):
        if p in self.coeffs:
            return self.coeffs[p]
        if like is not None:
            return np.zeros_like(like, dtype=float)
        if self.coeffs:
            return np.zeros_like(next(iter(self.coeffs.values())), dtype=float)
        return np.asarray(0.0)

    def leading_power(self, tol: float = 1e-14) -> int:
        for p in sorted(self.coeffs):
            a = self.coeffs[p]
            if np.nanmax(np.abs(a)) > tol:
                return p
        raise ValueError("zero Laurent series has no leading power")

    def evaluate(self, eps: float):
        if not self.coeffs:
            return np.asarray(0.0)
        out = np.zeros_like(next(iter(self.coeffs.values())), dtype=float)
        for p, a in self.coeffs.items():
            out = out + a * (eps ** p)
        return out

    def map_coefficients(self, fn: Callable[[np.ndarray], np.ndarray]) -> "LaurentSeries":
        return LaurentSeries({p: fn(a) for p, a in self.coeffs.items()}, self.pmin, self.pmax)

    def shift(self, dp: int) -> "LaurentSeries":
        return LaurentSeries({p + dp: a for p, a in self.coeffs.items()}, self.pmin, self.pmax)

    def scale(self, c) -> "LaurentSeries":
        c = np.asarray(c, dtype=float)
        return LaurentSeries({p: _mul_arrays(a, c) for p, a in self.coeffs.items()}, self.pmin, self.pmax)

    def __neg__(self):
        return self.scale(-1.0)

    def __add__(self, other):
        if not isinstance(other, LaurentSeries):
            other = LaurentSeries.const(other, pmin=self.pmin, pmax=self.pmax)
        pmin = min(self.pmin, other.pmin); pmax = max(self.pmax, other.pmax)
        out: Dict[int, np.ndarray] = {}
        for p in set(self.coeffs) | set(other.coeffs):
            a = self.coeffs.get(p)
            b = other.coeffs.get(p)
            if a is None: out[p] = np.array(b, copy=True)
            elif b is None: out[p] = np.array(a, copy=True)
            else: out[p] = a + b
        return LaurentSeries(out, pmin, pmax)

    __radd__ = __add__

    def __sub__(self, other):
        return self + (-other if isinstance(other, LaurentSeries) else -np.asarray(other, dtype=float))

    def __rsub__(self, other):
        return (-self) + other

    def __mul__(self, other):
        if not isinstance(other, LaurentSeries):
            return self.scale(other)
        pmin = max(self.pmin, other.pmin)
        pmax = min(self.pmax, other.pmax)
        out: Dict[int, np.ndarray] = {}
        for p, a in self.coeffs.items():
            for q, b in other.coeffs.items():
                k = p + q
                if pmin <= k <= pmax:
                    v = _mul_arrays(a, b)
                    out[k] = out[k] + v if k in out else np.array(v, copy=True)
        return LaurentSeries(out, pmin, pmax)

    __rmul__ = __mul__

    def reciprocal(self, *, tol: float = 1e-14) -> "LaurentSeries":
        """Formal pointwise reciprocal of a scalar-valued Laurent series.

        The leading coefficient must be nonzero at every radial grid point.
        """
        e0 = self.leading_power(tol=tol)
        a0 = self.coeffs[e0]
        if np.nanmin(np.abs(a0)) <= tol:
            raise ZeroDivisionError(f"Laurent leading coefficient nearly zero at power {e0}")
        # A = eps^e0 sum_{n>=0} a_n eps^n.  Compute B = eps^-e0 sum b_n eps^n.
        max_n = self.pmax + e0
        if max_n < 0:
            max_n = 0
        b: Dict[int, np.ndarray] = {0: 1.0 / a0}
        for n in range(1, max_n + 1):
            s = np.zeros_like(a0, dtype=float)
            for k in range(1, n + 1):
                ak = self.coeffs.get(e0 + k)
                if ak is not None and (n - k) in b:
                    s = s + ak * b[n - k]
            b[n] = -s / a0
        return LaurentSeries({n - e0: v for n, v in b.items()}, self.pmin, self.pmax)

    def __truediv__(self, other):
        if isinstance(other, LaurentSeries):
            return self * other.reciprocal()
        return self.scale(1.0 / np.asarray(other, dtype=float))

    def __rtruediv__(self, other):
        return LaurentSeries.const(other, pmin=self.pmin, pmax=self.pmax) / self


def series_sum(items: Iterable[LaurentSeries], *, pmin=-6, pmax=14) -> LaurentSeries:
    out = LaurentSeries({}, pmin=pmin, pmax=pmax)
    for x in items:
        out = out + x
    return out


def outer_product(a: LaurentSeries, b: LaurentSeries, c: LaurentSeries | None = None) -> LaurentSeries:
    """Series of c * outer(a,b) for arrays a,b shaped (...,m)."""
    if c is None:
        c = LaurentSeries.const(1.0, pmin=a.pmin, pmax=a.pmax)
    out: Dict[int, np.ndarray] = {}
    pmin = max(a.pmin, b.pmin, c.pmin); pmax = min(a.pmax, b.pmax, c.pmax)
    for pa, av in a.coeffs.items():
        for pb, bv in b.coeffs.items():
            for pc, cv in c.coeffs.items():
                p = pa + pb + pc
                if pmin <= p <= pmax:
                    v = av[..., :, None] * bv[..., None, :] * np.asarray(cv)[..., None, None]
                    out[p] = out[p] + v if p in out else np.array(v, copy=True)
    return LaurentSeries(out, pmin, pmax)


def mat_series_add(a: LaurentSeries, b: LaurentSeries) -> LaurentSeries:
    return a + b


def derivative(series: LaurentSeries, fn: Callable[[np.ndarray], np.ndarray]) -> LaurentSeries:
    return series.map_coefficients(fn)
