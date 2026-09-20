"""Light-ring diagnostics for the P5 Maxwell-Horndeski/SVT background.

The necessary on-shell identity is Eq. (85) of Kase & Tsujikawa (2023),
implemented in :mod:`ssz_p5.action.onshell`.  This helper deliberately keeps the
identity denominator-free, so a light ring (2 f - r f' = 0) is a regular point
of the diagnostic rather than a singular division.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class LightRingWitness:
    radius: float
    f: float
    h: float
    tensor_h: float
    a4: float
    geometric_factor: float
    zero_vector_required_tensor_f: float
    zero_vector_required_f_over_h: float
    electric_q_for_f_equals_h: float


def _linear_root(x0: float, x1: float, y0: float, y1: float) -> tuple[float, float]:
    if y1 == y0:
        raise ValueError("degenerate root bracket")
    t = -y0 / (y1 - y0)
    if not 0.0 <= t <= 1.0:
        raise ValueError("root is not bracketed")
    return x0 + t * (x1 - x0), t


def _lerp(a0: float, a1: float, t: float) -> float:
    return a0 + t * (a1 - a0)


def light_ring_witnesses(
    r: np.ndarray,
    f: np.ndarray,
    h: np.ndarray,
    tensor_h: np.ndarray,
    a4: np.ndarray,
    f_prime: np.ndarray,
    f_second: np.ndarray,
) -> list[LightRingWitness]:
    """Return denominator-free Eq.-85 witnesses at all resolved light rings.

    ``electric_q_for_f_equals_h`` is the product ``A0prime**2 * v8`` required
    at the light ring if one demands ``tensor_f = tensor_h``.  It is a local
    target only; it does not by itself construct a global action member.
    """
    arrays = [np.asarray(v, dtype=float) for v in (r, f, h, tensor_h, a4, f_prime, f_second)]
    if any(v.ndim != 1 for v in arrays) or len({len(v) for v in arrays}) != 1:
        raise ValueError("all inputs must be equal-length one-dimensional arrays")
    r, f, h, tensor_h, a4, f_prime, f_second = arrays
    if len(r) < 2 or not np.all(np.diff(r) > 0):
        raise ValueError("r must be strictly increasing")
    if not all(np.all(np.isfinite(v)) for v in arrays):
        raise ValueError("nonfinite input")
    if np.any(f <= 0) or np.any(h <= 0):
        raise ValueError("light-ring diagnostic requires f>0 and h>0")

    denominator = 2.0 * f - r * f_prime
    crossings = np.flatnonzero(denominator[:-1] * denominator[1:] <= 0.0)
    out: list[LightRingWitness] = []
    for i in crossings:
        # Avoid duplicate roots if an exact zero occupies two adjacent brackets.
        if i > 0 and denominator[i] == 0.0 and denominator[i - 1] == 0.0:
            continue
        radius, t = _linear_root(r[i], r[i + 1], denominator[i], denominator[i + 1])
        ff = _lerp(f[i], f[i + 1], t)
        hh = _lerp(h[i], h[i + 1], t)
        HH = _lerp(tensor_h[i], tensor_h[i + 1], t)
        aa = _lerp(a4[i], a4[i + 1], t)
        fp = _lerp(f_prime[i], f_prime[i + 1], t)
        fpp = _lerp(f_second[i], f_second[i + 1], t)
        geom = radius * fpp - radius * fp**2 / ff + 2.0 * fp - 2.0 * ff / radius

        # At 2f-rf'=0, Eq. 85 reduces to
        #   0 = geom*a4 + f^(3/2) F/(r sqrt(h)) - 2 r f h q,
        # q := A0prime^2 v8.
        required_f = -geom * aa * radius * np.sqrt(hh) / ff**1.5
        q_equal = (
            geom * aa + ff**1.5 * HH / (radius * np.sqrt(hh))
        ) / (2.0 * radius * ff * hh)
        out.append(
            LightRingWitness(
                radius=float(radius),
                f=float(ff),
                h=float(hh),
                tensor_h=float(HH),
                a4=float(aa),
                geometric_factor=float(geom),
                zero_vector_required_tensor_f=float(required_f),
                zero_vector_required_f_over_h=float(required_f / HH),
                electric_q_for_f_equals_h=float(q_equal),
            )
        )
    return out


def required_electric_q_profile(
    r: np.ndarray,
    f: np.ndarray,
    h: np.ndarray,
    tensor_f: np.ndarray,
    a4: np.ndarray,
    a4_prime: np.ndarray,
    f_prime: np.ndarray,
    f_second: np.ndarray,
) -> np.ndarray:
    """Solve Eq. 85 pointwise for q=A0prime**2*v8.

    This produces a background-equation target profile for a chosen
    ``tensor_f``/``a4`` pair.  It is not an action reconstruction.
    """
    r, f, h, tensor_f, a4, a4_prime, f_prime, f_second = [
        np.asarray(v, dtype=float)
        for v in (r, f, h, tensor_f, a4, a4_prime, f_prime, f_second)
    ]
    geom = r * f_second - r * f_prime**2 / f + 2.0 * f_prime - 2.0 * f / r
    lhs = (2.0 * f - r * f_prime) * a4_prime
    return (geom * a4 + f**1.5 * tensor_f / (r * np.sqrt(h)) - lhs) / (2.0 * r * f * h)
