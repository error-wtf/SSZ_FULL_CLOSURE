"""Physical units bridge: frozen mass scale -> observer frequencies.

Contract rule 6: the mass scale must come from an EXTERNAL, frozen prior —
never tuned mode-by-mode.  The conversion for the dimensionless SSZ
convention (omega_bar = omega * r_s / c) is

    f = c^3 / (4 pi G M) * Re(omega_bar)

with c, G fixed CODATA-style constants and M in kg.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from hashlib import sha256

# CODATA 2018 values (exact/frozen constants only; G is the 2018 measured value)
C_LIGHT = 299792458.0          # m/s (exact)
G_NEWTON = 6.67430e-11         # m^3 kg^-1 s^-2
M_SUN = 1.98892e30             # kg (frozen legacy solar mass parameter)


@dataclass(frozen=True)
class FrozenMassPrior:
    """An externally sourced, hash-frozen mass prior."""
    source: str          # e.g. "HEASARC/MAXI J1820+070 dynamical prior v1"
    mass_solar: float
    mass_solar_sigma: float
    frozen_note: str

    def to_json(self) -> dict:
        d = asdict(self)
        d["sha256"] = sha256(
            json.dumps(d, sort_keys=True).encode()).hexdigest()
        return d


def omega_bar_to_hz(omega_bar_re: float, mass_solar: float) -> float:
    """f = c^3/(4 pi G M) Re(omega_bar)."""
    m = mass_solar * M_SUN
    return C_LIGHT**3 / (4.0 * 3.141592653589793 * G_NEWTON * m) * omega_bar_re


def unit_conversion_factor_hz(mass_solar: float) -> float:
    """The single frozen number f/Re(omega_bar) for a given mass."""
    return omega_bar_to_hz(1.0, mass_solar)
