"""Direct zero-vector HSVT member: Horndeski P5 + epsilon_Y * Y.

For A0'=0 and constant epsilon_Y the deformation is exactly background-null and
changes the common 41-slot even action only through the pure vector coefficients
v1 and v10.  c3=e3=v5 remain zero on the pure Horndeski representative.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd

from ..jets.jet9d8 import profile_derivative

EPSILON_Y = 1.0e-2
REGION_SOURCES = {
    "exterior": "data/regression/ssz_p5_F2_exterior_horndeski_unreduced_39of41_2026-09-14.csv",
    "carrier": "data/authoritative/ssz_p5_F2_horndeski_carrier_unreduced_39of41_CORRECTED_2026-09-14.csv",
    "core": "data/authoritative/ssz_p5_F2_core_punctured_horndeski_unreduced_39of41_2026-09-14.csv",
}


def apply_eps_y(frame: pd.DataFrame, *, epsilon: float = EPSILON_Y, window: int = 9, degree: int = 8) -> pd.DataFrame:
    d = frame.copy().sort_values("x").reset_index(drop=True)
    r = d.x.to_numpy(float)
    if np.any(np.diff(r) <= 0):
        raise ValueError("strictly increasing x required")
    if "phi_r" not in d:
        if "phiprime" not in d:
            raise KeyError("phi_r/phiprime")
        d["phi_r"] = d.phiprime
    ph, f, h = d.phi_r.to_numpy(float), d.f.to_numpy(float), d.h.to_numpy(float)
    dr = lambda y, order=1: profile_derivative(r, np.asarray(y, float), order, window, degree)
    d["A0prime"] = 0.0
    d["phiprime"] = ph
    d["c3"] = 0.0
    d["e3"] = 0.0
    d["v5"] = 0.0
    # Locked same-action identities on the zero-vector branch.
    d["a5"] = dr(d.a2) - dr(d.a1, 2)
    d["v12"] = -d.v6 / (2.0 * d.h)
    # Delta f2 = epsilon Y, Appendix A at A0'=0, f3=f4=0.
    d["v1"] = d.v1.to_numpy(float) - r**2 * h**1.5 * ph**2 * epsilon / np.sqrt(f)
    d["v10"] = d.v10.to_numpy(float) + np.sqrt(f * h) * h * ph**2 * epsilon
    d["kappa"] = h * ph**2
    d["Z_A"] = 1.0 - 2.0 * epsilon * d.kappa
    slots = [*(f"a{i}" for i in range(1,10)),*(f"b{i}" for i in range(1,6)),*(f"c{i}" for i in range(1,7)),*(f"d{i}" for i in range(1,5)),*(f"e{i}" for i in range(1,5)),*(f"v{i}" for i in range(1,14))]
    if not np.isfinite(d[slots].to_numpy(float)).all():
        raise ValueError("nonfinite HSVT coefficient stream")
    return d


def build_region(root: Path, region: str, *, epsilon: float = EPSILON_Y) -> pd.DataFrame:
    if region not in REGION_SOURCES:
        raise KeyError(region)
    return apply_eps_y(pd.read_csv(root / REGION_SOURCES[region]), epsilon=epsilon)
