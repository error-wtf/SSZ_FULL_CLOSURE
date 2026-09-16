import numpy as np
import pandas as pd

from ..policy import numerical_policy
from ..provenance.blacklist import require_production_input
from ..types import P5Background


def validate_background(bg):
    n = len(bg.r)
    for key in ("r", "u", "phi", "f", "h", "phi_r", "X", "A0prime"):
        a = np.asarray(getattr(bg, key))
        if a.shape != (n,) or not np.all(np.isfinite(a)):
            raise ValueError(f"invalid background {key}")
    if not n or not np.all(np.diff(bg.r) > 0) or np.any(bg.r <= 0):
        raise ValueError(
            "r must be positive and strictly increasing; center is a separate analytic chart"
        )
    if np.any(bg.f <= 0) or np.any(bg.h <= 0):
        raise ValueError("nonpositive metric")
    if np.max(np.abs(bg.h * bg.phi_r**2 + 2 * bg.X)) > numerical_policy()["branch_identity_abs"]:
        raise ValueError("kappa != -2X")
    return bg


def from_frame(d):
    r = d["x"].to_numpy(float)
    ph = d["phi_r" if "phi_r" in d else "phiprime"].to_numpy(float)
    h = d.h.to_numpy(float)
    return validate_background(
        P5Background(
            r,
            d.u.to_numpy(float),
            d.phi.to_numpy(float),
            d.f.to_numpy(float),
            h,
            ph,
            d.X.to_numpy(float) if "X" in d else -h * ph**2 / 2,
            d.A0prime.to_numpy(float) if "A0prime" in d else np.zeros(len(d)),
        )
    )


def load_frozen_p5(path):
    require_production_input(path)
    return from_frame(pd.read_csv(path).sort_values("x").reset_index(drop=True))
