"""Action-jet -> Maxwell--Horndeski primitive bridge for the quartic G5=0 sector.

This module implements the subset of Kase--Tsujikawa Appendix A needed at the
P5 Inner -> punctured-Horndeski interface.  It does not import selected 41-slot
coefficients.  The primitive profiles F,G,H,a1,c4 are computed directly from
G3/G4 action jets on the fixed background.  It also exposes the exact linear
coefficient dc2/dG2XX, the background-null normal direction used to realize a
selected c2 profile without changing the background equations.

The deep subcore where G5 != 0 is intentionally outside this adapter and must
use a separate full-G5 bridge.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from ..jets.jet9d8 import derivative


def _col(d: pd.DataFrame, *names: str) -> np.ndarray:
    for name in names:
        if name in d:
            return d[name].to_numpy(float)
    raise ValueError(f"missing required column; expected one of {names!r}")


def quartic_g5zero_primitives(action: pd.DataFrame, *, window: int = 9, degree: int = 8) -> pd.DataFrame:
    """Derive F,G,H,a1,c4 and dc2/dG2XX from G5=0 action jets.

    Accepted geometry aliases are ``x``/``r_over_rs``, ``f``/``A_f`` and
    ``h``/``B_h``.  ``phi_r`` can be supplied directly; otherwise it is
    reconstructed from X=-h phi_r^2/2 using the negative P5 branch.
    """
    d = action.copy().sort_values("x" if "x" in action else "r_over_rs").reset_index(drop=True)
    r = _col(d, "x", "r_over_rs")
    f = _col(d, "f", "A_f")
    h = _col(d, "h", "B_h")
    if "phi_r" in d:
        ph = _col(d, "phi_r")
    else:
        X = _col(d, "X")
        ph = -np.sqrt(np.maximum(0.0, -2.0 * X / h))

    G4 = _col(d, "G4", "G4_background")
    G4X = _col(d, "G4X")
    G4XX = _col(d, "G4XX")
    G4phi = _col(d, "G4phi", "G4_phi")
    G4phiX = _col(d, "G4phiX")
    G3X = _col(d, "G3X")

    fp = derivative(r, f, 1, window=window, degree=degree)
    sqfh = np.sqrt(f * h)

    # Kase--Tsujikawa (2023), Eqs. (24)--(26), specialized to G5=0.
    F = 2.0 * G4
    H = 2.0 * G4 + 2.0 * h * ph**2 * G4X
    G = H.copy()

    # Appendix A, Eq. (A1), a1, specialized to G5=0.
    a1 = sqfh * (
        (G4phi + 0.5 * h * (G3X - 2.0 * G4phiX) * ph**2) * r**2
        + 2.0 * h * ph * (G4X - h * G4XX * ph**2) * r
    )

    # Appendix A, c4, specialized to G5=0.
    c4 = 0.25 * np.sqrt(f / h) * (
        (h * ph / f * (2.0 * G4X - 2.0 * h * G4XX * ph**2)) * fp
        + 4.0 * G4phi
        + 2.0 * h * (G3X - 2.0 * G4phiX) * ph**2
        + (4.0 * h * G4X * ph - 4.0 * h**2 * G4XX * ph**3) / r
    )

    # Exact coefficient of G2,XX in Appendix-A c2.  This is a transverse
    # normal jet: G2,XX is absent from the background equations on the curve.
    dc2_dG2XX = sqfh * (-0.5 * h * ph**3 * r**2)

    out = pd.DataFrame(
        {
            "x": r,
            "F_tensor_action": F,
            "G_tensor_action": G,
            "H_tensor_action": H,
            "a1_action": a1,
            "c4_action": c4,
            "dc2_dG2XX_action": dc2_dG2XX,
        }
    )
    if "u" in d:
        out.insert(0, "u", d["u"].to_numpy(float))
    return out
