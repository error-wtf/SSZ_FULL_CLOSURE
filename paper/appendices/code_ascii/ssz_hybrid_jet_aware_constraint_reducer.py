#!/usr/bin/env python3
"""
Jet-aware first-stage constraint reducer for the common Horndeski + U(1)-SVT
even-parity quadratic action.

Why this module exists
----------------------
For the constrained even system, freezing radial coefficient functions before
eliminating H0,H1,H2,h1,dA0,dA1 loses principal information.  The published
reduced kinetic matrix contains radial jets such as a4' and K2'.

This module therefore works on radial coefficient PROFILES and performs the
generalized-psi H0 reduction while retaining coefficient derivatives.

Common generalized dynamical variable
--------------------------------------
    psi = H2 + q h1 + p dphi'
    p = a1/a3
    q = L a4/a3

For EH + genuine SVT:
    a1 = 0, a3 = -r a4
so
    psi = H2 - (L/r) h1.

After substituting
    H2 = psi - q h1 - p dphi'
into the H0 constraint, the dphi'' and h1' terms cancel identically.

The remaining H0 constraint is algebraic in h1:
    0 = a3 psi'
        + B_eff dphi'
        + C_phi dphi
        + C_H2 psi
        + D_h1 h1
        + v2 V

with
    B_eff = B - a3 p' - C_H2 p
    B      = a2 - v2 v4/(2 v1)
    C_phi  = a5 + L a6 - v2 v5/(2 v1)
    C_H2   = a7 + L a8 - v2 v3/(2 v1)
    C_h1   = L (a9 - v2 v6/(2 v1))
    D_h1   = C_h1 - a3 q' - C_H2 q.

Thus, where D_h1 != 0,
    h1 = -(a3 psi' + B_eff dphi' + C_phi dphi + C_H2 psi + v2 V)/D_h1.

This is the first elimination step for a common same-action HSVT reducer.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Dict
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline

REQUIRED = (
    "a1","a2","a3","a4","a5","a6","a7","a8","a9",
    "v1","v2","v3","v4","v5","v6"
)

def radial_derivative(r: np.ndarray, y: np.ndarray) -> np.ndarray:
    r = np.asarray(r, float)
    y = np.asarray(y, float)
    order = np.argsort(r)
    rr = r[order]
    yy = y[order]
    if np.any(np.diff(rr) <= 0):
        raise ValueError("r values must be unique")
    dy = CubicSpline(rr, yy)(rr, 1)
    out = np.empty_like(dy)
    out[order] = dy
    return out

def generalized_psi_h0_reduction(df: pd.DataFrame, L: float,
                                 r_col: str = "x") -> pd.DataFrame:
    """Return the jet-aware H0 algebraic-constraint coefficients."""
    missing = [k for k in REQUIRED if k not in df.columns]
    if missing:
        raise KeyError(f"Missing coefficients: {missing}")

    r = df[r_col].to_numpy(float)
    qv = {k: df[k].to_numpy(float) for k in REQUIRED}

    if np.any(np.abs(qv["a3"]) < 1e-14):
        raise ZeroDivisionError("a3 crosses zero")
    if np.any(np.abs(qv["v1"]) < 1e-14):
        raise ZeroDivisionError("v1 crosses zero")

    p = qv["a1"]/qv["a3"]
    q = L*qv["a4"]/qv["a3"]
    pp = radial_derivative(r, p)
    qp = radial_derivative(r, q)

    B = qv["a2"] - qv["v2"]*qv["v4"]/(2*qv["v1"])
    Cphi = qv["a5"] + L*qv["a6"] - qv["v2"]*qv["v5"]/(2*qv["v1"])
    CH2 = qv["a7"] + L*qv["a8"] - qv["v2"]*qv["v3"]/(2*qv["v1"])
    Ch1 = L*(qv["a9"] - qv["v2"]*qv["v6"]/(2*qv["v1"]))

    Beff = B - qv["a3"]*pp - CH2*p
    Dh1 = Ch1 - qv["a3"]*qp - CH2*q

    # Explicit cancellation residuals before dropping the higher derivatives.
    dphi2_res = qv["a1"] - qv["a3"]*p
    h1prime_res = L*qv["a4"] - qv["a3"]*q

    result = pd.DataFrame({
        r_col: r,
        "p_a1_over_a3": p,
        "p_prime": pp,
        "q_La4_over_a3": q,
        "q_prime": qp,
        "B_raw": B,
        "B_eff": Beff,
        "C_phi": Cphi,
        "C_H2": CH2,
        "C_h1": Ch1,
        "D_h1_pivot": Dh1,
        "dphi_second_derivative_cancellation": dphi2_res,
        "h1_prime_cancellation": h1prime_res,
    })
    if "u" in df.columns:
        result.insert(0, "u", df["u"].to_numpy(float))
    return result

def h1_solution_coefficients(reduced: pd.DataFrame,
                             source_a3: np.ndarray,
                             source_v2: np.ndarray) -> pd.DataFrame:
    """Coefficients of h1 in terms of (psi', dphi', dphi, psi, V)."""
    D = reduced["D_h1_pivot"].to_numpy(float)
    return pd.DataFrame({
        "coeff_psi_prime": -np.asarray(source_a3,float)/D,
        "coeff_dphi_prime": -reduced["B_eff"].to_numpy(float)/D,
        "coeff_dphi": -reduced["C_phi"].to_numpy(float)/D,
        "coeff_psi": -reduced["C_H2"].to_numpy(float)/D,
        "coeff_V": -np.asarray(source_v2,float)/D,
    })
