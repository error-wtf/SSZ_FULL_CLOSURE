"""Necessary background identity; no reconstruction or coefficient back-fitting.

Kase & Tsujikawa, arXiv:2301.10362, Eq. (85), follows from E00 and E22.
The denominator-multiplied form remains meaningful at 2*f-r*f'=0.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from ..jets.jet9d8 import derivative
from ..provenance.manifest import sha256

CARRIER = (
    "data/production/"
    "ssz_p5_horndeski_carrier_through_light_rings_to_core_2026-09-12.csv"
)
ACTION = "SSZ_P5_HSVT_ACTION_MEMBER_2026-09-16.json"
SOURCE = "https://arxiv.org/html/2301.10362"


def identity_terms(r, f, h, tensor_f, a4, a4_prime, f_prime, f_second,
                   electric_prime=0.0, v8=0.0):
    """Return the two sides of Eq. (85), without dividing by its denominator."""
    lhs = (2 * f - r * f_prime) * a4_prime
    rhs = (
        (r * f_second - r * f_prime**2 / f + 2 * f_prime - 2 * f / r) * a4
        + f**1.5 * tensor_f / (r * np.sqrt(h))
        - 2 * r * f * h * electric_prime**2 * v8
    )
    return lhs, rhs


def luminal_e00_minus_e22(r, f, h, tensor_h, h_prime, tensor_h_prime, f_prime, f_second):
    """Independent subtraction of Eqs. (7),(9),(10), G4(phi), G5=0, A0'=0.

    C1 and C6 cancel; thus all G2 and G3 freedom cancels. tensor_h = 2*G4.
    """
    return (
        -tensor_h * h_prime / (2 * r)
        - h * tensor_h_prime / r
        - tensor_h * (h - 1) / r**2
        + h * tensor_h * f_second / (2 * f)
        - h * tensor_h * f_prime**2 / (4 * f**2)
        + tensor_h * f_prime * h_prime / (4 * f)
        + h * tensor_h_prime * f_prime / (2 * f)
        + tensor_h * h * f_prime / (2 * r * f)
    )


def audit_frozen_carrier(root: Path):
    """Recompute an interior witness from fixed production data, never old 41 slots.

    Decimation and action-chain-rule routes estimate numerical sensitivity.
    Their spread is not a rigorous interval enclosure. A conservative 100x
    spread is reported separately from the unchanged background policy.
    """
    root = Path(root)
    action = json.loads((root / ACTION).read_text())
    if action["background"]["vector_branch"] != "A0prime = 0":
        raise ValueError("This diagnostic requires the declared zero-vector production branch")
    d = pd.read_csv(root / CARRIER)
    r, f, h, phi, H, G4, G4phi, X = (
        d[key].to_numpy(float)
        for key in ("r_over_rs", "f", "h", "phi", "H_equals_F_equals_G",
                    "G4", "G4_phi", "X")
    )
    if not all(np.all(np.isfinite(v)) for v in (r, f, h, phi, H, G4, G4phi, X)):
        raise ValueError("nonfinite frozen carrier")
    if not (np.all(np.diff(r) > 0) and np.all(f > 0) and np.all(h > 0)
            and np.all(X < 0)):
        raise ValueError("invalid carrier chart")
    if not np.all(np.diff(phi) < 0):
        raise ValueError("scalar derivative sign must be taken from the frozen monotone chart")
    if not np.allclose(H, 2 * G4, rtol=0, atol=1e-12):
        raise ValueError("declared luminal carrier requires H=2G4")
    i = int(np.argmin(np.abs(r - 1.6)))
    if i < 20 or i >= len(r) - 20 or abs(r[i] - 1.6) > 0.001:
        raise ValueError("reference radius not resolved in the interior")
    a4 = 0.5 * np.sqrt(f * h) * H
    records = []
    for stride in (1, 2, 3, 4):
        indices = np.arange(i % stride, len(r), stride)
        x, ff, hh, aa, HH = (a[indices] for a in (r, f, h, a4, H))
        j = int(np.flatnonzero(indices == i)[0])
        lhs, rhs = identity_terms(
            x, ff, hh, HH, aa, derivative(x, aa), derivative(x, ff),
            derivative(x, ff, 2),
        )
        records.append({
            "method": "JET9D8", "stride": stride, "radius": float(x[j]),
            "lhs": float(lhs[j]), "rhs": float(rhs[j]),
            "residual": float(lhs[j] - rhs[j]),
        })
    # This route uses the frozen action jet and scalar-kinetic column, rather
    # than numerical second derivatives of f or numerical first derivatives of H.
    phiprime = -np.sqrt(-2 * X / h)
    fprime = -2 * phiprime / (1 + phi)**3
    fsecond = (6 * phiprime**2 / (1 + phi)**4
               - 2 * derivative(r, phiprime) / (1 + phi)**3)
    hprime = derivative(r, h)
    Hprime = 2 * G4phi * phiprime
    a4prime = 0.5 * np.sqrt(f * h) * (Hprime + 0.5 * H * (fprime / f + hprime / h))
    lhs, rhs = identity_terms(r, f, h, H, a4, a4prime, fprime, fsecond)
    difference = luminal_e00_minus_e22(r, f, h, H, hprime, Hprime, fprime, fsecond)
    records.append({
        "method": "action jets and scalar chain rule; remaining derivatives JET9D8",
        "radius": float(r[i]), "lhs": float(lhs[i]), "rhs": float(rhs[i]),
        "residual": float(lhs[i] - rhs[i]),
        "E00_minus_E22": float(difference[i]),
        "independent_subtraction_residual": float(
            (lhs - rhs + r * f * np.sqrt(f / h) * difference)[i]
        ),
    })
    values = np.array([item["residual"] for item in records])
    spread = float(np.ptp(values))
    policy = json.loads((root / "NUMERICAL_POLICY.json").read_text())
    tolerance = float(policy["background_residual_abs"])
    envelope = 100 * spread
    failed = bool(np.min(np.abs(values)) > tolerance + envelope)
    return {
        "status": "FAIL_ON_SHELL_IDENTITY" if failed else "NO_RESOLVED_CONTRADICTION",
        "scope": "Necessary frozen-carrier background identity, not a spectrum or theory no-go",
        "reference": {"url": SOURCE, "equations": [7, 9, 10, 24, 25, 85]},
        "inputs": [
            {"path": name, "sha256": sha256(root / name)}
            for name in (CARRIER, ACTION, "NUMERICAL_POLICY.json")
        ],
        "vector_branch": action["background"]["vector_branch"],
        "epsilon_Y": action["action"]["epsilon_Y"],
        "deformation_background_contribution": "zero: field strength is zero",
        "row_index_zero_based": i,
        "csv_line_one_based": i + 2,
        "radius": float(r[i]),
        "input_row": {key: float(value) for key, value in d.iloc[i].items()},
        "routes": records,
        "route_spread": spread,
        "diagnostic_sensitivity_envelope_100x_spread": envelope,
        "sensitivity_is_rigorous_interval_bound": False,
        "unchanged_background_policy_tolerance": tolerance,
        "consistency_checks": {
            "max_abs_f_minus_scalar_closure": float(np.max(np.abs(f - (1 + phi)**-2))),
            "max_abs_H_minus_2G4": float(np.max(np.abs(H - 2 * G4))),
            "witness_Hprime_action_minus_JET9D8": float(Hprime[i] - derivative(r, H)[i]),
            "witness_phiprime_X_minus_JET9D8": float(phiprime[i] - derivative(r, phi)[i]),
        },
        "interpretation": (
            "The prescribed f,h,F=H and A0'=0 violate a necessary on-shell identity. "
            "G2/G3 or core-only changes cannot remove this interior residual while "
            "these prescribed functions remain fixed. No replacement member is selected."
            if failed else "This test alone does not establish a globally on-shell action."
        ),
    }
