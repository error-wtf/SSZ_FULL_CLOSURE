"""Sagnac reference transport - closed form (Route 1) and edge contracts.

This module is part of the INDEPENDENT reference layer (SAG_REF_V1).
It contains NO SSZ assumptions and MUST NOT import any SSZ geometry or
postclosure module (enforced by tests/test_reference_transport_imports.py).

Ground truth (galilean kinematics, inertial readout):

    t_+ = L / (c - v)              co-rotating
    t_- = L / (c + v)              counter-rotating
    dt  = 2 L v / (c^2 - v^2)

Dimensionless form (units c = L = 1):

    t_+ = 1/(1 - beta),  t_- = 1/(1 + beta),  dt = 2 beta/(1 - beta^2)

Phase / JIF readout at carrier f = omega / (2 pi):

    dphi = omega * dt,     dJ = f * dt = dphi / (2 pi)
"""
from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class SagnacResult:
    t_plus: float
    t_minus: float
    delta_t: float
    delta_phi: float
    delta_J: float


def sagnac_times(beta: float, f_carrier: float = 1.0) -> SagnacResult:
    """Closed-form Sagnac transport in units c = L = 1.

    beta = v/c in (-1, 1).  f_carrier in cycles per time unit.
    """
    if not -1.0 < beta < 1.0:
        raise ValueError(f"beta must satisfy |beta| < 1, got {beta}")
    if f_carrier <= 0:
        raise ValueError("f_carrier must be positive")
    t_p = 1.0 / (1.0 - beta)
    t_m = 1.0 / (1.0 + beta)
    dt = t_p - t_m
    omega = 2.0 * math.pi * f_carrier
    return SagnacResult(t_plus=t_p, t_minus=t_m, delta_t=dt,
                        delta_phi=omega * dt, delta_J=f_carrier * dt)


def area_form(beta: float, f_carrier: float = 1.0, radius: float = 1.0) -> float:
    """Area-form Sagnac: dt = 4 Omega A / (c^2 - v^2) with Omega = beta/radius.

    In units c = 1: dt = 4 (beta/radius) (pi radius^2) / (1 - beta^2)
    = 4 pi beta radius / (1 - beta^2).  With L = 2 pi radius this equals
    the segment form 2 beta/(1-beta^2) exactly when radius = 1.
    """
    omega = beta / radius
    area = math.pi * radius * radius
    return 4.0 * omega * area / (1.0 - beta * beta)


def edge_contract_null(beta: float = 1e-15) -> float:
    """G105: v -> 0 must give dt -> 0."""
    return sagnac_times(beta).delta_t


def edge_contract_reversal(beta: float) -> tuple[float, float]:
    """G104: dt(-beta) must equal -dt(beta)."""
    return sagnac_times(beta).delta_t, sagnac_times(-beta).delta_t


def edge_contract_catch_up(beta: float = 0.999999) -> float:
    """t_+ diverges as beta -> 1 (co-rotator never catches the beacon)."""
    return sagnac_times(beta).t_plus
