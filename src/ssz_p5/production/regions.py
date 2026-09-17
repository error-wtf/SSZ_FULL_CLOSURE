"""Authoritative regional production classification.

This mapping is locked for the current Full-SVT production member.  Witness and
historical datasets may overlap these radii, but they do not override this map.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

import numpy as np
import pandas as pd

CENTRAL = "data/production/ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv"

OUTER_START = 0.5515230871346237
CENTRAL_START = 0.61
INNER_START = 0.71
CORE_START = 0.715


class ProductionRegion(StrEnum):
    WEAK_EXTERIOR_H = "weak_exterior_H"
    OUTER_SAME_ACTION_H_SVT = "outer_same_action_H_SVT"
    CENTRAL_EXACT_SVT = "central_exact_SVT"
    INNER_SAME_ACTION_SVT_H = "inner_same_action_SVT_H"
    PUNCTURED_H_CORE = "punctured_H_core"
    ANALYTIC_CENTER = "analytic_center"


@dataclass(frozen=True)
class RegionSpec:
    region: ProductionRegion
    u_min: float | None
    u_max: float | None
    notes: str


def classify_u(u: float) -> ProductionRegion:
    """Return the selected production region for u=r_s/r.

    The exact center is not represented by a finite u and is handled by the
    analytic-center chart separately.
    """
    u = float(u)
    if not math.isfinite(u) or u < 0:
        raise ValueError("u must be finite and non-negative; use classify_x(0) for center")
    if u < OUTER_START:
        return ProductionRegion.WEAK_EXTERIOR_H
    if u < CENTRAL_START:
        return ProductionRegion.OUTER_SAME_ACTION_H_SVT
    if u < INNER_START:
        return ProductionRegion.CENTRAL_EXACT_SVT
    if u < CORE_START:
        return ProductionRegion.INNER_SAME_ACTION_SVT_H
    return ProductionRegion.PUNCTURED_H_CORE


def classify_x(x: float) -> ProductionRegion:
    """Classify x=r/r_s by converting to u=1/x."""
    x = float(x)
    if not math.isfinite(x):
        raise ValueError("x must be finite")
    if x <= 0:
        if x == 0:
            return ProductionRegion.ANALYTIC_CENTER
        raise ValueError("x must be non-negative")
    return classify_u(1.0 / x)


def region_table() -> tuple[RegionSpec, ...]:
    return (
        RegionSpec(
            ProductionRegion.WEAK_EXTERIOR_H,
            None,
            OUTER_START,
            "weak/exterior Horndeski production member",
        ),
        RegionSpec(
            ProductionRegion.OUTER_SAME_ACTION_H_SVT,
            OUTER_START,
            CENTRAL_START,
            "re-solved C-infinity outer Horndeski+SVT handover",
        ),
        RegionSpec(
            ProductionRegion.CENTRAL_EXACT_SVT,
            CENTRAL_START,
            INNER_START,
            "explicit genuine-SVT production lobe; A0prime is generally nonzero",
        ),
        RegionSpec(
            ProductionRegion.INNER_SAME_ACTION_SVT_H,
            INNER_START,
            CORE_START,
            "re-solved C-infinity inner SVT+Horndeski handover",
        ),
        RegionSpec(
            ProductionRegion.PUNCTURED_H_CORE, CORE_START, None, "punctured general Horndeski core"
        ),
        RegionSpec(
            ProductionRegion.ANALYTIC_CENTER, None, None, "exact r=0 Taylor/analytic center chart"
        ),
    )


def audit_light_ring_regions(root: Path) -> dict:
    d = pd.read_csv(Path(root) / CENTRAL)
    out = []
    for x in (1.5999977950251312, 1.41616064):
        row = d.iloc[int(np.argmin(np.abs(d.x.to_numpy(float) - x)))]
        out.append(
            {
                "x_requested": x,
                "x": float(row.x),
                "u": float(row.u),
                "region": classify_u(float(row.u)),
                "A0prime": float(row.A0prime),
                "pure_H_null_vector_gate": "NOT_APPLICABLE",
            }
        )
    return {
        "cover": {
            "weak_exterior_H": "u < 0.5515230871346237",
            "outer_same_action_H_SVT": "0.5515230871346237 <= u < 0.61",
            "central_exact_SVT": "0.61 <= u < 0.71",
            "inner_same_action_SVT_H": "0.71 <= u < 0.715",
            "punctured_H_core": "u >= 0.715",
        },
        "witnesses": out,
        "strong_H_carrier_role": "principal/control witness; never global production background",
    }
