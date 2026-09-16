from __future__ import annotations

import numpy as np

from ..config import SLOT_NAMES
from ..types import Coefficients41, GateResult


def validate_41_schema(c: Coefficients41) -> list[GateResult]:
    out = []
    missing = [s for s in SLOT_NAMES if s not in c.slots]
    out.append(
        GateResult(
            "41_slots_present",
            "PASS" if not missing else "FAIL",
            missing,
            "no missing slots",
            "schema41",
        )
    )
    bad = [s for s in SLOT_NAMES if s in c.slots and not np.all(np.isfinite(c.slots[s]))]
    shapes = [
        s
        for s in SLOT_NAMES
        if s in c.slots and np.asarray(c.slots[s]).shape != np.asarray(c.background.r).shape
    ]
    out.append(
        GateResult("41_slot_shapes", "FAIL" if shapes else "PASS", shapes, "(N,)", "schema41")
    )
    out.append(
        GateResult("41_slots_finite", "PASS" if not bad else "FAIL", bad, "all finite", "schema41")
    )
    return out
