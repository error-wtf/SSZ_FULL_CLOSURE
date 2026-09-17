"""Bind production to the regional member rather than the historical null witness."""

import json
from pathlib import Path

import numpy as np

from .regions import classify_x, region_table

MEMBER_FILE = "SSZ_P5_REGIONAL_PRODUCTION_MEMBER_2026-09-17.json"


def load_regional_member(root: Path) -> dict:
    member = json.loads((root / MEMBER_FILE).read_text())
    expected = [(s.region.value, s.u_min, s.u_max) for s in region_table()]
    actual = [(s["region"], s["u_min"], s["u_max"]) for s in member["production_regions"]]
    if actual != expected:
        raise ValueError("regional member and resolver disagree")
    locked = member["locked_conventions"]
    if locked["v12"] != "-v6/(2h)" or locked["derivatives"] != "JET9D8 window=9 degree=8":
        raise ValueError("regional member changes locked conventions")
    return member


def validate_stream_regions(frame):
    """Check row coordinates and regional branches, allowing electric SVT rows."""
    x = frame["x"].to_numpy(float)
    if not np.isfinite(x).all() or np.any(x < 0):
        raise ValueError("invalid radial coordinates")
    expected = np.array([classify_x(v).value for v in x])
    if not np.array_equal(frame["region"].to_numpy(str), expected):
        raise ValueError("stream region does not match locked radial cover")
    if "u" in frame:
        positive = x > 0
        if not np.allclose(frame.loc[positive, "u"], 1 / x[positive], rtol=1e-12, atol=1e-13):
            raise ValueError("inconsistent x/u coordinates")
    ap = frame["A0prime"].to_numpy(float)
    if not np.isfinite(ap).all():
        raise ValueError("nonfinite electric background")
    pure = np.isin(expected, ["weak_exterior_H", "punctured_H_core", "analytic_center"])
    if np.any(np.abs(ap[pure]) > 1e-12):
        raise ValueError("electric background in a selected pure-H region")
    return expected
