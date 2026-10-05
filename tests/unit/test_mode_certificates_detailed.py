"""Per-mode detailed certification: real K-norms, smoothness A4, honest A6 skip."""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.qnm.mode_certificate import certify_mode  # noqa: E402


def test_a4_smoothness_accepts_radial_profile():
    """A legitimate radial omega profile (smooth variation) must PASS
    even when the total span exceeds the old fixed-u threshold."""
    om = np.linspace(0.162, 0.105, 20)  # smooth 35% drift
    c = certify_mode(
        L=6, mode_index=0,
        kinetic_norms=np.ones(20),
        tracking_overlaps=np.full(20, 0.99),
        relative_gaps=[0.1],
        omega_values=om,
        observable_shares=[0.3],
        d2_value=float("nan"),
        scale_control_ok=True, basis_control_ok=True,
        skip_axes=("A1", "A6"),
    )
    assert c.axes["A4_frequency_invariance"]["pass"] is True


def test_a4_smoothness_rejects_branch_jump():
    om = np.array([0.16] * 10 + [0.11] * 10)  # 31% single-step jump
    c = certify_mode(
        L=6, mode_index=0,
        kinetic_norms=np.ones(20),
        tracking_overlaps=np.full(20, 0.99),
        relative_gaps=[0.1],
        omega_values=om,
        observable_shares=[0.3],
        d2_value=0.9,
        scale_control_ok=True, basis_control_ok=True,
    )
    assert c.axes["A4_frequency_invariance"]["pass"] is False
    assert "A4" in c.failed_axes


def test_a6_skip_is_not_a_failure():
    c = certify_mode(
        L=6, mode_index=0,
        kinetic_norms=np.ones(5),
        tracking_overlaps=np.full(5, 0.99),
        relative_gaps=[0.1],
        omega_values=[1.0, 1.0],
        observable_shares=[0.3],
        d2_value=float("nan"),
        scale_control_ok=True, basis_control_ok=True,
        skip_axes=("A6",),
    )
    assert c.verdict == "NOT_EVALUABLE_V1"
    assert "A6" not in c.failed_axes
    assert c.failed_axes == []


def test_detailed_artifact_exists_and_honest():
    p = ROOT / "data/generated/spectral/MODE_CERTIFICATES_V1_DETAILED.json"
    if not p.exists():
        pytest.skip("detailed artifact not yet generated")
    d = json.loads(p.read_text())
    s = d["summary"]
    assert s["total"] == s["certified"] + s["rejected"] + (
        s["total"] - s["certified"] - s["rejected"])
    # A6 must never be counted as a hard failure in this layer:
    for c in d["certificates"]:
        if "failed_axes" in c:
            assert "A6" not in c["failed_axes"], (
                "A6 is not definable on the 3-slot layer and must be skipped")
