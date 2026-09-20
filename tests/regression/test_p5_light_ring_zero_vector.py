from __future__ import annotations

from pathlib import Path

import pandas as pd

from ssz_p5.action.light_ring_electric import light_ring_witnesses
from ssz_p5.jets.jet9d8 import derivative

ROOT = Path(__file__).resolve().parents[2]


def test_inner_light_ring_excludes_positive_zero_vector_tensor_branch():
    d = pd.read_csv(ROOT / "data/production/ssz_p5_horndeski_carrier_through_light_rings_to_core_2026-09-12.csv")
    r = d["r_over_rs"].to_numpy(float)
    f = d["f"].to_numpy(float)
    h = d["h"].to_numpy(float)
    H = d["H_equals_F_equals_G"].to_numpy(float)
    a4 = d["a4_Horndeski"].to_numpy(float)
    fp = derivative(r, f)
    rings = light_ring_witnesses(r, f, h, H, a4, fp, derivative(r, f, 2))
    assert len(rings) == 2
    inner, outer = sorted(rings, key=lambda x: x.radius)
    assert abs(inner.radius - 1.41616) < 5e-5
    assert inner.a4 > 0 and inner.tensor_h > 0
    assert inner.zero_vector_required_tensor_f < 0
    assert inner.zero_vector_required_f_over_h < -1.6
    assert inner.electric_q_for_f_equals_h > 1.5
    assert abs(outer.radius - 1.5) < 5e-5
    assert abs(outer.zero_vector_required_f_over_h - 1.0) < 1e-4
