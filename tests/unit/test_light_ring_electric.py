from __future__ import annotations

import numpy as np

from ssz_p5.action.light_ring_electric import light_ring_witnesses, required_electric_q_profile


def test_light_ring_identity_is_denominator_free_and_electric_control_is_finite():
    # Synthetic root at r=2 with f=1 and f'=1; use f''=1 and h=H=a4=1.
    r = np.array([1.9, 2.1])
    f = np.ones(2)
    h = np.ones(2)
    H = np.ones(2)
    a4 = np.ones(2)
    fp = np.array([0.95, 1.05])
    fpp = np.ones(2)
    w = light_ring_witnesses(r, f, h, H, a4, fp, fpp)
    assert len(w) == 1
    assert np.isfinite(w[0].zero_vector_required_tensor_f)
    assert np.isfinite(w[0].electric_q_for_f_equals_h)

    a4p = np.zeros(2)
    q = required_electric_q_profile(r, f, h, H, a4, a4p, fp, fpp)
    assert np.all(np.isfinite(q))
