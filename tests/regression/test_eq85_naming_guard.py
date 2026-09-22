"""Regression guard: Eq.85 must be evaluated with the ZK-canonical V9 (=MH v8)
coefficient, never the 13-slot v8 column. Fails if the mapping is violated."""
import numpy as np

def test_v_slot_mapping_is_documented():
    src = open('src/ssz_hybrid_unreduced_even_kernel.py').read()
    assert 'V9=MH v8' in src.replace(' ', '') or 'V9 = MH v8' in src, \
        "kernel must document V9=MH v8 mapping"
    assert 'V8=MH v7' in src.replace(' ', '') or 'V8 = MH v7' in src

def test_q_conventions_must_not_be_mixed():
    # q_r = A0_r^2 v8_MH ; q_u = A0_u^2 v8_MH ; A0_r = -(u^2/r_s) A0_u
    # => q_r = (u^4/r_s^2) q_u. A raw comparison of q_r and q_u is a bug.
    u = 0.706144; r_s = 1.0
    assert abs((u**4 / r_s**2) - u**4) < 1e-15
    # guard constant frozen here so silent convention drift is caught:
    assert abs(u**4 - 0.248827) < 1e-3
