"""Regression: operator symmetry table must match the golden producer.

The golden reducer (canonical_audit) exports K=-P20/2 and G=+P02/2 with
P20/P02 symmetric by construction, S the antisymmetric connection block,
and M mixed (symmetric P00 part + antisymmetric S' part).  The previous
sign table demanded K,G,M antisymmetric and made reduce_profile unusable
on the producer's own golden output.
"""

import numpy as np

from ssz_p5.reducer.canonical import validate_operator
from ssz_p5.types import ConstraintPivots, ReducedOperator

R = np.linspace(1.4, 1.6, 40)
N = len(R)


def _op(k_asym=0.0, g_asym=0.0, s_sym=0.0):
    def sym():
        return np.tile(np.eye(3), (N, 1, 1))

    K = sym().copy()
    K[:, 0, 1] += k_asym
    G = sym().copy()
    G[:, 1, 2] += g_asym
    S = np.zeros((N, 3, 3))
    S[:, 0, 1] += 1.0
    S[:, 1, 0] -= 1.0
    S[:, 0, 1] += s_sym
    M = sym().copy()  # mixed is allowed -> symmetric part passes finiteness
    return ReducedOperator(
        6, R, K, np.zeros((N, 3, 3)), G, S, M,
        ConstraintPivots(np.ones(N), np.ones(N), np.ones(N)), "test")


def test_symmetric_k_g_pass():
    validate_operator(_op())


def test_antisymmetric_k_rejected():
    try:
        validate_operator(_op(k_asym=1e-3))
    except ValueError as e:
        assert "matrix symmetry" in str(e)
    else:
        raise AssertionError("antisymmetric K must be rejected")


def test_symmetric_s_rejected():
    try:
        validate_operator(_op(s_sym=1e-3))
    except ValueError as e:
        assert "matrix symmetry" in str(e)
    else:
        raise AssertionError("symmetric S must be rejected")


def test_nonfinite_m_rejected():
    op = _op()
    op.M[0, 0, 0] = np.nan
    try:
        validate_operator(op)
    except ValueError as e:
        assert "M" in str(e)
    else:
        raise AssertionError("nonfinite M must be rejected")
