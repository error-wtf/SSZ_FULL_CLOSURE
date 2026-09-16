from dataclasses import replace

import numpy as np
import pytest

from ssz_p5.reducer.canonical import validate_operator
from ssz_p5.stability.finite_l import stability_diagnostics
from ssz_p5.types import ConstraintPivots, ReducedOperator


def operator():
    identity = np.tile(np.eye(3), (9, 1, 1))
    Z = np.zeros_like(identity)
    ones = np.ones(9)
    return ReducedOperator(
        6,
        np.arange(1.0, 10.0),
        identity,
        Z,
        identity,
        Z,
        Z,
        ConstraintPivots(ones, ones, ones),
        "synthetic-test",
    )


@pytest.mark.parametrize("kind", ["nan", "asymmetric", "symmetric_S", "zero_pivot", "nonzero_R"])
def test_reject_bad_operator(kind):
    op = operator()
    if kind == "nan":
        op.K[0, 0, 0] = np.nan
    if kind == "asymmetric":
        op.K[0, 0, 1] = 1
    if kind == "symmetric_S":
        op.S[0, 0, 0] = 1
    if kind == "zero_pivot":
        op.pivots.Dh1[0] = 0
    if kind == "nonzero_R":
        op = replace(op, R=np.ones_like(op.R))
    with pytest.raises(ValueError):
        validate_operator(op)


def test_symmetric_whitening_and_negative_kinetic():
    op = operator()
    assert stability_diagnostics(op)["min_radial"] == 1
    op.K[0, 0, 0] = -1
    with pytest.raises(ValueError):
        stability_diagnostics(op)
