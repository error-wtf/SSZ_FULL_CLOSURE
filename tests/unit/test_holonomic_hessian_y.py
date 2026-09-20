import numpy as np
import pandas as pd

from ssz_p5.production.holonomic_hessian_y import (
    TRANSVERSE,
    algebraic_response_matrix,
    chain_residual,
    complete_hessian,
    response_from_hessian,
)


def background(n=17):
    x = np.linspace(1.45, 1.35, n)
    f = 0.55 + 0.03 * (x - 1.4)
    h = 0.7 + 0.02 * (x - 1.4)
    ph = -0.8 - 0.05 * (x - 1.4)
    A = 0.3 + 0.04 * (x - 1.4)
    X = -0.5 * h * ph**2
    return pd.DataFrame(dict(x=x, f=f, h=h, phiprime=ph, A0prime=A, X=X))


def test_4d_completion_is_background_null():
    b = background()
    rng = np.random.default_rng(7)
    q = pd.DataFrame(rng.normal(size=(len(b), 6)), columns=TRANSVERSE)
    H, _ = complete_hessian(b, q)
    residual, normalized = chain_residual(b, H)
    assert np.max(np.abs(residual)) < 2e-10
    assert np.max(normalized) < 2e-12


def test_local_algebraic_matrix_matches_direct_response():
    b = background()
    rng = np.random.default_rng(11)
    q = rng.normal(size=(len(b), 6))
    H, _ = complete_hessian(b, q)
    emitted = response_from_hessian(b, H)
    M = algebraic_response_matrix(b)
    replay = np.einsum("nij,nj->ni", M, q)
    expected = emitted[["v5", "c3", "v1", "v4", "c2"]].to_numpy(float)
    assert np.max(np.abs(replay - expected)) < 2e-10


def test_y_sector_changes_response_without_breaking_holonomy():
    b = background()
    q = np.zeros((len(b), 6))
    q[:, 2] = np.linspace(-0.4, 0.7, len(b))  # XY
    q[:, 4] = 0.2  # FY
    q[:, 5] = -0.1  # YY
    H, _ = complete_hessian(b, q)
    residual, _ = chain_residual(b, H)
    emitted = response_from_hessian(b, H)
    assert np.max(np.abs(residual)) < 2e-10
    assert np.max(np.abs(emitted[["v1", "v4", "c2"]].to_numpy(float))) > 0
