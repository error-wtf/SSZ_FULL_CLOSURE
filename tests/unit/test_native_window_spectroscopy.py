import numpy as np

from ssz_p5.qnm.native_window import assemble_dirichlet_fem, solve_box_spectrum


def test_constant_scalar_box_recovers_dirichlet_laplacian():
    n = 96
    r = np.linspace(0.0, 1.0, n)
    u = r.copy()
    K = np.ones((n, 1, 1))
    G = np.ones((n, 1, 1))
    S = np.zeros((n, 1, 1))
    M = np.zeros((n, 1, 1))
    sp = solve_box_spectrum(r, u, K, G, S, M, modes=4)
    expected = (np.arange(1, 5) * np.pi) ** 2
    assert sp.negative_omega2_count == 0
    assert np.max(np.abs(sp.omega2 - expected) / expected) < 5e-3


def test_assembled_pencil_is_symmetric_with_antisymmetric_S():
    n = 40
    r = np.linspace(1.0, 2.0, n)
    K = np.repeat(np.eye(2)[None, :, :], n, axis=0)
    G = 2 * K.copy()
    M = 0.1 * K.copy()
    S = np.zeros_like(K)
    S[:, 0, 1] = 0.2
    S[:, 1, 0] = -0.2
    A, B, _, _ = assemble_dirichlet_fem(r, K, G, S, M)
    assert np.max(np.abs(A - A.T)) < 1e-12
    assert np.max(np.abs(B - B.T)) < 1e-12
