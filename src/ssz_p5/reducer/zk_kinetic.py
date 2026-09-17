"""Independent Zhang--Kase Eqs. (65)--(66), arXiv:2404.11910v3.

Only the Einstein + genuine-SVT branch (a1=b2=c1=0) is supported.
No background identity is substituted to force a positive kinetic matrix.
"""

import numpy as np

from ..jets.jet9d8 import profile_derivative


def zk_kinetic(d, L, *, window=9, degree=8):
    if not np.isfinite(L) or L < 6 or int(L) != L:
        raise ValueError("L must be an integer >=6")
    for name in ("a1", "b2", "c1"):
        if not np.isfinite(d[name]).all() or np.max(abs(d[name])) > 1e-12:
            raise ValueError("ZK kinetic formula requires the Einstein + SVT branch")
    r, f, h, a, A, v1, v6, v10 = (
        d[name].to_numpy(float) for name in ("x", "f", "h", "a4", "A0prime", "v1", "v6", "v10")
    )

    def dr(y):
        return profile_derivative(r, np.asarray(y), 1, window, degree)

    k1 = -2 * r**2 * a / f / (1 - f * v6**2 / (8 * a * v10))
    k2 = 4 * r**2 * a / f / (dr(f) / f + L / (r * h) - 2 / r + A * v6 / (2 * a))
    B = dr(a) / a + 1 / r + L / (2 * r * h)
    Q = f * d.b3.to_numpy() / (r * a)
    C = d.a5.to_numpy() + L * d.a6.to_numpy() - A * d.v5.to_numpy() / 2
    K = np.zeros((len(d), 3, 3))
    K[:, 0, 0] = (k1 + B * k2 - dr(k2) / 2) / L
    K[:, 0, 1] = K[:, 1, 0] = (Q * (k1 + B * k2 / 2) - C * k2 / (r * a) - dr(Q * k2) / 2) / (2 * L)
    K[:, 1, 1] = d.e1 + (
        f * d.b3 / (r**2 * a**2) * (f * d.b3 * k1 - 2 * C * k2) - dr(Q**2 * k2) / 2
    ) / (4 * L)
    K[:, 0, 2] = K[:, 2, 0] = v1 / (4 * L * r * a) * (f * v6 / v10 * k1 - 2 * A * k2)
    K[:, 1, 2] = K[:, 2, 1] = Q * K[:, 0, 2] / 2
    K[:, 2, 2] = f * v1**2 / (2 * L * r**2 * a * v10) * k1
    if not np.isfinite(K).all():
        raise ValueError("nonfinite ZK kinetic matrix")
    return K
