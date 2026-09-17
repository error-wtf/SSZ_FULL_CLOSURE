"""Independent kinetic-action Schur check in the fixed (psi, dphi, V) basis.

Uses neither constraint-map code nor the Euler-operator implementation. The
six columns represent (dot Y, dot Y'). H0 is eliminated by coefficient matching;
H1 and dA1 by a batched 2x2 solve. The symmetric mixed term is integrated by
parts using the same declared radial derivative service.
"""

import numpy as np

from ..config import SLOT_NAMES
from ..jets.jet9d8 import profile_derivative


def kinetic_schur(d, L, *, window=9, degree=8):
    if L < 6 or int(L) != L:
        raise ValueError("L must be an integer >=6")
    x = d.x.to_numpy(float)
    q = {c: d[c].to_numpy(float) for c in SLOT_NAMES}
    n = len(d)

    def dr(z):
        return profile_derivative(x, np.asarray(z), 1, window, degree)

    p = q["a1"] / q["a3"]
    t = L * q["a4"] / q["a3"]
    ch2 = q["a7"] + L * q["a8"] - q["v2"] * q["v3"] / (2 * q["v1"])
    b = q["a2"] - q["v2"] * q["v4"] / (2 * q["v1"]) - q["a3"] * dr(p) - ch2 * p
    pivot = L * (q["a9"] - q["v2"] * q["v6"] / (2 * q["v1"])) - q["a3"] * dr(t) - ch2 * t
    h = np.zeros((n, 6))
    h[:, 0] = -ch2 / pivot
    h[:, 1] = -(q["a5"] + L * q["a6"] - q["v2"] * q["v5"] / (2 * q["v1"])) / pivot
    h[:, 2] = -q["v2"] / pivot
    h[:, 3] = -q["a3"] / pivot
    h[:, 4] = -b / pivot
    H = -t[:, None] * h
    H[:, 0] += 1
    H[:, 4] -= p
    ef = np.zeros((n, 6))
    ef[:, 1] = 1
    ev = np.zeros((n, 6))
    ev[:, 2] = 1
    efp = np.zeros((n, 6))
    efp[:, 4] = 1

    def sym(a, b):
        return 0.5 * (a[:, :, None] * b[:, None, :] + b[:, :, None] * a[:, None, :])

    C = (
        q["e1"][:, None, None] * sym(ef, ef)
        + q["c1"][:, None, None] * sym(ef, H)
        + (L * q["d1"])[:, None, None] * sym(h, h)
    )
    Q = np.empty((n, 2, 2))
    Q[:, 0, 0] = L * q["b1"]
    Q[:, 1, 1] = L * q["v10"]
    Q[:, 0, 1] = Q[:, 1, 0] = L * q["v11"] / 2
    J = np.stack(
        [
            q["b2"][:, None] * efp
            + q["b3"][:, None] * ef
            + q["b4"][:, None] * H
            + (L * q["b5"])[:, None] * h,
            2 * q["v1"][:, None] * ev - (L * q["v6"] / 2)[:, None] * h,
        ],
        axis=1,
    )
    C -= 0.25 * np.swapaxes(J, 1, 2) @ np.linalg.solve(Q, J)
    cross = C[:, :3, 3:]
    K = C[:, :3, :3] - dr(0.5 * (cross + cross.swapaxes(1, 2)))
    return dict(
        K=K,
        quadratic_form=C,
        mixed_second=C[:, 3:, 3:],
        antisymmetric_cross=cross - cross.swapaxes(1, 2),
        Dh1=pivot,
        auxiliary_determinant=np.linalg.det(Q),
    )
