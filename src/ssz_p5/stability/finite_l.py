import numpy as np

from ..reducer.canonical import validate_operator


def stability_diagnostics(op):
    validate_operator(op)
    K = (op.K + op.K.swapaxes(1, 2)) / 2
    w, V = np.linalg.eigh(K)
    if np.min(w) <= 0:
        raise ValueError("nonpositive kinetic eigenvalue")
    inv = (V * (1 / np.sqrt(w))[:, None, :]) @ V.swapaxes(1, 2)
    radial = np.linalg.eigvalsh(inv @ ((op.G + op.G.swapaxes(1, 2)) / 2) @ inv)
    i = int(np.unravel_index(np.argmin(w), w.shape)[0])
    j = int(np.unravel_index(np.argmin(radial), radial.shape)[0])
    return {
        "L": op.L,
        "min_eig_K": float(w.min()),
        "radius_min_K": float(op.r[i]),
        "min_radial": float(radial.min()),
        "radius_min_radial": float(op.r[j]),
        "pass": bool(radial.min() > 0),
    }
