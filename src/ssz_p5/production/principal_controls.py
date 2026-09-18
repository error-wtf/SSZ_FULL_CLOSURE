"""Raw ZK f2-Hessian inversion for a pure SVT contribution.

Targets are v1,v4,c2 in the emitter's actual normalization, not v5,c3,e3.
This does not certify full hybrid assembly, holonomic completion or stability.
"""

import numpy as np

from ..numerics import module

CONTROLS = ("f2FF", "f2XF", "f2XX")
TARGETS = ("v1", "v4", "c2")


def raw_response(d):
    r, f, h, ph, A = (d[c].to_numpy(float) for c in ("x", "f", "h", "phiprime", "A0prime"))
    if not all(np.isfinite(v).all() for v in (r, f, h, ph, A)) or np.any(f <= 0) or np.any(h <= 0):
        raise ValueError("invalid control-map background")
    J = np.zeros((len(d), 3, 3))
    J[:, 0, 0] = r**2 * h**1.5 * A**2 / (2 * f**1.5)
    J[:, 1, 1] = -(h**1.5) * A * r**2 * ph / np.sqrt(f)
    J[:, 2, 1] = h**1.5 * A**2 * r**2 * ph / (2 * np.sqrt(f))
    J[:, 2, 2] = -0.5 * r**2 * np.sqrt(f * h) * h * ph**3
    return J


def restore_svt_principal(d, target, *, selected_lower, tolerance=1e-9):
    """Invert the triangular raw response and re-emit all affected coefficients.

    selected_lower explicitly supplies the independently chosen lower-order
    representation. At zero electric field the two vector channels must already
    equal their targets; no division or arbitrary endpoint filling is allowed.
    """
    if set(selected_lower) != {"v5", "c3", "e3"}:
        raise ValueError("explicit selected v5,c3,e3 required")
    target = np.asarray(target, float)
    if target.shape != (len(d), 3) or not np.isfinite(target).all():
        raise ValueError("invalid v1,v4,c2 targets")
    emitter = module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py")
    options = {"selected_" + k: v for k, v in selected_lower.items()}
    baseline = emitter.emit(d, **options)
    rhs = target - baseline[list(TARGETS)].to_numpy(float)
    J = raw_response(d)
    delta = np.zeros_like(rhs)
    for channel in range(3):
        residual = rhs[:, channel] - np.einsum("ni,ni->n", J[:, channel, :], delta)
        factor = J[:, channel, channel]
        nonzero = factor != 0
        impossible = (~nonzero) & (
            abs(residual) > tolerance * np.maximum(1, abs(target[:, channel]))
        )
        if impossible.any():
            row = int(np.flatnonzero(impossible)[0])
            raise ValueError(f"unreachable {TARGETS[channel]} at row {row}; zero raw response")
        delta[nonzero, channel] = residual[nonzero] / factor[nonzero]
    if not np.isfinite(delta).all():
        raise ValueError("nonfinite raw action controls")
    action = d.copy()
    for i, name in enumerate(CONTROLS):
        action[name] = (d[name].to_numpy(float) if name in d else 0.0) + delta[:, i]
    emitted = emitter.emit(action, **options)
    error = abs(emitted[list(TARGETS)].to_numpy(float) - target) / np.maximum(1, abs(target))
    if not np.isfinite(error).all() or np.max(error) > tolerance:
        raise ValueError("raw Hessian re-emission does not reproduce targets")
    return (
        action,
        emitted,
        {
            "principal_reemission": "PASS",
            "max_scaled_error": float(np.max(error)),
            "controls": list(CONTROLS),
            "targets": list(TARGETS),
            "zero_electric_rows": int(np.sum(d.A0prime.to_numpy() == 0)),
            "lower_order_action_realization": "NOT_CERTIFIED",
            "full_hybrid_assembly": "NOT_CERTIFIED",
        },
    )
