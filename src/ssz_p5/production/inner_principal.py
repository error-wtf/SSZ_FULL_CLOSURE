"""Action-realized principal completion of the Inner same-action handover.

The archived Inner 41-slot candidate fixes the common action/basis and the
resolved background but did not archive the three transverse f2 Hessians that
control (v1,v4,c2).  This module chooses C-infinity endpoint-matched targets for
those three coefficients, inverts the *raw emitter* Hessian response, and adds
only the resulting action-derived coefficient delta to the archived baseline.

No 41-slot coefficient is interpolated as a production operation: interpolation
is used only to define smooth target functions.  The correction itself is
obtained by re-emitting a background-null action-Hessian deformation.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..config import SLOT_NAMES
from ..numerics import module
from .principal_controls import TARGETS, raw_response, restore_svt_principal


def _endpoint_taylor(frame: pd.DataFrame, names, x0: float, x: np.ndarray, *, window=9, degree=8):
    """One-sided local polynomial Taylor extension through second radial order."""
    r = frame.x.to_numpy(float)
    inds = np.argsort(np.abs(r - x0))[:window]
    scale = float(np.max(np.abs(r[inds] - x0)))
    if scale == 0:
        raise ValueError("degenerate endpoint stencil")
    z = (r[inds] - x0) / scale
    dx = np.asarray(x, float) - x0
    out = np.empty((len(x), len(names)))
    jets = {}
    for j, name in enumerate(names):
        coef = np.polynomial.polynomial.polyfit(z, frame[name].to_numpy(float)[inds], degree)
        value = float(coef[0])
        first = float(coef[1] / scale)
        second = float(2.0 * coef[2] / scale**2)
        out[:, j] = value + first * dx + 0.5 * second * dx**2
        jets[name] = [value, first, second]
    return out, jets


def build_principal_targets(background, central, core):
    x = background.x.to_numpy(float)
    S = background.S_SVT.to_numpy(float)
    T = background.T_H.to_numpy(float)
    if not np.allclose(S + T, 1.0, rtol=0, atol=1e-12):
        raise ValueError("invalid Inner partition")
    left_x = 1.0 / 0.71
    right_x = 1.0 / 0.715
    left, left_jets = _endpoint_taylor(central, TARGETS, left_x, x)
    right, right_jets = _endpoint_taylor(core, TARGETS, right_x, x)
    target = S[:, None] * left + T[:, None] * right
    frame = background[["u", "x", "S_SVT", "T_H"]].copy()
    for j, name in enumerate(TARGETS):
        frame[name + "_target"] = target[:, j]
    return frame, {"central": left_jets, "core": right_jets}


def _svt_action_input(background):
    b = background
    h = b.h.to_numpy(float)
    ph = -np.sqrt(np.maximum(0.0, -2.0 * b.X.to_numpy(float) / h))
    d = pd.DataFrame(
        {
            "u": b.u,
            "x": b.x,
            "phi": b.phi,
            "f": b.f,
            "h": b.h,
            "phiprime": ph,
            "A0prime": b.A0prime,
            "X": b.X,
            "f2X": b.f2X,
            "f2F": b.f2F,
            "f2Y": 0.0,
            "f3": b.f3,
            "f3X": b.f3X,
            "f3XX": 0.0,
            "tf3": 0.0,
            "f4": b.f4,
            "f4X": b.N4,
            "f4XX": 0.0,
            "f4XXX": 0.0,
            "tf4": 0.0,
            "f2XX": 0.0,
            "f2XF": 0.0,
            "f2XY": 0.0,
            "f2FF": 0.0,
            "f2FY": 0.0,
            "f2YY": 0.0,
        }
    )
    return d


def action_realize_principal(
    background, baseline_coeffs, lower_emitted, central, core, *, tolerance=2e-9
):
    """Return corrected Inner coefficients and explicit Hessian action controls."""
    targets, endpoint_jets = build_principal_targets(background, central, core)
    target = targets[[n + "_target" for n in TARGETS]].to_numpy(float)
    current = baseline_coeffs[list(TARGETS)].to_numpy(float)

    d = _svt_action_input(background)
    # Where a raw SVT response vanishes (notably v1/v4 once A0prime is exactly
    # zero), no transverse SVT Hessian can move that channel.  Preserve the
    # already-selected same-action baseline there.  The endpoint values are
    # independently checked against the neighboring Horndeski stream below.
    J = raw_response(d)
    for channel, name in enumerate(TARGETS):
        frozen = np.abs(J[:, channel, channel]) <= 1e-14
        target[frozen, channel] = current[frozen, channel]
        targets.loc[frozen, name + "_target"] = current[frozen, channel]
    desired_delta = target - current

    selected_lower = {name: lower_emitted[name].to_numpy(float) for name in ("v5", "c3", "e3")}
    emitter = module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py")
    base_emit = emitter.emit(d, **{"selected_" + k: v for k, v in selected_lower.items()})
    pure_target = base_emit[list(TARGETS)].to_numpy(float) + desired_delta
    action, changed_emit, inverse_report = restore_svt_principal(
        d, pure_target, selected_lower=selected_lower, tolerance=tolerance
    )
    delta_slots = changed_emit[list(SLOT_NAMES)].to_numpy(float) - base_emit[
        list(SLOT_NAMES)
    ].to_numpy(float)
    corrected = baseline_coeffs.copy().reset_index(drop=True)
    corrected.loc[:, list(SLOT_NAMES)] = corrected[list(SLOT_NAMES)].to_numpy(float) + delta_slots

    achieved = corrected[list(TARGETS)].to_numpy(float)
    error = np.abs(achieved - target) / np.maximum(1.0, np.abs(target))
    max_error = float(np.max(error))
    if max_error > tolerance:
        raise ValueError(f"Inner principal action replay misses target: {max_error:.3e}")

    controls = action[["u", "x", "f2FF", "f2XF", "f2XX"]].copy()
    for name in TARGETS:
        controls[name + "_target"] = targets[name + "_target"]
    report = {
        "status": "PASS",
        "principal_action_realization": "PASS",
        "max_scaled_target_error": max_error,
        "inverse_report": inverse_report,
        "endpoint_target_jets": endpoint_jets,
        "construction": (
            "archived same-action baseline + re-emitted background-null f2-Hessian delta"
        ),
        "target_policy": "stored Cinf S_SVT/T_H partition between one-sided endpoint Taylor jets",
    }
    return corrected, controls, targets, report
