"""Inner lower-order targets from frozen central jets and the stored partition.

Targets are not injected into production without a realized action-control map.
The principal f2-Hessian map controls different slots and cannot fill that role.
"""

import json

import numpy as np
import pandas as pd

from ..provenance.manifest import sha256
from .sources import SOURCE_REGISTRY
from .lower_order_controls import invert_lower_order_targets

# v5/c3 multiply terms with at most one radial field derivative after constraints;
# e3 multiplies dphi^2 directly. Cphi cancels v5 through a5 and v2=A0prime*v1.
JET_ORDERS = {"v5": 1, "c3": 1, "e3": 0}


def endpoint_derivative_orders(name: str) -> tuple[int, ...]:
    """Return the *contractual* endpoint derivative orders for one target.

    This is the single source of truth for handover endpoint constraints.  A
    solver must not silently add higher derivative conditions: if ``e3`` is
    declared order 0, only its endpoint value is a hard patch condition.
    """
    if name not in JET_ORDERS:
        raise KeyError(f"unknown inner endpoint target: {name}")
    return tuple(range(JET_ORDERS[name] + 1))


def endpoint_constraint_count(names=None, sides: int = 2) -> int:
    """Number of scalar endpoint equations required by the frozen contract."""
    if names is None:
        names = tuple(JET_ORDERS)
    if sides <= 0:
        raise ValueError("sides must be positive")
    return sides * sum(len(endpoint_derivative_orders(name)) for name in names)


def left_endpoint_jets(r, values, endpoint, order):
    r, values = np.asarray(r, float), np.asarray(values, float)
    candidates = np.flatnonzero(r >= endpoint)
    if len(candidates) < 9:
        raise ValueError("nine central-side samples required")
    inds = candidates[np.argsort(abs(r[candidates] - endpoint))[:9]]
    scale = float(np.max(abs(r[inds] - endpoint)))
    coefficients = np.polynomial.polynomial.polyfit((r[inds] - endpoint) / scale, values[inds], 8)
    result = [float(coefficients[0])]
    if order:
        result.append(float(coefficients[1] / scale))
    return result


def partitioned_target(r, S, jets, endpoint):
    extension = np.full(len(r), jets[0], dtype=float)
    if len(jets) > 1:
        extension += jets[1] * (np.asarray(r) - endpoint)
    return np.asarray(S) * extension


def build_inner_targets(root, output):
    output.mkdir(parents=True, exist_ok=True)
    reg = SOURCE_REGISTRY
    paths = [
        reg["inner_same_action_SVT_H"]["background"],
        reg["inner_same_action_SVT_H"]["coeff_reference"],
        reg["central_exact_SVT"]["unreduced_even"],
        reg["central_exact_SVT"]["lower_jets"],
        reg["punctured_H_core"]["coeff_reference"],
    ]
    b, old, central, lower, core = [pd.read_csv(root / p) for p in paths]
    if not np.allclose(central.x, lower.x, rtol=0, atol=1e-13):
        raise ValueError("central lower-order source grid mismatch")
    if not np.allclose(b.S_SVT + b.T_H, 1, rtol=0, atol=1e-12):
        raise ValueError("invalid stored inner partition")
    if b.S_SVT.iloc[0] != 1 or b.S_SVT.iloc[-1] != 0:
        raise ValueError("wrong inner partition endpoints")
    if not np.allclose(old.v2, old.A0prime * old.v1, rtol=1e-12, atol=1e-10):
        raise ValueError("v5 cancellation prerequisite absent")
    if np.max(abs(core[["v5", "c3", "e3"]].to_numpy())) > 1e-12:
        raise ValueError("core lower-order endpoint is not zero")
    central["c3"], central["e3"] = lower.c3_selected, lower.e3_selected
    endpoint = 1 / 0.71
    jets = {
        name: left_endpoint_jets(central.x, central[name], endpoint, order)
        for name, order in JET_ORDERS.items()
    }
    targets = b[["x", "u", "S_SVT", "T_H"]].copy()
    for name, jet in jets.items():
        targets[name + "_target"] = partitioned_target(b.x, b.S_SVT, jet, endpoint)
    target_path = output / "INNER_LOWER_ORDER_TARGETS.csv"
    targets.to_csv(target_path, index=False)

    action_controls, emitted_lower, control_report = invert_lower_order_targets(b, targets)
    control_path = output / "INNER_LOWER_ORDER_ACTION_CONTROLS.csv"
    emitted_path = output / "INNER_LOWER_ORDER_REEMITTED.csv"
    action_controls.to_csv(control_path, index=False)
    emitted_lower.to_csv(emitted_path, index=False)

    report = dict(
        status="ACTION_REALIZED_PASS",
        production_changed=False,
        central_changed=False,
        central_endpoint_x=endpoint,
        central_endpoint_radial_jets=jets,
        required_radial_orders=JET_ORDERS,
        core_endpoint_radial_jets={k: [0.0] * (v + 1) for k, v in JET_ORDERS.items()},
        partition="stored S_SVT,T_H; no new transition function",
        derivative_policy="9 point degree 8 local polynomial, central side of full source grid",
        action_control_variables=control_report["controls"],
        control_map_residual=control_report["max_scaled_replay_error"],
        background_null_check="NOT_CERTIFIED_UNTIL_COMMON_HESSIAN_HOLONOMY",
        action_control_status="LOCAL_RESPONSE_REALIZED_AND_REEMITTED",
        common_action_warning=(
            "f2phiF/f2Xphi/f2phiphi cannot be promoted independently of f2FF/f2XF/f2XX; "
            "run the common Hessian holonomy audit"
        ),
        action_control_report=control_report,
        existing_emitter_interface=(
            "selected_v5/selected_c3/selected_e3 remain an output interface; production values are now supplied by the explicit lower-order action-control inverse"
        ),
        sources=[dict(path=p, sha256=sha256(root / p)) for p in paths],
        generator=dict(
            path="src/ssz_p5/production/inner_targets.py",
            sha256=sha256(root / "src/ssz_p5/production/inner_targets.py"),
        ),
        target_artifact=dict(path=str(target_path.relative_to(root)), sha256=sha256(target_path)),
        action_control_artifact=dict(path=str(control_path.relative_to(root)), sha256=sha256(control_path)),
        reemitted_lower_artifact=dict(path=str(emitted_path.relative_to(root)), sha256=sha256(emitted_path)),
    )
    (output / "INNER_LOWER_ORDER_CONTROL.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    return report
