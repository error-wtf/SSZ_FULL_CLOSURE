"""Normalize the selected inner representation and test both frozen interfaces."""

import json
from datetime import UTC, datetime

import numpy as np
import pandas as pd

from ..config import SLOT_NAMES
from ..jets.jet9d8 import profile_derivative
from ..provenance.manifest import sha256
from .central_export import slot_comparison
from .member import validate_stream_regions
from .regional_coefficients import select_lower
from .sources import SOURCE_REGISTRY
from .inner_targets import build_inner_targets
from .inner_principal import action_realize_principal
from .holonomic_hessian import audit_existing_split_controls

CENTRAL_EXPORT = "data/generated/central/central_exact_SVT_41of41.csv"


def endpoint_values(frame, x, order=0, window=9, degree=8):
    """Local JET-policy interpolation of existing radial jets at an interface.

    Derivatives are produced by the existing profile service before selecting
    the endpoint stencil. No endpoint coefficient is changed by this operation.
    """
    r = frame.x.to_numpy(float)
    values = frame[list(SLOT_NAMES)].to_numpy(float)
    if order:
        values = profile_derivative(r, values, order, window, degree)
    inds = np.sort(np.argsort(abs(r - x))[:window])
    scale = np.max(abs(r[inds] - x))
    return np.polynomial.polynomial.polyfit((r[inds] - x) / scale, values[inds], degree)[0]


def compare_interface(inner, neighbor, u):
    rows = []
    for order in (0, 1, 2):
        a = endpoint_values(inner, 1 / u, order)
        b = endpoint_values(neighbor, 1 / u, order)
        for i, name in enumerate(SLOT_NAMES):
            error = float(abs(a[i] - b[i]))
            rows.append(
                dict(
                    u=u,
                    x=1 / u,
                    slot=name,
                    radial_order=order,
                    inner_value=float(a[i]),
                    neighbor_value=float(b[i]),
                    abs_error=error,
                    relative_error=error / max(abs(float(b[i])), 1e-300),
                    scaled_error=error / max(1, abs(float(b[i]))),
                    finite=bool(np.isfinite(a[i]) and np.isfinite(b[i])),
                )
            )
    return pd.DataFrame(rows)


def export_inner(root, output):
    output.mkdir(parents=True, exist_ok=True)
    registry = SOURCE_REGISTRY["inner_same_action_SVT_H"]
    background = pd.read_csv(root / registry["background"])
    reference = pd.read_csv(root / registry["coeff_reference"])
    control_dir = root / "data/generated/inner_controls"
    control_report = build_inner_targets(root, control_dir)
    lower = pd.read_csv(control_dir / "INNER_LOWER_ORDER_REEMITTED.csv")
    if len(lower) != len(reference) or not np.allclose(lower.x, reference.x, rtol=0, atol=1e-12):
        raise ValueError("inner lower-order action-control grid mismatch")
    reference = reference.copy()
    for name in ("v5", "c3", "e3"):
        reference[name] = lower[name].to_numpy(float)
    for name in ("x", "u", "phi", "f", "h", "A0prime"):
        if not np.allclose(reference[name], background[name], rtol=0, atol=1e-12):
            raise ValueError(f"inner authoritative background mismatch: {name}")

    central = pd.read_csv(root / CENTRAL_EXPORT)
    core_path = SOURCE_REGISTRY["punctured_H_core"]["coeff_reference"]
    core = select_lower(pd.read_csv(root / core_path))
    principal, principal_controls, principal_targets, principal_report = action_realize_principal(
        background, reference, lower, central, core
    )
    principal_controls.to_csv(output / "INNER_PRINCIPAL_ACTION_CONTROLS.csv", index=False)
    principal_targets.to_csv(output / "INNER_PRINCIPAL_TARGETS.csv", index=False)

    # Historical 3D split-control diagnostic.  It is retained for comparison,
    # but is superseded as a certification gate by the 4D (phi,X,F,Y) audit.
    # Neither split inverse is a common-action certificate.
    lower_controls = pd.read_csv(control_dir / "INNER_LOWER_ORDER_ACTION_CONTROLS.csv")
    holonomic_frame, holonomic_report = audit_existing_split_controls(
        background, lower_controls, principal_controls
    )
    holonomic_frame.to_csv(output / "INNER_F2_HOLONOMIC_HESSIAN_AUDIT.csv", index=False)
    (output / "INNER_F2_HOLONOMIC_HESSIAN_AUDIT.json").write_text(
        json.dumps(holonomic_report, indent=2, allow_nan=False) + "\n"
    )
    reference = principal
    inner = select_lower(reference)
    # Stored residual columns are preserved, not recomputed from a sector subset.
    for name in ("JA", "metric_residual_14", "metric_residual_15"):
        inner[name] = background[name]
    inner["region"] = "inner_same_action_SVT_H"
    # Recompute the dependent holonomic/constrained slots after the action-level
    # lower-order realization. select_lower applies a5 and v12; v7 is the exact
    # auxiliary identity and must be refreshed on the same stream.
    inner["v7"] = inner.v2**2 / (4 * inner.v1)
    inner["selected_member"] = "action-realized Cinf v5/c3/e3 handover; current a5,v7,v12"
    mask = (inner.u >= 0.71) & (inner.u < 0.715)
    production = inner.loc[mask].reset_index(drop=True)
    validate_stream_regions(production)
    endpoints = pd.concat(
        [compare_interface(inner, central, 0.71), compare_interface(inner, core, 0.715)],
        ignore_index=True,
    )
    endpoints.to_csv(output / "INNER_ENDPOINT_COMPARISON.csv", index=False)
    comparisons = slot_comparison(inner, reference)
    allowed_changes = {"a5", "v5", "c3", "e3", "v7", "v12"}
    comparisons["change_policy"] = comparisons.slot.map(
        lambda name: (
            "action-realized lower-order/holonomic completion" if name in allowed_changes else "preserve"
        )
    )
    comparisons.to_csv(output / "INNER_SLOT_COMPARISON.csv", index=False)
    production.to_csv(output / "inner_same_action_SVT_H_41of41.csv", index=False)
    convergence = []
    for window, degree in ((9, 8), (7, 6), (11, 8)):
        normal = select_lower(reference, window=window, degree=degree)
        for u, other in ((0.71, central), (0.715, core)):
            a = endpoint_values(normal, 1 / u, window=window, degree=degree)
            b = endpoint_values(other, 1 / u, window=window, degree=degree)
            for name in ("a5", "v5", "c3", "e3"):
                i = SLOT_NAMES.index(name)
                convergence.append(
                    dict(
                        window=window,
                        degree=degree,
                        u=u,
                        slot=name,
                        inner_value=float(a[i]),
                        neighbor_value=float(b[i]),
                    )
                )
    pd.DataFrame(convergence).to_csv(output / "INNER_ENDPOINT_CONVERGENCE.csv", index=False)
    finite = bool(np.isfinite(production[list(SLOT_NAMES)].to_numpy()).all())
    # Value discontinuities must pass before any reduced operator can be certified.
    values = endpoints[endpoints.radial_order == 0]
    continuous = bool(values.finite.all() and (values.scaled_error < 1e-7).all())
    unchanged = comparisons[~comparisons.slot.isin(sorted(allowed_changes))]
    normalization = bool(finite and (unchanged.status == "PASS").all())
    paths = [
        registry["background"], registry["coeff_reference"], CENTRAL_EXPORT, core_path,
        "data/generated/inner_controls/INNER_LOWER_ORDER_TARGETS.csv",
        "data/generated/inner_controls/INNER_LOWER_ORDER_ACTION_CONTROLS.csv",
        "data/generated/inner_controls/INNER_LOWER_ORDER_REEMITTED.csv",
    ]
    code = [
        "src/ssz_p5/production/inner_export.py",
        "src/ssz_p5/production/inner_targets.py",
        "src/ssz_p5/production/lower_order_controls.py",
        "src/ssz_p5/production/inner_principal.py",
        "src/ssz_p5/production/principal_controls.py",
        "src/ssz_p5/production/holonomic_hessian.py",
        "src/ssz_p5/production/holonomic_hessian_y.py",
        "src/ssz_p5/production/regional_coefficients.py",
        "src/ssz_p5/jets/jet9d8.py",
        "tools/export_inner_production.py",
    ]
    cert = dict(
        INNER_DIRECT_41="PASS" if normalization and continuous else "FAIL",
        region="inner_same_action_SVT_H",
        domain="0.71 <= u < 0.715",
        rows=len(production),
        slot_count=41,
        finite=finite,
        DIRECT_ACTION_REPLAY_COMPLETE=False,
        AUTHORITATIVE_SELECTED_REPRESENTATION=True,
        LOWER_ORDER_ACTION_CONTROL=(
            "LOCAL_RESPONSE_PASS_NOT_COMMON_ACTION"
            if control_report["status"] == "ACTION_REALIZED_PASS" else "FAIL"
        ),
        lower_order_control_residual=float(control_report["control_map_residual"]),
        PRINCIPAL_ACTION_CONTROL=(
            "LOCAL_RESPONSE_PASS_NOT_COMMON_ACTION"
            if principal_report["status"] == "PASS" else "FAIL"
        ),
        principal_control_residual=float(principal_report["max_scaled_target_error"]),
        HOLONOMIC_F2_HESSIAN="3D_AUDIT_SUPERSEDED_BY_4D_Y_REACHABILITY",
        COMMON_ACTION_GATE="PENDING_FULL_ACTION_F3_F4",
        holonomic_f2_hessian_max_normalized=max(holonomic_report["max_normalized_chain_residual"].values()),
        lower_order_background_null=control_report["background_null_check"],
        SLOT_NORMALIZATION="PASS" if normalization else "FAIL",
        interface_values="PASS" if continuous else "FAIL",
        radial_jets="REPORTED; certification deferred until value continuity passes",
        reducer_compatibility="SCHEMA_ONLY; common single-action operator not certified",
        background_max_abs={
            name: float(abs(background[name]).max())
            for name in ("JA", "metric_residual_14", "metric_residual_15")
        },
        background_scope="stored actual resolved residuals; no background re-solve",
        convention_a5="a2_prime-a1_second-(A0prime*v4/2)_prime+A0prime*v5/2",
        convention_v12="-v6/(2h)",
        derivative_policy="JET9D8 before domain trimming",
        endpoint_method="local polynomial interpolation with the JET9D8 policy; no spline",
        generator="tools/export_inner_production.py",
        timestamp=datetime.now(UTC).isoformat(),
        sources=[dict(path=p, sha256=sha256(root / p)) for p in paths],
        implementation=[dict(path=p, sha256=sha256(root / p)) for p in code],
        outputs=[
            dict(path=str(p.relative_to(root)), sha256=sha256(p))
            for p in sorted(output.glob("*.csv"))
        ],
        central_changed=False,
        absolute_full_closure="NOT_CERTIFIED",
    )
    (output / "INNER_DIRECT_41_CERTIFICATE.json").write_text(
        json.dumps(cert, indent=2, allow_nan=False) + "\n"
    )
    ledger = root / "ABSOLUTE_CLOSURE_LEDGER.json"
    data = json.loads(ledger.read_text()) if ledger.exists() else {}
    data.update(
        baseline="0e13fb0", central_direct_41="PASS", inner_direct_41=cert["INNER_DIRECT_41"]
    )
    for key in ("core_direct_41", "analytic_center", "global_direct_41", "krgsm", "qnm"):
        data.setdefault(key, "PENDING")
    data["inner_evidence"] = str((output / "INNER_DIRECT_41_CERTIFICATE.json").relative_to(root))
    ledger.write_text(json.dumps(data, indent=2) + "\n")
    return cert
