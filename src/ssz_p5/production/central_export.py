"""Export the authoritative central 41-slot representation in common convention.

This is the coefficient import/normalization route explicitly allowed by the
regional execution contract, not an independent action-to-41 derivation.
"""

import json
import subprocess
from datetime import UTC, datetime

import numpy as np
import pandas as pd

from ..config import SLOT_NAMES
from ..provenance.manifest import sha256
from .member import MEMBER_FILE, validate_stream_regions
from .regional_coefficients import central_selected
from .sources import CENTRAL_PRODUCTION_SOURCES

DEPENDENCIES = {
    "c2": "selected f3XX versus the original principal coefficient/Hessian selection",
    "c6": "c2 propagated into c6",
    "e2": "c2 propagated into e2",
    "d3": "partial_phi(v6) versus the replay's total radial derivative/phiprime selector",
    "v9": "stored odd alpha4 and radial background jets",
    "v13": "higher mixed action jets and radial background derivatives",
    "e4": "v13, c4_prime, v6_prime and alpha6_prime",
    "a5": "corrected derivative-of-product convention",
}


def slot_comparison(actual, reference, tolerance=1e-7):
    if len(actual) != len(reference) or not np.allclose(actual.x, reference.x, rtol=0, atol=1e-13):
        raise ValueError("central slot comparison requires identical grids")
    rows = []
    for name in SLOT_NAMES:
        a, b = actual[name].to_numpy(float), reference[name].to_numpy(float)
        if not np.isfinite(a).all() or not np.isfinite(b).all():
            raise ValueError(f"nonfinite central coefficient: {name}")
        error = abs(a - b)
        scaled = error / np.maximum(1, abs(b))
        bad = np.flatnonzero(~np.isfinite(scaled) | (scaled > tolerance))
        i = int(bad[0]) if len(bad) else int(np.argmax(scaled))
        rows.append(
            dict(
                slot=name,
                direct_value=float(a[i]),
                reference_value=float(b[i]),
                abs_error=float(error[i]),
                relative_error=float(error[i] / max(abs(b[i]), 1e-300)),
                scaled_error=float(scaled[i]),
                first_bad_radius=float(actual.x.iloc[i]) if len(bad) else None,
                max_error=float(np.max(error)),
                max_scaled_error=float(np.max(scaled)),
                status="FAIL" if len(bad) else "PASS",
                probable_dependency=DEPENDENCIES.get(
                    name, "authoritative central coefficient and common convention"
                ),
            )
        )
    return pd.DataFrame(rows)


def export_central(root, output):
    output.mkdir(parents=True, exist_ok=True)
    # Normalize on the full guard grid; trim only after all radial derivatives.
    full = central_selected(root)
    full["region"] = "central_exact_SVT"
    full["source_member"] = MEMBER_FILE
    keep = (full.u >= 0.61) & (full.u < 0.71)
    central = full.loc[keep].reset_index(drop=True)
    validate_stream_regions(central)
    reference = pd.read_csv(root / CENTRAL_PRODUCTION_SOURCES["coeff_reference"])
    reference = reference.loc[keep].reset_index(drop=True)
    comparison = slot_comparison(central, reference)
    comparison.to_csv(output / "CENTRAL_SLOT_COMPARISON.csv", index=False)
    path = output / "central_exact_SVT_41of41.csv"
    central.to_csv(path, index=False)
    exact = pd.read_csv(root / CENTRAL_PRODUCTION_SOURCES["exact_regression"])
    if not np.allclose(exact.x, full.x, rtol=0, atol=1e-13):
        raise ValueError("exact central regression grid differs")
    # Only quantities actually present in both representations are comparable.
    exact_errors = {}
    for name in ("v1", "v4", "c2", "v13", "e4"):
        a, b = central[name].to_numpy(), exact.loc[keep, name].to_numpy()
        exact_errors[name] = float(np.max(abs(a - b) / np.maximum(1, abs(b))))
    odd_columns = [c for c in exact if c.startswith("alpha") or c.startswith("odd_") or c == "Codd"]
    odd = exact.loc[keep, ["u", "x", "A0prime"] + odd_columns].copy()
    odd["region"] = "central_exact_SVT"
    odd.to_csv(output / "central_exact_SVT_odd_reference.csv", index=False)
    passed = bool((comparison.status == "PASS").all() and max(exact_errors.values()) < 1e-7)
    code_paths = [
        "src/ssz_p5/production/central_export.py",
        "src/ssz_p5/production/regional_coefficients.py",
        "src/ssz_p5/production/sources.py",
        "src/ssz_p5/jets/jet9d8.py",
        MEMBER_FILE,
    ]
    sources = sorted({v for k, v in CENTRAL_PRODUCTION_SOURCES.items() if k != "role"})
    cert = dict(
        CENTRAL_DIRECT_41="PASS" if passed else "FAIL",
        certificate_scope=(
            "authoritative central coefficient import and common-convention normalization"
        ),
        independent_action_to_41_verified=False,
        finite_L_stability_certified=False,
        exact_ZK_comparable_slots="PASS" if max(exact_errors.values()) < 1e-7 else "FAIL",
        exact_ZK_slot_errors=exact_errors,
        stored_finite_L_matrix_available=False,
        region="central_exact_SVT",
        rows=len(central),
        slots=list(SLOT_NAMES),
        selection="preserve principal slots and encoded v5; selected c3/e3; corrected a5 and v12",
        numerical_policy="JET9D8 window=9 degree=8, derivatives before region trimming",
        tolerance=1e-7,
        generator="tools/export_central_production.py",
        git_commit=subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=root, text=True
        ).strip(),
        timestamp=datetime.now(UTC).isoformat(),
        sources=[dict(path=p, sha256=sha256(root / p)) for p in sources],
        implementation=[dict(path=p, sha256=sha256(root / p)) for p in code_paths],
        outputs=[
            dict(path=str(p.relative_to(root)), sha256=sha256(p))
            for p in (
                path,
                output / "CENTRAL_SLOT_COMPARISON.csv",
                output / "central_exact_SVT_odd_reference.csv",
            )
        ],
        background_recalculation="DIAGNOSTIC_ONLY",
        absolute_full_closure="NOT_CERTIFIED",
    )
    (output / "CENTRAL_DIRECT_41_CERTIFICATE.json").write_text(
        json.dumps(cert, indent=2, allow_nan=False) + "\n"
    )
    return cert
