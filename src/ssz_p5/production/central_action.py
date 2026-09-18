"""Reconstruct central action inputs and test the on-shell prerequisites.

References: Zhang--Kase, arXiv:2404.11910v3, Eqs.12,13,18,131.
The original selected member remains immutable. c3/e3 are supplied lower-order
controls; this replay is not an independently varied 41-slot certificate.
"""

import numpy as np
import pandas as pd

from ..jets.jet9d8 import profile_derivative
from ..numerics import module
from .regional_coefficients import HESSIANS
from .sources import SOURCE_REGISTRY


def central_action_inputs(root):
    registry = SOURCE_REGISTRY["central_exact_SVT"]
    background = pd.read_csv(root / registry["unreduced_even"])
    source = pd.read_csv(root / registry["exact_regression"])
    jets = pd.read_csv(root / registry["lower_jets"])
    for other in (source, jets):
        if len(other) != len(background) or not np.allclose(
            background.x, other.x, rtol=0, atol=1e-13
        ):
            raise ValueError("central action sources require identical radial grids")
    d = background[["u", "x", "phi", "f", "h", "phiprime", "A0prime", "X"]].copy()
    for name in ("f2", "f2X", "f3", "f4"):
        d[name] = source[name]
    d["f2F"], d["f2Y"], d["tf4"] = 1.0, 0.0, 0.0
    d["f4X"], d["f3X"], d["tf3"] = source.N4, source.f3X_integrated, source.tilde_f3
    for name in (
        "f3phi",
        "f3phiX",
        "f4phi",
        "f4phiX",
        "f4phiXX",
        "f2phi",
        "f2phiX",
        "f2phiF",
        "f2phiY",
    ):
        d[name] = jets[name]
    d["f3XX"] = jets.f3XX_selected
    d["f4XX"] = jets.f4XX_recovered
    d["f4XXX"] = jets.f4XXX_selected
    for src, dst in HESSIANS.items():
        d[dst] = source[src]
    # Explicitly labelled controls, not independently regenerated coefficients.
    d["selected_c3"], d["selected_e3"] = jets.c3_selected, jets.e3_selected
    if not np.isfinite(d.to_numpy(float)).all():
        raise ValueError("nonfinite central action input")
    return d


def replay_central_action(d, *, window=9, degree=8):
    r, f, h, ph, A = (d[c].to_numpy() for c in ("x", "f", "h", "phiprime", "A0prime"))
    v5 = (
        np.sqrt(h / f)
        * A
        * (
            2 * h * ph**2 * (h * d.f4phiX - r**2 * d.f2phiY)
            + 4 * r * h * ph * d.f3phi
            + 8 * (1 - h) * d.f4phi
            + r**2 * d.f2phiF
        )
    )
    return module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py").emit(
        d,
        window=window,
        degree=degree,
        selected_v5=v5,
        selected_c3=d.selected_c3,
        selected_e3=d.selected_e3,
        v6_phi_selector="action",
    )


def background_residuals(d, *, window=9, degree=8):
    """Evaluate the electric central action, without importing stability targets."""
    r, f, h, ph, A = (d[c].to_numpy() for c in ("x", "f", "h", "phiprime", "A0prime"))

    def dr(y, order=1):
        return profile_derivative(r, np.asarray(y), order, window, degree)

    f2, f2x, f2f, f2y, f3, f3x, f4, f4x, f4xx, tf4 = (
        d[c].to_numpy()
        for c in ("f2", "f2X", "f2F", "f2Y", "f3", "f3X", "f4", "f4X", "f4XX", "tf4")
    )
    fp, hp = dr(f), dr(h)
    e00 = r * f * hp - (
        f * (1 - h)
        + r**2 * (f * f2 - h * A**2 * (f2f - 2 * h * ph**2 * f2y))
        - 2 * r * h**2 * ph * A**2 * f3
        + h * A**2 * (4 * (h - 1) * f4 - h**2 * ph**2 * (f4x + 2 * tf4))
    )
    e11 = r * h * fp - (
        f * (1 - h)
        + r**2 * (f * f2 + f * h * ph**2 * f2x - h * A**2 * (f2f - 4 * h * ph**2 * f2y))
        - 2 * r * h**2 * ph * A**2 * (3 * f3 - h * ph**2 * f3x)
        + h
        * A**2
        * (
            4 * (3 * h - 1) * f4
            - h * (9 * h - 4) * ph**2 * f4x
            + h**3 * ph**4 * f4xx
            - 10 * h**2 * ph**2 * tf4
        )
    )
    ja = (
        np.sqrt(h / f)
        * A
        * (
            r**2 * (f2f - 2 * h * ph**2 * f2y)
            + 4 * r * h * ph * f3
            + 8 * (1 - h) * f4
            + 2 * h**2 * ph**2 * (f4x + 2 * tf4)
        )
    )
    a4 = np.sqrt(f * h) / 2
    v6 = 2 * h**1.5 * A / (r * np.sqrt(f)) * (r * ph * f3 - 4 * f4 + h * ph**2 * (f4x + 2 * tf4))
    v10 = (
        -np.sqrt(f * h)
        / (2 * r)
        * (
            r * (f2f - 2 * h * ph**2 * f2y)
            + 2 * h * ph * f3
            + h * fp / f * (r * ph * f3 - 4 * f4 + h * ph**2 * (f4x + 2 * tf4))
        )
    )
    alpha7 = (1 - 4 * h * A**2 * f4 / f) / (4 * r**2 * np.sqrt(f * h))
    # Undivided Eq.131 remains meaningful also when A'=0. Never divide by A'.
    lhs = f**2 * h * r**2 * A * dr(v6)
    rhs = a4 * (
        2 * f**2 * (r * hp + 2 * h)
        + h * r**2 * fp**2
        - f * r * (r * fp * hp + 2 * h * (r * dr(f, 2) + fp))
    )
    rhs -= f * h * r**2 * (A * (4 * v10 * A + v6 * fp) + f * v6 * dr(A) + 8 * alpha7 * f**2)
    return dict(
        E00=e00,
        E11=e11,
        JA=ja,
        JA_prime=dr(ja),
        eq131_lhs=lhs,
        eq131_rhs=rhs,
        eq131_residual=lhs - rhs,
    )


def audit_central_action(root):
    from ..config import SLOT_NAMES
    from ..policy import numerical_policy
    from ..provenance.manifest import sha256
    from ..reducer.kinetic_schur import kinetic_schur
    from ..reducer.zk_kinetic import zk_kinetic
    from .member import MEMBER_FILE
    from .regional_coefficients import central_selected

    inputs = central_action_inputs(root).sort_values("x").reset_index(drop=True)
    reference = central_selected(root).sort_values("x").reset_index(drop=True)
    emitted = replay_central_action(inputs)
    production = (inputs.u >= 0.61) & (inputs.u < 0.71)
    interior = (inputs.u > 0.62) & (inputs.u < 0.70)
    i = int(np.argmin(abs(inputs.u.to_numpy() - 0.69)))
    residuals = background_residuals(inputs)
    slot_errors = {
        name: float(
            np.max(
                abs(emitted.loc[production, name] - reference.loc[production, name])
                / np.maximum(1, abs(reference.loc[production, name]))
            )
        )
        for name in SLOT_NAMES
    }
    K = zk_kinetic(reference, 6)
    schur = kinetic_schur(reference, 6)["K"]
    convergence = []
    for stride, window, degree in ((1, 9, 8), (2, 9, 8), (4, 9, 8), (1, 7, 6), (1, 11, 8)):
        data = inputs.iloc[::stride].reset_index(drop=True)
        j = int(np.argmin(abs(data.u.to_numpy() - 0.69)))
        values = background_residuals(data, window=window, degree=degree)
        convergence.append(
            dict(
                stride=stride,
                window=window,
                degree=degree,
                x=float(data.x.iloc[j]),
                u=float(data.u.iloc[j]),
                **{name: float(v[j]) for name, v in values.items()},
            )
        )
    policy = numerical_policy()
    validation = {
        "metric_background_equations": all(
            np.max(abs(residuals[name][production])) < policy["background_residual_abs"]
            for name in ("E00", "E11")
        ),
        "published_kinetic_matches_schur": bool(
            np.max(abs(K[interior] - schur[interior]) / np.maximum(1, abs(schur[interior]))) < 2e-7
        ),
        "published_kinetic_positive_L6": bool(np.all(np.linalg.eigvalsh(K[production]) > 0)),
    }
    paths = [MEMBER_FILE] + [
        SOURCE_REGISTRY["central_exact_SVT"][name]
        for name in ("unreduced_even", "exact_regression", "lower_jets", "coeff_reference")
    ]
    code = [
        "src/ssz_p5/production/central_action.py",
        "src/ssz_p5/reducer/zk_kinetic.py",
        "src/ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py",
        "src/ssz_p5/jets/jet9d8.py",
        "src/ssz_p5/production/regional_coefficients.py",
        "src/ssz_p5/reducer/kinetic_schur.py",
        "NUMERICAL_POLICY.json",
    ]
    return dict(
        status="DIAGNOSTIC_ONLY",
        full_selected_background_verified=False,
        scope="central_exact_SVT action consistency and independent finite-L kinetic replay",
        formula_source="https://arxiv.org/html/2404.11910v3",
        equations=[12, 13, 18, 65, 66, 131],
        validation=validation,
        sources=[dict(path=p, sha256=sha256(root / p)) for p in paths],
        implementation=[dict(path=p, sha256=sha256(root / p)) for p in code],
        production_rows=int(production.sum()),
        background_max_abs={
            name: float(np.max(abs(v[production]))) for name, v in residuals.items()
        },
        metric_tolerance=policy["background_residual_abs"],
        published_schur_max_scaled=float(
            np.max(abs(K[interior] - schur[interior]) / np.maximum(1, abs(schur[interior])))
        ),
        witness=dict(
            x=float(inputs.x.iloc[i]),
            u=float(inputs.u.iloc[i]),
            A0prime=float(inputs.A0prime.iloc[i]),
            background={name: float(v[i]) for name, v in residuals.items()},
            K_published=K[i].tolist(),
            eigenvalues=np.linalg.eigvalsh(K[i]).tolist(),
        ),
        derivative_convergence=convergence,
        central_action_replay_slot_errors=slot_errors,
        replay_scope=(
            "v5 regenerated; c3/e3 retain encoded controls; no direct 41 certificate issued"
        ),
        global_direct_41="NOT_CERTIFIED",
        absolute_closure="NOT_CERTIFIED",
    )
