"""Reproducible principal-feasibility controls for the electric-hybrid search.

This module deliberately does *not* claim background/on-shell closure.  It takes
an already emitted electric ZK coefficient stream and adds Horndeski primitive
responses generated through the common Appendix-A Maxwell--Horndeski emitter.
The purpose is to test whether the finite-L kinetic obstruction is controllable
without violating the algebraic constraint pivots.

The stored recipe is a research-search checkpoint, not a production member.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from ..coefficients.mh_general_primitives import SLOTS, emit_from_primitives
from ..reducer.kinetic_schur import kinetic_schur

CONTROL_COLUMNS = {
    "a1": "a1_control",
    "F_tensor": "F_tensor_control",
    "H_tensor": "H_tensor_control",
}


def load_search_inputs(root: Path):
    root = Path(root)
    recipe_path = (
        root
        / "data/generated/absolute_attempt_2026-09-19/ELECTRIC_HYBRID_PRINCIPAL_RECIPE.json"
    )
    recipe = json.loads(recipe_path.read_text())
    baseline = pd.read_csv(root / recipe["baseline"]).sort_values("x").reset_index(drop=True)
    basis = pd.read_csv(root / recipe["primitive_basis"]).sort_values("x").reset_index(drop=True)
    if len(baseline) != len(basis) or not np.allclose(baseline.x, basis.x, rtol=0, atol=1e-13):
        raise ValueError("electric-hybrid baseline/control basis grid mismatch")
    return baseline, basis, recipe


def _reference_frame(baseline: pd.DataFrame) -> pd.DataFrame:
    d = pd.DataFrame(
        {
            "u": baseline.u,
            "x": baseline.x,
            "phi": baseline.phi,
            "f": baseline.f,
            "h": baseline.h,
            "phiprime": baseline.phiprime,
            "A0prime": baseline.A0prime,
            "a1": 0.0,
            "c2": 0.0,
            "c4": 0.0,
            "F_tensor": 0.0,
            "G_tensor": 0.0,
            # Keep a nonzero reference H so the emitter's diagnostic scalar oracle
            # remains finite; all 41 coefficient responses are linear and the
            # reference contribution is subtracted below.
            "H_tensor": 1.0,
        }
    )
    return d


def primitive_response(
    baseline: pd.DataFrame,
    basis: pd.DataFrame,
    control: str,
    shape: np.ndarray,
) -> dict[str, np.ndarray]:
    if control not in CONTROL_COLUMNS:
        raise ValueError(f"unsupported control {control!r}")
    shape = np.asarray(shape, float)
    if shape.shape != (len(baseline),):
        raise ValueError("control shape length mismatch")
    ref = _reference_frame(baseline)
    e0 = emit_from_primitives(ref, regularize_photon_root=True)
    d = ref.copy()
    d[control] = d[control].to_numpy(float) + basis[CONTROL_COLUMNS[control]].to_numpy(float) * shape
    e1 = emit_from_primitives(d, regularize_photon_root=True)
    return {slot: e1[slot].to_numpy(float) - e0[slot].to_numpy(float) for slot in SLOTS}


def build_recipe_stream(root: Path) -> pd.DataFrame:
    baseline, basis, recipe = load_search_inputs(root)
    u = baseline.u.to_numpy(float)
    q = (u - float(recipe["polynomial_center_u"])) / float(recipe["polynomial_scale_u"])

    total = {control: np.zeros_like(u) for control in CONTROL_COLUMNS}
    for control, amps in recipe["polynomial_amplitudes"].items():
        total[control] += (
            float(amps["const"])
            + float(amps["linear"]) * q
            + float(amps["quad"]) * q * q
        ) * basis[CONTROL_COLUMNS[control]].to_numpy(float)

    for key in ("ramp_1", "ramp_2", "late_ramp"):
        spec = recipe[key]
        ramp = 0.5 * (1.0 + np.tanh((u - float(spec["center_u"])) / float(spec["width_u"])))
        shape = ramp ** int(spec["power"])
        for control, amplitude in spec["amplitudes"].items():
            total[control] += (
                float(amplitude) * shape * basis[CONTROL_COLUMNS[control]].to_numpy(float)
            )

    ref = _reference_frame(baseline)
    e0 = emit_from_primitives(ref, regularize_photon_root=True)
    changed = ref.copy()
    for control, delta in total.items():
        changed[control] = changed[control].to_numpy(float) + delta
    e1 = emit_from_primitives(changed, regularize_photon_root=True)

    out = baseline.copy()
    for slot in SLOTS:
        out[slot] = baseline[slot].to_numpy(float) + (
            e1[slot].to_numpy(float) - e0[slot].to_numpy(float)
        )
    return out


def kinetic_feasibility_audit(root: Path, *, umin: float = 0.62, umax: float = 0.70):
    d = build_recipe_stream(root)
    u = d.u.to_numpy(float)
    production = (u > umin) & (u < umax)
    inds = np.flatnonzero(production)
    rows = []
    for L in (6, 12, 20, 42, 110, 420, 1000):
        result = kinetic_schur(d, L)
        K = result["K"]
        eig = np.linalg.eigvalsh((K + K.swapaxes(1, 2)) / 2)
        j = int(inds[np.argmin(eig[inds, 0])])
        rows.append(
            {
                "L": L,
                "min_eig_K": float(eig[j, 0]),
                "u": float(u[j]),
                "x": float(d.x.iloc[j]),
                "min_abs_Dh1": float(np.min(np.abs(result["Dh1"][inds]))),
                "min_abs_auxiliary_determinant": float(
                    np.min(np.abs(result["auxiliary_determinant"][inds]))
                ),
                "pass_kinetic": bool(np.all(eig[inds, 0] > 0)),
            }
        )
    return d, rows
