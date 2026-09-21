#!/usr/bin/env python3
"""Project the raw Horndeski a1 control onto the local background-EOM tangent.

This is the exact anti-shortcut test requested by the strong-field transition
contract.  The primitive a1 sensitivity found earlier is off-shell until it is
shown to survive projection onto a genuine action-jet nullspace.

We use the quartic G5=0 local action coordinates

    p = (Delta G4_XX, Delta G3_X, Delta G2, Delta G2_X),

construct the full local metric-background Jacobian B=d(E00,E11,E22)/dp, and
project the action-derived primitive gradient d a1/dp onto ker(B).  No selected
41-slot coefficient is treated as an independent control.

A second diagnostic constructs, row by row around the inner light ring, the
unit background-null direction that maximizes positive Delta a1.  Its induced
(Delta a1, Delta c4) profiles are passed through the common Maxwell-Horndeski
primitive emitter and the same kinetic Schur reducer used elsewhere.  This
checks whether the a1-maximizing *on-shell* Horndeski tangent improves all
finite-L kinetic modes simultaneously.

The result is a local tangent/rank diagnostic only, never a finite transition
certificate.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from ssz_p5.coefficients.mh_general_primitives import SLOTS, emit_from_primitives
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central
from ssz_p5.production.strong_field_transition import U_INNER_LIGHT_RING
from ssz_p5.reducer.kinetic_schur import kinetic_schur
from ssz_p5.jets.jet9d8 import derivative

OUTDIR = ROOT / "data/generated/strong_field_transition_2026-09-20"
REPORT = OUTDIR / "TRANSITION_HORNDESKI_A1_PROJECTION.json"
TABLE = OUTDIR / "TRANSITION_HORNDESKI_A1_PROJECTION.csv"
LS = (6, 20, 42, 1000)
SAMPLE_POINTS = (0.7013, 0.7030, U_INNER_LIGHT_RING, 0.7080, 0.7090)
REL_SVD_TOL = 1.0e-10
BUNDLE_HALF_WIDTH_U = 1.25e-3
FD_AMPLITUDE = 1.0e-4


def _sym(a: np.ndarray) -> np.ndarray:
    return (a + a.swapaxes(1, 2)) / 2.0


def _bump(u: np.ndarray, center: float, half_width: float) -> np.ndarray:
    z = (np.asarray(u, float) - center) / half_width
    out = np.zeros_like(z)
    m = np.abs(z) < 1.0
    out[m] = np.exp(-1.0 / (1.0 - z[m] ** 2)) / np.exp(-1.0)
    return out


def _local_responses(frame: pd.DataFrame) -> dict[str, np.ndarray]:
    """Return B rows and primitive gradients for q,p,g0,gx on each grid row.

    Coordinates:
      q  = Delta G4_XX,
      p  = Delta G3_X,
      g0 = Delta G2,
      gx = Delta G2_X.

    The formulas are the unfactored version of the exact quartic-Horndeski
    tubular background construction in production/horndeski_tubular.py.  We do
    not invert a preferred 2x2 pivot; the full 3x4 Jacobian is retained and its
    numerical nullspace is obtained by SVD.
    """
    d = frame.copy().sort_values("x").reset_index(drop=True)
    r = d.x.to_numpy(float)
    f = d.f.to_numpy(float)
    h = d.h.to_numpy(float)
    X = d.X.to_numpy(float)
    ph = d.phiprime.to_numpy(float)
    fp = derivative(r, f, 1)
    hp = derivative(r, h, 1)
    fpp = derivative(r, f, 2)
    phipp = derivative(r, ph, 1)
    Xp = derivative(r, X, 1)
    Xphi = Xp / ph

    def eoms(C1, C2, C4, C6, C7, C9):
        E00 = (C1 + C2 / r) * phipp + (ph * C1 / (2 * h) + C4 / r) * hp + C6 + C7 / r
        E11 = -(ph * C1 / (2 * h) + C4 / r) * h * fp / f + C9 - 2 * ph * C1 / r - (ph * C2 / (2 * h) + (h - 1) * C4) / r**2
        E22 = (
            (C2 * fp / (4 * f) + C1 + C2 / (2 * r)) * phipp
            + (2 * h * C4 - ph * C2) / (4 * f) * (fpp - fp**2 / (2 * f))
            + C4 * fp * hp / (4 * f)
            + C7 * fp / (4 * f)
            + (ph * C1 / h + C4 / r) * hp / 2
            + C6
            + C7 / (2 * r)
        )
        return np.stack([E00, E11, E22], axis=1)

    z = np.zeros_like(r)
    o = np.ones_like(r)
    # Unit q=G4XX response. Holonomy on X=X_b(phi) gives G4_phiX=-q X_phi.
    Bq = eoms(
        -2 * h**2 * Xphi * ph**2,
        4 * h**3 * ph**3,
        2 * h**2 * ph**4,
        -2 * h * Xphi**2 * ph**2,
        4 * h**2 * Xphi * ph**3,
        z,
    )
    # Unit p=G3X response, with G3_phi=-p X_phi on the curve.
    Bp = eoms(-h**2 * ph**2, z, z, -h * Xphi * ph**2, z, -h * Xphi * ph**2)
    # Unit G2 value and G2X responses.
    Bg0 = eoms(z, z, z, o, z, -o)
    Bgx = eoms(z, z, z, z, z, -h * ph**2)
    B = np.stack([Bq, Bp, Bg0, Bgx], axis=2)  # n x 3 x 4

    sqfh = np.sqrt(f * h)
    da1_q = sqfh * (h * Xphi * ph**2 * r**2 - 2 * h**2 * ph**3 * r)
    da1_p = sqfh * (0.5 * h * ph**2 * r**2)
    dc4_q = 0.25 * np.sqrt(f / h) * (
        (h * ph / f * (-2 * h * ph**2)) * fp
        + 4 * h * Xphi * ph**2
        - 4 * h**2 * ph**3 / r
    )
    dc4_p = 0.25 * np.sqrt(f / h) * (2 * h * ph**2)
    ga1 = np.stack([da1_q, da1_p, z, z], axis=1)
    gc4 = np.stack([dc4_q, dc4_p, z, z], axis=1)
    return {"B": B, "ga1": ga1, "gc4": gc4}


def _null_projection(B: np.ndarray, g: np.ndarray) -> dict:
    U, s, Vt = np.linalg.svd(B, full_matrices=True)
    tol = REL_SVD_TOL * (s[0] if len(s) else 1.0)
    rank = int(np.sum(s > tol))
    N = Vt[rank:].T
    if N.shape[1]:
        proj = N @ (N.T @ g)
    else:
        proj = np.zeros_like(g)
    normg = float(np.linalg.norm(g))
    normp = float(np.linalg.norm(proj))
    if normp > 0:
        direction = proj / normp
        if float(g @ direction) < 0:
            direction = -direction
    else:
        direction = np.zeros_like(g)
    return {
        "rank": rank,
        "null_dim": int(B.shape[1] - rank),
        "singular_values": [float(x) for x in s],
        "relative_svd_tolerance": REL_SVD_TOL,
        "projected_norm": normp,
        "raw_norm": normg,
        "projected_fraction": normp / max(normg, 1e-300),
        "best_unit_null_direction": direction,
        "a1_gain_unit_direction": float(g @ direction),
        "background_residual_norm": float(np.linalg.norm(B @ direction)),
    }


def _reference_primitive_frame(base: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame({
        "u": base.u.to_numpy(float),
        "x": base.x.to_numpy(float),
        "phi": base.phi.to_numpy(float),
        "f": base.f.to_numpy(float),
        "h": base.h.to_numpy(float),
        "phiprime": base.phiprime.to_numpy(float),
        "A0prime": base.A0prime.to_numpy(float),
        "a1": np.zeros(len(base)),
        "c2": np.zeros(len(base)),
        "c4": np.zeros(len(base)),
        "F_tensor": np.zeros(len(base)),
        "G_tensor": np.zeros(len(base)),
        "H_tensor": np.ones(len(base)),
    })


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    build = build_onshell_central(ROOT)
    action = build.action.sort_values("x").reset_index(drop=True)
    stream = build.direct41.sort_values("x").reset_index(drop=True)
    resp = _local_responses(action)
    u = action.u.to_numpy(float)

    sampled = []
    for requested in SAMPLE_POINTS:
        i = int(np.argmin(np.abs(u - requested)))
        pr = _null_projection(resp["B"][i], resp["ga1"][i])
        direction = pr.pop("best_unit_null_direction")
        sampled.append({
            "u_requested": float(requested),
            "u_grid": float(u[i]),
            **pr,
            "best_unit_null_direction": {
                "G4XX": float(direction[0]),
                "G3X": float(direction[1]),
                "G2": float(direction[2]),
                "G2X": float(direction[3]),
            },
            "c4_gain_unit_direction": float(resp["gc4"][i] @ direction),
        })

    # Build the pointwise background-null, a1-maximizing Horndeski bundle near LR.
    # A guarded local grid is mathematically identical at the sampled LR row
    # under JET9D8 but avoids recomputing thousands of irrelevant points.
    local_mask_a = np.abs(action.u.to_numpy(float) - U_INNER_LIGHT_RING) < 3.0e-3
    local_mask_s = np.abs(stream.u.to_numpy(float) - U_INNER_LIGHT_RING) < 3.0e-3
    action_h = action.loc[local_mask_a].sort_values("x").reset_index(drop=True)
    stream_h = stream.loc[local_mask_s].sort_values("x").reset_index(drop=True)
    resp_h = _local_responses(action_h)
    u_h = action_h.u.to_numpy(float)
    shape = _bump(u_h, U_INNER_LIGHT_RING, BUNDLE_HALF_WIDTH_U)
    dirs = np.zeros((len(action_h), 4))
    da1 = np.zeros(len(action_h))
    dc4 = np.zeros(len(action_h))
    bg_residual = np.zeros((len(action_h), 3))
    null_ranks = np.zeros(len(action_h), dtype=int)
    for i in np.flatnonzero(shape > 0):
        pr = _null_projection(resp_h["B"][i], resp_h["ga1"][i])
        v = pr["best_unit_null_direction"]
        dirs[i] = shape[i] * v
        da1[i] = float(resp_h["ga1"][i] @ dirs[i])
        dc4[i] = float(resp_h["gc4"][i] @ dirs[i])
        bg_residual[i] = resp_h["B"][i] @ dirs[i]
        null_ranks[i] = pr["rank"]

    ref = _reference_primitive_frame(action_h)
    e0 = emit_from_primitives(ref, regularize_photon_root=True)
    changed = ref.copy()
    changed["a1"] = da1
    changed["c4"] = dc4
    e1 = emit_from_primitives(changed, regularize_photon_root=True)
    delta = {slot: e1[slot].to_numpy(float) - e0[slot].to_numpy(float) for slot in SLOTS}

    idx = int(np.argmin(np.abs(stream_h.u.to_numpy(float) - U_INNER_LIGHT_RING)))
    kinetic = {}
    all_positive = True
    for L in LS:
        K0 = _sym(kinetic_schur(stream_h, L)["K"])
        w0 = np.linalg.eigvalsh(K0[idx])[0]
        vals = []
        for eps in (-FD_AMPLITUDE, FD_AMPLITUDE):
            d = stream_h.copy()
            for slot in SLOTS:
                d[slot] = d[slot].to_numpy(float) + eps * delta[slot]
            K = _sym(kinetic_schur(d, L)["K"])
            vals.append(float(np.linalg.eigvalsh(K[idx])[0]))
        slope = (vals[1] - vals[0]) / (2 * FD_AMPLITUDE)
        kinetic[str(L)] = {
            "baseline_lambda_min": float(w0),
            "lambda_minus_eps": vals[0],
            "lambda_plus_eps": vals[1],
            "d_lambda_min_d_bundle_amplitude": float(slope),
        }
        all_positive = all_positive and slope > 0

    active = shape > 0
    max_bg = float(np.max(np.abs(bg_residual[active]))) if np.any(active) else 0.0
    lr_i = int(np.argmin(np.abs(u_h - U_INNER_LIGHT_RING)))
    lr_proj = next(x for x in sampled if abs(x["u_requested"] - U_INNER_LIGHT_RING) < 1e-8)

    table = pd.DataFrame({
        "u": u_h,
        "shape": shape,
        "delta_G4XX": dirs[:, 0],
        "delta_G3X": dirs[:, 1],
        "delta_G2": dirs[:, 2],
        "delta_G2X": dirs[:, 3],
        "delta_a1": da1,
        "delta_c4": dc4,
        "delta_E00_linear": bg_residual[:, 0],
        "delta_E11_linear": bg_residual[:, 1],
        "delta_E22_linear": bg_residual[:, 2],
        "local_background_rank": null_ranks,
    })
    table.to_csv(TABLE, index=False)

    report = {
        "status": "PASS_A1_SURVIVES_BACKGROUND_NULL_PROJECTION" if lr_proj["projected_fraction"] > 1e-6 else "OPEN_A1_BACKGROUND_NULL_PROJECTION",
        "scope": "local quartic-G5=0 Horndeski tangent projection and kinetic sensitivity; not a finite transition member",
        "action_coordinates": ["G4XX", "G3X", "G2", "G2X"],
        "background_rows": ["E00", "E11", "E22"],
        "rank_policy": {
            "relative_svd_tolerance": REL_SVD_TOL,
            "reason": "the third singular scale can sit at the numerical/identity floor; use the full SVD nullspace rather than the historically singular preferred 2x2 pivot",
        },
        "sample_points": sampled,
        "light_ring_bundle": {
            "u_light_ring_grid": float(u_h[lr_i]),
            "compact_half_width_u": BUNDLE_HALF_WIDTH_U,
            "max_abs_linear_background_residual": max_bg,
            "delta_a1_at_light_ring": float(da1[lr_i]),
            "delta_c4_at_light_ring": float(dc4[lr_i]),
            "kinetic_sensitivity": kinetic,
            "all_tested_L_improve": bool(all_positive),
            "interpretation": (
                "The raw a1 clue survives projection onto a genuine local Horndeski background-null tangent. "
                "However the a1-maximizing on-shell Horndeski bundle is not a common all-L stabilizer if any reported slope is non-positive. "
                "The next predictor must therefore solve the joint H+SVT null-bundle feasibility problem rather than patching a1 alone."
            ),
        },
        "promotion": {
            "a1_offshell_clue_upgraded_to_onshell_tangent_clue": bool(lr_proj["projected_fraction"] > 1e-6),
            "a1_bundle_is_finite_transition_solution": False,
            "a1_bundle_is_common_all_L_stabilizer": bool(all_positive),
            "joint_H_SVT_continuation_required": True,
        },
    }
    REPORT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({
        "status": report["status"],
        "lr_projected_fraction": lr_proj["projected_fraction"],
        "lr_a1_gain": lr_proj["a1_gain_unit_direction"],
        "max_background_residual": max_bg,
        "kinetic_slopes": {L: kinetic[str(L)]["d_lambda_min_d_bundle_amplitude"] for L in LS},
        "all_L_positive": all_positive,
    }, indent=2))
    print(f"WROTE {REPORT}")
    print(f"WROTE {TABLE}")
    return 0 if report["status"].startswith("PASS") else 2


if __name__ == "__main__":
    raise SystemExit(main())
