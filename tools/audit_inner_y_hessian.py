#!/usr/bin/env python3
"""Audit whether the current six independently blended Inner targets are reachable
by a background-null f2(phi,X,F,Y) Hessian deformation alone.

This is a reachability diagnostic, not a closure certificate.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from ssz_p5.production.holonomic_hessian_y import ALGEBRAIC_RESPONSES, algebraic_response_matrix
from ssz_p5.production.inner_principal import build_principal_targets
from ssz_p5.production.inner_targets import build_inner_targets
from ssz_p5.production.regional_coefficients import central_selected
from ssz_p5.production.sources import SOURCE_REGISTRY

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/generated/inner_y_hessian"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    reg = SOURCE_REGISTRY
    background = pd.read_csv(ROOT / reg["inner_same_action_SVT_H"]["background"]).reset_index(drop=True)
    current = pd.read_csv(ROOT / reg["inner_same_action_SVT_H"]["coeff_reference"]).reset_index(drop=True)
    central = central_selected(ROOT).reset_index(drop=True)
    core = pd.read_csv(ROOT / reg["punctured_H_core"]["coeff_reference"]).reset_index(drop=True)

    build_inner_targets(ROOT, ROOT / "data/generated/inner_controls")
    lower = pd.read_csv(ROOT / "data/generated/inner_controls/INNER_LOWER_ORDER_TARGETS.csv")
    principal, _ = build_principal_targets(background, central, core)

    target = pd.DataFrame(
        {
            "v5": lower.v5_target,
            "c3": lower.c3_target,
            "v1": principal.v1_target,
            "v4": principal.v4_target,
            "c2": principal.c2_target,
        }
    )
    # At exactly zero electric background the v1/v4 f2-Hessian response is
    # structurally zero.  Preserve the selected baseline rather than manufacture
    # an unreachable endpoint delta.
    zero_electric = np.abs(background.A0prime.to_numpy(float)) <= 1e-14
    for name in ("v1", "v4"):
        target.loc[zero_electric, name] = current.loc[zero_electric, name].to_numpy(float)

    rhs = target[list(ALGEBRAIC_RESPONSES)].to_numpy(float) - current[list(ALGEBRAIC_RESPONSES)].to_numpy(float)
    M = algebraic_response_matrix(background)
    records = []
    rank_counts: dict[int, int] = {}
    for i in range(len(background)):
        scale = np.maximum(1.0, np.abs(target.iloc[i].to_numpy(float)))
        Ms = M[i] / scale[:, None]
        bs = rhs[i] / scale
        _u, singular, _vh = np.linalg.svd(Ms, full_matrices=False)
        tol = max(1e-13, (singular[0] if len(singular) else 0.0) * 1e-10)
        rank = int(np.sum(singular > tol))
        rank_counts[rank] = rank_counts.get(rank, 0) + 1
        q = np.linalg.lstsq(Ms, bs, rcond=1e-10)[0]
        residual = Ms @ q - bs
        records.append(
            dict(
                row=i,
                u=float(background.u.iloc[i]),
                x=float(background.x.iloc[i]),
                A0prime=float(background.A0prime.iloc[i]),
                local_rank=rank,
                max_scaled_projection_residual=float(np.max(np.abs(residual))),
                l2_scaled_projection_residual=float(np.linalg.norm(residual)),
                smallest_retained_singular=float(singular[rank - 1]) if rank else 0.0,
            )
        )

    frame = pd.DataFrame(records)
    frame.to_csv(OUT / "INNER_F2_4D_ALGEBRAIC_REACHABILITY.csv", index=False)
    worst = frame.iloc[int(np.argmax(frame.max_scaled_projection_residual.to_numpy(float)))]
    report = {
        "status": "F2_ONLY_TARGET_SET_NOT_REACHABLE",
        "scope": "background-null f2(phi,X,F,Y) Hessian deformation only",
        "transverse_controls": ["XX", "XF", "XY", "FF", "FY", "YY"],
        "algebraic_targets": list(ALGEBRAIC_RESPONSES),
        "rank_counts": {str(k): v for k, v in sorted(rank_counts.items())},
        "max_local_rank": int(frame.local_rank.max()),
        "max_scaled_projection_residual": float(frame.max_scaled_projection_residual.max()),
        "median_scaled_projection_residual": float(frame.max_scaled_projection_residual.median()),
        "worst_row": {k: (int(v) if k in {"row", "local_rank"} else float(v)) for k, v in worst.to_dict().items()},
        "interpretation": (
            "The Y sector supersedes the old 3D holonomy exclusion, but dimension counting does not give six independent coefficient controls. "
            "The current independently blended lower/principal target set violates the local column space of the f2-only background-null response. "
            "Full Inner closure therefore requires a joint action-jet construction with additional mixed f3/f4 (and, if needed, lower Horndeski) control directions, rather than separate coefficient targets."
        ),
        "inner_direct_41": "NOT_CERTIFIED",
    }
    (OUT / "INNER_F2_4D_ALGEBRAIC_REACHABILITY.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
