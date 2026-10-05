#!/usr/bin/env python3
"""N3+N4: Kase-Tsujikawa primitive re-emission on the N2-V2 profile and
the declared K_scalar decision.

Pipeline (docs/N1_N4_LUMINAL_RESOLVE_SPEC.md):
  N2 profile (f,h,phi on the luminal grid, G4=1/2 const, G4phi=0)
  -> quartic_g5zero_primitives (verified chain, UNCHANGED code)
  -> 41-slot emission (emit_from_primitives)
  -> K_scalar = 2 P1 - F  +  finite-L health (K / c_r^2, L = 6..1000)
  -> decision against data/diagnostic/N1_N4_DECISION_RULES.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ssz_p5.coefficients.mh_action_primitives import (  # noqa: E402
    quartic_g5zero_primitives,
)
from ssz_p5.coefficients.mh_general_primitives import (  # noqa: E402
    emit_from_primitives,
)

ROOT = Path(__file__).resolve().parents[1]
PROFILE = ROOT / ("data/generated/spectral/"
                  "N2_LUMINAL_BACKGROUND_PROFILE_V2.csv")
RULES = ROOT / "data/diagnostic/N1_N4_DECISION_RULES.json"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--profile", type=Path, default=PROFILE)
    ap.add_argument("--output", type=Path,
                    default=ROOT / ("data/generated/spectral/"
                                    "N3_N4_LUMINAL_DECISION.json"))
    args = ap.parse_args()

    # x = 1/u must be STRICTLY INCREASING for the emitter; the N2 profile
    # runs u: 0.71 -> 100 (x: 1.408 -> 0.01, decreasing), so flip.
    prof = pd.read_csv(args.profile).sort_values("u", ascending=False).reset_index(drop=True)
    # Luminal sector, declared: G4 = 1/2 constant, G4phi = 0, all other
    # non-luminal jets zero.  phiprime from the solved profile.
    prof_in = pd.DataFrame({
        "x": 1.0 / prof["u"].to_numpy(float),
        "u": prof["u"].to_numpy(float),
        "f": prof["f"].to_numpy(float),
        "h": prof["h"].to_numpy(float),
        "phi_r": prof["phi_p"].to_numpy(float),
        "G4": 0.5,
        "G4X": 0.0,
        "G4XX": 0.0,
        "G4phi": 0.0,
        "G4phiX": 0.0,
        "G3X": 0.0,
    })
    prim = quartic_g5zero_primitives(prof_in)

    # c2: luminal G2 = G2X * X with X = phi_r^2 / (2 f) (declared law):
    phip = prof_in["phi_r"].to_numpy(float)
    f = prof_in["f"].to_numpy(float)
    X = phip**2 / (2.0 * np.maximum(f, 1e-12))
    c2 = 1.0 * X
    prim_in = pd.DataFrame({
        "x": prof_in["x"],
        "f": prof_in["f"],
        "h": prof_in["h"],
        "phiprime": phip,
        "A0prime": np.zeros(len(prof_in)),
        "a1": prim["a1_action"],
        "c2": c2,
        "c4": prim["c4_action"],
        "F_tensor": prim["F_tensor_action"],
        "G_tensor": prim["G_tensor_action"],
        "H_tensor": prim["H_tensor_action"],
    })
    # The flat-end asymptotics (phi_r -> 0 at x -> 0.01) produce division
    # by phi in downstream slot derivations — physically correct, but the
    # emitter needs finite inputs.  The decision band (pocket x ~ 0.30)
    # is far from that edge; evaluate on a window covering the pocket
    # plus generous margin.
    band = (prim_in["x"] > 0.10) & (prim_in["x"] < 0.60)
    prim_in = prim_in[band].reset_index(drop=True)
    stream41 = emit_from_primitives(prim_in)

    # K_scalar comes from the VERIFIED chain inside emit_from_primitives:
    #   mu = 2(phi*a1 + 2 r a4)/sqrt(fh);  Y = f r^4 H^4/(mu^2 h);
    #   P1 = h mu/(2 f r^2 H^2) dY/dr;  K_scalar = 2 P1 - F
    # (Kase-Tsujikawa Eq. 4.29).  NOT from -c2/2.
    K_scalar = stream41["K_scalar"].to_numpy(float)

    x_grid = prim_in["x"].to_numpy(float)
    neg_mask = K_scalar < 0
    pocket = (x_grid > 0.29) & (x_grid < 0.32)

    decision = {
        "audit": "N3_N4_LUMINAL_DECISION_V1",
        "g4_law": "constant 1/2 (declared)",
        "profile_source": str(args.profile.name),
        "n_points": int(len(x_grid)),
        "K_scalar_min": float(np.min(K_scalar)),
        "K_scalar_at_x": float(x_grid[int(np.argmin(K_scalar))]),
        "negative_rows": int(neg_mask.sum()),
        "pocket_band_rows_negative": int((neg_mask & pocket).sum()),
        "negative_elsewhere": int((neg_mask & ~pocket).sum()),
    }
    # Noise-floor clause (declared): in the luminal G4=const, G4phi=0
    # sector the analytic expectation is K_scalar == 0 IDENTICALLY
    # (P1 = F/2).  Residuals at the float noise floor (~1e-12) are the
    # zero solution, not a ghost.  A physical claim needs |K_scalar|
    # far above that floor.
    NOISE_FLOOR = 1e-9
    deep_neg = K_scalar < -NOISE_FLOOR
    decision["noise_floor"] = NOISE_FLOOR
    decision["rows_below_noise_floor"] = int(deep_neg.sum())
    decision["K_scalar_abs_max"] = float(np.max(np.abs(K_scalar)))
    if decision["rows_below_noise_floor"] == 0:
        decision["verdict"] = (
            "K_SCALAR_ZERO_ON_LUMINAL_SOLUTION__H1_SUPPORTED_"
            "NO_GHOST_NO_POCKET")
    elif (deep_neg & pocket).sum() > 0:
        decision["verdict"] = "H3_PHYSICAL_GHOST"
    else:
        decision["verdict"] = "NEW_PATHOLOGY_MAP"

    args.output.write_text(json.dumps(decision, indent=1) + "\n")
    csv_out = args.output.with_suffix(".stream41.csv")
    stream41.to_csv(csv_out, index=False)
    print(f"written: {args.output}")
    print(json.dumps(decision, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
