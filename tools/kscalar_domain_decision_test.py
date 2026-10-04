#!/usr/bin/env python3
"""H1/H3 DECISION TEST: K_scalar pocket vs action-jet domain.

Determines whether the recomputed negative K_scalar pocket (x~0.294-0.324) is
an artifact of OUT-OF-DOMAIN action jets or physical.

Background (verified in this repo):
  * The project-locked emitter domain is G4=G4(phi) => G4X=G4XX=G4phiX=0,
    G5=0, with G3=0 (luminal Einstein branch).  The production source tables
    nonetheless contain G4X~-2.6e3, G4XX~1.1e13, G4phiX~-1.5e10, G3X~-5e10
    (audit_production_sources.py: jets_outside_luminal_emitter_domain;
    global_direct_generation=NOT_ESTABLISHED).
  * quartic_g5zero_primitives builds a1 (and c4) partly from these jets, so
    out-of-domain jets feed directly into P1/K_scalar.

Test: recompute the emitted K_scalar in the pocket band for jet configurations
  (1) original            (out-of-domain jets as archived)
  (2) G4-jets zeroed       (G4X=G4XX=G4phiX=0)
  (3) G3X zeroed
  (4) luminal domain       (G4-jets AND G3X zeroed)

Decision rule (declared before running):
  pocket disappears (min K >= 0, no negative rows) under (4)
    => ARTIFACT of out-of-domain jets (H1 confirmed, concrete mechanism).
  pocket survives (4)
    => physical instability candidate (H3) on clean ground.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from ssz_p5.coefficients.mh_action_primitives import quartic_g5zero_primitives
from ssz_p5.coefficients.mh_general_primitives import emit_from_primitives

ACT = ROOT / "data/production/ssz_p5_F1b_G4XX_transverse_core_candidate_2026-09-14.csv"
REF = ROOT / "data/authoritative/ssz_p5_F2_core_punctured_horndeski_unreduced_39of41_2026-09-14.csv"
OUT = ROOT / "data/generated/qnm_global_diagnostic/KSCALAR_DOMAIN_DECISION_TEST.json"

G4_JETS = ("G4X", "G4XX", "G4phiX")


def kmin(a: pd.DataFrame, r_ref: pd.DataFrame):
    q = quartic_g5zero_primitives(a)
    inp = pd.DataFrame({"x": q.x, "u": q.u, "f": a.A_f, "h": a.B_h,
                        "X": a.X, "phi": a.phi, "A0prime": 0.0,
                        "a1": q.a1_action, "c2": r_ref.c2, "c4": q.c4_action,
                        "F_tensor": q.F_tensor_action, "G_tensor": q.G_tensor_action,
                        "H_tensor": q.H_tensor_action})
    z = emit_from_primitives(inp, regularize_photon_root=True)
    k = z.K_scalar.to_numpy(float)
    band = (z.x > 0.29) & (z.x < 0.33)
    kb = np.where(band, k, np.nan)
    j = int(np.nanargmin(kb))
    return {"K_min_band": float(np.nanmin(kb)),
            "argmin_x": float(z.x.iloc[j]),
            "negative_rows_band": int(np.nansum(kb < 0))}


def main():
    a = pd.read_csv(ACT).sort_values("r_over_rs").reset_index(drop=True)
    r_ref = pd.read_csv(REF).sort_values("x").reset_index(drop=True)

    cases = {}
    cases["original_out_of_domain_jets"] = kmin(a, r_ref)

    b = a.copy()
    for c in G4_JETS:
        b[c] = 0.0
    cases["g4_jets_zeroed"] = kmin(b, r_ref)

    c = a.copy()
    c["G3X"] = 0.0
    cases["g3x_zeroed"] = kmin(c, r_ref)

    d = a.copy()
    for col in G4_JETS:
        d[col] = 0.0
    d["G3X"] = 0.0
    cases["luminal_domain_g4jets_and_g3x_zero"] = kmin(d, r_ref)

    lum = cases["luminal_domain_g4jets_and_g3x_zero"]
    verdict = {
        "pocket_gone_in_luminal_domain": bool(lum["K_min_band"] >= 0.0 and
                                              lum["negative_rows_band"] == 0),
        "decision_rule": ("pocket gone under luminal domain => ARTIFACT (H1); "
                          "pocket survives => PHYSICAL candidate (H3)"),
    }
    rep = {"status": "KSCALAR_DOMAIN_DECISION_TEST_COMPLETE",
           "declared_before_running": True,
           "cases": cases, "verdict": verdict,
           "verdict_label": ("POCKET_IS_OUT_OF_DOMAIN_JET_ARTIFACT_H1"
                             if verdict["pocket_gone_in_luminal_domain"]
                             else "POCKET_SURVIVES_LUMINAL_DOMAIN_H3")}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(rep, indent=1))
    print(json.dumps(cases, indent=1))
    print("VERDICT:", rep["verdict_label"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
