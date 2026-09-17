#!/usr/bin/env python3
"""Run implemented regional steps with explicit scope, without global promotion."""

import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from ssz_p5.config import SLOT_NAMES
from ssz_p5.production.regional_coefficients import (
    central_selected,
    regenerate_outer_svt_sector,
    OUTER_HESSIANS,
)
from ssz_p5.production.sources import SOURCE_REGISTRY
from ssz_p5.provenance.manifest import sha256

if __name__ == "__main__":
    out = ROOT / "build/regional"
    out.mkdir(parents=True, exist_ok=True)
    inputs, outer = regenerate_outer_svt_sector(ROOT)
    refpath = ROOT / SOURCE_REGISTRY["outer_same_action_H_SVT"]["coeff_reference"]
    ref = pd.read_csv(refpath)
    np.testing.assert_allclose(ref.x, outer.x, rtol=0, atol=1e-14, err_msg="outer grid mismatch")
    errors = {
        c: float(np.max(abs(outer[c] - ref[c]) / np.maximum(1, abs(ref[c])))) for c in SLOT_NAMES
    }
    inputs.to_csv(out / "outer_svt_action_jets.csv", index=False)
    outer.to_csv(out / "outer_svt_sector_41.csv", index=False)
    central_selected(ROOT).to_csv(out / "central_selected_41.csv", index=False)
    report = dict(
        outer_svt_raw_regression="PASS" if max(errors.values()) < 1e-7 else "FAIL",
        max_scaled_error=max(errors.values()),
        slot_errors=errors,
        scope="Outer ZK/SVT sector replay plus central existing-member normalization; not complete H+SVT assembly or global direct export",
        sources=[
            dict(path=p, sha256=sha256(ROOT / p))
            for p in [
                SOURCE_REGISTRY["outer_same_action_H_SVT"]["background"],
                OUTER_HESSIANS,
                SOURCE_REGISTRY["outer_same_action_H_SVT"]["coeff_reference"],
                SOURCE_REGISTRY["central_exact_SVT"]["unreduced_even"],
            ]
        ],
        global_direct_41="NOT_CERTIFIED",
    )
    (out / "REGIONAL_REGENERATION.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    raise SystemExit(0 if report["outer_svt_raw_regression"] == "PASS" else 2)
