#!/usr/bin/env python3
"""Inspect the specifically supplied action profiles without substituting old 41 tables."""

import json
from pathlib import Path

import numpy as np
import pandas as pd

from ssz_p5.provenance.manifest import sha256

ROOT = Path(__file__).resolve().parents[1]
SOURCES = (
    "ssz_p5_horndeski_carrier_through_light_rings_to_core_2026-09-12.csv",
    "ssz_p5_F1b_G4XX_transverse_core_candidate_2026-09-14.csv",
    "ssz_p5_F1b_FULL_G5_subcore_onshell_candidate_2026-09-14.csv",
    "ssz_p5_F1b_FINAL_Cinf_center_to_punctured_handover_2026-09-14.csv",
)


def audit(root=ROOT):
    records = []
    for name in SOURCES:
        path = root / "data/production" / name
        d = pd.read_csv(path)
        r = d["r_over_rs" if "r_over_rs" in d else "r"].to_numpy(float)
        unsupported = {}
        for key in ("G4X", "G4XX", "G4phiX", "G5X", "G5phi"):
            if key in d and np.max(np.abs(d[key])) > 1e-12:
                i = int(np.argmax(np.abs(d[key])))
                unsupported[key] = {"value": float(d[key].iloc[i]), "r": float(r[i])}
        records.append(
            {
                "path": str(path.relative_to(root)),
                "sha256": sha256(path),
                "rows": len(d),
                "radius_min": float(r.min()),
                "radius_max": float(r.max()),
                "columns": list(d.columns),
                "jets_outside_luminal_emitter_domain": unsupported,
            }
        )
    return {
        "scope": "source-to-emitter contract; not a claim of physical instability",
        "sources": records,
        "emitter": "src/ssz_p5_mh_luminal_g4phi_emitter_JET9D8_2026-09-16.py",
        "emitter_domain": "G4=G4(phi), G4X=0, G5=0; c2 profile required",
        "global_direct_generation": "NOT_ESTABLISHED",
        "notes": [
            "Finite sampled intervals do not themselves specify the analytic endpoints.",
            "Core action jets require their quartic/quintic coefficient terms.",
            "No electric handover slots or old reduced matrices are production inputs.",
        ],
    }


if __name__ == "__main__":
    payload = audit()
    output = ROOT / "data/diagnostic/PRODUCTION_SOURCE_CONTRACT.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2) + "\n")
    print(output)
