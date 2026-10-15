#!/usr/bin/env python3
"""Fail-closed gate for the actual SSZ spectral-weight question."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/generated/spectral_selection/SPECTRAL_WEIGHT_GATE.json"


def main():
    direct = ROOT / "data/certificates/DIRECT_GLOBAL_KRGM_CERTIFICATE.json"
    direct_alt = ROOT / "data/certificates/SSZ_P5_DIRECT_GLOBAL_KRGM_CERTIFICATE.json"
    if not direct.is_file() and not direct_alt.is_file():
        payload = {
            "status": "SPECTRAL_SELECTION_NOT_YET_EVALUABLE",
            "blocker": "MISSING_CURRENT_MEMBER_CENTER_TO_INFINITY_DIRECT_KRGM",
            "scope": "observable modal residues Z_n(r) / rho(r,omega)",
            "principal_local_null_is_not_spectral_weight_null": True,
        }
        OUT.parent.mkdir(parents=True, exist_ok=True)
        OUT.write_text(json.dumps(payload, indent=2) + "\n")
        print(json.dumps(payload))
        return 2

    payload = {
        "status": "DIRECT_KRGM_PRESENT_REQUIRES_RESIDUE_SOLVER",
        "blocker": "GLOBAL_GREEN_RESIDUE_SOLVER_NOT_RUN",
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload))
    return 3


if __name__ == "__main__":
    raise SystemExit(main())
