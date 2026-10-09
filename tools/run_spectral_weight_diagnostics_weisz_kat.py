#!/usr/bin/env python3
"""SSZ_SPECTRAL_WEIGHT_DIAGNOSTICS — Step 1: validate diagnostics on the KNOWN
Weisz 1978 FK reproduction (Known-Answer-Test) before touching the SSZ operator.

Diagnostics (frozen in SSZ_SPECTRAL_WEIGHT_DIAGNOSTICS_PREREGISTRATION_V1):
  - selectivity C = max_i p_i, p_i = w_i / sum(w)
  - spectral entropy S = -sum p_i ln p_i
  - participation ratio PR_i (basis-invariant)

KAT acceptance (frozen): computed C for case (13,10,0.3) must reproduce the
archived 0.9712 (dominant q0_residue) within 1%.
"""
import csv
import json
import math
from pathlib import Path

ROOT = Path("/home/error/physics/clones/SSZ_FULL_CLOSURE")
CSV = ROOT / "data/generated/spectral/weisz1978/WEISZ1978_FK_MODES.csv"
REF = ROOT / "data/generated/spectral/weisz1978/WEISZ1978_FK_REPRODUCTION.json"
OUT = ROOT / "artifacts/SSZ_SPECTRAL_WEIGHT_DIAGNOSTICS_WEISZ_KAT_V1.json"


def diagnostics(weights):
    total = sum(weights)
    p = [w / total for w in weights]
    C = max(p)
    S = -sum(pi * math.log(pi) for pi in p if pi > 0)
    # PR on the "mode axis" (1D mode-index space; true spatial PR needs eigenvectors)
    pr = (sum(p) ** 2) / (len(p) * sum(pi * pi for pi in p))
    return {"selectivity_C": C, "spectral_entropy_S": S,
            "participation_ratio_PR": pr, "n_modes": len(p)}


def main():
    rows = list(csv.DictReader(open(CSV)))
    cases = {}
    for r in rows:
        key = (int(r["N"]), int(r["M"]), float(r["beta"]))
        cases.setdefault(key, []).append(float(r["q0_residue"]))

    ref = json.load(open(REF))
    results = {}
    kat_ok = None
    for key, w in sorted(cases.items()):
        res = diagnostics(w)
        results[str(key)] = res
        if key == (13, 10, 0.3):
            kat_ok = abs(res["selectivity_C"] - 0.9712499743350262) < 0.01 * 0.9712499743350262

    report = {
        "artifact": "SSZ_SPECTRAL_WEIGHT_DIAGNOSTICS_WEISZ_KAT_V1",
        "purpose": "Known-Answer-Test of the frozen diagnostics on the archived Weisz reproduction",
        "preregistration": "SSZ_SPECTRAL_WEIGHT_DIAGNOSTICS_PREREGISTRATION_V1",
        "diagnostics_per_case": results,
        "KAT_case_13_10_03": {
            "archived_selectivity": 0.9712499743350262,
            "computed": results.get("(13, 10, 0.3)") or results.get("13, 10, 0.3"),
            "pass_1pct": bool(kat_ok),
        },
    }
    OUT.write_text(json.dumps(report, indent=1) + "\n")
    print(json.dumps(report["KAT_case_13_10_03"], indent=1))
    print("wrote", OUT)
    return 0 if kat_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
