#!/usr/bin/env python3
"""Fast deterministic verification of the self-contained Codex handoff.

This validates handoff completeness and regional source placement. It does not
mislabel remaining scientific direct-closure work as completed.
"""

from pathlib import Path
import json
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from ssz_p5.production.regions import ProductionRegion, classify_x, classify_u
from ssz_p5.production.sources import SOURCE_REGISTRY

SLOTS = [
    *[f"a{i}" for i in range(1, 10)],
    *[f"b{i}" for i in range(1, 6)],
    *[f"c{i}" for i in range(1, 7)],
    *[f"d{i}" for i in range(1, 5)],
    *[f"e{i}" for i in range(1, 5)],
    *[f"v{i}" for i in range(1, 14)],
]

checks = []


def check(name, ok, detail=""):
    checks.append((name, bool(ok), detail))
    print(f"{'PASS' if ok else 'FAIL':4s}  {name}: {detail}")


for region, items in SOURCE_REGISTRY.items():
    for k, v in items.items():
        if k == "role":
            continue
        p = ROOT / v
        check(f"source:{region}:{k}", p.exists(), v)

# Important pre-staged schemas and domains.
cases = [
    (
        "outer",
        ROOT / "data/prestaged/direct41/outer_resolved_raw_41of41.csv",
        41,
        0.5515230871346237,
        0.61,
    ),
    (
        "central",
        ROOT / "data/prestaged/direct41/central_selected_corrected_41of41.csv",
        41,
        0.61,
        0.715,
    ),
    (
        "inner",
        ROOT / "data/prestaged/direct41/inner_selected_candidate_41of41.csv",
        41,
        0.71,
        0.715,
    ),
    ("core", ROOT / "data/prestaged/direct41/core_selected_41of41.csv", 41, 0.71, 100.0),
    (
        "weak",
        ROOT / "data/prestaged/direct41/weak_exterior_39of41.csv",
        39,
        0.020065706,
        0.5515230871346237,
    ),
]
for name, p, nslots, lo, hi in cases:
    d = pd.read_csv(p)
    present = sum(c in d.columns for c in SLOTS)
    check(f"schema:{name}", present >= nslots, f"{present}/41 slots; rows={len(d)}")
    if "u" in d:
        umin = float(d.u.min())
        umax = float(d.u.max())
        check(f"domain:{name}", umin <= hi and umax >= lo, f"u=[{umin:.12g},{umax:.12g}]")

odd = ROOT / "data/prestaged/odd/central_svt_odd_profile.csv"
d = pd.read_csv(odd)
check(
    "central_odd_profile",
    {"Codd", "odd_angular_margin_lb", "u"}.issubset(d.columns),
    f"rows={len(d)}",
)

# The two previously misclassified points must remain central genuine-SVT.
check(
    "x=1.599997795 region",
    classify_x(1.5999977950251312) == ProductionRegion.CENTRAL_EXACT_SVT,
    classify_x(1.5999977950251312).value,
)
check(
    "stable-light-ring region",
    classify_x(1.41616064) == ProductionRegion.CENTRAL_EXACT_SVT,
    classify_x(1.41616064).value,
)

for p in ["START_HERE_CODEX.md", "CODEX_FINAL_EXECUTION_CONTRACT.md", "LICENSE", "CITATION.cff"]:
    check(f"repo:{p}", (ROOT / p).exists(), p)

failed = [x for x in checks if not x[1]]
payload = {
    "checks": len(checks),
    "passed": len(checks) - len(failed),
    "failed": len(failed),
    "handoff_ready": not failed,
    "note": "Handoff readiness is not the same as an Absolute Full Closure certificate.",
}
(ROOT / "CODEX_HANDOFF_VERIFICATION.json").write_text(json.dumps(payload, indent=2) + "\n")
print("\n", json.dumps(payload, indent=2))
raise SystemExit(0 if not failed else 2)
