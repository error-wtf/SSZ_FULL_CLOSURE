#!/usr/bin/env python3
"""Freeze the 19-Sep action-handover falsification/rank checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ssz_p5.production.action_handover import (
    compact_horndeski_control_rank,
    endpoint_action_diagnostics,
    one_control_a2_probe,
)
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument(
        "--out",
        type=Path,
        default=Path("data/generated/absolute_attempt_2026-09-19/ACTION_HANDOVER_CONTROL_SPACE_AUDIT.json"),
    )
    args = ap.parse_args()
    root = args.root.resolve()
    build = build_onshell_central(root)
    report = {
        "status": "PENDING_JOINT_BACKGROUND_CONSTRAINED_SOLVE",
        "central_member": "ELECTRIC_HYBRID_ONSHELL_CENTRAL",
        "central_bulk_locked_u": [0.62, 0.70],
        "buffers": [[0.61, 0.62], [0.70, 0.71]],
        "endpoint_action_diagnostics": endpoint_action_diagnostics(root, build),
        "one_control_a2_probe": {
            side: one_control_a2_probe(root, side, build) for side in ("outer", "inner")
        },
        "compact_six_primitive_rank": {
            side: compact_horndeski_control_rank(root, side, build) for side in ("outer", "inner")
        },
        "certification": {
            "cinf_partition_available": True,
            "one_control_f4_bridge": False,
            "six_primitive_compact_control_rank_tested": True,
            "joint_background_E00_E11_Ephi_JA_resolve": False,
            "direct41_handover_certified": False,
            "global_action_replay_complete": False,
        },
        "next_gate": "solve joint SVT + (a1,c2,c4,F,G,H) compact controls against full background EOM, then re-emit Direct-41",
    }
    out = args.out if args.out.is_absolute() else root / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    print(f"WROTE {out}")
    # This audit passes as a *checkpoint* only if the dangerous one-control
    # branch is rejected and the intended compact primitive bundle has rank 6.
    ok = (
        not report["one_control_a2_probe"]["outer"]["accepted"]
        and not report["one_control_a2_probe"]["inner"]["accepted"]
        and report["compact_six_primitive_rank"]["outer"]["rank"] == 6
        and report["compact_six_primitive_rank"]["inner"]["rank"] == 6
    )
    return 0 if ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
