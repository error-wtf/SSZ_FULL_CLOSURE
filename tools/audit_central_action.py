#!/usr/bin/env python3
"""Replay the central action and its independent on-shell/kinetic checks."""

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from ssz_p5.production.central_action import audit_central_action, central_action_inputs, replay_central_action

if __name__ == "__main__":
    output = ROOT / "build/regional"
    output.mkdir(parents=True, exist_ok=True)
    inputs = central_action_inputs(ROOT)
    inputs.to_csv(output / "central_action_jets.csv", index=False)
    replay_central_action(inputs).to_csv(output / "central_action_replay_41.csv", index=False)
    report = audit_central_action(ROOT)
    (output / "CENTRAL_ACTION_CONSISTENCY.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({key: report[key] for key in ("status", "validation", "witness", "background_max_abs")}, indent=2))
    raise SystemExit(0)  # Diagnostic only; never a production gate.
