#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
registry = json.loads((ROOT / "REPO_EVIDENCE_REGISTRY.json").read_text())
for entry in [*registry["audits"], *registry.get("full_checks", [])]:
    prefix = "python "
    command = " ".join(entry["command"])
    print(f"{entry['id']}: {prefix}{command}  [expected exit {entry.get('expected_exit', 0)}]")
