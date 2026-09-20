#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

from ssz_p5.observability import build_inventory

ROOT = Path(__file__).resolve().parents[1]

if __name__ == "__main__":
    print(json.dumps(build_inventory(ROOT), indent=2))
