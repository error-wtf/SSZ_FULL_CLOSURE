#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
python tools/generate_global_action_cover.py
python tools/verify_codex_handoff.py
pytest -q
python ssz_p5_full_pipeline.py --strict
printf '\nBaseline verified. Now read CODEX_FINAL_EXECUTION_CONTRACT.md and continue forward only.\n'
