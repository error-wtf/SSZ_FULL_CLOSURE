# REPRODUCE HANDOFF (relative to ZIP root)

Env: CPython 3.14, Linux. `python3.14 -m venv .venv && . .venv/bin/activate && python -m pip install -r requirements.lock && python -m pip install -e . --no-deps`

1. Manifest verify: `python tools/release_manifest.py --check` (may differ from shipped MANIFEST if files added; use SHA256SUMS.txt for the ZIP inventory)
2. Historical generator replay: collect the inputs listed in HISTORICAL_DATA_LINEAGE.json from data/production|authoritative|diagnostic + archive/full_working_snapshot into one flat dir, copy `historical_sources/build_ssz_p5_selected_41stream_v2_2026-09-16.py` there, set its `B=` to that dir, run it; compare `ssz_p5_SELECTED_41STREAM_V2_2026-09-16.csv` against `data/diagnostic/...` (expect 19411 rows, 40/41 slots <1e-10, a5 ~1.7e-8)
3. Instrumentation no-op: re-run with the PRE/POST export hooks (see Q2_INSTRUMENTATION_NOOP_AUDIT.json); require final equality 0.0
4. Eq85 replay: IMPLEMENT FIRST per CHATGPT_HANDOFF section I (independent from paper; do NOT copy src/ssz_p5/action/light_ring_electric.py); evaluate both rings on the carrier CSV; targets inner F_req=-1.6386 (vs H=+1.0), outer +0.9721
5. Context guard: `python tools/context_guard.py --claim "healthy K does not exist"` (must BLOCK); `python tools/context_guard.py --apply-relation SVT_SPLIT_OUTER_HANDOVER --target LOBE_CENTRAL` (must BLOCK)
6. Q2 continuation: follow CHATGPT_HANDOFF section B decision tree
7. Tests: `OPENBLAS_NUM_THREADS=1 PYTHONPATH=src python -m pytest -q tests/unit/test_eq47_projection.py tests/unit/test_angular_laurent_series.py tests/unit/test_angular_universal_provenance.py` (expect 7 passed)
