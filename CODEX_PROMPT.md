# Codex task — repository engineering only; do not redo SSZ research

You are receiving a frozen SSZ P5 research repository. Your first instruction is **do not redo completed physics or re-open old branches**.

Read in this exact order:
1. `RESEARCH_STATE.md`
2. `DO_NOT_RECOMPUTE.md`
3. `ARTIFACT_CLASSIFICATION.json`
4. `README.md`
5. existing tests and source modules

Use `data/production/` as current production data. Use `data/qnm/` for already-performed spectral work. Use `archive/full_working_snapshot/` only for provenance/recovery.

Do not infer authority from historical filenames. Do not change a scientific PASS/FAIL by weakening tolerances. Do not replace an archived result with a fresh derivation unless the repository itself contains a reproducibility test that explicitly requires regeneration.

Primary software goal: turn the frozen scientific assets into a clean, maintainable, tested package while preserving all numerical results and provenance. Refactor code only with regression tests against the archived outputs. Never silently reinterpret test-field/eikonal QNM results as a coupled HSVT QNM spectrum.
