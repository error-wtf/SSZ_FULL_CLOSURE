# Same-action checkpoint integration — 2026-09-20

Source: `SSZ_P5_2026-09-20_FULL_CLOSURE_CHECKPOINT_SAME_ACTION.zip`.

SHA-256: `6276934f1b359d874814aa8cde7f20ad05f36848121f37f2b215fc91aa4a1fd0`.

All 780 archive files match the supplied unpacked directory byte-for-byte.
Additional older reports, visualization code and research files are retained;
runtime caches are excluded from Git. Source snapshot files were not modified.

## Fresh verification

- Python 3.14.4 with the supplied `requirements.lock`.
- Full pytest: **118 passed, 1 failed**, exit 1, 218.81 seconds.
- Failure: `test_electric_hybrid_recipe_replays_frozen_search_checkpoint`,
  maximum scaled difference `4.2670605451605024e-7`, tolerance `5e-8`.
- Single-thread BLAS reproduces the failure. CSV round-trip parsing gives
  `4.2708693289730614e-7`, slot `a5`, row 3740; it does not repair the mismatch.
- `verify_absolute_closure`: exit 2; direct-global-KRGM certificate absent.
- The supplied 119/119 report is retained as source evidence, not presented as
  the result of this fresh run. No emitter, derivative policy or test tolerance
  was changed to hide this regression.

## Inner continuation

The imported Inner solver exported `historical_reference + emitted_delta`.
The corrected continuation targets `desired - emit(baseline_action)` and exports
`emit(completed_action)` directly. It also verifies grid alignment, checks finite
output and reuses one sparse factorization during refinement. Imported numerical
results remain unchanged in `data/generated/inner_free_e3`.

The direct-total-action run currently stops at an unreachable structurally zero
control row: maximum absolute target residual 1.1671147268319149. All 72 rows are
`v1` at the Inner-Core end. At the exact endpoint `S_SVT=f2F=A0prime=0`, so the
supplied SVT action alone does not contain the required Maxwell/Horndeski baseline.
The historical coefficient offset had concealed this assembly gap. This is not
a physical no-go theorem. No new action is certified and no
new coefficient stream is published as a successful production result.

```bash
PYTHONPATH=src python tools/solve_inner_free_e3_common_action.py --output data/generated/inner_total_action
pytest -q
python ssz_p5_full_pipeline.py --strict
```

Next: supply the explicit common Horndeski/Maxwell baseline action, perform
unreduced H + SVT - shared assembly and bind its action jets before solving the
remaining SVT controls. Background, interfaces, global KRGSM and coupled
QNM gates remain mandatory downstream.
