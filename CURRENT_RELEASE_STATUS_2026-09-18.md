# Current release status — 2026-09-18

This file records the exact working state packaged into `SSZ_P5_BEST_AVAILABLE_WORKING_2026-09-18.zip` after the latest Inner lower-order-control integration.

## Reproducible software / archive checks

- Focused lower-order-control tests: **2 / 2 PASS**.
- `python -m pytest -q`: **85 / 85 PASS**.
- The full release pipeline is regenerated after manifest refresh before packaging.

## New Inner implementation

The Inner handover now contains an explicit background-null lower-order action-control inverse for
`q1=f2_{phi F}`, `q2=f2_{X phi}`, `q3=f2_{phi phi}` -> `v5,c3,e3`.

The control replay is PASS, finite, background-null by transverse-jet construction, and preserves the frozen Central member and the principal-symbol Hessian controls. The generated artifacts are:

- `data/generated/inner_controls/INNER_LOWER_ORDER_TARGETS.csv`
- `data/generated/inner_controls/INNER_LOWER_ORDER_ACTION_CONTROLS.csv`
- `data/generated/inner_controls/INNER_LOWER_ORDER_REEMITTED.csv`
- `data/generated/inner_controls/INNER_LOWER_ORDER_CONTROL.json`

`SLOT_NORMALIZATION = PASS` and the 1,400-row Inner stream is finite. `INNER_DIRECT_41` remains **FAIL** because the current interface certification still finds dependent-slot / derivative continuity residuals, dominated by `a5`; that is the next implementation target.

## Strict closure

- `FULL_CONSTRUCTIVE_CLOSURE = PASS`
- `DIRECT_GLOBAL_KRGM_EXPORT = OPEN`
- `PERFECT_PASS = false`

No tolerance or rejection gate has been weakened.
