# Inner production export

Baseline: `0e13fb0`. The central export is unchanged.

```bash
python tools/export_inner_production.py
```

The importer consumes the registered resolved inner background and selected
inner coefficient reference. It preserves all selected slots except the two
explicit convention updates (`a5`, `v12`), computes derivatives on the full
1,401-row guard grid, then exports the 1,400 production rows on
`0.71 <= u < 0.715`. All 41 slots are finite. The actual stored background
residuals remain `8.881784e-16`, `4.662937e-15`, `3.072771e-15`; no background
was re-solved and no diagnostic sector sum was used.

## Result

- `SLOT_NORMALIZATION = PASS`
- `INNER_DIRECT_41 = FAIL`: endpoint continuity is not satisfied.
- `DIRECT_ACTION_REPLAY_COMPLETE = false`
- `AUTHORITATIVE_SELECTED_REPRESENTATION = true`

The selected inner reference explicitly encodes `v5=c3=e3=0`. At the central
interface, the frozen central export has estimated boundary values:

| Slot | Inner at u=0.71 | Frozen central limit |
| --- | ---: | ---: |
| v5 | 0 | -3873.246323 |
| c3 | 0 | 337.204402 |
| e3 | 0 | 323.087359 |
| a5 | 696.606674 | -2030.771487 |

Boundary values are local polynomial evaluations under the existing JET9D8
policy, using the frozen central end stencil; central was not reimported or
modified. Window variations 7/6, 9/8 and 11/8 retain the zero-versus-nonzero
lower-order discrepancy. The `a5` difference therefore cannot be attributed
solely to differentiation at a cut edge. The core interface and radial orders
0, 1 and 2 for every slot are also recorded in the endpoint CSV.

The next correction must reconcile the **inner selected lower-order completion**
with the frozen central endpoint. Merely setting the central values to zero,
relabeling a failed interface as PASS, or inventing a coefficient interpolation
would change the selected representation without its required action/provenance
justification. This is a concrete interface failure, not a no-go statement about
the geometry or the full theory.

## Evidence and tests

`data/generated/inner/` contains the normalized stream, all-slot comparison,
endpoint comparison, derivative-window comparison and hashed certificate.
`ABSOLUTE_CLOSURE_LEDGER.json` preserves central PASS and records inner FAIL;
core/global/operator/spectral stages remain pending.

Two focused software tests passed: polynomial endpoint jets and detection of an
isolated lower-order slot jump. No central, outer, full-suite or spectral retest
was run for this step.
