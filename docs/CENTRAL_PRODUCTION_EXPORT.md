# Central production representation

The central export reproduces the authoritative selected coefficient member on
`0.61 <= u < 0.71`, with its nonzero electric background. Run:

```bash
python tools/export_central_production.py
```

Outputs in `data/generated/central/` contain 3,809 rows, all 41 slots, a complete
slot comparison and a scoped certificate binding source and output SHA256 hashes.
The stored common ZK quantities `v1`, `v4`, `c2`, `v13`, `e4` agree exactly.
The odd quantities are deterministically extracted from the same exact source;
this extraction alone is not an odd-stability certificate.

## Correct source selection

`CENTRAL_PRODUCTION_SOURCES` registers the unreduced coefficient table, corrected
profile, exact ZK table, selected lower jets, normalized reference and odd profile.
The production path preserves authoritative principal coefficients and encoded
`v5`, takes the explicit `c3/e3` completion, then applies the locked `a5/v12`
convention. Derivatives are evaluated on the guard grid before production trimming.

The attempted partial action replay used a different mix of principal and later
transverse selections. At `u=0.6900037509377344`, `x=1.4492674838955182`, the supplied
`f3XX_selected=-38.03585966937742` contributes approximately `13.04545295` to `c2`.
That replay gives `49.28084462`; the authoritative coefficient is `36.23539167`.
A diagnostic replay with that jet zero reproduces the original `c2` to a maximum
scaled error of `5.49e-9`. This isolates the mismatch; it does not authorize changing
the selected action jet to zero. The production export instead uses the complete
stored coefficient representation, as explicitly permitted by the execution
contract. `d3` also depends on the partial-phi versus total-radial selector; the
outer emitter's selected derivative prescription cannot silently replace its
stored central value. Remaining partial-replay differences are recorded in
`PARTIAL_ACTION_REPLAY_DIAGNOSTIC.csv`.

## Certificate scope

`CENTRAL_DIRECT_41=PASS` in this local certificate means **authoritative coefficient
import and common-convention normalization**. It does not mean independent action
to 41 derivation, positive finite-L kinetic eigenvalues, or global closure. These
limitations are machine-readable in the certificate. The existing kinetic gate
and absolute-closure requirements are unchanged. No partial background diagnostic
has been added to production gates. Background recalculation remains
`DIAGNOSTIC_ONLY`; the full selected background has not been independently verified
by this export.

The prior five focused tests passed. Two additional slot-report tests passed,
including detection of changed coefficients, nonfinite values and a mismatched
grid. The published kinetic-formula cross-check remains frozen; this export does
not rerun it. The next production segment is the inner same-action handover.
