# Inner action-control update — 2026-09-18

The current working snapshot now contains executable action-level control maps for
both classes that were previously missing from the Inner handover.

## 1. Lower-order normal-jet inverse

Implemented in `src/ssz_p5/production/lower_order_controls.py` using the stored
triangular control map

- `q1 = f2_{phi F}` → `v5`,
- `q2 = f2_{X phi}` → `c3`,
- `q3 = f2_{phi phi}` → `e3`.

The C-infinity targets are generated from the frozen central endpoint jets and the
stored `S_SVT,T_H` partition.  The control replay has a maximum scaled residual
of approximately `1.2e-11`; all controls are finite.  When `A0prime=0`, the
`v5` target is identically zero and the reduced `(c3,e3)` block is used without
singular division.  The background remains unchanged because these are transverse
normal jets.

## 2. Principal Hessian inverse

Implemented in `src/ssz_p5/production/inner_principal.py` using the already-tested
raw Appendix-A response

`(f2FF,f2XF,f2XX) -> (v1,v4,c2)`.

The target functions use the stored C-infinity partition between one-sided
second-order endpoint Taylor jets.  Production coefficients are not interpolated:
the target is inverted to action Hessians and the complete emitter delta is added
to the same-action baseline.  Maximum scaled target replay residual is below
`4e-16`.

At the central interface the controlled values now match the frozen central
member at machine precision:

- `v5`: ~`1.6e-15` scaled error,
- `c3`: ~`3.4e-15`,
- `e3`: ~`7.1e-13`,
- `v1`: ~`1.4e-15`,
- `v4`: exact,
- `c2`: ~`3.8e-15`.

At the core interface `v5,c3,e3,v4` are exact, `c2` is roundoff-level, and
`v1` differs by only ~`3.4e-9` where the raw SVT vector response has already
vanished (`A0prime=0`).

## 3. Remaining Inner Direct-41 failure

`INNER_DIRECT_41` remains `FAIL`, but the reason has changed materially.  It is
no longer a lower-order or principal-control failure.  The remaining value-level
interface residuals are confined to derivative-sensitive higher slots, dominated
by

- `a5` (especially the core interface),
- `d3`,
- `e4`,
- `v13`,
- a small central `v9` mismatch.

These require the remaining higher mixed/holonomic action-jet re-emission across
the Inner handover.  No tolerance was relaxed and no coefficient endpoint was
manually overwritten.

Software status after this update: **86 tests pass**.
