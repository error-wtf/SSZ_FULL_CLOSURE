# Current release status — 2026-09-18

This is the current reproducible intermediate state. It is not an Absolute Full Closure release.

## Software status

- `pytest -q`: **96 / 96 PASS**.
- New 4D `f2(phi,X,F,Y)` holonomic completion tests: PASS.
- New Inner `f2`-only algebraic reachability audit: executed and preserved.

## Corrected Inner finding

The previous 3D `(phi,X,F)` common-Hessian audit was incomplete because the SVT action
contains the `Y` direction. The 4D background-null completion with transverse controls
`(XX,XF,XY,FF,FY,YY)` is now implemented in
`src/ssz_p5/production/holonomic_hessian_y.py`.

However, the six transverse Hessian entries do **not** become six independent coefficient
controls on the static electric background. For the five algebraic channels
`(v5,c3,v1,v4,c2)`, the audited local map has maximum rank **4**. The present lower/principal
targets were built independently and have a maximum scaled projection residual of about
`9.996e-1`; therefore they cannot be promoted as a single background-null `f2` action.

This does **not** imply that the Inner handover or P5 geometry is impossible. It means the
control space must be enlarged to the full action jets already present in Appendix A, most
importantly the mixed `f3/f4` directions, and the target should be constructed at action level
rather than by independently blending six coefficient functions.

Evidence:

- `data/generated/inner_y_hessian/INNER_F2_4D_ALGEBRAIC_REACHABILITY.json`
- `data/generated/inner_y_hessian/INNER_F2_4D_ALGEBRAIC_REACHABILITY.csv`

## Gates that remain unchanged

- `central_direct_41`: scoped normalization PASS.
- Frozen selected Central finite-L kinetic gate: **FAIL** for `L=6,12,20,42`.
- `inner_direct_41`: **NOT CERTIFIED**.
- Core Direct-41 / analytic center / global Direct-41 / KRGSM / coupled QNM: pending.
- `v1.0.0`: not released.

No tolerance or rejection gate was weakened.
