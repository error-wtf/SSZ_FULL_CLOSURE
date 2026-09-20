# SSZ P5 intermediate status — full-action gate

Date: 2026-09-18

This snapshot advances the Inner closure audit beyond the earlier 3D Hessian test.

## Implemented

- 4D background-null `f2(phi,X,F,Y)` Hessian completion.
- Exact `Y=4XF` electric-background invariant handling.
- Appendix-A `f2`-only response for `v5,c3,e3,v1,v4,c2`.
- Sparse/JET-compatible derivative support for the derivative-sensitive `e3` channel.
- Reproducible local reachability audit on the 1,401-row Inner guard grid.
- Full software suite: **96/96 PASS**.

## Result

The old 3D holonomy FAIL is superseded; omitting `Y` was too restrictive. But the stronger
4D test does not yield six independent coefficient controls. For the current static electric
background, the local algebraic map from `(XX,XF,XY,FF,FY,YY)` to
`(v5,c3,v1,v4,c2)` has ranks 1--4 across the guard grid, never 5.

The current independently blended lower/principal target set is outside that local response
image (maximum scaled projection residual about `0.9996427`). This falsifies the proposed
`f2`-only common inverse, **not** the existence of a full SVT action handover.

## Next hard gate

Use the additional action-jet directions already present in Zhang--Kase Appendix A, especially
mixed `f3/f4` jets, construct the Inner handover directly in action-jet space, then derive
`v5,c3,e3` from the background equations and re-emit the complete guarded 41-slot stream.
Only that stream can be used for interface, kinetic, KRGSM and QNM certification.
