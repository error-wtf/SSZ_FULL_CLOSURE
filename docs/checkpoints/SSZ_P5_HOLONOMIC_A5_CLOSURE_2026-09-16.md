# SSZ P5 — Holonomic Maxwell–Horndeski `a5` closure
**Date:** 2026-09-16

The exact Maxwell–Horndeski Appendix-A identities were combined so that the selected pure-Horndeski `a5` no longer requires a second numerical differentiation of the archived `a1` column.  For `A0'=0`,

```math
a_2=a_1'-P a_1+R,\qquad a_5=-P'a_1-Pa_1'+R',
```
with `P=phi''/phi'+h'/(2h)` and `R=r/phi'(f'/f-h'/h)a4`.

## Independent regression
On the smooth strong-H carrier the regenerated `a5` has median scaled relative residual **8.111757e-07**, 95% **5.433474e-06**, and maximum **6.636212e-04**.

## Punctured-core selected member
The old pointwise core `a5` column is derivative-noise contaminated.  A deterministic selected member was therefore generated from the first-order `a2` identity.  Comparing zero smoothing with the very small `s=1e-8` `a2` regularization gives median `a5` change **1.935829e-07**, 95% **3.416096e-05**, and maximum **1.495892e-04**.  The corresponding selected `a1` remains within a maximum scaled distance **1.877440e-04** of the archived full-rank control target.

All three gates pass.  This fixes the old missing-`-a1''` implementation error without inventing a physical instability or differentiating a rounded coefficient table twice.
