# Inner action-control update — superseded split-inverse interpretation

The split lower-order and principal inverses implemented earlier on 2026-09-18 remain valid
as **local response inverses**. They are not a common single-action certificate.

The former 3D holonomy test omitted `Y`. The corrected 4D implementation is
`src/ssz_p5/production/holonomic_hessian_y.py`. It shows that the static electric
background collapses the `f2`-only coefficient response: the five algebraic channels
`(v5,c3,v1,v4,c2)` have local rank no greater than 4. The current independently blended
coefficient targets therefore cannot all arise from one background-null `f2` Hessian.

See `data/generated/inner_y_hessian/INNER_F2_4D_ALGEBRAIC_REACHABILITY.json`.

The correct continuation is a full action-jet construction using mixed `f3/f4` controls and
direct Appendix-A/background-equation re-emission of the lower slots.
