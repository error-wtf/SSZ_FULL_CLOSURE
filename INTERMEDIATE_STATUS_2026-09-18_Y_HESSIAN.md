# SSZ P5 intermediate status — 2026-09-18

This is an honest intermediate research snapshot, not a final/full-closure release.

## Reproducible software state

- `pytest -q`: **93 passed** on 2026-09-18.
- `central_direct_41`: PASS.
- `inner_direct_41`: still OPEN / not certified PASS.
- Core/global/KRGSM/QNM closure remains downstream of the inner action closure.

## Important correction to the previous holonomy audit

The previous common-Hessian audit treated `f2` only as a function of `(phi, X, F)` and therefore imposed a symmetric 3x3 Hessian. That audit is now **superseded as a certification gate** because the SVT action/emitter contains the `Y` sector as well.

The relevant even-parity action dependence is at least

`f2(phi, X, F, Y)`

with a symmetric 4x4 Hessian. The on-curve background-null chain rule is therefore

`H_f2 @ (phi', X', F', Y') = 0`.

The six transverse Hessian directions may be represented by

`(XX, XF, XY, FF, FY, YY)`

with the mixed/phi rows reconstructed by the chain rules. This is the next production inverse to implement and test against the coupled lower/principal response channels.

Accordingly, the prior 3D result `INNER_F2_HOLONOMIC_HESSIAN_AUDIT=FAIL` must **not** be interpreted as proof that a single smooth SVT action is impossible. It is retained as historical diagnostic evidence only.

## Included post-ZIP artifacts

This snapshot also includes the later generated artifacts that were not present in the earlier packaged ZIP:

- `data/generated/inner_intermediate/inner_candidatefree_guarded.csv`
- `data/generated/inner_intermediate/core_action_oracle_41.csv`
- `data/generated/inner_intermediate/inner_after_guard_coreoracle.csv`

and the supporting research note:

- `docs/research_notes/KAPPA_Y_TRANSITION_NOTE_2026-09-18.md`

## Next hard gate

Implement the coupled 4D `(phi,X,F,Y)` holonomic Hessian response/inverse, including the reduced endpoint treatment where vector-background responses degenerate, then re-emit the full Inner Direct-41 stream from action jets and test both interfaces before promoting `INNER_DIRECT_41`.
