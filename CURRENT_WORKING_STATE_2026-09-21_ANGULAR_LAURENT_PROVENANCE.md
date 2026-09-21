# Current working state — Angular Laurent provenance pipeline — 2026-09-21

HISTORICAL CHECKPOINT: superseded by [CURRENT_WORKING_STATE_2026-09-21_EQ47.md](CURRENT_WORKING_STATE_2026-09-21_EQ47.md).
The M23-first-divergence interpretation below is retained only as provenance.

## Scope freeze

No Dense-branch angular verdict, no continuation beyond approximately `u=0.708`, no Inner work, and no QNM work were performed in this checkpoint.

## New executable infrastructure

A formal Laurent reducer in `eps=1/L` now propagates the common 41-slot action through the accepted generalized-psi constraint maps, radial product rules, adjoints, and Euler reduction.  The pipeline exports both reduced `K/M` Laurent coefficients and the raw cubic characteristic Laurent polynomial.

Key files:

- `src/ssz_p5/stability/laurent.py`
- `src/ssz_p5/stability/angular_universal_laurent.py`
- `src/ssz_p5/stability/angular_oracles.py`
- `tools/audit_angular_m5_symbolic.py`
- `tools/audit_angular_universal_provenance.py`
- `tools/run_angular_universal_pipeline.py`

## Results

- Formal Laurent vs accepted finite-L common-action reducer: PASS at the method-regression tolerance on `0.62 <= u < 0.70`.
- Genuine-SVT `delta A0` source provenance: symbolic PASS for the action-derived minus branch.
- Published reduced `m5` plus shortcut: retained as `PUBLISHED_REDUCED_SHORTCUT_CONFLICT`; the action is not edited to match it.
- Project-to-paper mass mapping uses the explicit convention bridge `M_paper=-M_project`.
- `M11_0`, `M12_0`, `M13_0`, `M22_0`, and `M33_0` match the action-derived minus structure tightly; the first growing subleading discrepancy in this safety lineage is `M23_0` toward the upper end of the genuine-SVT interval.
- The normalized raw characteristic Laurent polynomial still disagrees strongly with the published coupled-angular reconstruction made from the same extracted coefficient set.  This confirms that the remaining issue is before root sorting and after the already-correct leading structure.

## Maxwell-Horndeski audit B

The conversation records a later exact GM-GHS/EMD formal-Laurent PASS with `(1-z)^3`.  The exact oracle input from that run is not present in the attached safety checkpoint/library search used for this build.  Therefore this repository checkpoint contains the oracle hook and documents the upstream result, but classifies the exact replay here as `NOT_REPLAYED_FROM_ATTACHED_INPUT` rather than inventing a machine PASS.

## Current gate

`ANGULAR = OPEN`.

The next algebraic task is the componentwise mapping of common-action Laurent coefficients into the published `tilde K_ij^(0,1), tilde M_ij^(0,1)` convention, with `M23_0` as the first current divergence marker.  Dense-branch evaluation and `0.708 -> 0.715` continuation remain held.

`ABSOLUTE_FULL_CLOSURE = false`.
