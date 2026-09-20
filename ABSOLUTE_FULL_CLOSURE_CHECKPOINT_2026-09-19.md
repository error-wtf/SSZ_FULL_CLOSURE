# SSZ P5 — Absolute Full Closure working checkpoint

**Date:** 2026-09-19  
**Status:** advanced working checkpoint; **ABSOLUTE FULL CLOSURE NOT YET CERTIFIED**.

This checkpoint supersedes informal chat-only progress claims. Only artifacts physically present in this archive count as reproducible evidence.

## Established / preserved

- Central electric-hybrid single-action radial K/G scan: 7/7 finite-L gates PASS in the packaged production audit.
- Compact interface control-space rank: 6/6 outer and 6/6 inner in the packaged handover audit.
- The current experimental angular evaluator is action-derived/high-L and uses the current unreduced-action convention rather than the superseded printed-shortcut oracle.

## New angular work in this checkpoint

A functional 42-parameter lower-order/mass-slot continuation over
`(c4,c5,d3,e4,v13,v5,a9)` found a **numerical slot-space angular existence witness** with all five monitored angular quantities positive:

- min `cOmega_minus^2` = 0.9147506333592901
- min `cOmega_plus^2` = 54.535957134605965
- min `M1/cond55` = 69.26309766130805
- min `M2/cond56` = 803.1587044795431
- min angular discriminant = 1584.7418797217242

This is **not yet an accepted action-level closure certificate**. It demonstrates that the lower-order angular obstruction is controllable in coefficient space.

The witness reduces to four essential target profiles (`c4,c5,e4,a9`), but the first inversion back to pure Horndeski primitives reveals large correlated primitive excursions and changes additional 41-slot coefficients. Therefore the slot witness must not be relabeled as an action PASS.

## Genuine-SVT A2 family route

The action-level A2 on-shell family was scanned by changing the outer initial value of `v6` and rebuilding the central member through the actual A2/metric/Hessian/G2XX chain.

Materialized members included here:

- `Delta v6(0.61) = -1.5`
- `Delta v6(0.61) = -1.8`

For the `-1.8` member the direct angular audit gave approximately:

- min `cOmega_minus^2` = -118.27836843648466
- min `cOmega_plus^2` = -118.27836843648466
- min `cond55` = -236.5567368729693
- min `cond56` = -820.4690378513294
- min discriminant = -74271.59374866784
- min vector angular branch = 0.045645783004934326
- min alpha7 = 3.0778571121438407

A later `-2.0` scan attempt did not complete within the checkpoint run, so no claim is made for that point here.

## Remaining hard gates before a FINAL release

1. Find an **action-level** angular FULL PASS, most plausibly by continuing the genuine-SVT A2 family and/or adding holonomic background-null mass controls.
2. Replay the full background equations and all seven finite-L radial K/G gates on that same member.
3. Solve and export both joint HSVT handovers with the same action.
4. Re-emit one global Direct-41 stream.
5. Build and certify global `K,R,G,S,M` from that same stream.
6. Run same-operator QNM / spectral convergence and boundary-condition audit.
7. Only then set `ABSOLUTE_FULL_CLOSURE=true` and create a FINAL archive.

## Test note

A fresh full pytest run was started while preparing this checkpoint but exceeded the available packaging time before completion. The previously packaged base reported 118/118 PASS. This checkpoint does **not** claim a new full-suite count for the added experimental scripts.
