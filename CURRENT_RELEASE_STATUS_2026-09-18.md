# Current release status — 2026-09-18

**Absolute Full Closure is not yet certified.**  The current best production
candidate is now the explicit zero-vector HSVT member
`P5 Horndeski global + epsilon_Y Y`, `epsilon_Y=0.01`.

## Software

All **103/103 tests PASS** when executed in bounded groups (52 unit, 2 regression,
1 integration, 42 negative/guard, 1 release-auditor).

## Production-member decision

The 2026-09-17 regional electric member is retained as a local research witness,
but is rejected for an Absolute-Closure claim because:

- the direct finite-L Central kinetic gate fails for L=6,12,20,42; and
- direct Appendix-A replay from one holonomic `f2+f3+f4` action reproduces `v5`
  but not the frozen selected `c3/e3` pair.

The replacement candidate is `SSZ_P5_HSVT_PRODUCTION_CANDIDATE_2026-09-18.json`.
Its full constructive audit has **56 PASS + 7 PASS_EXACT, zero failures**.

## Direct implementation status of the healthy candidate

- Exterior direct finite-L KRGM: **PASS** at L=6,12,20,42,110,420,1000.
- Strong-H carrier direct finite-L KRGM: **PASS** at all required L.
- Punctured-Core exact Horndeski invariants: **PASS**.
- Analytic-center / C-infinity handover: **PASS constructive**.
- Punctured-Core direct 41-slot -> KRGM: **FAIL implementation**.  Historical
  rounded/selected Core coefficients suffer the known deep-Core cancellation;
  direct reducer eigenvalues cannot override the simultaneously positive exact
  Horndeski invariants.  A covariant G4XX/G5 action emitter is still missing.
- Global direct KRGM certificate: **OPEN, blocked by Core**.
- Coupled QNM: **BLOCKED by design** until the same-operator global certificate exists.
- v1.0.0 / Absolute Full Closure: **not released**.

See `ABSOLUTE_CLOSURE_ATTEMPT_2026-09-18.md` and
`data/generated/absolute_attempt_2026-09-18/HSVT_DIRECT_REGION_AUDIT.json`.

## Generalized-psi Core descriptor checkpoint

The former Core label `FAIL implementation deep-core cancellation` has been
refined.  The direct 3-field Schur-reduced `M` block remains numerically
ill-conditioned, but the constraint itself is no longer an unresolved blocker.
A profile-aware eight-field descriptor now performs only the safe differential
change of variables

```text
H2 = psi - (L a4/a3) h1 - (a1/a3) dphi'
```

and keeps `H0` and `h1` explicit instead of dividing by the deep-Core `D_h1`
pivot.  On both the strong-H carrier and the punctured Core, for
`L=6,12,20,42,110,420,1000`, the descriptor reproduces the closed JET9D8 H0
constraint and cancels the forbidden `dphi''` and `h1'` terms.  The audit is
`CORE_GENERALIZED_PSI_DESCRIPTOR = PASS`.

The remaining implementation task is now narrower: linearize this validated
constraint-preserving descriptor into the same-operator radial spectral/QNM
problem without reintroducing an explicit deep-Core Schur complement.  QNM
remains blocked until that descriptor-to-spectral equivalence is certified.
