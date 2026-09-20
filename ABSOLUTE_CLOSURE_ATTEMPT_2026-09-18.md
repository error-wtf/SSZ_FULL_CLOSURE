# SSZ P5 Absolute Closure attempt — 2026-09-18

## Result

**Absolute Full Closure is not yet certified.**  This run executed or audited every
remaining gate that can be reached from the present repository and replaced the
failed regional-electric production choice by an explicit healthy candidate:

`P5 Horndeski global + epsilon_Y Y`, `epsilon_Y=0.01`, `A0prime=0`.

The healthy candidate passes the complete constructive auditor and all seven
finite-L direct tests on the exterior and strong-H carrier.  The one remaining
implementation blocker is the deep punctured Core direct coefficient/reducer
stream.  QNM is therefore correctly kept closed.

## What changed in this run

1. Added `production/full_action_lower.py`: direct Appendix-A `v5,c3,e3` from one
   holonomic `f2+f3+f4` action, including total-action chain rules.
2. Audited the 2026-09-17 regional electric Central member.  `v5` replays, but
   the selected `c3/e3` do not simultaneously replay from the frozen common
   action jets.  That member is not used for an Absolute-Closure claim.
3. Added `production/hsvt_eps_y.py`: direct zero-vector `Delta f2=epsilon_Y Y`
   deformation generated from the pure Horndeski carrier, not copied from a
   prepared HSVT coefficient table.
4. Re-ran the historical full constructive HSVT auditor: 56 PASS + 7 PASS_EXACT,
   zero constructive/numerical failures; only `DIRECT_GLOBAL_KRGM_EXPORT` is OPEN.
5. Direct finite-L profile-reducer audit of the new HSVT candidate:
   exterior PASS and carrier PASS at L=6,12,20,42,110,420,1000.
6. Direct punctured-Core coefficient/reducer audit FAILS at every scanned L.
   This is the already documented deep-core cancellation problem: the exact
   Maxwell-Horndeski Core invariants remain positive (`F,G,H,K_scalar > 0`).
   A raw rounded/core-selected 41 table is therefore not accepted as the physical
   no-ghost verdict.  A direct action-level G4XX/G5 Core emitter is still required.
7. Analytic-center / C-infinity handover remains constructive PASS.
8. Coupled QNM is intentionally BLOCKED until a valid direct-global KRGM
   certificate exists.  No archived test-field/Jost root is promoted.

## Software verification

103 tests collected and all passed when run in bounded groups:

- 57 unit
- 2 regression
- 1 integration
- 42 negative/guard
- 1 release-auditor

The monolithic pytest command exceeds the execution-shell time window; the same
98 tests are fully executed by the grouped commands above.

## Current closure chain

```text
P5 geometry                        PASS constructive
HSVT epsilon_Y action member       PASS constructive
Exterior direct finite-L KRGM      PASS all required L
Strong-H carrier direct finite-L   PASS all required L
Punctured Core exact invariants    PASS constructive/exact theorem
Punctured Core direct 41 -> KRGM   FAIL implementation
Analytic center                    PASS constructive
Global direct KRGM certificate     OPEN (Core blocker)
Coupled QNM convergence            BLOCKED
Absolute Full Closure / v1.0.0     NOT CERTIFIED
```

## Exact next engineering task

Do **not** tune the old Core coefficient dump.  Emit the punctured Core directly
from its covariant `G4XX` and subcore `G5` action tables, with all lower-order
Appendix coefficients and JET9D8 derivatives generated on the common guard grid.
The regenerated Core must reproduce the positive exact invariant theorem and
make the direct finite-L reducer positive before the global certificate or QNM
gate can be opened.

## Follow-up: constraint-preserving Deep-Core descriptor

The Core investigation was pushed one level earlier than the unstable 3-field
`M` block.  `src/ssz_p5/reducer/unreduced_descriptor.py` now builds the complete
profile-aware eight-field Euler operator and a generalized-psi descriptor which
keeps the H0 multiplier and h1 constraint variable explicit.

This removes the numerically dangerous `1/D_h1` step from the Deep Core.  Using
the same JET9D8 derivative service as production, the descriptor reproduces the
independently closed H0 constraint for every required multipole.  The worst
scaled closed-form difference is about `4.84e-9` on the carrier and `1.98e-9`
in the punctured Core; forbidden higher H0 terms are at or below about
`2.84e-14`.  The H0-square identity is zero on the audited candidate.

Accordingly the Core status is refined to:

```text
Core exact principal/kinetic invariants         PASS_EXACT_CONDITIONED
Core generalized-psi descriptor constraint      PASS
Core 3-field explicit Schur M block              SUPERSEDED AS DEEP-CORE PATH
Descriptor -> radial spectral/QNM linearization OPEN
```

Evidence:
`data/generated/absolute_attempt_2026-09-18/CORE_GENERALIZED_PSI_DESCRIPTOR_AUDIT.json`.
