# Current Working State — Strong-Field Action Continuation — 2026-09-20

Absolute Full Closure is **not certified**.  This checkpoint intentionally
separates historical coefficient reproduction from physical member closure.

## Current physical architecture

- Historical Central Direct-41 normalization: **scoped PASS only**.
- Electric-Hybrid bulk `0.62<u<0.70`: **production candidate bulk PASS** for
  all seven finite-L kinetic/radial gates.
- Physical right strong-field transition: **OPEN**, solved as one action problem
  on `0.70<=u<=0.715`; the historical `u=0.71` seam is not physical.
- Stable inner light ring: `u ~= 0.7061345733`, inside that transition.
- Historical Inner coefficient member: **rejected as final physical member**.
- Global absolute Direct-41 / KRGSM / same-operator QNM: **not run** until the
  strong-field transition is action-closed.

## New hard findings

1. The unmodified Hybrid extension loses even-sector kinetic positivity near
   `u ~= 0.7013`, before the inner light ring.
2. The electric constitutive pivot remains finite and far from zero in the
   electric branch (`min |chi_E| > 357` in the audited right buffer), so this is
   not merely a multiple-root vector-constraint degeneracy.
3. The existing background-null `G2_XX` lift is an exact kinetic null direction
   for the failing right-transition mode in the current reducer.
4. A standalone `G4_XX` background-null quartic tubular repair is rejected by a
   rank-deficient background compensation pivot.
5. Mode-projected Horndeski primitive sensitivity identifies `a1` as a strong
   positive kinetic lever; it must be realized inside a jointly compensated
   Horndeski+SVT action bundle before promotion.
6. The historical Electric-Hybrid regression mismatch was isolated entirely to
   the confirmed `c6` formula correction.  The legacy snapshot is preserved,
   a corrected reference is versioned separately, and no tolerance was relaxed.

## Solver rules

- action first; coefficients are outputs;
- no historical interior coefficient targets;
- only contract-declared endpoint jet orders may be imposed;
- prefer background-null/background-preserving normal directions;
- continuation must track physical modes and constraint rank;
- finite L is not an analytic high-L proof;
- angular evaluation must first reproduce the exact historical SVT oracle and
  the Maxwell-Horndeski limit;
- final certification requires absolute `Action -> C41` re-emission;
- QNM uses the same final operator with regular-center + outgoing-infinity
  boundary conditions.

## Software status

The repaired checkpoint passes the complete test matrix when run in the two
explicit groups used for release verification:

- Unit: `73/73 PASS`.
- Negative + regression + integration + release: `52/52 PASS`.
- Total grouped matrix: `125/125 PASS`.

A monolithic `pytest -q` invocation hit the execution timeout in the packaging
environment without reporting a failure before timeout; no test is omitted from
the two successful groups above.

## Primary evidence

- `FULL_CLOSURE_WORKING_STATUS_2026-09-20_STRONG_FIELD.json`
- `docs/STRONG_FIELD_TRANSITION_CONTRACT_2026-09-20.md`
- `data/generated/strong_field_transition_2026-09-20/STRONG_FIELD_TRANSITION_AUDIT.json`
- `data/generated/strong_field_transition_2026-09-20/HYBRID_RIGHT_TRANSITION_MODE_SCAN.csv`
- `data/generated/strong_field_transition_2026-09-20/TRANSITION_PRIMITIVE_SENSITIVITY.json`
- `data/generated/absolute_attempt_2026-09-20/ELECTRIC_HYBRID_C6_REFERENCE_MIGRATION.json`
