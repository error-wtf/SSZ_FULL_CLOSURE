# Working-snapshot integration — 2026-09-18

## Executed checks

- Software suite: **83 passed in 30.21 seconds** in the active repository.
- Ruff: package, tests and both new CLI tools pass.
- Import provenance: **63/63** original hashes verified.
- MH primitive emitter: **3,299 rows × 41 slots**, finite and bitwise equal to
  the archived executable on identical input.
- Stored CSV comparison: maximum scaled discrepancy **8.780694998334242e-05**
  in `a5`; this is a separate diagnostic, not covered by the bitwise port claim.
- Raw Hessian map: three tests cover emitter Jacobian responses, inverse
  roundtrip, and unreachable targets at the zero-electric endpoint.
- Inner target construction: two tests cover one-sided frozen-central jets and
  partitioned endpoint matching. No production CSV is overwritten.

## What was implemented

`src/ssz_p5/production/principal_controls.py` implements the raw triangular
`f2FF/f2XF/f2XX -> v1/v4/c2` map, re-emits the coefficients and rejects impossible
zero-response requests. It explicitly does not certify full hybrid assembly or
holonomic/background-null completion.

`src/ssz_p5/production/inner_targets.py` builds lower-order targets from the frozen
central sources and existing `S_SVT,T_H`. The targets and source hashes are in
`data/generated/inner_controls/`. The central coefficient export remains frozen.

`src/ssz_p5/coefficients/mh_general_primitives.py` is the imported primitive-level
emitter connected to the existing JET9D8 derivative service. The byte-original
implementation and working experiments are retained under `archive/working_2026-09-18`.
Archived scripts with `/mnt/data` paths remain historical records, not portable CLI tools.

## Exact remaining implementation boundary

The imported `solve_central_A2_member_restored.py` solves the response in
`v1/v4/c2`, then assigns archived `v5/c3/e3` values. It does not implement the
lower-order action-control inverse required by the Inner handover contract.
A coefficient-level target is not proof that a fixed-background action realizes it.

Accordingly, the existing Inner interface certificate remains FAIL and global
Direct-41, KRGSM, same-operator coupled QNM and Absolute Full Closure remain
uncertified. No tolerances, production selection or strict rejection gates have
been relaxed. The known central kinetic result is retained, not recomputed as
part of this integration. A later strict release verification executes the
existing required checks and records its actual exit code.
