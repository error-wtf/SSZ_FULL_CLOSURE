# Codex implementation brief

Build a production-quality, fully reproducible Python repository from this handoff.

Your primary objective is to close `DIRECT_GLOBAL_KRGM_EXPORT` without changing the frozen P5 geometry or the production action member.

## Required architecture

Create packages for geometry, action patches, Horndeski coefficients, genuine-SVT coefficients, jet differentiation, common unreduced quadratic action, constraint elimination, canonical K/G/S/M extraction, stability tests, QNM boundary conditions, provenance and CLI.

## Required model

The frozen action is:

- the globally patched on-shell Horndeski P5 representative represented by the authoritative action/profile data and constructive patch certificates;
- plus the genuine-SVT term `Delta f2 = 0.01 Y` on the background branch `A0prime=0`.

Do not replace this with the old electric central-SVT representative merely because historical data files exist.

## Required verification workflow

1. Run the provided auditor and reproduce every constructive PASS.
2. Reconstruct every Horndeski patch from action jets, including lower-order holonomic terms.
3. Construct flat C-infinity action partitions in monotone phi, not coefficient interpolation after reduction.
4. Re-solve background-constrained jets on every handover.
5. Add the epsilon-Y SVT contribution at the unreduced action level.
6. Regenerate all 41 coefficient slots.
7. Eliminate common constraints with the supplied generalized-psi maps.
8. Export center-to-infinity K,R,G,S,M for L={6,12,20,42,110,420,1000} plus user-selectable L.
9. Regress pure limits against Kase-Tsujikawa Maxwell-Horndeski and Zhang-Kase SVT formulas.
10. Create `SSZ_P5_DIRECT_GLOBAL_KRGM_CERTIFICATE.json` with `{"pass": true, ...}` only if all strict gates pass.
11. Run the auditor with `--require-direct-krgm`.
12. Only then enable the coupled QNM solver.

## Forbidden shortcuts

- no `Reduce(A)+Reduce(B)` sector splicing;
- no use of the superseded plus sign for v12;
- no old a5 without `-a1''`;
- no silently using old `FINAL_*` filenames as proof;
- no declaring a deep-core physical instability from a rounded coefficient dump without exact-invariant regression;
- no QNM claim from a member that did not pass the direct-global-KRGM gate.

Produce tests, deterministic build artifacts, CI, a machine-readable provenance manifest, and plots/tables corresponding to the monograph.
