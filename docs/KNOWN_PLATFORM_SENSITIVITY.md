# Known platform sensitivity — RESOLVED (2026-10-05)

## Status: RESOLVED
Both previously failing tests now pass on Kali (numpy 2.4.6 system BLAS)
AND on GitHub CI (ubuntu, numpy 2.5.3). The tolerances were relaxed with
measured justification; the frozen-stream identity remains enforced by
the exact member-hash gate, which is platform-independent.

## Affected tests and frozen tolerances
1. `tests/regression/test_g05_member_roundtrip.py`
   - old: `rtol=0, atol=1e-10` (held only on the CI BLAS)
   - new: `rtol=5e-4, atol=1e-8`
   - measured cross-platform drift: worst 1.66e-4 relative on
     `f2phiphi` (BLAS summation-order sensitivity)
   - identity gate: `action_member_sha256` hash check unchanged
2. `tests/unit/test_electric_hybrid_controls.py`
   - old: replay max err `< 5e-8`
   - new: `< 5e-7`; the x-grid stays pinned at `atol=1e-13`
   - measured: 3.2e-7 cross-platform

## Why this is not tolerance-shopping
- The G05 test's purpose is "no silent member switching". That identity
  is enforced by the SHA-256 member hash (unchanged, exact).
- The numerical gate now catches formula-level changes (which move
  values by >> 1e-3) while tolerating BLAS dispatch noise (<< 1e-3).
- Both relaxations are annotated in the test bodies with measured numbers.
