# Known platform sensitivity of two regression tests (2026-10-05)

## Status
`TRUE_FULL_CLOSURE_PASS` on GitHub Actions CI (ubuntu, python 3.14, numpy 2.5.3,
`reproduce: success`). Locally on Kali (numpy 2.4.6 system / 2.5.3 venv) two
regression tests fail by tiny numeric margins.

## Affected tests
1. `tests/regression/test_g05_member_roundtrip.py::test_member_builder_reproduces_locked_stream`
   - asserts `atol=1e-10` on rebuilt 3809-row action stream
   - local max abs diff (numpy 2.5.3 venv): 4.94e-08 (494x atol)
   - the values agree to ~9 significant digits; the assertion is far tighter
     than cross-BLAS reproducibility can guarantee.
2. `tests/unit/test_electric_hybrid_controls.py::test_electric_hybrid_recipe_replays_corrected_frozen_checkpoint`
   - asserts max err < 5e-8 on a 4000x41 replay
   - local: 3.2e-7 (6.4x tolerance)

## Cause
Dense linear algebra (BLAS/LAPACK dispatch) differs between the CI image and
Kali; numpy 2.4.6 (system) vs 2.5.3 (CI) additionally changes summation order
in a few kernels. The tests were written against the CI platform's exact
bit patterns.

## What this does NOT mean
- Not an SSZ formula error: values agree to 8-9 significant digits.
- Not a member corruption: member hash is platform-independent and gated
  elsewhere.

## Correct remediation (separate work order, not done here)
Relax these two tolerances to cross-platform values (atol 1e-6 for G05,
1e-5 for electric hybrid) WITH justification comments, OR pin the BLAS
backend in CI and locally. Deliberately not done in this pass to avoid
masking a genuine regression behind a tolerance bump: first re-run both on
the CI runner locally (docker with ubuntu+numpy 2.5.3) and confirm the
diff pattern is pure BLAS dispatch.
