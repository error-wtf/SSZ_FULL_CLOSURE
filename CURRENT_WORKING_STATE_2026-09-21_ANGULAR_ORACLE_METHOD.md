# Current working state — Angular oracle method checkpoint (2026-09-21)

This checkpoint extends the dense strong-field K/R continuation checkpoint without claiming an angular PASS for the mixed H+SVT branch.

## Preserved finite branch result

The previously accepted dense strong-field continuation remains the current principal-sector result: the finite action-derived branch crosses the stable inner light ring with positive Even kinetic matrix and positive radial characteristic speeds on the tested dense grid and finite-L set. That result is preserved unchanged in the repository evidence.

## Angular methodology now frozen

The angular gate is now treated as an independent falsification layer.

1. The genuine U(1)-SVT reference is evaluated in the published Zhang–Kase Eq. (86–95) convention. On the certified witness interval `0.62 <= u < 0.70`, the stored reconstruction has `cV = 1` to machine precision and both coupled angular branches positive. The minimum stored values are approximately `cminus = 20.95426095` and `cplus = 639.41926830`.
2. The historical Sep13 angular dataset is retained as provenance but classified `LEGACY_CONVENTION_MISMATCH`; it is not an oracle for the currently adopted published convention.
3. Applying the SVT-specific simplified formula directly to the mixed H+SVT dense branch is retained only as a diagnostic artifact. Its negative values are NOT a physical angular FAIL for the mixed branch.
4. The mixed H+SVT angular gate remains OPEN until a universal large-L reducer, derived from the common 41-slot action with asymptotic constraint elimination performed before taking the angular limit, reproduces both the published genuine-SVT and Maxwell–Horndeski limits.

## Release semantics

- `DENSE_KINETIC_GATE = PASS` on the previously certified dense interval.
- `DENSE_RADIAL_GATE = PASS` on the previously certified dense interval.
- `MIXED_H_SVT_ANGULAR_GATE = OPEN`.
- `CONTINUATION_0p708_TO_0p715 = HELD` until angular reducer validation.
- `ABSOLUTE_FULL_CLOSURE = false`.

## Reproducible angular artifacts

See `data/generated/angular_oracle_2026-09-21/` and run:

```bash
python tools/audit_angular_oracle_checkpoint.py
```

The audit intentionally distinguishes published-oracle validation, legacy-convention provenance, and mixed-branch diagnostics.

## Universal reducer contract

The mandatory derivation/validation contract is frozen in `docs/ANGULAR_UNIVERSAL_REDUCER_CONTRACT_2026-09-21.md`. In particular, constraints are eliminated order-by-order in large L before the angular limit is taken; continuous mixed-to-SVT and mixed-to-MH epsilon regressions and eigenvector-overlap mode tracking are required before production use.
