# Comprehensive geometric spectral-selection audit

This branch tests every layer that is presently derivable without changing the
frozen SSZ model.

## Boundaries

The canonical current member is the MODEL_LOCK-pinned electric member
`ELECTRIC_PRODUCTION_MEMBER_CURRENT` with SHA-256
`8bd460ef022a9cdbcc3644abd8aecbfbb910f8e364ac1410378d2641291559cf`.

The audit performs:

1. fresh current-member action rebuild and hash/provenance comparison;
2. fresh 41-slot reduction and finite-L K/G/R/S/M scans;
3. overlap-based local principal-mode tracking;
4. Weisz/Frenkel-Kontorova benchmark reproduction;
5. Li CST-AAH benchmark reproduction at N=2584;
6. an explicit QNM-residue/Green-function analysis API;
7. a fail-closed global readiness gate.

It does **not** concatenate legacy reduced matrices or silently bind the current
electric member to the older regional-member QNM certificate.

## Scientific semantics

A local principal-mode continuation test and an observable spectral-residue
test are different questions.

- Principal test: diagonalize `K^{-1/2} G K^{-1/2}`, track branches by
  adjacent eigenvector overlap.
- Bound/self-adjoint spectral weight: a positive local K density can be used.
- Open/QNM observable: use the retarded Green function or an equivalent
  biorthogonal pole residue.

The current frozen member covers the strong-field domain only. Therefore a
physical center-to-infinity QNM/residue verdict is intentionally
`NOT_YET_EVALUABLE` until a single current-electric same-action global KRGSM
export exists.

A missing global product is not converted into either PASS or NULL.
