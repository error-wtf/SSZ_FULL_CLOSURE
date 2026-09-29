# Radial spectral-selection probe

This probe implements the pre-registered observable

```text
Z_n(r) = psi_n(r)^dagger K(r) psi_n(r)
pi_Z(r) = argsort_n[-Z_n(r)]
```

after global K-normalization of the exported eigenfunctions.

It is deliberately fail-closed.  It returns `NOT_YET_EVALUABLE` until both
the hash-bound `DIRECT_GLOBAL_KRGM_CERTIFICATE.json` and a matching
`GLOBAL_CANONICAL_KRGM_SPECTRAL_EXPORT.npz` exist.  It never substitutes the
historical principal KRG stream or the local strong-field production member for
the missing global coupled eigenoperator.

When the physical export exists, the only terminal scientific outcomes are:

- `SPECTRAL_SELECTION_PASS`: at least one pair of local mode weights robustly inverts ordering with radius.
- `SPECTRAL_SELECTION_NULL`: no pairwise ordering inversion is detected on the exported radial grid.

The probe does not claim that synthetic unit tests constitute SSZ evidence.
