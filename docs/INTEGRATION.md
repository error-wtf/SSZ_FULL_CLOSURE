# Integration changes and remaining boundary

The complete prepared snapshot was chosen over the smaller prepared package.
Content-hash comparison found no additional scientific data unique to the
smaller package; unique items were operational scaffolding and caches.
Original manifests for both are preserved under `provenance/source_manifests/`.

Repairs: portable legacy module/data paths; shared scalar and cached profile
JET implementations; forwarded console arguments; actual production blacklist;
strict source inventory; nonzero audit exits on numerical/hash failures; stale
report removal; validated certificate guard; clean-install CI and dependency lock.

The cached derivative now solves the transposed Vandermonde system directly
instead of explicitly inverting the matrix. A new nonuniform-grid third-derivative
comparison exposed a boundary roundoff discrepancy; the solve removes it without
relaxing the regression tolerance. The scalar derivative remains byte-numerically
identical to the archived local-polynomial implementation in the golden test.

`ssz_p5.reducer.canonical.reduce_profile` wraps the existing profile-aware
operator. It validates schema, radial order, branch identities, auxiliary identity,
constraint pivots and matrix structure. `stability_diagnostics` uses symmetric
whitening and records radii at the extrema. These APIs do not create a missing
global action patch generator.

## Direct export and QNM

The supplied dedicated certificate ledger still reports direct global export
OPEN. The prepared CLI itself gates D2-D5, and the canonical spec explicitly
prohibits the archived selected global stream as a production replacement.
Targeted filename/certificate inspection found principal global K/R/G and local
outer KRGM diagnostics, but no separately certified global action-derived
K/R/G/S/M product. No new existence verdict is inferred from this implementation
gap. No coupled spectrum or full-closure v1 release is issued.

Future direct certificates must satisfy `schemas/direct_global_krgm_certificate.schema.json`
and the required gate names defined in `ssz_p5.qnm.gate.REQUIRED_GATES`, with all
seven default L values and every output bound to its SHA-256. The validator
checks declared evidence/provenance; mathematical derivation is not certified
merely by checking JSON.
