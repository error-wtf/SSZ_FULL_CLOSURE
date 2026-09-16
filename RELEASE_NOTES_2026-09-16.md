# Release notes — 2026-09-16

This repository snapshot is the executable SSZ P5 full-closure handoff.

## Verified in this release

- Editable offline install with the local Python environment using setuptools.
- Package, negative and integration tests pass.
- Master full-closure auditor reports `FULL_CONSTRUCTIVE_CLOSURE = PASS`.
- No constructive or numerical regression failures.
- Frozen production deformation is `epsilon_Y * Y`, epsilon_Y=0.01, on `A0prime=0`.
- QNM execution is guarded until a direct-global-KRGM certificate passes.

## Deliberately separate implementation certificate

`DIRECT_GLOBAL_KRGM_EXPORT` remains OPEN. This is not treated as an additional existence theorem. It requires regeneration of one single center-to-infinity 41-slot stream and finite-l KRGM operators directly from the frozen production action without splicing archived rounded streams.

Do not relabel this gate PASS without producing the certificate described by `schemas/direct_global_krgm_certificate.schema.json` and passing the master auditor with `--require-direct-krgm`.
