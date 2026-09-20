# SSZ P5 — Repaired Full Working Release Status — 2026-09-20

This archive is the repaired and complete 2026-09-20 working repository snapshot.
It is **not** labeled ABSOLUTE_FULL_CLOSURE, because the strict Inner Direct-41 common-action replay remains open.

## Repository verification

- Negative tests: **42/42 PASS**
- Unit tests: **67/67 PASS**
- Regression + integration + release-auditor tests: **10/10 PASS**
- Total: **119/119 PASS**
- `tools` import/package defect repaired by adding the repository-local `tools/__init__.py` marker.
- Missing files from the damaged uploaded archive were restored only when their SHA-256 matched the prior 2026-09-20 manifest/baseline source.

## Current regional Direct-41 status

- `WEAK_DIRECT_41 = PASS`
- `OUTER_DIRECT_41 = PASS`
- `CENTRAL_DIRECT_41 = PASS`
- `CORE_DIRECT_41 = PASS`
- `INNER_DIRECT_41 = FAIL`
  - slot normalization: PASS
  - direct common-action replay: incomplete
  - lower-order/principal local response controls: available, but not yet one common holonomic action

Therefore the strict center-to-infinity Direct-Global-KRGM certificate and same-operator QNM certificate are **not yet issued**.

## Scope of this ZIP

The archive contains the complete repaired repository tree, restored infrastructure/tools, generated regional Direct-41 artifacts and certificates, current-session research artifacts, manifest, and SHA256SUMS.
