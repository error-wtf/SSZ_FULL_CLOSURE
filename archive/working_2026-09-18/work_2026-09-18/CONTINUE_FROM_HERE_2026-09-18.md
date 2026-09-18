# SSZ P5 — Continue From Here (2026-09-18)

This bundle preserves the best current state after the 2026-09-18 direct-closure work.

## Reproduced release state

- `python -m pytest -q`: **11/11 PASS**
- `python ssz_p5_full_pipeline.py`: **43/43 PASS**
- constructive full-closure auditor: **PASS**
- `DIRECT_GLOBAL_KRGM_EXPORT`: **OPEN**
- `python ssz_p5_full_closure_auditor.py --require-direct-krgm`: exits nonzero because the dedicated regenerated-action certificate does not yet exist.

## Frozen production target

Use `SSZ_P5_REGIONAL_PRODUCTION_MEMBER_2026-09-17.json` as the production definition:

`weak H -> outer same-action H/SVT -> central exact SVT -> inner same-action SVT/H -> punctured H core -> analytic center`.

Locked conventions remain:

- `v12 = -v6/(2h)`
- holonomic `a5 = a2' - a1'' - (A0' v4/2)' + A0' v5/2`
- derivative service: JET9D8, window 9, degree 8
- assemble unreduced first, eliminate common constraints once, retain antisymmetric `S`.

Do not replace this regional member with the historical 2026-09-16 global `Horndeski + epsilon Y` simplification.

## 2026-09-18 direct work preserved here

`central_direct/` contains the direct Central background/A2 regeneration experiments, including metric-on-shell audits and regenerated 41-slot candidates.

`outer_control/` contains the current Outer same-action/control experiments and sensitivity/optimization code.

`reports/` contains the strict auditor output for this exact bundle state.

## Remaining strict release task

Regenerate one center-to-infinity regional 41-slot action stream directly from the covariant/action-jet sources and reduce it to one canonical finite-L `K,R,G,S,M` stream without relying on archived rounded/split lower-order tables. Only after this reproduces all regional exact/oracle gates should `SSZ_P5_DIRECT_GLOBAL_KRGM_CERTIFICATE.json` be created with `pass=true`.

No such certificate is included in this bundle.
