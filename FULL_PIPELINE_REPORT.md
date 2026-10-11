# FULL PIPELINE REPORT — SSZ_FULL_CLOSURE

**Ausgeführt:** 2026-10-11T14:52:35+02:00
**git HEAD:** `261ac3a93504087b2742cc579e154c2535c0c0f1`
**Gesamtergebnis: PASS — alle Stages grün**

## Stage: Evidence-Registry (run_all_evidence) — PASS

> ==> angular_m5_symbolic_provenance
>     expected=0 actual=0 met=True
> ==> angular_universal_provenance
>     expected=0 actual=0 met=True
> ==> angular_universal_pipeline
>     expected=0 actual=0 met=True
> data/generated/repo_observability/ALL_EVIDENCE_RUN.json
> data/generated/repo_observability/ALL_EVIDENCE_RUN.md

## Stage: Zertifizierte Pfeiler-Artefakte — PASS

| SSZ_FK_COUPLED_MODEL_V1.json | PASS | sha256 7f3f306bfac87879… | 2,517 B |
| SSZ_WEISZ_RESPONSE_BRIDGE_V1.json | PASS | sha256 3df042a18fd247f6… | 3,699 B |
| SSZ_SPECTRAL_RESPONSE_FRAMEWORK_V1.json | PASS | sha256 87dd328a25df17ac… | 3,054 B |

## Stage: Release-Metadaten + Demos — PASS

| git HEAD | PASS | 261ac3a93504 |
| MODEL_LOCK.json | PASS | 261ac3a93504 |
| EVIDENCE_INDEX.json | PASS | 261ac3a93504 |
| Demo-GIFs | PASS | 01_fk_coupled_gates_ladder.gif, 02_swdiag_resolution_drop.gif, 03_bridge_projections.gif, 04_srf_blind_recovery.gif, 05_chain_overview.gif |

---

## Interpretation (ehrlich)

- Alle grünen Ergebnisse belegen das **implementierte Modell und die
  synthetische Auswertekette** — sie sind KEINE neue astrophysikalische
  Entdeckung. Zertifiziertes Urteil der SW-Diagnostik bleibt
  **PEAKS_NOT_ROBUST**; 56,4 % bzw. comparable Kennzahlen sind
  numerische Kandidaten, keine gesicherte neue Physik.
- Reproduktion: `PYTHONPATH=src python tools/run_full_pipeline.py`
  (Laufzeit ~10–25 min, pytest dominiert). Einzelnachweise in
  `full-report.md` § 4.5 und `REPRODUCE.md`.
