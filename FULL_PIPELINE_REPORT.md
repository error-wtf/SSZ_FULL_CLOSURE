# FULL PIPELINE REPORT — SSZ_FULL_CLOSURE

**Ausgeführt:** 2026-10-11T17:09:03+02:00
**git HEAD:** `5245a0fe5ca15e864bf1dc95021437622cd495be`
**Gesamtergebnis: PASS — alle Stages grün**

## Stage: pytest Suite — PASS

>     cls = validator_for(schema)
> 
> -- Docs: https://docs.pytest.org/en/stable/how-to/capture-warnings.html
> =========================== short test summary info ============================
> SKIPPED [1] tests/unit/test_bridge_catalog_schemas.py:22: smoke catalog is a placeholder string, not a document
> 318 passed, 1 skipped, 5 warnings in 164.94s (0:02:44)

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

| git HEAD | PASS | 5245a0fe5ca1 |
| MODEL_LOCK.json | PASS | 5245a0fe5ca1 |
| EVIDENCE_INDEX.json | PASS | 5245a0fe5ca1 |
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
