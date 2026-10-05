# Branch status: spectroscopy-real-data-20261004 (Stand 2026-10-05)

## Was dieser Branch ist
Der aktuelle Forschungsstand (42+ Commits vor main auf Basis
qnm-global-final-20261002@5c78ea3).  main ist für Spektroskopie-Arbeit
NICHT aktuell; das ZIP d8de7e7e ist eingefrorener Regression-Snapshot.

## Verifikationszustand (alles CI- oder artefaktbelegt)
- pytest: 269/269 PASSED (lokal numpy 2.5.3 venv UND CI ubuntu)
- ruff gate: All checks passed
- release manifest: 1286 files verified
- Reproduktionsvertrag: exakt 2 registrierte Blocker werden reproduziert
  (finite_l.central_regional_kinetic = K_scalar-Pocket;
   absolute_closure = physical_qnm_claim_allowed=false)
- Li-2023 CST-AAH: 6/6 Gates PASS (jc=1937=Paper-Fig1)
- Weisz-1978: Pflichtregression grün
- MODE_CERTIFICATE_V1: 5/7 NOT_EVALUABLE (brauchen Per-Mode-Run), 2/7
  REJECTED aus echter Physik (L=420/1000 Near-Degeneracy)
- REAL_SPECTROSCOPY_BRIDGE_V1: alle 8 Gates implementiert
  (observations/ + spectroscopy/ + compare-Tool + Schemas)
- 55-Hz-GTI-Negativkontrolle: SPURIOUS_GTIL_GAP_ARTIFACT_REJECTED
- Branch-Protection: force-push/delete verboten, required checks
  audit+reproduce

## Die zwei echten wissenschaftlichen Blocker (unverändert)
1. K_scalar-Pocket (x≈0.29–0.32, stencil-/routen-robust negativ):
   Entscheidung über N1–N4 (echtes Background-Re-Solve im luminalen
   Domain) steht aus.  Bis dahin: keine globalen QNM-Claims.
2. Per-Mode-Zertifizierung: A1/A4/A6 brauchen einen detaillierten
   K-Normen+D2-Lauf pro Modus.

## Ehrlichkeitsregeln dieses Branches
- Wir faken nie, wir fitten nie (kein mode-weise Mass-Fitting;
  FrozenMassPrior ist extern + gehasht).
- "Perfekte Mode" = konvergiert + basisinvariant + physikalisch
  admissibel + beobachtbar — NICHT Im(ω)=0.
- Jeder CI-Fail wird entweder behoben oder als registrierter Blocker
  reproduziert; nichts wird weggeworfen oder versteckt.
