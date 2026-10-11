# PEER_REVIEW — SSZ_FULL_CLOSURE

Leitfaden für Peer Review: was die Kernbehauptungen sind, welche
Dateien sie tragen, welche Kontrollen existieren, was explizit **nicht**
behauptet wird.

## 1. Kernbehauptungen und ihre tragenden Dateien

| # | Behauptung | tragende Datei(en) |
|---|------------|--------------------|
| K1 | Das Validierungssystem ist intern geschlossen, reproduzierbar, cross-gecheckt, falsifizierbar (`TRUE_FULL_CLOSURE_PASS`) | `artifacts/true_full_closure/physics_dependency_graph.json` (15 Knoten), `tools/evaluate_true_closure.py` — Verdict wird AUSSCHLIESSLICH aus dem Graphen berechnet, nichts von Hand deklariert |
| K2 | Quellefrei: Materie/Licht bewegen sich geodätisch in der einen Geometrie; eine injizierte Kraft würde gemessen | README § Key results 1 (ΔE-Diskriminierung 10⁶×), `tests/scientific/` Transport-Tests |
| K3 | Alle 43 Gates (G00–G140) passen | `docs/TRUE_CLOSURE_EXPLAINED.md` (Formeln + Routen je Gate), `tests/` |
| K4 | FK Coupled V1: 6/6 Gates PASS, G-E `NEW_STRUCTURE_NOT_COLLAPSIBLE` (SSZ-FK ≠ FK × Redshift) | `artifacts/SSZ_FK_COUPLED_MODEL_V1.json`, full-report.md § 4.5 |
| K5 | SW-Diag: `PEAKS_NOT_ROBUST` (ehrliches negatives Urteil, kein Resonanz-Claim) | full-report.md § 4.5 |
| K6 | Weisz Response Bridge Run 6: Stage 1 + Stage 2 PASS | `artifacts/SSZ_WEISZ_RESPONSE_BRIDGE_V1.json` |
| K7 | SRF Run 6 Stage 3: 5/5 Gates PASS, blinde Recovery 3/3, Null median max-z = 0 | `artifacts/SSZ_SPECTRAL_RESPONSE_FRAMEWORK_V1.json` |

## 2. Wo ein Reviewer anfangen sollte

1. `tools/run_full_pipeline.py` lesen (~230 Zeilen, keine Dependencys
   außer stdlib) — er zeigt die gesamte Verifikationslogik.
2. `PYTHONPATH=src python tools/run_full_pipeline.py` ausführen
   (~10–25 min) und `FULL_PIPELINE_REPORT.md` gegen die Erwartungen in
   `REPRODUCE.md` § 5 vergleichen.
3. `tools/evaluate_true_closure.py` + den Dependency-Graphen prüfen:
   Sind die 15 Knoten-Abhängigkeiten zirkelfrei? Trägt jedes Gate seine
   Toleranz und Messzahl im Graphen?
4. Artefakte manuell gegen die README-Zahlen prüfen (Stichprobe):
   z. B. `python -c "import json; a=json.load(open('artifacts/SSZ_FK_COUPLED_MODEL_V1.json')); print({k: g.get('pass') for k,g in a['gates'].items()})"`.

## 3. Eingebaute Kontrollen (negative controls)

* **8-Korruptions-Batterie** (G119): künstliche Viererkkraft,
  Metrik-Perturbation, falsches Christoffel-Vorzeichen, falsche
  Energieformel, verschobener Ring, falsche Amplituden-/Phasen-/
  Rotverschiebungsgesetze — alle 8 werden zuverlässig erkannt
  (README-Tabelle mit clean/corrupted-Werten).
* **2 registrierte Rejections** mit erwartetem Exit 2:
  `regional_kinetic_rejected_member`, `frozen_onshell_rejected_zero_vector`
  (`REPO_EVIDENCE_REGISTRY.json`) — das System lehnt dokumentierte
  schlechte Mitglieder bewusst ab, der Exit-Code ist Teil des Vertrags.
* **SRF-Null-Test**: Detektor-Rauschen ohne Injektion → median max-z = 0
  (false-alarm control, Gate S3_G4).
* **Immutables Mitglied**: ein einziger Member-Hash (`8bd460ef…`)
  überall; `MODEL_LOCK.json` = Manifest = CSV-sha256 (G01).
* **Metadata-Stempel**: `MODEL_LOCK.json`/`EVIDENCE_INDEX.json` müssen
  den aktuellen HEAD tragen (`test_release_metadata_current`) — alte
  Stamps an einem sauberen Baum sind nur als Ancestor erlaubt.

## 4. Bekannte Grenzen / offene Punkte (bewusst nicht versteckt)

* `PEAKS_NOT_ROBUST` bleibt stehen: die Spektralgewicht-Peaks sind
  auflösungsabhängig (maxG 0.05696 → 0.03417 von n=600 → n=1000).
  **Daraus wird keine Resonanz behauptet.**
* G-E `NEW_STRUCTURE_NOT_COLLAPSIBLE` ist ein **Befund**, keine
  Schwäche — aber seine physikalische Interpretation ist offen.
* Appendix B (Segmentgitter) ist **nicht herleitbar** — dokumentiert,
  nicht wegoptimiert.
* Der historische Forschungszweig (`spectroscopy-real-data-20261004`,
  README-Block oben): B6/QNM-Status und die V4-Branch-Arbeit sind dort
  ehrlich als offen markiert; global physikalische QNM-Claims bleiben
  auf diesem Zweig geschlossen.
* B6 (QNM-Randbedingungen) im Bridge-Repo ist OFFEN — Box-Modes erlauben
  KEINE QNM-Claims.

## 5. Was NICHT behauptet wird

1. **Keine neue astrophysikalische Entdeckung.** Alle grünen Ergebnisse
   belegen das implementierte Modell und die synthetische Auswertekette.
2. Keine empirische/observationale Bestätigung von SSZ. Externe Claims:
   keine (registriert).
3. Die 56,4 %-Kennzahl und Vergleichbares = numerischer Kandidat, keine
   gesicherte physikalische Signatur.
4. Materie-Kopplung ist eine dokumentierte **coupling assumption**, kein
   abgeleitetes Gesetz (Alternativen: kinetisch m·s²/D statt m/D²;
   Substratphase φ(r) = Q·ℓ(r) mit ℓ = ∫s dρ).
5. FK-Gitter-Mapping (r_n = R+n·Δr vs. x_n = na+u_n) ist ein zulässiges
   Mehrskalenmodell — **keine** direkte radiale Materie-Gitterinterpretation.
6. `TRUE_FULL_CLOSURE_PASS` gilt für den registrierten frozen corpus;
   er macht keine Aussage über tieferliegende SSZ-Dynamik jenseits
   davon (das ist der nächste Forschungsschritt).

## 6. Reproduktionsvertrag

`REPRODUCE.md` § Full Pipeline (2026-10-11) ist der verbindliche
Ablauf; ein Reviewer braucht nur python ≥ 3.11 + `requirements.lock`
und ~25 min. Jede Zahl im README/Report ist entweder (a) live aus einem
Artefakt gelesen (Demos, Pipeline) oder (b) im Dependency-Graphen mit
Test-Node-IDs verankert. Abweichungen = Bug, bitte melden.
