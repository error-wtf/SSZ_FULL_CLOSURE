# FULL PIPELINE REPORT — 2026-10-11 (Dokumentations-Lauf)

**Repo:** `/home/error/physics/clones/SSZ_FULL_CLOSURE`
**Branch:** `spectroscopy-real-data-20261004` @ HEAD `5245a0fe5ca15e864bf1dc95021437622cd495be`
**Ausführend:** lokaler Dokumentations-Agent (Subagent), Arbeitsauftrag vom 2026-10-11.
**Disziplin:** KEIN push, KEIN commit — nur Working-Tree-Änderungen; Owner reviewt.

---

## 1. Auftrag

| # | Aufgabe | Status |
|---|---------|--------|
| A | Repo-Inventur (README, docs/, artifacts/, tools/, tests/) | ✅ erledigt |
| B | Full Pipeline REAL ausführen (venv, `OMP_NUM_THREADS=2`) | ✅ 2 Läufe, Final: **GESAMT PASS** |
| C | `figures/` + `tools/make_figures_v1.py` (PNGs + GIFs aus echten Artefakten) | ✅ 6 PNG + 3 GIF |
| D | `demos/` — 4 interaktive Demo-Skripte | ✅ 4 Skripte, alle laufen |
| E | README-Überarbeitung (Kernergebnisse, 6a-Abschnitt, Figuren/Demos, Peer-Review) | ✅ erledigt |
| F | Dieser Bericht | ✅ erledigt |

## 2. Umgebung

- Python: `/home/error/physics/venv/bin/python` — **3.14.7**
- matplotlib 3.11.2, numpy 2.5.3, scipy 1.18.1
- `OMP_NUM_THREADS=2`, `PYTHONPATH=src` für alle Läufe
- OS: Linux (Kali), Kernel 7.1.5+kali-amd64
- ffmpeg vorhanden (GIFs wurden mit PillowWriter gerendert)

## 3. Pipeline-Läufe (aufgabe B)

### 3.1 Lauf 1 — GESAMT: **FAIL** (erwarteter, aufgeklärter Fehler)

```
PYTHONPATH=src OMP_NUM_THREADS=2 /home/error/physics/venv/bin/python \
    tools/run_full_pipeline.py
```

- Stage 1 (pytest): **317 passed, 1 failed** in 164.94 s
- FAIL: `tests/unit/test_import_provenance.py::test_release_metadata_current`
- Ursache: `MODEL_LOCK.json` und `EVIDENCE_INDEX.json` trugen noch
  `git_commit: 261ac3a…` ( alter Release-Stand), HEAD ist aber `5245a0f…`.
  Der Owner hat beim 6a-Commit die Release-Metadaten nicht aktualisiert.
  Der Test ist als Provenanz-Preflight genau dafür da — er hat korrekt versagt.
- Stage 4 bricht daraufhin ab → GESAMT FAIL. Alle anderen Stages/Tests grün.

### 3.2 Fix (legitimes Repo-Werkzeug, Working-Tree-only)

```
PYTHONPATH=src /home/error/physics/venv/bin/python tools/refresh_release_metadata.py
# -> {'MODEL_LOCK.json': '5245a0fe…', 'EVIDENCE_INDEX.json': '5245a0fe…'}
```

`tools/refresh_release_metadata.py` ist das vom Repo selbst vorgesehene
Release-Stamping-Werkzeug (Docstring: „Run BEFORE committing a release").
Verifiziert: `pytest tests/unit/test_import_provenance.py` → **4 passed**.

### 3.3 Lauf 2 — GESAMT: **PASS** ✅

```
PYTHONPATH=src OMP_NUM_THREADS=2 /home/error/physics/venv/bin/python \
    tools/run_full_pipeline.py
```

- **Stage 1 pytest: 318 passed, 1 skipped, 5 warnings in 164.94 s** — PASS
  - Skip: `tests/unit/test_bridge_catalog_schemas.py:22` — „smoke catalog is a
    placeholder string, not a document" (dokumentiertes Platzhalter-Skip, kein Fehler)
- **Stage 2 Evidence-Registry (run_all_evidence): PASS**
  - `angular_m5_symbolic_provenance` expected=0 actual=0 met=True
  - `angular_universal_provenance` expected=0 actual=0 met=True
  - `angular_universal_pipeline` expected=0 actual=0 met=True
  - Output: `data/generated/repo_observability/ALL_EVIDENCE_RUN.json` + `.md`
- **Stage 3 Zertifizierte Pfeiler-Artefakte: PASS**

  | Artefakt | Status | Größe |
  |---|---|---|
  | `artifacts/SSZ_FK_COUPLED_MODEL_V1.json` | PASS | 2,517 B |
  | `artifacts/SSZ_WEISZ_RESPONSE_BRIDGE_V1.json` | PASS | 3,699 B |
  | `artifacts/SSZ_SPECTRAL_RESPONSE_FRAMEWORK_V1.json` | PASS | 3,054 B |

- **Stage 4 Release-Metadaten + Demos: PASS**
  - git HEAD = MODEL_LOCK = EVIDENCE_INDEX = `5245a0fe5ca1` ✅
  - Demo-GIFs: `01_fk_coupled_gates_ladder.gif` … `05_chain_overview.gif` (5 Stück, via `demos/animate_core_results.py`)
- Report: `FULL_PIPELINE_REPORT.md` (automatisch geschrieben), **GESAMT: PASS**
- Gesamtdauer: ~9 min (pytest dominiert)
- Log: `/root/.hermes/cache/scratch/pipeline_20261011_run2.log` (außerhalb des Repos)
- Interpretations-Zeile des Repo-Reports (wörtlich übernommen): „Alle grünen
  Ergebnisse belegen das implementierte Modell und die synthetische
  Auswertekette — sie sind KEINE neue astrophysikalische Entdeckung."

## 4. Figuren (Aufgabe C)

Erzeugt von `tools/make_figures_v1.py` (478 Zeilen) am 2026-10-11 mit dem
venv-Python. **Alle geplotteten Zahlen werden zur Laufzeit aus `artifacts/*.json`
bzw. der zertifizierten NPZ geladen — nichts ist synthetisch oder hardcodiert.**

| Datei | sha256 | Datenquelle(n) |
|---|---|---|
| `figures/F1_geometry_metric_functions.png` | `16fc04071f48b095924af123c264c98faadc79f61ca494eb3da28b999969fbe8` | FK-Modell: f(r), h(r), φ(r)-Profil aus `artifacts/SSZ_FK_COUPLED_MODEL_V1.json` |
| `figures/F2_fk_coupled_gates.png` | `b7702f0cece7596bd88cae21cbd47ee4f05709e2451e500af9542e205e0e341d` | FK: Gate-Summary, G-C-Stabilitätsleiter (6 rungs, min eig > 0), G-E-Residuum aus `SSZ_FK_COUPLED_MODEL_V1.json` |
| `figures/F3_proper_length_probe_6a.png` | `4eafa7e5bbdf51e66f0d0ec5a41ba1127da6bac00df6b9256a718e5e4f967318` | 6a-Prüfstand: 4 Varianten × 3 Matchings, Residuum 0.5642 → 0.0043, aus `artifacts/SSZ_ACTION_DERIVED_MATTER_PROBE_V1.json` |
| `figures/F4_swdiag_resolution_drop.png` | `2d40e8edf8d4355fcc5d59176964728811fcebf4fed2744b02ef5860ad9d1793` | SW-Diagnostik: maxG 0.05696 (n=600) → 0.03417 (n=1000), 4 Grids (Dirichlet/Robin), aus `artifacts/SSZ_SPECTRAL_WEIGHT_DIAGNOSTICS_GW_V1.json` |
| `figures/F5_weisz_bridge_projections.png` | `21027101a726bded427d69e1c8bd1afb064342eaedea58754be3ef2913b9ee70` | Bridge Run 6: Projektionen A/B/C + FK-Referenz aus `artifacts/SSZ_WEISZ_RESPONSE_BRIDGE_V1.json` |
| `figures/F6_srf_blind_recovery.png` | `b374671e93e6e5e3952042b6183def180ebd9330c1eaf8e2bbe8b595be3772de` | SRF: injizierte vs. blind rekonstruierte Moden (3/3, Fehler < 1 bin) aus `artifacts/SSZ_SPECTRAL_RESPONSE_FRAMEWORK_V1.json` |
| `figures/A1_fk_stability_ladder.gif` | `a84cead98e87a8db9cdf2f5435c7aa314995245781ec52f5cbfcd36f80f54ebd` | FK G-C-Leiter, animiert rung für rung (matplotlib FuncAnimation → PillowWriter) |
| `figures/A2_proper_length_6a_sweep.gif` | `27f95a8dab10494057065f11bd16909045bf2817bef1eac332f1728ee584fb83` | 6a: Residuums-Abfall über die 12 Kombinationen, animiert |
| `figures/A3_srf_blind_recovery.gif` | `be069a0020c78f250f6be94771096b8018c2b651241425a3636a013c66eaa9ae` | SRF: Spektrum mit Injektion + Recovery-Frame, animiert |

Reproduktion: `OMP_NUM_THREADS=2 /home/error/physics/venv/bin/python tools/make_figures_v1.py`
(~2 min inkl. GIF-Encoding).

**Technischer Vorfall (dokumentationspflichtig):** Erster Lauf erzeugte nur die
6 PNGs; `A1_fk_stability_ladder.gif` schlug mit `IndexError` fehl, dessen echte
Ursache vom `FuncAnimation.finish()`-Pfad maskiert wurde: `ax.lines.clear()`
existiert in matplotlib 3.11.2 nicht. Fix: `upd(i)` auf das robuste
`ax.clear()`-Muster umgebaut (pro Frame Achsen neu konfigurieren), wie es A2/A3
bereits verwendeten. Keine Daten-/Logikänderung — nur Achsen-Handling.
Danach alle 9 Figuren erfolgreich.

## 5. Demos (Aufgabe D)

Alle Demos lesen live aus den Artefakten und laufen in < 1 s (nur I/O):

| Datei | Zeilen | sha256 | zeigt |
|---|---|---|---|
| `demos/demo_1_ge_not_collapsible.py` | 80 | `fe1657ad2d88ba81068e380a710ad959fb7f0555103fd6ae2ee19cb3567b9f66` | G-E nicht kollabierbar: bestmögliches uniform-α Residuum 0.5642 (56 %), 6a-Prüfstand → NEW_STRUCTURE_NOT_COLLAPSIBLE |
| `demos/demo_2_proper_length_probe.py` | 86 | `00c15f1ae88e2f0d9d2f6f28342b5e8abf458647c29a329e876e785dc58ebe1c` | Proper-Length φ=Q·l: Residuum 0.5642→0.0043 (Faktor ~131), Urteil bleibt NEW_STRUCTURE_NOT_COLLAPSIBLE 12/12 |
| `demos/demo_3_peaks_not_robust.py` | 67 | `493da21417c41e0c34b68343f04c556c7b421387e217c0c163b334134ec51095` | Ehrliches Negativ: maxG 0.05696→0.03417 mit Auflösung → PEAKS_NOT_ROBUST |
| `demos/demo_4_recovery_and_selection_rule.py` | 76 | `d8738e93943d0efd6cef34e08db5440657a7864e71bada122f04806e934a12a8` | Kalibrierte Maschine: blind Recovery 3/3 + Null-Test (median max-z = 0.0) + exakte Selektionsregel B = 0 |

Verifiziert am 2026-10-11: alle 4 Skripte exit 0 (je ~0.02 s).
Bestand: `demos/animate_core_results.py` (5 animierte GIFs, `demos/output/`)
bleibt unverändert und wird von der Pipeline (Stage 4) geprüft.

## 6.README-Überarbeitung (Aufgabe E)

`README.md` komplett überarbeitet (gleiche Datei, Working-Tree-Änderung):
- 6a Proper-Length-Prüfstand-Abschnitt ergänzt (fehlte bisher komplett):
  0.5642 → 0.0043, 4 Varianten × 3 Matchings, P0/P1 bit-genau, P2 = FAIL,
  Urteil NEW_STRUCTURE_NOT_COLLAPSIBLE 12/12
- Alle 9 Figuren mit relativen Markdown-Links eingebettet + Captions
  (Daten, Aussage, Artefakt-Quelle); GIF-Dateipfade zusätzlich angegeben
- Kern-Metriken jeweils mit Zahl + Repo-Pfad belegt
- Reproduktions-Anleitung, Demo-Liste, Peer-Review-Abschnitt
- Keine Code-/Artefakt-Änderungen dafür

## 7. Bekannte Fehler, Lücken, Einschränkungen (ehrlich)

1. **Lauf 1 FAIL (aufgeklärt):** stale Release-Metadaten (261ac3a vs. 5245a0f)
   — vom Provenanz-Test korrekt abgefangen; gefixt mit dem Repo-eigenen
   `tools/refresh_release_metadata.py`. Der Fix ist eine Working-Tree-Änderung
   an `MODEL_LOCK.json`/`EVIDENCE_INDEX.json`, die der Owner mitcommitten muss.
2. **1 Skip in pytest:** Bridge-Katalog-Smoke-Test (Platzhalter-String statt
   Dokument) — dokumentiertes, gewolltes Skip.
3. **6a P2 = FAIL:** der Prüfstand selbst weist P2 als dokumentierten FAIL aus
   (Teil des Artefakts, kein Pipeline-Fehler). Urteil NEW_STRUCTURE_NOT_COLLAPSIBLE
   bleibt über alle 12 Kombinationen stabil.
4. **SW-Diagnostik: PEAKS_NOT_ROBUST:** Resonanzkandidaten fallen mit der
   Grid-Auflösung (0.057 → 0.034). Das ist das zertifizierte, ehrliche Ergebnis —
   KEIN GW-Resonanz-Claim. Alle Pipeline-Pässe belegen die Software-Kette,
   nicht neue Astrophysik.
5. **GIF-Vorfall A1** (matplotlib-API, siehe § 4) — behoben, nur Achsen-Handling
   betroffen.
6. **Keine Push-Disziplin eingehalten:** keine Commits, keine Pushes; git status
   zeigt alle neuen Dateien als untracked, HEAD unverändert bei `5245a0f`.
7. **Scope:** alles ist Effective-Model-Verifikation (Software-Charakter);
   „Full Closure" im physikalischen Sinne ist damit NICHT beansprucht.

## 8. Peer-Review-Leitfaden (was wäre zu prüfen)

1. **Metadaten-Fix nachvollziehen:** `git diff MODEL_LOCK.json EVIDENCE_INDEX.json`
   — es darf sich nur `git_commit` (261ac3a… → 5245a0f…) geändert haben.
2. **Figuren gegen Artefakte prüfen:** jede Zahl in `figures/` muss in der
   genannten Artefakt-JSON auffindbar sein (z. B. maxG-Werte in
   `SSZ_SPECTRAL_WEIGHT_DIAGNOSTICS_GW_V1.json`, Residuen in
   `SSZ_ACTION_DERIVED_MATTER_PROBE_V1.json`).
3. **Demos als schnelle Review-Hooks:** jede Demo printet ihre Zahlen live aus
   den Artefakten — Abweichungs-Kontrolle in < 1 s pro Demo.
4. **Pipeline selbst laufen lassen:** `REPRODUCE.md` folgen (Laufzeit ~10–25 min);
   erwartetes Ergebnis 318 passed / 1 skipped / GESAMT PASS auf `5245a0f`.
5. **Schwächste Stellen des Repros:** (a) die 3 Evidence-Registry-Checks sind
   notwendige, nicht hinreichende Kriterien; (b) SW-Diag „PEAKS_NOT_ROBUST"
   limitiert jede spekulative Interpretation der Spektren; (c) 6a-P2-FAIL ist
   offen und im Artefakt dokumentiert.
6. **Freeze prüfen:** `git status -sb` — HEAD muss `5245a0f` bleiben, keine
   neuen Commits.

## 9. Verwandte Dateien (alle Working-Tree-Änderungen des 2026-10-11)

```
 M EVIDENCE_INDEX.json            (git_commit-Stempel 5245a0f)
 M FULL_PIPELINE_REPORT.md        (Lauf-2-Output der Pipeline)
 M MODEL_LOCK.json                (git_commit-Stempel 5245a0f)
 M README.md                      (Überarbeitung, Aufgabe E)
?? demos/demo_1_ge_not_collapsible.py
?? demos/demo_2_proper_length_probe.py
?? demos/demo_3_peaks_not_robust.py
?? demos/demo_4_recovery_and_selection_rule.py
?? figures/                       (9 Dateien: 6 PNG + 3 GIF)
?? tools/make_figures_v1.py
?? docs/FULL_PIPELINE_REPORT_2026-10-11.md   (dieser Bericht)
```
