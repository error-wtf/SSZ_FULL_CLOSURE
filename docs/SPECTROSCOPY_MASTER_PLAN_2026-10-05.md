# MASTER PLAN: Von der Pipeline zur Falsifikation (2026-10-05)

Normative Linie (vereinbart): Pipeline beweisen → einfrieren → Theorie-
Entscheid → zwei Solvers → Prediction-Katalog → GW-blind-vergleich.

## Phase B — Transport Bridge (NEU, 2026-10-05)

Eigenes Repo: error-wtf/SSZ-Transport-Bridge (typisierte Verträge,
zwei read-only Backends, NULL gemeinsame Formeln).  Die gleichen
strukturellen Prüfungen (forward/parity/convergence/inverse) laufen auf
beiden Backends mit demselben Prüf-Code:
  - Sagnac: validated reference (SAGNAC_REFERENCE_CLOSURE_PASS)
  - SSZ: synthetic inversion gegen KNOWN frozen geometry (done first —
    vor echten Astronomiedaten!)
Synthetic-Inversion-Pfad: {ΔΦ, Δt, z} → {f(r), h(r)} — später
gesunde-Operator-Modi.  Semantische Disziplin übernommen: Muster
dürfen mathematisch inspirieren, erst der Operator entscheidet, was
physikalisch eine Mode ist.

## Phase P — Pipeline-Beweis (Stage A härtet)

P1  Negativkontrolle (DONE): NICER 1200120107 / 55 Hz GTI-Kamm-Artefakt
    wird erzeugt und verworfen (test_nicer_gap_artifact_control.py).
P2  Positivkontrollen (NEXT): synthetische Strain-injection mit
    BEKANNTEN Parametern (f, gamma, A) in realistische NICER-GTI-Struktur;
    der Finder muss f/gamma innerhalb der deklarierten Fehlermaße
    REKONSTRUIEREN — sonst ist die Pipeline nicht vertrauenswürdig.
P3  Pipeline-Freeze: analysis_config + Code-SHA + GTI-Liste werden als
    FROZEN_PIPELINE_RECORD eingefroren. Danach keine Änderung ohne neuen
    Record. (Schutz vor unbewusstem SSZ-Anpassen.)

## Phase T — Theorie-Entscheid (parallel, ein Blocker)

S_covariant → P1 → K_scalar über N1–N4 (echtes Re-Solve im luminalen
Domain, docs/N1_N4_LUMINAL_RESOLVE_SPEC.md).
  T-A: Tasche bestätigt  → Member raus/ändern; kein globaler QNM.
  T-B: positiv gefunden   → Implementierungs-/Oracle-Fehler gefunden.
Erst danach: globaler QNM-Schalter.

## Phase S — Zwei unabhängige Produktionssolver (nach T)

S1  Compactified Jost / spectral determinant
S2  Exterior complex scaling (ECS)
Akzeptanz: beide finden denselben komplexen Pol ω_R - i·γ, stabil
unter radialer Auflösung, äußerer Domain, ECS-Winkel, Basiswechsel,
Root-Continuation. → MODE_CERTIFICATE_V2 (globale Achsen).

## Phase C — Erster SSZ-Prediction-Katalog

{ℓ, n, ω_R, ω_I, Q, Z_q, IPR, mode shape} — Frequenz + räumliche
Modenform + Sichtbarkeit + Lokalisation (Li-Diagnostik).  Daraus der
observable-mode catalogue (Weisz-Selektivität): wenige scharfe Linien,
nicht 500 Eigenwerte auf einen Plot.

## Phase G — GW blind comparison (erster direkter physikalischer Angriff)

Ringdown statt NICER: geometry → h(t) = Σ A_n e^{-t/τ} cos(2πf_n t+φ),
f = ω_R/2π, τ = 1/|ω_I|.  GR prediction vs. SSZ prediction:
  * dieselben Daten (GWOSC public strain),
  * dieselbe Noise-Behandlung,
  * dieselben externen Mass-/Spin-Priors,
  * KEINE nachträgliche Frequenzverschiebung, KEINE Mode-Auswahl nach
    Datensichtung, KEINE freie SSZ-Skalierung pro Peak.
Ausgänge:
  A  identisch innerhalb Fehler → Nulltest, keine Trennkraft
  B  SSZ schlechter            → Member/Perturbationsdynamik unter Druck
  C  SSZ besser (≥2 unabhängige Moden, f UND τ, ohne freie Parameter)
     → erstmals wissenschaftlich interessant

## Phase O — Out-of-sample

Ein Ereignis = Entwicklung.  Danach freeze und weitere Ereignisse blind.
Endzustand: SSZ macht eine numerische Spektralvorhersage, die die Natur
annehmen oder töten kann.  Danach ggf. SSZ Spectroscopy Toolkit.
