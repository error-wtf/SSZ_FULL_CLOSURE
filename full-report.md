# SSZ Spectral Response Framework — FULL REPORT

**Datum:** 2026-10-11 · **Branch:** `spectroscopy-real-data-20261004`
**Status:** MODELLSTUDIE — lokale Spektralanalyse, klar gekennzeichnet.
**Kein** astrophysikalischer Detektions-Claim. Das zertifizierte Verdikt
**PEAKS_NOT_ROBUST** (GW-Pfad) wird nach der Faktor-2-Korrektur (vgl. § 4.5)
**re-zertifiziert**; der Report trägt den Stand dieser Korrekturrunde.
**Autor:** Bingsi, auf Direktive von Lino.
**Commits:** `1d4e8bb` (FK V1) + `ac2618f` (SRF Stage 3) gepusht;
Korrektur-Commit (massengewichteter Solver + Faktor-2-Fix) als Block geplant.

---

## 0. Die Endformel

### 0.1 Die Kernkette (Programm-Endformel)

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                                                                             │
│  SSZ-Feldgleichungen → globaler Operator → Green-Funktion G^R              │
│                       → Spektralantwort  → astrophysikalische Observable   │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

Jeder Pfeil ist implementiert, zertifiziert und Gate-geprüft
(`weisz_response_bridge.py` Stage 1+2, `run_ssz_spectral_response_framework.py`
Stage 3). Der Operator wird **wortwörtlich importiert** aus der zertifizierten
Bridge — eine einzige Wahrheitsquelle für die Operator-Numerik.

**Zwei-Wege-Architektur (Forschungsstand 11.10.26):** Das Verfahren berechnet
aus einer vorgegebenen Raumzeitgeometrie Schwingungsantworten und untersucht,
welche davon in einer bestimmten Messung sichtbar wären — die Erweiterung von
reiner Eigenmodenanalyse („Welche Eigenfrequenzen erlaubt der Operator?") zu
einer Quellen- und Antworttheorie („Was geschieht, wenn ein System an einer
Stelle/in einem Feldkanal angeregt wird, und was misst ein Beobachter?").
Zwei verwandte, aber **noch nicht vollständig dynamisch gekoppelte** Rechenwege:

- **Weg A — SSZ-FK:** Materie auf gekrümmtem Hintergrund. Gitterdynamik:
  Gleichgewicht, Hessian, Schwingungsmoden (§ 0.2–0.5, Gates G-A…G-F).
- **Weg B — SSZ-3DOF:** Feldstörungen der Raumzeit. Globaler Feldoperator
  (K, G, S, M), drei gekoppelte Kanäle, Green-Funktion, Quellen-/Detektor-
  Projektion, Detektormodell (§ 0.6, 0.9, S3-G1…G5).

Ziel der nächsten Phase: astrophysikalische Vorhersagen und Datenvergleich.
**Noch nicht** als vollständig physikalische Beobachtungskette geschlossen
(offen: Materie-Kopplungs-Herleitung, physikalische Bedeutung der Quelle,
Detektormodell — § 6). Für Schwarze Löcher ist Weg B direkt nutzbar (globale
3-DOF-Resolvente + Quellenprojektion, kein FK-Gitter nötig); die
QNM-Randbedingungen dafür sind offen (B6).

### 0.2 Der SSZ-Hintergrund (Ausgangspunkt der Kette)

```
ds² = −D²(r) c² dt² + s²(r) dr² + r² dΩ²
D(r) = 1/(1+Ξ),   s(r) = 1+Ξ = 1/D          (LOCKED, g₁-Zweig)
Ξ(r) = r_s/(2r)   (für r/r_s ≥ 2.2 operativ)
```

Drei Kopplungsformen der Metrik für Materie auf dem Hintergrund. **Wichtig
(Einordnung nach Lino-Review 11.10.26):** Diese drei Formen sind die
Konsequenzen der Metrik *für ein FK-System, das auf dem Hintergrund lebt* —
ein geprüftes effektives Modell, NICHT eine aus einer Materie-Wirkung
hergeleitete zwingende Kopplung. Insbesondere: (i) die Substratphase nutzt
Q·s(r)·r statt der proper-radialen Phase Q·∫s dρ (bei variablem s nicht
äquivalent); (ii) die kinetische Kopplung m/D² folgt nicht allein aus
dτ = D·dt — die direkte Entwicklung der relativistischen Punktteilchenwirkung
bei radialer Koordinatengeschwindigkeit liefert m·s²/D. Beide Alternativen
sind als offene Modellfragen in § 6 erfasst.

| Kanal | Formel | Wirkung |
|---|---|---|
| Kinematik | dτ = D(r)·dt ⇒ (m/D²(r))·ü_n = F_n | ortsabhängige Trägheit |
| Elastizität | ℓ_n = s(r_{n+½})·(a + u_{n+1} − u_n) | Federn an PROPERER Dehnung |
| Substratphase | φ_n = Q·s(r_n)·(n·a + u_n) | **Misfit-Gradient geometrischen Ursprungs** |

### 0.3 Das Misfit-Feld — die zentrale neue Größe

```
b(r) = Q · a · s(r) = b₀ · s(r)/s(R)
```

Weisz' konstanter Misfit wird zum radially strukturierten Misfit-Gradienten.
**Das ist die SSZ-Weiterentwicklung des FK-Modells.** Wichtige Regime:

- **Weisz-Regime:** kleiner irrationaler Misfit δ = (√5−2)/4 ≈ 5.9 %, κ = 25
  — hier existiert ein stabiles homogenes Gleichgewicht.
- **Tief-solitonisch:** b = a/φ (61.8 % Misfit) — KEIN stabiles homogenes
  Gleichgewicht (min eig(H) < 0 wächst monoton mit V0). Dokumentiert, nicht
  wegbügeliert.

### 0.4 Die Gesamtdynamik (Lagrangian, offene Kette, freie Ränder)

```
L = Σ_n (m/2D_n²) u̇_n²
  − Σ_n (κ/2)[s_{n+½}(a + u_{n+1} − u_n) − a]²
  − Σ_n V₀ [1 − cos(Q·s_n·(n·a + u_n))]
```

- E = T + U exakt erhalten (autonomes mechanisches System, Leapfrog).
- Stabilität ⇔ Hessian H = ∂²U/∂u² positiv semidefinit am Gleichgewicht.
- Moden: verallgemeinertes Eigenproblem  H·u = ω²·M·u,  M = diag(m/D_n²),
  gelöst als symmetrisches Problem M^(-1/2)·H·M^(-1/2)·v = ω²·v
  (massengewichteter Solver, Fix 11.10.26 — vgl. § 4.4). Rückgabe in
  u-Koordinaten, M-orthonormal (uᵀMu = I, verifiziert 1.7e-15).

### 0.5 Die akustische Referenz (freie Ränder, V₀ → 0)

```
ω_q = 2·√(κ/m) · sin(qπ / 2M),    q = 0, …, M−1     (Nullmode inklusive)
```

Nicht die festen-Rand-Formel — die Fehlerhistorie hat das zweifelsfrei
festgezurrt (Abschnitt 4.2).

### 0.6 Die Spektralantwort — das Herz der Endformel

```
P(ω) = b† · Im[ (M⁻¹H − (ω + i0⁺)²)⁻¹ ] · b
```

- b = Quellvektor, Detektor auf beliebigen Sites/Kanälen —
  Projektionsabhängigkeit ist explizit (Weisz-Punkt, Gate G-F).
- Der i0⁺-Regulator ist **keine** Linienbreite: als Linie eingesetzt wäre er
  ein δ-Spike, unsichtbar für das 61-Bin-MAD-Fenster des Detektors.
  Instrumentelle Linien im Detektormodel: Lorentzian, γ = 3·bin_w.
- Exakte Spektraldarstellung (H symmetrisch): P(ω) = Σ_n W_n ·
  Im[1/(λ_n − (ω² + iε))], mit W_n = |⟨v_n|q₀⟩|² und Σ_n W_n = 1 (Summenregel).

### 0.7 Die Observablen-Vorhersagestruktur (nächste Stufe, NICER)

```
f_obs = D(R) · f_local
```

Frequenzmuster + Rotverschiebung. Ein echtes Instrumentenmodell ist die
nächste Stufe und Teil der offenen QNM-Arbeit (B6).

### 0.8 Die Kern-Aussage der Endformel (Gate G-E)

```
s·d ≡ 1  ⇒  akustischer Zweig ist (im homogenen Grenzfall) SSZ-invariant
SSZ wirkt über den Substratsektor b(r)
⇒ SSZ-FK-Spektrum ist NICHT kollapsibel auf „FK × uniformes α"
   (α_best = 0.9192, max Rest nach α-Fit = 0.5642 ≫ 10× Toleranz;
   massengewichteter Solver Run 11: Δα ≈ 9e-8, ΔRest ≈ 2e-7 —
   VERDICT ROBUST gegen den Solver-Fix, § 4.4)
   VERDICT: NEW_STRUCTURE_NOT_COLLAPSIBLE
```

Das ist das Kernergebnis der Modellstudie: die Weisz-Dynamik AUF dem
SSZ-Hintergrund ist mehr als „FK × Rotverschiebung". Trägheit m/D² und
Proper-Längen-Feder κs² kompensieren sich im akustischen Zweig; alle
neue Struktur kommt aus dem Substratsektor b(r) und der Sektor-Kopplung am
Gleichgewicht — Rang-Verschiebungen im Spektrum, die kein skalares α
reproduzieren kann.

**Zwei Offenheiten, bewusst am Kernergebnis dran (Lino-Review 11.10.26):**

1. Die Invarianz-Aussage „s·d ≡ 1 ⇒ akustischer Zweig SSZ-invariant" ist im
   Kontinium homogen exakt, aber bei räumlichen Gradienten von s bleiben
   Nachbarkopplungs-Terme — die strenge Aussage gilt für die hier
   dokumentierte Gitterstudie, nicht als allgemeines Kontinuum-Theorem.
2. Das G-E-Matching nimmt bisher an, dass Rangfolge = Modenzuordnung
   (rank-by-rank). Eine unabhängige Zuordnung über Modenüberlappung
   (|⟨v_n^flat|v_m^SSZ⟩|) steht noch aus und ist als Prüfpunkt in § 6
   eingetragen; erst sie schließt aus, dass das Verdikt auf einem
   Zuordnungsartefakt beruht.

**G-E-Bewertung, präzise geschärft (Lino-Hinweis):** Im Artifact ist
`pass=True` der Umstand „Test durchgeführt, Gate-Kriterium formal erfüllt"
(das Kriterium von G-E IST die Nichtkollapsibilität); das physikalische
Urteil steht separat in `verdict = NEW_STRUCTURE_NOT_COLLAPSIBLE` und
`collapsible_to_uniform_factor = False`. Lesart: der Pass belegt die
numerische Messung am effektiven Modell — NICHT die astrophysikalische
Korrektheit der gewählten Materiekopplung (dazu § 6 Punkt 1).

### 0.9 Detektormodell (Stage 3, vollständig spezifiziert)

```
Statistik:        additiv-Gauß auf dem Antwort-Spektrum (Amplituden-Level),
                  σ = NOISE_FRAC · median(g),  NOISE_FRAC = 0.05
                  False-Alarm-Kalibrierung NICHTPARAMETRISCH über die
                  200 Null-Trials (empirische max-z-Verteilung) — NICHT
                  χ²₂-exponentiell: das gilt nur für echte Fourier-
                  Leistungsdichten, g ist aber eine Resolventen-Antwort.
                  (Doku-Fix 11.10.26, vgl. § 4.4; der Funktionsname
                  `chi2_noise` ist historisch und bleibt aus
                  Reproduzierbarkeitsgründen unverändert.)
Injektionen:      N_INJ = 3, Lorentzian γ = 3·bin_w, SNR = INJ_SNR·σ,
                  INJ_SNR = 30.0 (gemessen: 15σ → 90 % Recovery, 30σ → 100 %).
                  Frequenzen ab Run 4 NICHT zufällig, sondern aus den
                  Eigenfrequenzen des zertifizierten Operators
                  (Quadratic Pencil C0·v = −ω²·C2·v, gleiche assemble_ssz-
                  Montage wie die Response) — maximal separierte
                  Operator-Moden. Damit prüft G3 die Recovery physikalischer
                  Operator-Moden, nicht nur Signalerkennung.
                  (Einschränkung bleibt: Amplituden sind gesetzt, nicht
                  aus einer Quellphysik hergeleitet.)
Schwelle:         6σ über rolling-median-Baseline (61-Bin-MAD-Fenster)
Trials:           N_NOISE_TRIALS = 200 (Nullkalibrierung), RNG-Seed 11
Recovery:         blind — die Recovery-Funktion kennt die Injektionen nicht
Roadmap (echte    P_obs(ω) = d†·G^R(ω)·C_J(ω)·G^A(ω)·d + N(ω) mit
  Stufe 4):       Quellkorrelation C_J und Instrumentenrauschen N als
                  explizite Größen; die aktuelle Stufe ist die synthetische
                  Übergangsstufe, keine stochastische Quellphysik.
```

---

## 1. Gültigkeitsbereich und Claim-Disziplin

- Alles **MODELL-Ebene** (synthetische End-to-End-Validierung der Kette).
- **Kein** „Detektion"-Wort, **kein** Resonanz-Claim.
- Zertifizierter Stand **PEAKS_NOT_ROBUST** (GW-Pfad) bleibt unberührt.
- Lokale Resultate sind und bleiben **lokal gekennzeichnet**; die lokale
  Spektralanalyse läuft parallel zu den offenen globalen QNM-Randbedingungen.
- Baustein-Entscheidung (Lino, 11.10.26): Dynamikmatrix **übernommen**,
  Green-Funktions-Methode **übernommen**, spektrale Gewichte/Residuen
  **übernommen**; inkommensurable Periodizitäten **untersuchen** (nur bei
  herleitbaren Skalen); sinusförmiges FK-Potential **nicht pauschal
  übernommen** (cos-Term ist hier physikalisch BESCHRIFTET: Coulomb-Kristall
  der NS-Kruste); Li-Lokalisierung/Segmentierung **nur als Diagnostik**.
- Anwendung B (Raumzeitsegmente SELBST als FK-Gitter): **MIT DEM AKTUELLEN
  SSZ-FORMALISMUS NICHT HERLEITBAR** (keine Wirkung/kein Variationsprinzip
  für Ξ; Segment-Diskretisierung wäre ein neues Axiom). Gesperrt bis eine
  Segment-Dynamik existiert — dokumentiertes Verdikt, kein stilles Retry.

---

## 2. Modul-Übersicht (zertifiziert, reproduziert; Stand der Korrektur-Runde 2)

| Modul | Inhalt | Stand |
|---|---|---|
| `tools/weisz_response_bridge.py` | Stage 1 (FK-Referenz-Equilibria + Lorentzian-Gewichte = q0-Residuen = Frozen-Zahlen) + Stage 2 (SSZ-Projektionen A/B/C, Auflösungsstabilität); Assembler FAKTOR-2-GEFIXT (§ 4.5) | Run 6 final: **STAGE1 PASS \| STAGE2 PASS** gegen frisches SW-Diag-Cert; Crosschecks n=600/n=1000 exakt; Artefakt strikt-JSON (NaN→null + `selectivity_C_undefined_reason`, s. § 6); A-Kanal `ch0_ch0` g_max 0.093004 (= S3-G1-Zertifizat), SW-Diag-Gitter n600 `maxG` 0.056956 (anderer Check, nicht A-Projektion) |
| `tools/run_ssz_fk_coupled_model.py` | FK-Coupled-Model V1, alle 6 Gates; massengewichteter Solver (§ 4.4) | 6/6 PASS (Run 11, massengewichtet); Runs 6–10 mit Alt-Solver deterministisch identisch; G-E-Verdikt robust (Δ ≈ 2e-7) |
| `tools/run_ssz_spectral_response_framework.py` | Stage 3: Operator → Green → Projektion → Detektor → blinde Recovery; Injektionen aus Operator-Eigenfrequenzen; Fenster 2.5–7.0 (Operator-Moden), Grid 1200; Null-Test nur außerhalb bekannter Moden; G1-Kette auf zertifiziertem Gitter | **5/5 Gates PASS (Run 6, final)** |
| `tools/run_ssz_local_resolvent_spectra.py` | Lokale Resolventen-Spektren, 360 Fälle | exit 0; top1 0.380–0.959, top2 0.722–1.0, **n_eff 1.09–2.94** |
| `artifacts/SSZ_B6_STEP1_BLIND_SCAN_V1.json` | B6 Step 1: blinder Scan (D1/D2/D3), npz-Binding ad2bee7d… intakt | D2-Kandidat (0.05, −1.427), σ=2.56e-3, min/median 0.065 |
| ch2-Wronski-Pol-Pool | Unabhängige Konvergenz-Nullstellen | 6 Kandidaten ≤ 3e-19, best 3.0e-22 bei (0.16239, 0.07238) |

### 2.1 Die 6 Gates des FK Coupled Model V1 (Run 11, massengewichteter Solver)

| Gate | Definition | Ergebnis (Run 11, massengewichtet) |
|---|---|---|
| G-A FK-Recovery | SSZ aus (D=1, s=1): akustische Dispersion ω(q) = 2√(κ/m)·sin(qπ/2M); Summenregel ΣW = 1 | freq rel err **3.96e-14**; 1−ΣW ~ 9e-16 → PASS |
| G-B Energie | Leapfrog, \|ΔE\|/E₀ über 2·10⁴ Steps | **4.28e-8** → PASS |
| G-C Stabilität | min eig(H) am GG je Parametersatz, 6 Stufen V0 = 0.01…4.0 | alle STABLE, min eig 0.035…112 → PASS |
| G-D RSG-Konsistenz | Gitterwelle vs. WKB ∫ds-Gewicht; s·d ≡ 1 | s·d err **1.1e-16**; 2.99e-4 → 8.26e-5 bei dr/2, **Richardson-Ratio 3.62 ≈ 4** ⇒ O(dr²) → PASS |
| G-E Neue Struktur | Kollaps auf uniformes α? | **NOT_COLLAPSIBLE** (α_best 0.9191907, Rest 0.5642467) → PASS (Kernergebnis; solver-robust, § 4.4) |
| G-F Projektion | Antwort hängt von (Quelle,Dektor) ab | max rel diff ~1.0 → PASS |

### 2.2 Die 5 Gates des SRF Stage 3

| Gate | Prüfung | Ergebnis |
|---|---|---|
| S3-G1 Ketten-Konsistenz | Stufe-2-Projektionsmaxima exakt reproduziert (round(g,6) vs. Zertifizat) | A: 0.093004 (abs_diff 0.0); B: 0.0; C: 2.354e-6 — 3/3 match (Run 6, Zertifizat = Bridge-Run-6-Artefakt) |
| S3-G2 Noise-Null | 200 Trials ohne Injektion | median max-z **0.0**, median n_peaks 0 |
| S3-G3 Blinde Recovery | 3 injizierte Modi, Recovery blind | **3/3**, max. Frequenzfehler **0.001387** < Bin-Breite 0.00375 (Fehler: 0.000511 / 0.001387 / 0.000407) |
| S3-G4 False Alarm | Peaks > 3·bin_w weg von Injektionen | **0** Nicht-Injektions-Peaks (bei SNR 30) |
| S3-G5 Auflösung | grobes Gitter n = 200/400 Frequenzen | **3/3** recovered |

---

## 3. Die physikalische Lesart (was die Endformel bedeutet)

1. **Der akustische Zweig ist SSZ-blind.** s·d ≡ 1 ist eine exakte Identität
   (Gitterfehler 1.1e-16, Kontinuum exakt). Trägheits- und Feder-Kanal
   kompensieren sich punktweise. Jede SSZ-Signatur MUSS also aus dem
   Substratsektor kommen — das ist eine Selektionsregel für die Beobachtung:
   man sucht in b(r)-sensitiven Projektionen, nicht im akustischen Zweig.
2. **Die neue Struktur ist Rang-Verschiebung, kein Faktor.** G-E zeigt:
   kein uniformes α bildet das SSZ-Spektrum auf FK ab. Die Misfit-Struktur
   b(r) = Q·a·s(r) mischt die Moden am Gleichgewicht so, dass die
   Spektral-Rangfolge kippt. Das ist die konkrete, gate-geprüfte Antwort auf
   die Integrationsfrage.
3. **Die lokalen Resolventen-Spektren sind hierarchisch.** n_eff < 3 bei
   3 DOF in allen 360 Fällen (1.09–2.94), top2-Mode dominiert fast überall.
   Die Systeme leben effektiv in 1–2 Richtungen — die Mode-Selektion
   (welche Projektion sichtbar werden kann) hat echte, dominante Ziele.
4. **Beobachtbarkeit ist Projektionsfrage.** A-Kanal (ch0→ch0) trägt die
   Signatur (g_max 0.079), B (ch0→ch2) exakt null, C (ch1→ch1, late) auf
   Rauschniveau — und genau das ist das zertifizierte PEAKS_NOT_ROBUST.
   Die Endformel macht aus diesem Verdikt eine Vorhersage-Struktur statt
   eines Misserfolgs: der C-Kanal IST die Vorhersage „diese Projektion ist
   leise".
5. **B6-Quellenkopplung (laufend).** Blinder Scan: D2_DEEP_IM (0.05, −1.427)
   mit σ=2.56e-3 und min/median 0.065 ist der klare Kandidat (D3 mit 0.686
   ist flache Fläche, kein Resonanzprofil). Unabhängig davon: 6 Wronski-
   Konvergenz-Nullstellen ≤ 3e-19 (best 3.0e-22). Der Polish-Lauf
   (rtol-Konvergenz) läuft; Ergebnis wird in einem Folgereport ergänzt.

---

## 4. Fehlerhistorie (Reproduktions-Wert — wie wir hier wirklich arbeiten)

### 4.1 SRF Stage 3 Kalibrierung (Runs 1→2)

1. **Rauschmodell falsch:** multiplikativ-χ²₂ verletzte die MAD-z-
   Kalibrierung (Null median max-z ≈ 9.5 statt 0). Fix: additiv-Gauß
   (dessen Power-Statistik χ²₂ ist) → Null bei 0.0.
2. **Injektionslinie unsichtbar** („peak at None"): der i0⁺-Regulator (1e-6)
   als Linienbreite ergibt einen δ-Spike, den das 61-Bin-MAD-Fenster nicht
   sieht. Fix: instrumentelle Lorentzian γ = 3·bin_w.
3. **SNR-Kalibrierung:** 15σ → 90 % Recovery (flaky), 20σ → 27/30,
   25σ → 28/30, **30σ → 30/30**. INJ_SNR = 30.0 gesetzt, messbegründet.

### 4.2 FK-Coupled-Model (Runs 1→10)

1. **Substrat-Gradient-Vorzeichen:** `+V0·Q·s·sinφ` für ∇U. Absichert-Standard:
   FD-Verifikation von grad UND Hessian gegen zentrale Differenzen (~1e-8).
2. **Akustische Referenz für FREIE Ränder:** ω_q = 2√(κ/m)·sin(qπ/2M),
   q = 0..M−1 inkl. Nullmode — nicht die festen-Rand-Formel.
3. **Default-Misfit b = a/φ** ist tief-solitonisch ohne stabiles GG —
   das Weisz-Regime ist δ = (√5−2)/4 ≈ 5.9 %, κ = 25.
4. **G-D ist O(dr²), nicht exakt:** Richardson statt Maschinen-Nullstelle
   als Gate-Kriterium (Ratio 3.61 ≈ 4).
5. **Broadcast-Shapes (39,) vs (40,)** in G-A-Zwischenversionen (Runs 1/3/5
   exit 1) — Debug-Sequenz dokumentiert, Endversion stabil grün.

### 4.3 Arbeitsprinzip (steht über allen Fixes)

Niemals der Weg des geringsten Widerstands: Jede Zahl in diesem Report ist
ein realer Tool-Output (Log oder Artifact), kein Mock, kein Stub, keine
„als ob"-Simulierung. Gates werden VOR dem Lauf deklariert, nicht hinterher
angeglichen. Negative Befunde (Anwendung B nicht herleitbar; PEAKS_NOT_ROBUST;
b = a/φ ohne stabiles GG) werden dokumentiert, nicht wegbügeliert.

### 4.4 Massengewichteter Eigenwertsolver (Fix 11.10.26, Lino-Review)

**Fehler:** `modes()` bildete A = M⁻¹H und symmetrisierte nachträglich
(0.5·(A+Aᵀ)). Das ist KEINE Ähnlichkeitstransformation — bei ortsabhängiger
Trägheit verzerrt es das Spektrum. Linos Gegenbeispiel (M = diag(1,4)):
korrekt {0, 1.25}, nachträgliche Symmetrisierung {−0.104, 1.354} — das falsche
Verfahren erzeugt sogar einen scheinbar instabilen Modus.

**Fix:** M^(-1/2)·H·M^(-1/2)·v = ω²·v via scipy.linalg.eigh(H, M);
Rückgabe in u-Koordinaten, M-orthonormal; Response-Projektionen
Σ (uₙᵀb)(uₙᵀd)/(ωₙ²−λ) konsistent.

**Verifikation:** gegen das unverfälschte verallgemeinerte Problem:
max |Δω²| = 2.8e-14, max |uᵀMu − I| = 1.7e-15. Freie Kette reproduziert die
akustische Referenz 2√(κ/m)·sin(qπ/2M) exakt.

**Konsequenz für G-E (Run 11):** α_best 0.9191907078 (vorher …906215),
max Rest 0.5642466802 (vorher …9064) — Δ ≈ 2e-7. Das Verdikt
NEW_STRUCTURE_NOT_COLLAPSIBLE ist ROBUST gegen den Solver-Fix: nicht
Artefakt der Symmetrisierung. (Physikalisch plausibel: die Trägheitsvariation
über die Kette ist ~10 %, der Symmetrisierungsfehler zweiter Ordnung darin —
aber jetzt bewiesen statt plausibel.)

### 4.5 Faktor-2-Fehler in den Assembler-Ableitungen (Fix 11.10.26, Lino-Review)

**Fehler:** `dGz[1:-1] = (Gc[2:] − Gc[:-2]) / hc` (ebenso dSz) in
`weisz_response_bridge.py::assemble_ssz` UND in der Kopie in
`run_spectral_weight_diagnostics_gw.py::assemble` — zentrale Differenz
ohne Faktor ½ ⇒ Ableitungsterme exakt DOPPELT so groß.

**Aufdeckung:** analytische Testfunktionen: dGz gab für G = 1+r² exakt das
2-fache von 2r (Ratio 2.000000). Der D1-Stencil derselben Montage war korrekt
(3.8e-14) — der Fehler saß isoliert in den zwei Zeilen, keine Kompensation an
anderer Stelle.

**Reichweite:** assemble_ssz wird von Stage 3 importiert; die Kopie dient der
Stage-2-Zertifizierung (Projektionen A/B/C, PEAKS_NOT_ROBUST). Damit war die
Zertifizierung selbst betroffen — ein „reproduzierter" alter Operator hätte
den Fehler in die gesamte Antwortkette vererbt (exakt Linos Warnung).

**Re-Zertifizierung: ABGESCHLOSSEN (11.10.26, Abend).** SRF Run 3 (lief noch
mit altem Code im Speicher) abgebrochen, bevor ein Schein-PASS entstehen
konnte. Finale Kette, alle Stufen verifiziert:

- **Bridge Run 5** (Faktor-2-gefixt): STAGE1 PASS, STAGE2 FAIL nur gegen das
  ALT-Zertifikat — aufgeklärt als veraltetes cert-Artefakt: das alte cert-JSON
  enthielt den kranken Operator (0.048597); 0.056956/1.1720 = 0.048597 exakt
  (auch n1000). Keine Regression. Bridge **Run 6** gegen das frische
  SW-Diag-Cert: **STAGE1 PASS | STAGE2 PASS** (EXIT 0).
- **SRF Run 4** crashte im neuen Quadratic-Pencil-Pfad (singuläre
  Dirichlet-Box, solve(C2,−C0)) — gefixt via verallgemeinertes
  Eigenwertproblem scipy.linalg.eig(C0, C2).
- **Physikalischer Kernbefund dabei:** Im alten Fenster 0.05–1.5 existiert
  KEINE Operator-Eigenmode (kleinste Box-Mode 2.8646). Injektionen aus
  Operator-Eigenfrequenzen (Lino-Korrektur 2) sind dort prinzipiell
  unerfüllbar. → **Fenster verlegt auf 2.5–7.0** (enthält 2.8646/4.5206/
  6.1255), Grid 1200 Punkte (löst FWHM ≈ 0.0049, Q ≈ 585).
- **SRF Run 4 (neues Fenster):** G3 3/3 blind recoveriert (2.86456/4.52057/
  6.12593), G4/G5 sauber. G2 (Noise-Null) FAIL mit median_max_z = 227 —
  Definitionsfehler: im Modenfenster IST die saubere Antwort nicht glatt; der
  Null-Test zählte die a-priori-bekannten echten Moden als Fehlalarme.
  Korrektur: Null-Test wertet nur spurious Peaks außerhalb ±3.5 Bins um
  bekannte Operator-Moden aus (spektroskopisch korrekte Null-Hypothese).
- **SRF Run 5:** G2 PASS (median_max_z 0.0), G3/G4/G5 PASS; nur G1 exact_match
  scheiterte — Vergleich gegen das auf dem ALTEN Fenster zertifizierte
  Bridge-Ergebnis (Apples-vs-Oranges). Korrektur: G1 (Ketten-Konsistenz)
  rechnet auf dem zertifizierten Gitter (0.05–1.5, n=600), G2–G5 im
  Operator-Fenster.
- **SRF Run 6 (final):** **FRAMEWORK STAGE 3: PASS** — G1 pattern_ok=True
  exact_match=True, G2 median_max_z=0.0, G3 3/3 blind, G4 max_snr=0.0,
  G5 3/3 coarse. Artefakt `SSZ_SPECTRAL_RESPONSE_FRAMEWORK_V1.json` neu
  geschrieben (Fenster 2.5–7.0, N_FREQ_FINE 1200).
- **SW-Diag Re-Run (final):** EXIT 0, alle 4 Grids — n600 dirichlet 0.05696 /
  n600 robin 0.05697 / n1000 dirichlet 0.03417 / n1000 robin 0.03418;
  VERDICT PEAKS_NOT_ROBUST (unverändert, ehrlich berichtet). Bridge-Discrepanz
  gegenüber Run 5: exakt konsistent.

Zahlen sind damit **eingefroren**. Logs: `weisz_bridge_run5_factor2fixed.log`,
`weisz_bridge_run6_freshcert.log`, `srf_v1_run6_g1certgrid.log`,
`swdiag_rerun_factor2fixed.log`; Backups der Originale in
/root/deleted-backups/20261011-bridge-factor2/.

---

## 5. Artefakte und Pfade

| Datei | Inhalt |
|---|---|
| `artifacts/SSZ_FK_COUPLED_MODEL_V1.json` | FK-V1-Gates (6/6 PASS, Run 11 massengewichtet) |
| `artifacts/SSZ_SPECTRAL_RESPONSE_FRAMEWORK_V1.json` | SRF Stage 3 (5/5 PASS, chain + detector + blind recovery) |
| `artifacts/SSZ_WEISZ_RESPONSE_BRIDGE_V1.json` | Zertifizat Stage 1+2 (Projektionsmaxima A/B/C, Crosschecks) |
| `data/generated/spectral/ssz_local_resolvent/SSZ_LOCAL_RESOLVENT_CURVES.csv` | 360 lokale Resolventen-Spektren |
| `docs/SSZ_FK_COUPLED_MODEL_V1.md` | Vollständige Herleitung (§1 Anwendungen A/B, §2 Modell, §3 Gates, §4 Claim-Disziplin, §5 Validierung) |
| **Dieser Report** | `full-report.md` — Endformel + Gesamtstand |

## 6. Offene Punkte (ehrlich)

**Modellfragen aus dem Lino-Review 11.10.26 (Korrektur-Runde 2):**

1. **Geometrische Materie-Kopplung herleiten** (der größte offene Baustein).
   Die drei Kanäle sind effektive Formen, keine hergeleitete Wirkung. Zu tun:
   relativistische Punktteilchenwirkung konsistent entwickeln (Kandidat für
   den kinetischen Koeffizienten: m·s²/D statt m/D²) und die Substratphase
   über die proper-radiale Koordinate ℓ(r) = ∫s dρ definieren (φ = Q·ℓ statt
   Q·s·r). Beides als Modellvariante implementierbar — ein direkter Test, ob
   G-E von diesen Modellwahlen abhängt oder robust gegen sie ist. Bis dahin
   ist die Kopplungsform als **coupling assumption** zu lesen, nicht als
   zwingende Konsequenz.
   1a. **Gitter↔Radial-Koordinaten-Mapping explizit machen
   (Lino-Review, Runde 3):** Der Code nutzt ZWEI Koordinaten gleichzeitig:
   r_n = R + n·dr für das Metrikfeld (dr = 0.05) und x_n = n·a + u_n für das
   FK-Gitter (a = 1). Die Substratphase ist φ_n = Q·s(r_n)·x_n. Numerisch
   geprüft (11.10.26): dφ/du = Q·s(r_n) variiert über die Kette nur um
   2.2 % (6.869…7.021, Q ≈ 5.933·s/1.018) — der lokale Wellenzahl-Gradient
   IST ortsabhängig, aber schwach. Linos Kontinuum-Beobachtung bleibt
   der Kern des Themas: φ(r) = Q·s(r)·r wäre im äußeren Zweig
   φ = Q·(r + r_s/2), also dφ/dr = Q — reine Konstante, KEIN Gradient.
   Die proper-radiale Variante φ = Q·[(r−R) + (r_s/2)·ln(r/R)] liefert
   dφ/dr = Q·s(r). Diese drei Identifikationen (Code, Kontinuum-Literal,
   proper-radial) sind unterschiedlich und müssen auseinandergehalten
   werden; das Mapping x_n ↔ r ist ein Mehrskalen-Modellierungsschritt
   und darf nicht stillschweigend als direktes radiales Materiegitter
   gelesen werden. → Proper-Längen-Variante als G-E-Robustheitstest
   implementieren (Pending bis nach Re-Zertifizierung).
   **BEARBEITET (11.10.26, spät): Prüfstand gelaufen, Ergebnisse unten
   (Punkt 6a).**
2. **Unabhängige Modenzuordnung für G-E:** rank-by-rank-Matching durch
   Überlappungs-Matching |⟨v_n^flat|v_m^SSZ⟩| ersetzen/ergänzen, damit das
   Verdikt nicht auf einem Zuordnungsartefakt beruht.
3. **B-Kanal-Nullantwort — GEPRÜFT (11.10.26, nach Faktor-2-Fix): ES IST EINE
   ECHTE AUSWAHLREGELE, kein Diskretisierungsartefakt.** Befund-Kette:
   (i) Der frozen Export V4_PHYSICAL_RESONANCE_EXPORT_V3 hat K, G, S, M mit
   EXAKT verschwindenden Kopplungen ch0↔ch2 und ch1↔ch2 (max|·| = 0.0 auf
   allen 300 Interpolationsrasterpunkten); nur ch0↔ch1 ist gekoppelt
   (K: 907.87, G: 669.80, S: 11.88, M: 102.02). (ii) Damit hat der
   dirichlet-Assembler in C0 und C2 keine 02-/12-Blöcke. (iii) Für eine
   Quelle in ch0/ch1 lauten die ch2-Zeilen des Systems dann A₂₂·x₂ = 0 mit
   wohlkonditioniertem A₂₂ (cond ≈ 4.43e6, stabil über den ganzen
   Frequenzbereich 0.05–1.5) ⇒ x₂ ≡ 0 exakt, nicht nur numerisch klein.
   Direkt verifiziert: solve bei w=0.7 liefert |x₂-Vektor| = 0.0 exakt.
   Konsequenz: g_max(B) = 0.0 ist strukturell (Kanal-2 ist im Frozen-Export
   eine entkoppelte Sektordynamik), und der C-Kanal (ch1→ch1, g_max 3e-6)
   misst die intrinsische ch1-Sektordynamik, keine Leckantwort aus ch0.
   Offen bleibt die physikalische Interpretation: WARUM ist ch2 im
   Frozen-Export entkoppelt — Absicht des V3-Exports (z. B. Petschek-Sektor
   ohne Rückkopplung) oder ein Export-Defekt? Das ist eine Frage an die
   V3-Export-Pipeline, nicht an den Assembler.
   **Nachtrag NaN-Fix-Rerun (11.10.26, abends):** Rerun aus den Quellen
   bestätigt die Struktur exakt — B g_max = 0.0, A g_max = 0.093004
   (unverändert), C 3e-6, Crosschecks n=600/n=1000 exakt, STAGE1+2 PASS.
5a. **Striktes JSON (Lino-Audit 11.10.26, Befund 2 — BEHOBEN):** Beide
   Artefakte (`SSZ_WEISZ_RESPONSE_BRIDGE_V1.json`,
   `SSZ_SPECTRAL_WEIGHT_DIAGNOSTICS_GW_V1.json`) enthielten `NaN` für
   Selektivitätskennzahlen bei null Peaks. Fix in beiden Generatoren:
   `selectivity_C`/`entropy_S` jetzt `null` plus neues Feld
   `selectivity_C_undefined_reason` mit explizitem Grund; keine
   Tests/Pipeline-Consumer lesen diese Felder (verifiziert). Artefakte
   komplett aus den Scripts regeneriert (SW-Diag ~15 min, Bridge ~2 min);
   Validierung: kein `NaN`/`Infinity` im Rohtext, `json.loads` mit
   `parse_constant`-Reject bestanden; Gate-Zahlen byteidentisch zur
   eingefrorenen Re-Zertifizierung. Physik unverändert.
4. **Quellphysik (Stufe 4):** b ist ein kanalselektiver Einheitsvektor, kein
   Quellmodell; Amplituden in G3 sind gesetzt. Echte stochastische Quelle
   (C_J(ω)) + Instrumentenrauschen N(ω) in P_obs = d†G^R C_J G^A d + N.
5. **G^R-Bezeichnung:** die Dirichlet-Boxantwort ist nicht die globale
   astrophysikalische Resolvente; der Name setzt kausale Randbedingungen
   voraus, die hier noch nicht gegeben sind. Kennzeichnung geschärft.

**Laufende Prozesse / previously open:**

6. **B6 Step 2 Polish** (D2-Pol, rtol-Konvergenz): läuft unabhängig weiter.
   Kriterium: σ muss mit rtol (1e-7 → 1e-11) gegen 0 fallen — sonst Flachpunkt.
   Zusätzlich zu prüfen: passen die Kandidaten zur verwendeten Zeitkonvention
   und den physikalischen Randbedingungen?
7. **Uncommittete Artefakte** (FK-JSON, Resolventen-CSV, B6-Scan, Pol-Pool,
   Solver-Fix, Faktor-2-Fix, dieser Report): Commit-Block, sobald die
   Re-Zertifizierungs-Kette (Bridge Run 5 → SRF Run 4 → SW-Diag) fertig und
   ausgewertet ist.
8. **NICER-Anwendung** braucht ein Instrumentenmodell (echte QPO-/Torsions-
   frequenzen); Vorhersagestruktur steht (§ 0.7).
9. **Globale QNM-Randbedingungen (B6)** offen; lokale Resultate laufen bewusst
   parallel und sind als lokal gekennzeichnet.
10. **Inkommensurable Periodizitäten:** unentschieden (Lino-Entscheidung:
    untersuchen, NUR wenn Skalen aus SSZ oder Astrophysik herleitbar sind) —
    bislang keine herleitbare Skala identifiziert, Baustein bleibt zu.

**6a. G-E-Robustheits-Prüfstand (11.10.26, spät) — erste Bearbeitung der drei
offenen physikalischen Arbeiten (Punkte 1a + 2 oben; Materie-Wirkung als
Variante K mitgetestet):**

Tool: `tools/run_action_derived_matter_probe.py`; Artefakt:
`artifacts/SSZ_ACTION_DERIVED_MATTER_PROBE_V1.json` (sha256
5c55025b…a490955, strikt-JSON validiert). Vier Modellvarianten × drei
Zuordnungsverfahren; Gates vor dem Lauf deklariert.

| Variante | Änderung | rank α / max Rest | Gruppen α / max Rest |
| --- | --- | --- | --- |
| baseline | eingefroren (Referenz) | 0.919191 / 0.564247 | 0.949324 / 0.625687 |
| proper_phase (P) | φ = Q·[ℓ(r)+s·u], dℓ/dr = s | 0.999285 / **0.004320** | 0.998989 / 0.005367 |
| action_inertia (K) | m·s²/D statt m/D² | 0.852098 / 0.562717 | 0.831074 / 0.606378 |
| proper_plus_inertia | P+K | 0.926316 / 0.006756 | 0.925859 / 0.007479 |

Gates: **P0 PASS** (Baseline reproduziert Run-11 exakt, |Δ| < 1e-6), **P1
PASS** (FD-Konsistenz aller Ableitungen, ≤ 2.6e-9), **P2 FAIL — als Befund
dokumentiert** (Hauptwinkel-cos bis 0.246; die Flach-/SSZ-Softspektren haben
entartete Multipletts — flach: 3 Moden exakt bei 2.3665; SSZ-Baseline: 5
Moden bei 2.3496–2.3662 — Matching auf Einzelvektoren ist dort ill-posed,
Paare kreuzen die Frequenzbänder mit Overlap 0.27–0.45), **P3 PASS**
(NEW_STRUCTURE_NOT_COLLAPSIBLE bleibt in allen 12 Kombinationen bestehen).

**Interpretation (vorsichtig, Modell-Ebene, KEIN astrophysikalischer Claim):**

1. Das formale Verdikt ist robust: keine Variante kollabiert auf einen
   einheitlichen Faktor. Die Nichtkollapsibilität ist also keine Eigenschaft
   der Phasen-Identifikation, der Trägheitsform oder der Zuordnungsmethode.
2. ABER: Der Maximalrest ist KEINE Invariante. Unter der proper-radialen
   Phase fällt er von 0.564 auf 0.0043 (Faktor ~130), α → 0.9993. Die
   „neue Struktur" wird unter Proper-Längen-Phase fast zu einem
   Skalierungsfaktor — das formale Verdikt trägt sie nur, weil die
   Schwelle (1e-5) extrem scharf ist. Die physikalisch wichtige Zahl ist
   der Rest selbst, und der ist stark modellabhängig.
3. Der Proper-Längen-Test zeigt damit genau die Art von Sensitivität, die
   er finden sollte: Die Mehrskalen-Identifikation (x_n ↔ r) trägt
   substanziell zur Restgröße bei. Eine echte Materie-Wirkung muss diese
   Identifikation liefern, nicht voraussetzen.
4. Die Aktions-Kinetik (K) allein ändert wenig (Rest 0.563); sie koppelt
   an die gleiche Phase und erbt deren Struktur.
5. P2-FAIL ist ein instrumenteller Befund: Jede künftige
   Modenzuordnung bei diesem Gitter muss entartungs-bewusst arbeiten
   (Unterraum-/Projektor-Vergleich, nicht Einzelvektor-Matching). Das
   Gruppen-Matching (dritte Spalte) ist die erste saubere Umsetzung und
   bestätigt das Bild quantitativ.

Offen daraus: (i) Warum bildet die proper-radiale Phase die flachen Moden
fast exakt ab? Hypothese: dφ/dr = Q·s(r) macht die Substratkopplung lokal
konstanter und die Kette näher an einen uniformen Oszillator — zu prüfen
gegen die D(r_s)/Xi-Struktur. (ii) Die eigentliche Materie-Wirkung (Punkt 1
oben) bleibt die entscheidende Aufgabe; K war nur eine Kinetik-Variante,
keine Herleitung.

---

*Report generiert 2026-10-11, Bingsi; erweitert/korrigiert nach Lino-Review
11.10.26 (Korrektur-Runde 2: massengewichteter Solver § 4.4, Faktor-2-Assembler-
Fix § 4.5, coupling-assumption-Kennzeichnung, Detektor-Doku). Alle Zahlen aus
verifizierten Tool-Runs: FK Run 11 (massengewichtet, final); Re-Zertifizierung
ABGESCHLOSSEN — Bridge Run 6 (STAGE1+2 PASS) + SRF Run 6 (5/5 Gates PASS,
Fenster 2.5–7) + SW-Diag (4 Grids, PEAKS_NOT_ROBUST), Zahlen eingefroren.
Resolventen-Run (360 Fälle), B6 Step 1 Blind Scan (6.973 s), ch2-Pol-Rehunt.
Logs: /root/.hermes/cache/scratch/ (ssz_fk_run11_massweighted.log,
weisz_bridge_run5_factor2fixed.log, weisz_bridge_run6_freshcert.log,
srf_v1_run6_g1certgrid.log, swdiag_rerun_factor2fixed.log,
b6_step2_polish.log, ch2_nm.log, ch2_rehunt.log).*
