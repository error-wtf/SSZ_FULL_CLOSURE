# SSZ Operator-Assembly Fehlerklärung — M3/M4 GÜLTIGKEITSENTSCHEIDUNG

Datum: 2026-10-10 · Referenz-Commit: 9648f6e2e23c21ff80710f69b8e338351e2a0226
Repo: github.com/error-wtf/SSZ_FULL_CLOSURE, Branch spectroscopy-real-data-20261004
Kontrolle: /root/.hermes/cache/scratch/m3_variable_kgsm_control.py — 9/9 PASS (Operator-Niveau)

---

## 1. Bestätigte Fehler (zertifiziert durch M3-Kontrolle)

### BUG-1: ML≡0 in den ECS-Adaptern (Faktor: Nulloperator)
Dateien (identischer kopierter Assembly-Block):
- tools/run_ssz_domain_expansion_v1.py (D1–D3-Scan, Zeile ~80: `ML = sp.lil_matrix(...)`, niemals `+=`)
- tools/run_ssz_domain_expansion_v2_l12.py
- tools/run_ssz_domain_expansion_v3_l20.py
- tools/run_ssz_domain_expansion_v4_l42.py

Shift-Invert-Operator OP = (A−σL)⁻¹L ≡ 0 (M3 C3a/C3b: max|OPx| = 0.0 exakt).
Jeder "Kandidat" aus diesen Pfaden ist ein Artefakt von eigs() auf dem
Nulloperator; die Residual-Gates (res ≤ 1e-6 gegen A v = λ L v mit L=0)
sind trivial nicht erfüllbar → n=0 in ALLEN ECS-Kandidatendateien ist
KEIN physikalisches Nullresultat, sondern ein Konstruktionsartefakt.

Zusatzbefund v3-Box-Run (run_ssz_domain_v3.py, ECS_SSZ_CANDIDATES_V3):
ML ist hier NICHT null, sondern mit M_phys befüllt. Das ist die FALSCHE
Gewichtsmatrix: Jost (B0 = w²K − (M+S'/2)) und ECS-Discovery
(C2 = blk(Kc)) benutzen K_phys als w²-Koeffizient. Der v3-Run löst
damit ein anderes Eigenwertproblem. → ebenfalls invalid.

### BUG-2: Faktor-2-Stencil in den zentralen Ableitungen
Dateien:
- tools/run_spectral_weight_diagnostics_gw.py (Zeilen 50/52: /hc statt /2hc)
- tools/run_ecs_discovery_v2.py (Zeilen 111/113: /h[1:] statt /2h)

M3 C1a: ratio exakt 2.000000 für quadratische Profile.
M3 C1b: AUCH lineare Profile sind exakt verdoppelt (zentrale Differenz
ist für lineare Funktionen exakt, der Nennersherr macht daraus 2×).
M3 C2b/c: für die tatsächlichen SSZ-Profile G_phys, S_phys exakt 2×
gegenüber np.gradient. Weil G_phys um Faktor ~55 variiert, sind G′ und
S′ physikalisch bedeutsame Terme — der Fehler verändert den Operator.

Konsequenz für die KAT-Zertifizierung: Die ECS_V2_2-Zertifizierung
(ECS_V2_2_CERTIFICATE.json, status ECS_V2_2_CERTIFIED) lief auf
run_ecs_square_well_v22.py — dort ist ML korrekt befüllt (Zeile 104
`ML[...] += ml[i,j]`) und das Profil konstant (V0=−2.5, dG=dS=0), der
Stencil-Bug ist dort also UNSICHTBAR (0×2 = 0). Die Solver-Zertifizierung
selbst bleibt gültig. Der Fehler liegt ausschließlich in den
SSZ-Adaptern, die das zertifizierte Verfahren auf die variablen
SSZ-Profile übertragen.

---

## 2. Gültigkeitsentscheidung je Befund

### BLEIBT GÜLTIG (unverändert, keine Neuberechnung nötig)
| Artefakt | Grund |
|---|---|
| ECS_V2_2_CERTIFICATE, JOST_V2_2_CERTIFICATE, COUPLED_RESONANCE_SOLVER_V2_2_CERTIFICATE | Zertifizierung auf konstantem Referenzproblem; run_ecs_square_well_v22.py Assembly sauber (ML befüllt, Galerkin). Kein ML- und kein Stencil-Pfad betroffen. |
| JOST_SSZ_*_CANDIDATES (alle Tiles/L) | Jost-Bein nutzt der()/deriv_arrays() mit korrekter Sekantenformel (M3-verifiziert) und keinen ML-Pfad. |
| NO_CERTIFIED_RESONANCE_IN_SCANNED_DOMAIN (Katalog V2, Original-Box) — EINSCHRÄNKUNG | Jost-Bein gültig; das ECS-V2-Bein (ECS_V2_CANDIDATES via run_ecs_discovery_v2.py) ist vom Stencil-Bug betroffen → der DOPPEL-Null bleibt bestehen, aber das ECS-Bein muss mit korrektem Stencil bestätigt werden, bevor der Befund als unabhängig doppelt getragen gilt. Status: GÜLTIG MIT EINSEITIGER TRAGKRAFT (nur Jost). |
| ECS_SSZ_CANDIDATES_V3.json | WIRD NEU GERECHNET (siehe unten) — der aktuelle Inhalt (n=255 "Kandidaten" via falsche M-Matrix) ist invalid, aber auch ohne ihn bleibt der Jost-Null der Box bestehen. |

### INVALID — UNDER_REVIEW (Neuberechnung erforderlich)
| Artefakt | Bug | Neuberechnung |
|---|---|---|
| ECS_SSZ_D{1,2,3}_CANDIDATES_V4.json | ML≡0 | läuft (M4-Rerun A) |
| ECS_SSZ_L12/L20/L42_D{1,2,3}_CANDIDATES_V1.json | ML≡0 | ausstehend (nach M4-A-Muster) |
| SSZ_DOMAIN_EXPANSION_V1_MATCHING.json + V2_MATCHING_L12/L20/L42 | "common" war leer, weil ECS-Seite konstruktiv leer | nach ECS-Reruns |
| NO_CERTIFIED_RESONANCE_IN_EXPANDED_DOMAIN_V1.json, NO_CERTIFIED_RESONANCE_L12/L20/L42_V1.json | Doppel-Null beruht hälftig auf Nulloperator | nach ECS-Reruns; die Jost-Seite bleibt gültig |
| artifacts/SSZ_SPECTRAL_WEIGHT_DIAGNOSTICS_GW_V1.json | Stencil-Faktor 2 in dGz/dSz | läuft (M4-Rerun B) |
| ECS_SSZ_CANDIDATES_V3.json (Box-Run ECS-Bein) | ML=M statt K (falsche Gewichtsmatrix) | ausstehend |

Alle 22 betroffenen Artefakte sind SHA-256-verzeichnet in
/root/.hermes/cache/scratch/m4_under_review_hashes.json (Stand vor
Neuberechnung; die Originale bleiben unverändert im Repo).

### NICHT ABGELEITET, ABER ABHÄNGIG — Hinweis
Die Zusammenfassungen in Memory/Sitzungslog ("NUMERICALLY_VERIFIED NULL,
0 Peaks über 6σ, max|G| 0.029-0.049") erben den Status der GW_V1-Datei
und sind bis zum Abschluss von M4-B als UNDER_REVIEW zu führen.

---

## 3. Was die Neuberechnungen zeigen (ABGESCHLOSSEN 10.10., alle 4 Läufe verifiziert)

### M4-RERUN A — ECS D1–D3 mit gefüllter K-Massenmatrix (ML = P1-Masse von K_phys)
Wall: 346.5 s, h=0.02, θ=45°, identische σ-Lattice + Refinement-Regel wie V1.
Ergebnis (ECS_SSZ_{D1,D2,D3}_CANDIDATES_V1_FIXED_ML.json in Scratch):
- D1: 135 Rohkandidaten (res ≤ 1e-6), Re ∈ [1.50, 1.84], Im ∈ [−0.34, −0.05],
  hiervon 5, die beim halbierten Refinement-Lattice wiedergefunden werden.
- D2: 0 Kandidaten. D3: 0 Kandidaten. (Der Produktions-Null in D2/D3
  reproduziert sich am korrigierten Operator.)
- ABER: Theta-Stabilität der 5 D1-Kandidaten (Gate G3-Logik, θ=40/45/50°,
  Datei M4_THETA_STABILITY_D1.json): 0/5 stabil. Vier finden sich bei
  anderen Thetas nicht wieder (>0.02 Abstand), einer hat Spread 5.5e-3
  > Gate 2e-3. → KEIN einziger D1-ECS-Kandidat besteht das
  Zertifikats-Gate. Der ECS-Null der erweiterten Domäne bleibt nach
  Korrektur ALSO BEstanden — jetzt erstmals auf einem validen Operator.
  (Einschränkung: G3-Spread allein ist notwendig, nicht hinreichend;
  eine vollständige C-R-Kette inkl. Jost-Gegenprobe steht noch aus.)

### M4-RERUN B — GW-Diagnostik mit 2h-Stencil (n=400 Zertifizierungspass)
Dateien: g_n400_{dirichlet,robin}_{prod,fixed}.npy + zwei Summary-JSONs.
- max|G|: prod 0.0729 → fixed 0.0855 (+17 %), median +15 %.
- Der in der GW_V1 beobachtete 40-%-Abfall 600→1000 Knoten bleibt
  unexplained durch den Stencil allein; die gitterkonsistente
  Punktquellen-Normierung ist weiterhin offen (siehe Abschnitt 5).
- Peak-Suche (identischer lokaler 6σ-MAD-Filter wie Produktion):
  prod 0 Peaks (max_sig 1.32σ), fixed 0 Peaks (max_sig 1.33σ),
  je für Dirichlet und Robin.
- → Das GW-Nullresultat übersteht die Stencil-Korrektur: die
  Antwortkurve ändert sich um ~15 % in der Amplitude, aber es
  entstehen KEINE Peaks. Der Befund "keine auflösbaren Peaks im
  Fenster [0.05, 1.5] bei n=400, beide BCs" ist am korrigierten
  Operator reproduziert. (Die 6σ-Aussage bleibt trotzdem ein
  Peakfilter-Statement, kein Rauschmodell mit FAR-Kalibrierung.)

### FEHLEND FÜR DEN ABSCHLUSS (offen, nicht Teil dieser Entscheidung)
- L12/L20/L42-ECS-Reruns (gleiches ML-Fix-Muster wie M4-A).
- v3-Box-ECS-Rerun mit ML=K statt ML=M.
- Jost-Gegenprobe der 5 D1-Kandidaten (C-R-Kette).
- Gitterkonsistente Quell-Normierung + echte retarded-Green-Observable
  A_b(w) = −(1/π) Im[b†G^R(w)b] statt |G_00|; QNM-Randbedingungen.

---

## 3b. M5-WEAK-FORM-KONTROLLE (10.10., zweite Nutzer-Prüfung) — TEILWIDERRUF

Kontrolle: /root/.hermes/cache/scratch/m5_weak_form_control.py — 8/8 PASS.
Anlass: Nutzer-Review des M4-Runners fand drei weitere Operatordefekte.
Alle drei sind am Code VERIFIZIERT und quantifiziert:

### BUG-3: Phantom-Term −G′ in der schwachen Form (schwerwiegend)
Starke Form (identisch in Jost V2.2 `coeffs()` A1 = Gp − S und GW-Runner
C0 = blk(G)D2 + blk(dG−S)D1): G Ψ″ + (G′−S) Ψ′ + (ω²K − M − S′/2) Ψ = 0.
Partiell-Integration (Randterm vernachlässigt) ⇒ Divergenzform:
  ∫ G Ψ′v′ + ∫ S Ψ′v + ∫ (M + S′/2) Ψv = ω² ∫ K Ψv.
Der D1-Koeffizient ist **+S**. Produktion (run_ssz_domain_expansion_v1.py:108)
UND der M4-Fix (m4_rerun_ecs_fixed.py:95) assemblieren stattdessen
`s_ab = S − G′` ⇒ zusätzlicher Phantom-Term −∫ G′ Ψ′ v′, der in der
schwachen Form nicht existiert. Größenordnung auf dem SSZ-Profil
(Innenbereich r ≤ 30): ‖G′‖ = 4.18e5 gegen ‖S‖ = 7.25e1 — der
Phantom-Term dominiert den korrekten Term um **Faktor 5.8e3**.
Hinweis: Die Ableitung selbst ist im M4-Fix korrekt (2h-Form, mediane
rel. Abweichung np.gradient vs. 2h = 2.9e-14 auf dem nichtuniformen
Gitter) — der Fehler ist die schwache Form selbst, nicht der Stencil.
Betroffen: ALLE ECS-SSZ-Pfade inkl. M4-RERUN A. **Nicht** betroffen:
GW-Pfad (starker Operator, FD, korrekt dG − S) und Jost-Beine.

### BUG-4: Kein ω²-Term im ECS-Außengebiet
m4_rerun_ecs_fixed.py:84-87 (und Pendant in Produktion v1): Außengebiet
mit Km = 0, Gm = I3, Sm = Mm = 0 ⇒ die äußere Gleichung enthält keinen
ω²K-Term und damit keine frequenzabhängige Wellenpropagation. Der
Konturabschnitt repräsentiert nicht die auslaufende Wellenbedingung.
(Vgl. ECS_V2_2-Quadrattopf-Zertifikat: dort korrektes Außenverhalten —
Defekt nur im SSZ-Adapter.)

### BUG-5: Zentrumsbasis geladen, nie angewendet
m4_rerun_ecs_fixed.py:49 lädt V4_CENTER_REGULAR_BASIS_V2.npz
(r_ref = 0.05, Exponenten [0.539, 2.000, 1.002]); es existiert KEINE
weitere Referenz auf `cb`. Der Zentrumsknoten behält seine natürliche
(Neumann-artige) RB. Jost V2.2 startet dagegen MIT der regulären Basis
(Y0, dY0) ⇒ ECS und Jost lösen bislang NACHWEISLICH nicht dasselbe
Randwertproblem; eine Jost-Gegenprobe der ECS-Kandidaten wäre so
ohnehin bedeutungslos gewesen.

### TEILWIDERRUF der M4-Bilanz (Abschnitt 3, Netto-Bilanz)
1. D1-ECS (135/5, 0/5 theta-stabil): Kandidaten entstanden auf einer
   Assembly mit BUG-3 (Phantom-Term, Faktor 5.8e3 dominant), BUG-4 und
   BUG-5 → als EVIDENZ INVALID, nicht nur als Interpretation fraglich.
   Die Aussage "jetzt erstmals auf einem validen Operator" (Abschnitt 3)
   wird HIERMIT ZURÜCKGEZOGEN. M4_THETA_STABILITY_D1.json (gleiches
   build_fem) erbt den Status.
2. D2/D3-ECS-Null "reproduziert am validen Operator": ebenfalls
   ZURÜCKGEZOGEN — reproduziert zwar, aber auf demselben defekten
   Operator (BUG-3/4/5). Status zurück auf UNDER_REVIEW.
3. GW-Null (0 Peaks): BLEIBT — starker FD-Operator mit korrekter
   Formulierung (dG − S), vom M4-B-Rerun mit korrigiertem Stencil
   getragen; Normierungsfrage separat offen (Abschnitt 3).
4. Gültig bleiben unverändert: ECS/JOST/COUPLED_V2_2-Zertifikate
   (Konstantprofil, saubere Assembly), alle Jost-SSZ-Nulls.

### NEUE REIHENFOLGE (Nutzer-Vorgabe 10.10., übernommen)
NICHT sofort L12/L20/L42. Zuerst eine FEM-Neuassembly mit:
  (a) korrekter schwacher Form (D1-Koeffizient +S, kein G′-Term),
  (b) variablen K, G, S, M im gesamten Innenbereich,
  (c) ω²K-Term + analytisch fortgesetzte Koeffizienten im Außengebiet
      (Muster: run_ecs_discovery_v2.py F3, 1/r-Fit),
  (d) Zentrums-RB aus V4_CENTER_REGULAR_BASIS_V2 (Muster: F4,
      Basis-Replace an den ersten beiden Gitterzeilen),
  verifiziert GEGEN den Jost-Operator (C-R-Kette am selben npz).
Erst danach Neubewertung der D1-Kandidaten und L-Scans.

### UNDER_REVIEW-ERWEITERUNG
Zusätzlich zu den 22 Artefakten des SHA-Manifests (Stand 9648f6e) sind
jetzt auch die M4-Rerun-A-Artefakte UNDER_REVIEW:
ECS_SSZ_{D1,D2,D3}_CANDIDATES_V1_FIXED_ML.json, M4_THETA_STABILITY_D1.json,
M4_RERUN_ECS_SUMMARY.json sowie m4_rerun_ecs_fixed.py /
m4_theta_stability.py als dokumentiert-defekte Runner (Bug-3/4/5).
Unberührt: M4-RERUN-B (GW) samt M4_RERUN_GW_SUMMARY.json,
M4_RERUN_GW_ROBIN_ONLY.json, m4_rerun_gw_stencil.py.

### NETTO-GÜLTIGKEITSBILANZ (nach Korrektur) — ÜBERHOLT durch Abschnitt 3b
> Status 10.10. (zweite Prüfung): Punkte 1 und 2 dieser Bilanz sind durch
> BUG-3/4/5 WIDERRUFEN; maßgeblich ist die Bilanz in Abschnitt 3b.
> Dieser Abschnitt bleibt aus Dokumentationsgründen unverändert stehen.
1. D2/D3-ECS-Null: REPRODUZIERT am validen Operator (stark — vorher
   wertlos, weil Nulloperator).
2. D1-ECS: es existieren echte Pencil-Eigenwerte nahe Re≈1.5–1.55,
   aber keiner ist theta-stabil → kein zertifizierbarer Kandidat.
   Der Nullbefund steht, mit neuem, jetzt validem Beweis.
3. GW-Null (0 Peaks): REPRODUZIERT am korrigierten Stencil.
4. Der DOPPEL-Null der erweiterten Domäne (NO_CERTIFIED_RESONANCE_*)
  ist damit qualitativ GESICHERT, aber die publizierten Artefakte
  tragen den Beweis bislang auf invaliden Operatoren → die
  UNDER_REVIEW-Markierung bleibt bestehen; die Neuberechnung
  (dieser Rerun) LIEFERT den validen Ersatz.


---

## 3c. M6: FEM-Neuassembly + Jost-Gegenprobe am selben npz (2026-10-10, Commit 3fa959f + M6-Dateien)

Runner: `m6_fem_jost_crosscheck.py` · Ergebnis: `M6_FEM_JOST_CROSSCHECK.json`
· Log: `M6_RUN.log` · Ausführung: 5/5 Gates PASS, 0 FAIL, wall 39,7 s.

Implementiert alle drei Nutzer-Punkte gegen den M4-Rerun:
1. **Schwache Form**: D1-Koeffizient **+S** (Element-Term `de[i,j]*s_ab` mit
   `s_ab = S`, KEIN `S − G′`), Testschleifen korrekt transponiert
   (`A[test=jj, trial=ii]`) — der S-Term war im ersten Entwurf dieses
   Runners selbst noch transponiert (antisymmetrisch → Vorzeichenflip);
   vom Smoke-Test G1b abgefangen, VOR dem Lauf.
2. **Außenzone**: volles ω²K — K, G, S, M via dokumentierter 1/r-Fortsetzung
   der V3-Profile auf dem rotierten Konturabschnitt (fit_rel_max 2,2e−2,
   Interpolation, kein Fit).
3. **Zentrums-RB**: `Psi'(r0) = Trel·Psi(r0)`, `Trel = dY_reg·Y_reg⁻¹`,
   einseitiger 2.-Ordnungs-Stencil in Mesh-Zeile 0 (bc_row_norm 1,5e2,
   Zeile aktiv).

Gates (vor dem Lauf deklariert):
- **G1** String-Kontrolle der Weak-Assembly: rel_max 3,3e−4 (tol 6,2e−6). PASS
- **G1b** FEM-weak vs. GW-validierter FD-strong-Operator, IDENTISCHES BVP
  (Realteil [r0,30], Dirichlet): rel 0,0–0,4 % (tol 2 %). PASS — die
  korrigierte schwache Form löst denselben Operator.
- **G2a/G2b** Struktur: PASS.
- **G3** Band-Scan Re ω 0,6–4,4 / Im ω −0,55..−0,02, 112 Shifts,
  Shift-Invert-Residuum-Gate 1e−8: **0 Kandidaten**. (D1-Band in Re 1,2–2,2
  eingeschlossen.) Kein Pol des korrigierten Operators in der Domäne.
- **G4** Jost-Gegenprobe (zertifizierte Metrik: rohes σ_min der 6×6-
  Matching-Matrix, Kalibration JOST_ECS_MATCHING_V2 „smooth, NO isolated
  minimum"): σ_min glatt, min 6,1e−4 / max 1,8e−3, **Kontrast 3,1** (kein
  isoliertes 2D-Minimum). Alle 5 M4-Feinkandidaten: σ_min 8,0–8,9e−4,
  Spanne < 12 %, auf Umgebungsniveau — kein Kandidat ist als Resonanz
  hervorgehoben. r_match-Leiter (20/30/40): monoton fallend an Gitter-Minimum
  UND Kandidat (Faktor ~2,8), kein struktureller Unterschied. PASS.

### Konsequenz
- Der korrigierte ECS/FEM-Operator und der Jost-Operator tragen jetzt
  DENSELBEN Nullbefund im D1-Fenster: **kein theta-stabiler Pol, kein
  Jost-Minimum, keine FEM-Kandidatin mit res ≤ 1e−8.**
- Die fünf D1-Feinkandidaten aus M4 (ECS_SSZ_D1_CANDIDATES_V1_FIXED_ML)
  sind damit endgültig widerlegt: nicht theta-stabil (M4), nicht im
  FEM-Spektrum des korrigierten Operators (M6 G3), nicht als Jost-
  Struktur sichtbar (M6 G4). Status der Datei: INVALID bleibt, jetzt
  mit positiver Gegenprobe statt nur Konstruktionskritik.
- M4-Rerun-B (GW-Null) bleibt unberührt und wird durch M6 G1b zusätzlich
  getragen (FD-strong und FEM-weak stimmen am echten Profil überein).
- **D2/D3 und L12/L20/L42**: die Neuberechnung auf M6-Assembly ist jetzt
  freigegeben (Nutzer-Reihenfolge eingehalten: erst Operator-Level-
  Verifikation, dann L-Scans). Bis dahin bleibt die UNDER_REVIEW-
  Markierung der ECS-Kandidatendateien bestehen.

## 4. Vorgehen (Nutzer-Vorgaben, eingehalten)

- Zuerst variable-K/G/S/M-Kontrolle (M3, 9/9 PASS), dann Reproduktion. ✓
- Zertifizierte Solver-Kontrollen unangetastet. ✓ (ECS/JOST_V2_2-Zertifikate bleiben stehen; M3 belegt, warum sie nicht vom Bug berührt sind)
- Keine Änderung von SSZ-Parametern, Branches, Beobachtungsdaten. ✓
- Keine Fits. ✓
- Alte Befunde erhalten (Originale unverändert), Fehler versioniert via SHA-Manifest + dieser Report. ✓
- Produktionsdateien nicht geändert; Reruns als exakte Kopien in Scratch. ✓

## 5. M7: Jost/FEM-Konsistenznachweis auf korrigiertem Integrator (10.10., Nutzer-Vorgabe nach M6-Code-Review)

Auslöser (Nutzer-Code-Review von M6, alles bestätigt):
1. `m6_fem_jost_crosscheck.py::_rk4` hatte `h = d*(x1-x0)/n_steps` mit
   `d = ±1` — bei x1<x0 ist d=−1 UND (x1−x0)<0, also h>0: die INWARD-Leg
   integrierte 45→70 nach AUSSEN (M7-A3 bewies: 0 von 1500 Samples im
   Zielgebiet). M6-G4 mass sigma_min an einer fiktiven Stelle.
2. `outgoing_basis()` ersetzte bei Re(k)<0 nur den Eigenwert (k→−k) ohne
   den Eigenvektor — M7-B mass die Residuen: geflippte Paare haben
   Residuum 4.4e−1 statt 5e−13. Der M6-Flip war kein gültiges Eigenpaar.
3. M6-G3 zertifizierte 0 Kandidaten bei 0 gescannten (Band-Scan ohne
   Fund); M6-G2 nutzte Toleranz 3e−2 statt dokumentierter 1e−6
   (gemessen 2.19e−2).

M7 (m7_jost_fem_consistency.py, SHA256 338b0c7f…, RUN12 = offizieller
Lauf, M7_RUN.log + M7_JOST_FEM_CONSISTENCY.json, 7/7 PASS):

- A  RK4-Richtungskontrolle (y''+k²y=0): outward 3.0e−10, inward 3.0e−9
  (tol 1e−8) + Bug-Reproduktion (Walk 8→15, nie im Zielgebiet). PASS.
- B  Eigenpaar-Erhaltung am echten SSZ-Profil: alle 3 gewählten
  outgoing-Kanäle res ≤ 5e−13 (tol 1e−10), 0 Flips nötig; Flip-Residuen
  4.4e−1/4.4e−1/2.0e−2 dokumentieren den M6-Defekt. PASS.
- C  M6-G4-Kriterien auf korrigiertem Apparat (frozen M6-Gates):
  contrast=4.5 (vorher 3.1), alle 5 D1-Kandidaten unauffällig
  (1.70–1.75e−2 vs Umgebung). Das M6-G4-URTEIL (keine Jost-Struktur,
  D1-Kandidaten widerlegt) hält dem korrigierten Apparat stand. PASS.
- D  Referenzpol-Recovery (attraktives PT V=−l(l+1)sech², l=2.5; exakte
  gebundene Pole 2.5i/1.5i/0.5i): wiedergefunden mit Fehler 0, 2.2e−16,
  2.4e−5 (Gate 5e−3). Die Jost-Nullstellensuche ist damit an exakter
  Referenz kalibriert — inkl. der beidseitig INWARD integrierten Beine
  (exakt die reparierte Operation). PASS.
- E  FD-strong (n=1351) vs korrigierter Jost, Box [r0,30], identische
  Zentrums-RB (Trel) + Dirichlet: 6/6 tiefste Moden innerhalb 5e−3
  (max 4.2e−3), 52 Jost-Moden gegen 3186 FD-Moden durchgehend ≤ 2.7e−3
  für w>3. Keine Diskrepanz. PASS.
  Entwurfskorrekturen (alle vor dem jeweiligen Lauf deklariert):
  (i) FEM-weak-Box-Pencil junk-dominiert (983 Spurmoden < w=0.1; bestes
  Eigenpaar-Residuum 4.8e−2) → FD-strong-Operator aus M6-G1b verwendet;
  (ii) Jost-Modensuche mit Becken-Verfeinerung (ein echtes Becken lag
  zwischen den 0.015-Gridpunkten; Konvergenzleiter FD n=451→901→1351
  belegt O(h²)-Verschiebungen bis 7e−3 bei w≈3.5).

Entscheidung:
- Die M6-G4-Absage der fünf D1-Kandidaten ist JETZT auf belastbarem
  Apparat reproduziert (M7-C) und die Jost-Seite ist an Referenzpolen
  kalibriert (M7-D). Status der D1-Datei: INVALID bleibt endgültig.
- Die FEM-Seite von M6 ist als Zertifizierungsinstrument im Box-Fall
  NICHT belastbar (junk-dominiertes Pencil, M7-E-Diagnose) — die
  G1b-Grenzkontrolle bleibt gültig, aber M6-G3-artige Kandidaten-
  Zertifizierungen auf diesem Pencil sind ohne Residuen-Gate pro Paar
  unzulässig. Für D2/D3 + L12/L20/L42 gilt: Freigabe NUR mit dem
  M7-Apparat (FD-strong + korrigierter Jost + M7-E-Methodik).
- Manifest/CI: MANIFEST.json + SHA256SUMS müssen mit den Audit-Dateien
  regeneriert und committet werden (P5-reproducibility-Workflow war
  deshalb rot, lokales --check war exit 1).

CI-Verifikation (Commit 8ee40752, live per GitHub-API abgerufen):
- Spectroscopy source-of-truth audit = SUCCESS
- SSZ P5 reproducibility = SUCCESS (Fix bestätigt; vorheriger Stand
  9596047 war dort ROT: "Manifest differs from current files")

Artefakte (SHA256):
- m7_jost_fem_consistency.py  338b0c7fcf660c232439f55e1b340dc689846004290d6f736ecdb8bd905eb67a
- M7_JOST_FEM_CONSISTENCY.json 8dc12bd418f9520b55382a02991d747735c97040cc04a2ae7ee09dd096991a22
- M7_RUN.log  7b8ba89d14c5229cbb9374dfebe4e24fdb1ff77ada8e222140eff2541bc45490
