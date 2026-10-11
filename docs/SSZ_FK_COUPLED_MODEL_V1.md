# SSZ-FK Coupled Model V1 — Herleitung, Freiheitsgrade, Energiebedingungen

**Datum:** 2026-10-11 · **Branch:** spectroscopy-real-data-20261004
**Status:** MODELLSTUDIE — keine astrophysikalische Detektions- oder Resonanz-Claim
**Autor:** Bingsi, auf Direktive von Lino (Integrationsfrage statt weiterer
isolierter Weisz-Reproduktionen)

## 0. Forschungsfrage (Lino, wörtlich wiederhergestellt)

> Kann die eigentliche Weisz-/Frenkel-Kontorova-Dynamik als sinnvoller
> physikalischer Bestandteil eines SSZ-Modells verwendet werden, um damit
> neue astrophysikalische Vorhersagen zu berechnen?

Nicht: vergleiche zwei fertige Spektren. Sondern: leite κ_SSZ(r) und V_SSZ(r)
aus den SSZ-Gleichungen oder einer begründeten Materiekopplung her und
untersuche das gesamte Modell.

## 1. Zwei Anwendungen — getrennt behandelt

### Anwendung A: strukturierte Materie im SSZ-Hintergrund (NEUTRONENSTERNKRUSTE)

**Physikalische Begründung:** Eine NS-Kruste ist materiell ein
FK-System: Coulomb-Kristall von Kernen (Substratperiodizität V0·cos) +
elastische Scherung (Federkonstante κ). Beide Terme existieren unabhängig
von SSZ (Standard-Physik der Kruste). SSZ tritt über die LOCKED-Metrik ein:

    ds² = −D²(r) c² dt² + s²(r) dr² + r² dΩ²,   D = 1/(1+Ξ), s = 1+Ξ = 1/D

Für die Kruste (R_NS ≈ 10–12 km, r_s ≈ 4.4 km für 1.5 M☉ → r/r_s ≈ 2.3–2.7)
ist der Regime-Katalog eindeutig: g₁-Zweig operativ, Ξ(r) = r_s/(2r).

**Kopplung (3 zwingende Konsequenzen der Metrik, keine Ansatze):**

1. **Kinematik (Zeitdilatation):** Eigenzeit dτ = D(r) dt. Die Newton-
   Gleichung ü ∝ F gilt in Proper-Zeit; in Koordinatenzeit t:
       (m / D(r_n)²) · ü_n = F_n .
   Die Trägheit ist ortsabhängig: M_n = m/D_n².

2. **Elastizität (Proper-Länge):** Federenergie hängt von der PROPEREN
   Dehnung ab. Properer Abstand benachbarter Sites:
       ℓ_n = s(r_{n+½}) · (a + u_{n+1} − u_n)
   mit a = natürliche (proper) Gitterkonstante. Federterm:
       U_el = Σ_n (κ/2) (ℓ_n − a)².

3. **Substratphase (Segmente pro Wellenlänge):** Die Kristall-Periodizität
   2π/Q ist eine PROPERE Länge. Die Phase am Site n:
       φ_n = Q · s(r_n) · (n·a + u_n).
   Damit ist der **Inkommensurabilitäts-Parameter radially strukturiert:**
       b(r) = Q · a · s(r) = b₀ · (1 + Ξ(r))/(1 + Ξ(R))·…   (b₀ = Q·a·s(R))
   Weisz' konstanter Misfit wird zum **Misfit-Gradienten geometrischen
   Ursprungs** — genau das ist die SSZ-Weiterentwicklung des FK-Modells.

**Gesamtpotential + kinetischer Term (offene Kette, freie Ränder):**

    L = Σ_n (m/2D_n²) u̇_n²  −  Σ_n (κ/2)[s_{n+½}(a+u_{n+1}−u_n) − a]²
        − Σ_n V0 [1 − cos(Q s_n (n a + u_n))]

**Freiheitsgrade:** M reelle Skalare u_n (M = Sitezahl). Kein Gauge, keine
Reduktion nötig — alle drei SSZ-Kanäle (ψ, dφ, V) des kontinuierlichen
Operators sind hier NICHT die Freiheitsgrade; die Kruste koppelt als
Effektiv-Materiefreiheitsgrad an den Hintergrund (f, h, φ → D, s), nicht an
die radialen Feldmoden direkt. Diese Trennung ist bewusst und wird als
Näherung dokumentiert (keine Rückwirkung der Gitterschwingung auf D, s —
linearisierte Materie auf festem Hintergrund, konsistent mit Weisz, der
ebenfalls auf festem Substrat linearisiert).

**Energiebedingungen:**
- E = T + U ist exakt erhalten (autonomes mechanisches System).
- Stabilität ⇔ Hessian H = ∂²U/∂u² positiv semidefinit am Gleichgewicht.
  Die Substratterme tragen −V0 Q² s_n² cos(φ_n) zur Diagonale bei —
  bei großem β kann der Misfit-Gradient Regions mit cos(φ_n) < 0 erzeugen
  (lokale Instabilität). Der Code misst min eig(H) pro Parametersatz und
  meldet STABLE / UNSTABLE — das Ergebnis wird nicht weggebügelt.
- Tilt/Gravitation (Depinning) ist bewusst NICHT im Basismodell (wäre
  weiterer, separat zu begründender Term; hier nur der s-Gradient als
  neuer Misfit — minimales Modell).

### Anwendung B: Raumzeitsegmente SELBST als FK-Gitter

**Verdict: MIT DEM AKTUELLEN SSZ-FORMALISMUS NICHT HERLEITBAR.**

Begründung aus dem Kanon (ssz-complete-documentation, Known Limitations):
- Ξ(r) ist axiomatisch definiert; es existiert keine Wirkung/kein
  Variationsprinzip für Ξ. Segmente sind geometrische Zählgrößen
  (N₀ = 4 pro Periode), keine dynamischen Variablen mit Lagrange-
  Funktion — es gibt nichts, woraus κ und V0 für ein Segment-Gitter
  abgeleitet werden könnten.
- Der operative Kontinuumsoperator (K/G/S/M auf dem r-Gitter) ist ein
  radiales Feldsystem; eine Diskretisierung zu Segment-DOFs mit
  Federkopplung wäre ein NEUES Axiom („Segmente haben Scherelasticität"),
  nicht eine Herleitung.
- Falsifizierbare Konsequenz dieser Entscheidung: Anwendung B bleibt
  gesperrt, bis eine Segment-Dynamik (Wirkung für Ξ) existiert. Wir
  beanspruchen NICHT, dass sie unmöglich ist — nur, dass sie aus dem
  bestehenden Kanon nicht folgt.

## 2. Numerisches Modell (tools/run_ssz_fk_coupled_model.py)

- Offene Kette, M Sites, r_n = R + n·Δ (Krusten-Column), a = κ = m = 1.
- SSZ-Hintergrund: D_n = 1/(1 + r_s/(2 r_n)), s_n = 1/D_n (g₁-Zweig,
  LOCKED für r/r_s ≥ 2.2 — Parameterwahl respektiert den Regime-Katalog).
- Gleichgewicht: BFGS auf U mit analytischem Gradient; Hessian analytisch.
- Moden: verallgemeinertes Eigenproblem ω² v = M⁻¹ H v (M = diag(m/D_n²)).
- Response (Weisz-Maschinerie, explizite Projektionen):
      P(ω) = b† Im[(M⁻¹H − (ω+i0)²)⁻¹] b
  mit (Quelle, Detektor)-Paaren auf verschiedenen Sites/Kanälen.
- Nichtlineare Energieerhaltung: symplektischer Leapfrog auf der vollen
  nichtlinearen Dynamik.

## 3. Gates (VOR dem Lauf deklariert)

| Gate | Definition | Schwelle |
|---|---|---|
| G-A FK-Recovery | SSZ aus (D=1, s=1): akustische Dispersion offener Kette ω(q)=2√(κ/m)·sin(qπ/2(M+1)) für V0→0; Summenregel Σ_n W_n = 1 (q0-Probe) | rel err < 1e-6; ΣW=1 < 1e-10 |
| G-B Energie | Leapfrog, |ΔE|/E₀ über 2·10⁴ Steps | < 1e-6 |
| G-C Stabilität | min eig(H) am Gleichgewicht je Parametersatz | > −1e-8 ⇒ STABLE, sonst UNSTABLE (wird berichtet, nicht repariert) |
| G-D RSG-Konsistenz | Gitterwelle (V0=0, SSZ an): ω numerisch vs. WKB ∝ ∫ds-Gewicht — Gitterversion der locked Formel Φ = k∫s(r)dr | rel err < 1e-2 |
| G-E Neue Struktur | ν_SSZ vs ν_FK (gleiche Bare-Parameter): bestes uniformes Redshift-α anpassen; Struktur ist NEU ⇔ kollabierbar zu EINEM Faktor? Kriterium: max relative Restabweichung nach α-Fit | > 10 × G-A-Toleranz ⇒ NOT_COLLAPSIBLE (echte neue Struktur) |
| G-F Projektion | Antwort hängt von (Quelle,Dektor)-Projektion ab (Weisz-Punkt) | max|P_A−P_B|/max P > 10⁻³ |

## 4. Claim-Disziplin

- Alles MODELL-Ebene. Kein „Detektion"-Wort, kein Resonanz-Claim, der
  zertifizierte Stand PEAKS_NOT_ROBUST des GW-Pfads bleibt unberührt.
- NICER-Anwendung (echte QPO-/Torsionsmoden-Frequenzen) ist die nächste
  Stufe und braucht ein Instrumentenmodell — hier nur die Vorhersage-
  struktur (Frequenzmuster + Rotverschiebung f_obs = D(R)·f_local).

## 5. Validierungsstand (Run 10, 2026-10-11 — alle 6 Gates PASS)

| Gate | Ergebnis |
|---|---|
| G-A FK-Recovery | freq rel err **3.96e-14**, Summenregel 1−1e-15 → PASS |
| G-B Energieerhaltung | max rel ΔE **4.28e-8** über 2·10⁴ Leapfrog-Steps → PASS |
| G-C Stabilitätsleiter | alle 6 Stufen (V0=0.01…4.0) STABLE, min eig 0.035…112, grad_inf ≤ 1.6e-6 → PASS |
| G-D RSG-Konsistenz | s·d ≡ 1 (err 1.1e-16); Gitterabweichung 3.23e-4 bei dr=0.05 → 8.95e-5 bei dr=0.025, **Richardson-Ratio 3.61 ≈ 4** ⇒ O(dr²)-Diskretisierung, kein Modelleffekt → PASS |
| G-E Neue Struktur | α_best=0.919, max Rest nach α-Fit 0.564 ≫ 10×G-A-Tol ⇒ **NEW_STRUCTURE_NOT_COLLAPSIBLE** — das SSZ-FK-Spektrum kollabiert NICHT auf einen uniformen Redshift-Faktor |
| G-F Projektion | max rel Projektionsdifferenz ~1.0 → PASS |

### 5.1 Physikalische Lesart von G-E (Kern-Ergebnis der Modellstudie)

Das SSZ-Hintergrundfeld erzeugt im FK-Sektor KEINE bloße uniforme
Rotverschiebung, sondern strukturell neue Modenmischung:
- Trägheit m/D(r)² und Proper-Längen-Feder κ·s(r)² kompensieren sich im
  akustischen Zweig exakt (s·d≡1, G-D) — der akustische Sektor ist
  SSZ-invariant.
- Der Substratsektor trägt die Ortsabhängigkeit b(r)=Q·a·s(r), und die
  Wechselwirkung beider Sektoren am Gleichgewicht erzeugt Rang-verschie-
  bungen im Spektrum, die kein skalares α reproduzieren kann.
- Dies ist die erste konkrete, gecheckte Antwort auf die Integrations-
  frage: die Weisz-Dynamik AUF dem SSZ-Hintergrund ist mehr als
  „FK × Rotverschiebung".

### 5.2 Fehlerhistorie (Reproduktions-Wert)

Runs 1–9 deckten iterativ auf: (i) Substrat-Gradient-Vorzeichen
(`+V0·Q·s·sin φ` für ∇U) — FD-Verifikation von grad UND Hessian gegen
zentrale Differenzen (Fehler ~1e-8) ist der Absichert-Standard;
(ii) akustische Referenz für FREIE Ränder: ω_q = 2√(κ/m)·sin(qπ/2M),
q=0..M−1 (Nullmode inklusive) — nicht die festen-Rand-Formel;
(iii) Default-Misfit: b = a/φ (61.8% Misfit) ist tief-solitonisch und
besitzt KEIN stabiles homogenes Gleichgewicht (min eig < 0 wächst
monoton mit V0) — das Weisz-Regime ist kleiner irrationaler Misfit
δ=(√5−2)/4 ≈ 5.9% bei κ=25; (iv) G-D auf dem Gitter ist O(dr²), nicht
exakt — Richardson statt Maschinen-Nullstelle als Gate.
