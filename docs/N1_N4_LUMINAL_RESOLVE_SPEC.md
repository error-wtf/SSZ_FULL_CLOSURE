# N1-N4 SPEC: Luminal-domain background re-solve (der saubere H1/H3-Test)

Status: SPEC, nicht implementiert. Aufwand: eigenes Arbeitspaket (Multi-Day).
Ziel: entscheidet H1 (Pocket = Jet-Artefakt) gegen H3 (echter Ghost) SAUBER.

## Hintergrund (warum noetig, Stand 05.10.26)

Der Decision-Test (kscalar_domain_decision_test.py) zeigt:
- Das K_scalar-Pocket (x 0.294-0.324) koppelt an die out-of-domain Jets
  (G4X~ -2.6e3, G4XX~1.1e13, G4phiX~ -1.5e10, G3X bis -5e10).
- Jet-NULLIERUNG auf der vorhandenen Kette ist KEINE Reparatur: sie erzeugt
  eine andere Kette mit 888/3299 negativen K-Zeilen (x 0.337-1.011,
  K_min -0.0739). Beweis: CORE_QUARTIC_LUMINAL_CLEAN_41.csv.
- Grund: Das archivierte (f,h,phi,G2)-Profil ist NICHT als Loesung der
  luminalen Hintergrundgleichungen konstruiert. Setzt man die Jets null,
  ist das Profil off-shell -> die daraus emittierte 41er-Kette kann nicht
  gesund sein. E00-Residuum-Check (05.10.26): median |E00| ~ 0.13 bei
  G2-Spalte vs. G2 aus E00=0 ~ 0.09 — das Profil loest E00 nicht.

## Was zu implementieren ist

### N1: Luminaler Aktionssektor (symbolisch, fertig vorbereitet)
Luminaler Domain: G4 = G4(phi), G4X = G4XX = G4phiX = 0, G3 = G3X = G3phi = 0,
G5 = 0, G2 = G2(X,F) mit G2F = 1 (Maxwell-Normalisierung).
Luminal-spezialisierte C-Jets (bereits per sympy abgeleitet, kt_mh_background):
  C1 = -2 G4phi h,  C2 = C3 = 0,  C4 = -2 G4,  C6 = G2 - 2 G4phi h ph^2,
  C7 = -4 G4phi h ph,  C8 = 2 G4 (1-h),  C10 = -G4 h.
Damit E00, E11, E22 (existieren schon als kt_e00/11/22 mit Substitution) plus
Skalargleichung (Kombination bzw. eigene Variation nach phi).

### N2: Hintergrund-Solver
Unbekannte: f(r), h(r), phi(r)  [A0' = 0 im elektrolosen Sektor oder Maxwell-
current-Gleichung (G2F sqrt(h/f) r^2 A0')'=0 falls geladen].
Da das Projekt-Asymptot f -> 1/4, h -> 1, phi -> 1 ist (KEIN asymptotisch
flacher RH), sind die Randbedingungen:
  am äusseren Rand (u klein / r gross): f -> 1/4, h -> 1, phi -> 1, phi' -> 0
  am inneren Rand (u = 0.71 bzw. r ~ 1.4): Regularität / match an den
  vorhandenen Kern-Regime-Übergang (genaues RB noch festzulegen gegen die
  F1b-Konvention; aus der bestehenden Kette ablesen und DOKUMENTIEREN).
Methode: Relaxation/Collocation auf dem vorhandenen x-Grid (3299 Punkte)
oder solve_ivp mit Shooting; scipy steht bereit.
Konvention "u = 1/r" beachten (Repo nutzt u-Spalte).

### N3: Primitive + Emission
Auf der NEUEN Loesung: quartic_g5zero_primitives (unveraendert! dieselbe
Kase-Tsujikawa-Kette, die seitenexakt verifiziert ist) -> 41-slot emission
(emit_from_primitives) mit dem GLEICHEN c2-Zielprofil wie bisher
(sha 253cfe04...) bzw. — sauberer — c2 aus dem Domain neu bestimmt, WO die
G2XX-Normalsteuerung erreichbar ist (reachable_fraction dokumentieren).

### N4: Bewertung (Decision-Regeln VOR dem Lauf deklarieren)
1. K_scalar auf der N3-Kette (Pocket-Band UND global).
2. finite-L Health (K/c_r^2) fuer L = 6..1000 (existierendes Reducer-Modul).
3. Vergleiche mit (a) Original-Kette, (b) jet-nullierter Kette.
Decision:
  K >= 0 ueberall + Health pass  => H1 BESTAETIGT (Pocket = Artefakt der
     out-of-domain Jets). H2 wird zur Dokumentationsluecke. Freigabe-Pfad
     fuer F.4 (globaler Health-Rerun) offen.
  K < 0 im Pocket-Band           => H3 LEBT auf sauberem Boden: echter
     Even-Parity-Geist im luminalen Sektor bei x~0.30 — physikalisches
     Resultat erster Guete.
  K < 0 an ANDEREM Ort           => neue Pathologie-Map = Ausgangspunkt
     fuer den naechsten Lokalisierungsschritt.

## Abgrenzungen
- Nicht Teil von N1-N4: deep full-G5 subcore (u >= 50) — bleibt separater
  Bridge (bestehende Abgrenzung).
- Nicht Teil: QNM-Solver-Freischaltung. Das folgt erst NACH N4-H1 + F.4.
- Anti-Circularity: N4-Ergebnis wird VOR dem ersten Lauf hier deklariert
  (geschehen), nie nachtraeglich angepasst.

## RESOLUTION V3 (05.10.26): Solve konvergiert — System hat nur die triviale Loesung

Status: N2-V3 ausgefuehrt (`tools/run_luminal_background_solve_v3.py`).
Alle offiziellen N2-Blocker (Commit e8a9584) behoben. Artefakte:
`data/generated/spectral/N2_LUMINAL_SOLVE_RESULT_V3.json`,
`N2_LUMINAL_BACKGROUND_PROFILE_V3.csv` (3299 Punkte, VOLLE Domain
u in [0.71, 100]), `N2_LUMINAL_J0_FAMILY_V3.csv`,
`N3_N4_LUMINAL_DECISION_V3.json`.

### (a) Was V2 tatsaechlich auferlegte (Forensik, maschinell verifiziert)
- V2 differentiierte mit `numpy.gradient(..., u)` — das ist d/du, NICHT
  d/dr. Korrekt: d/dr = -u^2 d/du bei r = 1/u.
- V2 nutzte X = +phi_u^2/(2f); statische radiale Konvention des Repo
  (`geometry/p5.py`): X = -h phi_r^2/2 (Branch-Identitaet h phi_r^2+2X=0).
- V2's E11 fehlte der C9-Rest -h*phi'^2 exakt (sympy: E11_vollstaendig -
  E11_V2 = -h*ph**2). V2's cost 2.4e-15 war ein Artefakt des thereby
  veraenderten Systems.
- E22 wurde nie aufuerlegt oder zertifiziert.

### (b) E22-Klaerung (symbolisch + numerisch, nicht per Annahme)
Im konstanten Sektor (G4=1/2, G4phi=0, G2=X, A0'=0):
- E00 = 0 loest algebraisch nach h', E11 = 0 nach f' (jeweils erster Ordnung).
- Die Skalar-Konsistenz schliesst das System exakt:
  (r^2 sqrt(fh) phi')' = 0  <=>  J := r^2 sqrt(fh) phi' = J0 (konstant).
- E22 == 0 IDENTISCH auf der Loesungsmannigfaltigkeit {E00, E11, J=J0}
  (sympy-Zertifikat; f''-Kanal durch Differenzieren der E11-Loesung).
  E22 ist KEIN fake: als Gleichung ist es unabhaengig (f''-Kanal, Koeffizient
  -h/(2f) != 0), als Bedingung ist es eine Bianchi-Schatten — es wird
  aufuerlegt UND numerisch mit-zertifiziert (trivial branch: max|E22|
  1.4e-7, median 7.1e-11).

### (c) Exakte Reduktion (alles sympy-verifiziert)
h' = (1-h)/r + rX,  f' = f(1-h)/(hr) - fXr/h,  phi' = J0/(r^2 sqrt(fh)).
On-shell K_scalar = r^2 phi'^2/2 = J0^2/(2 r^2 f h) >= 0 (analytisch).
Die Familie ist ein 1-Parameter-Shooting in J0.

### (d) No-Go-Theorem (rigoros, alle Zahlen aus dem Archiv)
Auf jedem Flat-End-verankerten Zweig (f(r_min)=1/4, h(r_min)=1):
h <= 1 und f' >= J0^2/(2 r^3 h)  =>
f(r_core) >= 1/4 + (J0^2/4)(1/r_min^2 - 1/r_core^2).
Das Archiv erzwingt am Core-Rand (u=0.71) J0^2 = -2 r_core^4 f X = 0.7315,
also f(r_core) >= 1829 gegenueber Archiv 0.2927 — um Faktor 6.25e3
inkompatibel. Die FAMILIE (76 J0-Werte, dicht) bestaetigt: h(r_core) in
[0.9913, 1.0] fuer alle regulaeren Zweige vs. Archiv 0.3891.
Ergebnis: die EINZIGE regulaere Loesung des korrigierten Systems mit den
Projekt-Randbedingungen ist der triviale Zweig (f,h,phi) = (1/4, 1, 1).
Konvergenz-Ziele auf der vollen 3299-Punkt-Domain erreicht:
max|E00| = max|E11| = 9.9e-10 (Ziel <= 1e-6), median 1.3e-14 (Ziel <= 1e-10);
max|dJ/dr| = 5.0e-9; max|E22| = 1.4e-7. Collocation-Kreuzcheck (200 Punkte,
sparse trf): Warm-Start vom Epsilon-Zweig max|res| 2.1e-6 — konsistent.

### (e) N3/N4-Entscheid strikt nach deklarierten Regeln (fail-closed)
Auf dem trivialen Zweig gilt K_scalar == 0 IDENTISCH -> unter dem
deklarierten 1e-12 Noise-Floor -> `noise_floor_degenerate`: N4 = BLOCKED /
NOT_EVALUABLE. Kein H1/H3/NEU-Verdikt: H1 benoetigt ein nichttriviales
on-shell Hintergrundprofil; die 41-Slot-Emission ist im Chart phi' = 0
singuulaer (a6, d3, e4, a2-Teile dividieren durch phi_r; ausgefuehrter
Emissionsversuch dokumentiert), c2 = sqrt(fh) phi'(G2X/2 - h G2XX phi'^2/2)
r^2 = 0. Keine neuen Verdikt-Namen. Entscheidungsregeln unveraendert.

### (f) Nächster deklarierter Schritt
Luminaler Solve mit VARIERENDEM G4(phi)-Gesetz (a1 != 0) — der einzige Weg
zu einem nichttrivialen on-shell Hintergrund im luminalen Sektor. Erst dann
sind H1/H3 auf sauberem Boden entscheidbar. Der Pocket-Befund (K < 0 bei
x~0.30) bleibt damit OFFEN, nicht widerlegt.

### (g) RESOLUTION V4 (06.10.26): G4(phi)-Solve ausgefuehrt — H3 GEFEUERT

Status: `tools/run_luminal_background_solve_v4.py` (Artikel:
`N2_LUMINAL_SOLVE_RESULT_V4.json`, `N2_LUMINAL_V4_EPS_FAMILY.csv`,
`N3_N4_LUMINAL_DECISION_V4.json`, Stream41-CSVs fuer zertifizierte Zweige).

Gesetz (DEKLARIERT im Tool, erste Implementation): G4 = 1/2 + a1(phi-1),
G4phi = a1, G4phiphi = 0; G4X/G3/G5-Familien = 0; G2 = G2X*X (G2X=1),
G2F = 1, A0' = 0.

Architektur = V3-Struktur, verallgemeinert (sympy-verifiziert):
- E00 -> h', E11 -> f' algebraisch.  Mit a1 != 0 enthaelt E00 zusaetzlich
  phi'' (Gesetz-Kopplung): hp = alpha*phpp + beta, numerisch komponiert.
- E22 kanzeliert nach Substitution von h', f' und f'' (= totale Ableitung
  der f'-Loesung, EXPLIZITE metrische Kettenregel — dr_total traegt nur
  den Jet-Chain) zu einer phpp-LINEAREN Gleichung: die V3-Bianchi-Schatten-
  Struktur generalisiert auf das variierende Gesetz.
- PITFALL: sp.expand auf diesen verschachtelten rationalen Ausdruecken
  haengt (>13 min); cancel/together + subs bleiben sub-sekunden.
- Fail-closed: coef_phpp proportional phiphi' -> nur degenerierter Punkt
  ist phiphi'=0 exakt (trivialer Zweig, dort ist phi''=0 Loesung).

Anker (alle runtime, fail-closed):
- a1=0 reproduziert V3 exact_rhs: 5.551e-17 ueber Zufallsstates.
- eps=0 reproduziert den trivialen Zweig exakt.
- Alle Struktur-Fakten als Assertions im Lauf.

Familie: Shooting vom Flat-End (u=100) mit Haar-Amplitude eps;
164/164 Zweige erreichen den Core-Rand REGULAER (V3: keiner).

Zertifizierung: Interior-Fenster (Rand-Margin 8 Punkte, jet9d8-Stencils
einseitig am Rand — deklariert); voll-Domain-Stats separat berichtet.
Convergierte Zweige: E00/E11 int <= 1e-6, E22 int <= 1e-4.

N4-ENTSCHEID (Regeln UNVERAENDERT aus N1_N4_DECISION_RULES.json):
Auf zertifizierten nichttrivialen Zweigen ist K_scalar GLOBAL negativ
(inkl. Pocket-Band): sign(K) = -sign(a1*eps), |K| skaliert mit |eps|
(z.B. a1=-0.5, eps=-0.3: K in [-1.79e-3, -6.8e-4], Stencil-Diff 3.2e-12).
Verdikt: **H3_PHYSICAL_GHOST** — ein Even-Parity-Geist des variierenden
Gesetzes, aber NICHT pocket-lokalisiert: er ist ein Ueberall-Geist des
Zweigs, skaliert mit der Haar-Amplitude. Der historische Pocket-Befund
des Archiv-Members (K<0 NUR bei x~0.30) wird damit NICHT als artefakt
bestaetigt noch widerlegt — das Archiv-Member bleibt ein anderes Objekt
(sein K-Band ist lokalisiert; unsere Zweige sind global negativ).

Stream41-Emission auf dem nichttrivialen Hintergrund: ERFOLGREICH
(V3-Blocker `phi'=0 -> singular chart` ist auf Haar-Zweigen weg).

Ehrlichkeitsgrenzen:
- K wird numerisch ueber die V3-zertifizierte Identitaet
  K = F r/2 (f'/f - h'/h), F = 2 G4(phi), gemessen (keine geschlossene
  Form gefunden — Versuch dokumentiert, max|diff| ~ Skala).
- Fail-closed-Klausel: |K| unter max(Stencil-Diff, Residuen-Skala)
  -> NOT_EVALUABLE (getriggert fuer kleine eps).
- finite-L Health (B6/B8) und BC-Schalt offen; QNM-Claims bleiben
  geblockt bis zur Bridge-Entscheid — aber NUR noch auf der
  ghost-freien Seite (a1*eps < 0) sinnvoll fortsetzbar.
