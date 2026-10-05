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
