# KORREKTUR zur H1-Interpretation (2026-10-05, nachfolgende Pruefung)

## Der Faktbestand des Decision-Tests bleibt unveraendert gueltig

tools/kscalar_domain_decision_test.py Ergebnisse (reproduzierbar):
  original jets:            K_min(band) = -0.055866, 54 negative Zeilen im Band
  G4-Jets nulliert:         K_min(band) = -0.006822, 36
  G3X nulliert:             K_min(band) = -0.062290, 49
  luminal (alle genullert): K_min(band) = +0.000504,  0

## Was die Weiterpruefung aber zeigt (GLOBAL, nicht nur Band)

Die luminal-nullierte Kette CORE_QUARTIC_LUMINAL_CLEAN_41.csv ist GLOBAL
NICHT gesund:
  - 888 von 3299 Zeilen K_scalar < 0
  - negatives Gebiet x = 0.337 .. 1.011 (breitbandig, K_min_global = -0.0739)
  - das POCKET-Band (0.29-0.33) ist sauber (0 negative Zeilen)

Die ORIGINAL-Kette ist dagegen NUR im Pocket-Band negativ (54 Zeilen,
x 0.294-0.324, sonst positiv).

## Korrekte Interpretation

Das Nullsetzen der Jets ist KEINE "Befreiung derselben Theorie von einem
Fehler", sondern die Wahl einer ANDEREN Theorie (G3X=0 aendert die
dynamischen Terme massiv). Diese jet-freie Theorie hat ein VOLLIG ANDERES,
breitbandiges negatives Gebiet jenseits x>0.337.

Der Pocket-Fund bleibt ein Jet-Effekt (faktisch korrekt). ABER die Aussage
"H1 bestaetigt, Tasche ist Artefakt und die Kette kann sauber freigeschaltet
werden" ist ZU STARK und wird HIERMIT ZURUECKGENOMMEN.

Neue, ehrliche Lage:
  1. Das spezifische Pocket bei x~0.30 ist an out-of-domain Jets gekoppelt
     (faktisch stabil).
  2. Weder die Original-Jet-Kette (Pocket) noch die jet-nullierte Kette
     (breitbandiges negatives Gebiet) ist global gesund.
  3. Es ist BISHER NICHT gezeigt, dass eine im gelockten luminalen Domain
     VOLLSTAENDIG regenerierte Action-Kette (Emitter + c2-Ziel + konsistente
     Background-Loesung, nicht "Spalten nullieren") K ueberall positiv macht.
     Der Unterschied: Echte Regeneration loest die gekoppelten
     Hintergrundgleichungen IM Domain; einfache Spalten-Nullierung erzeugt
     eine inkonsistente (u.U. off-shell) Kette.

## Was davon unangetastet bleibt

- Formel-Kette Kase-Tsujikawa: korrekt implementiert (separat verifiziert).
- Stencil-/Routen-Robustheit des recomputierten K: unangetastet.
- Der archivierte positive Oracle: weiterhin unerklärt (H2 bleibt offen).
- physical_qnm_claim_allowed = false bleibt.

## Konkreter Weg (F.3, praezisiert)

Die echte Pruefung ist NICHT "Jets nullieren", sondern:
  N1. Emitter (luminal, G4=G4(phi)) auf das c2-Ziel ansetzen;
  N2. Hintergrundgleichungen im luminalen Domain NEU loesen (f, h, phi
      konsistent mit G3=0, G4X=0);
  N3. Darauf quartische Primitive + 41-Slot-Emission;
  N4. Erst dann K_scalar/finite-L Health bewerten.
Erst N4-Ergebnis entscheidet H1-gegen-H3 sauber.
