# H1-Test Fortschritt: Emitter-Domain-Verletzung im Pocket-Band (2026-10-05)

## Was geprueft wurde (H1 = Input-Primitiven-Provenienz)

Frage: Sind die G-Spalten (G4X, G4XX, G3X, G2X, ...) des eingegebenen
Quartic-Core-Action-Profils konsistent mit dem dokumentierten
Emitter-Domainvertrag des Projekts?

## Befund (aus dem Repo selbst, nicht neu interpretiert)

1. tools/audit_production_sources.py (vorhanden, ausgefuehrt) meldet fuer
   ALLE drei F1b-Quelltabellen `jets_outside_luminal_emitter_domain`:

   - F1b_G4XX_transverse_core_candidate: G4X ~ -2.6e3, G4XX ~ 1.1e13,
     G4phiX ~ -1.5e10
   - FULL_G5_subcore: zusaetzlich G5X ~ 3.4e6
   - FINAL_Cinf_handover: gleiche Groessenordnung

2. Der dokumentierte Emitter (src/ssz_p5_mh_luminal_g4phi_emitter_JET9D8)
   ist explizit auf
       G4 = G4(phi), G4X = 0, G5 = 0
   festgelegt und VERWEIGERT per ValueError das Emiten mit G4X != 0
   (numerical_policy exact_zero_abs).

3. d.h.: Die Produktionstabellen enthalten Spalten, die der projektgelockte
   Emitter-Domain NICHT zulaesst. Der Emitter wirft bei diesen Eingaben eine
   Exception; die Tabellen sind also NICHT aus diesem Emitter entstanden.

4. Im POCKET-BAND (x 0.294-0.324) ist G4X = -2384..-2529, G4XX = 2e11..4e11,
   G3X = -3e8..-5e8, G2X = -1e5..-1e6. G4_background variiert (54 unique
   Werte) um 0.5001 - also nicht strikt konstant.

5. audit_production_sources.json selbst dokumentiert:
       global_direct_generation: "NOT_ESTABLISHED"

## Logische Konsequenz

Die negative K_scalar-Tasche kann JETZT mit einer konkreten, repo-internen
Provenienzluecke erklaert werden KANDIDAT: Die quartischen Primitive
(a1, c4, F, H) werden in mh_action_primitives.quartic_g5zero_primitives aus
G4X/G4XX/G4phiX-Termen mit aufgebaut (z.B. a1 enthaelt den Term
2 h phi (G4X - h G4XX phi^2) r). Wenn diese G4X/G4XX-Spalten
a) nicht zum projektgelockten luminalen Emitter-Domain gehoeren und
b) global_direct_generation NOT_ESTABLISHED ist,
dann ist die Tasche hoechstwahrscheinlich ein Artefakt von Action-Jets, deren
Herkunft ausserhalb des gelockten Domains liegt und die den per Contract
geforderten Term G4X=0 verletzen.

Die drei Hypothesen (H1 Primitive-Herkunft / H2 Oracle-Herkunft /
H3 echte Physik) sind damit KONKRETISIERT:
  H1 ist jetzt die bevorzugte Arbeitshypothese mit konkretem Mechanismus.
  H3 (echter Ghost) waere erst NACH einer sauberen, im gelockten Domain
  regenerierten Action-Kette wieder zulaessig.

## Naechster Schritt (abgeleitet)

Regeneriere den quartischen Core KOMPLETT im luminalen Domain
(G4=G4(phi), G4X identisch 0, G5=0) mit dem vorhandenen Emitter und dem
vorhandenen c2-Ziel, und berechne K_scalar auf dieser regenerierten Kette:
  - K verschwindet / wird positiv  => Tasche war Domain-Verletzungs-Artefakt
    (H1 bestaetigt, H2 als Dokumentationsluecke).
  - K bleibt negativ                => H3 lebt und ist jetzt sauber begruendet.

Das ist genau tools/build_core_quartic_action_direct41.py +
tools/regenerate_core_quartic_action41.py mit erzwungenem Domain-Check.

## Status

physical_qnm_claim_allowed=false bleibt.
Aber: Die Beweislast hat sich VERSCHOBEN. Der empfohlene naechste
Rechenlauf ist die Domain-konforme Regeneration, nicht neue Physik.
