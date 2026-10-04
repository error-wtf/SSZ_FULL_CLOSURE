# P1/K_scalar Action-Level Verification (2026-10-05, Bingsi)

## Was verifiziert wurde (Fall A vs. Fall B Vorbereitung)

Referenzquelle identifiziert, geladen und SEITENGENAU gelesen:
Kase & Tsujikawa, "Black hole perturbations in Maxwell-Horndeski theories",
Phys. Rev. D 107, 104045 (2023), arXiv:2301.10362.
Lokal: /home/error/physics/papers/Kase_Tsujikawa_2023_Maxwell_Horndeski_PRD107.104045.pdf

## Formelkette: Code vs. Paper (alle geprueft)

| Groesse | Paper | Code (repo) | Ergebnis |
|---------|-------|-------------|----------|
| F (G5=0) | Eq.(3.10): F=2G4 | F=2*G4 | IDENTISCH |
| H (G5=0) | Eq.(3.8): H=2G4+2h phi'^2 G4,X | H=2G4+2h ph^2 G4X | IDENTISCH |
| G (G5=0) | Eq.(3.9): G=H | G=H | IDENTISCH |
| a1 (G5=0) | Appendix A, S.22 (VISUELL verifiziert): sqrt(fh)[(G4phi+0.5h(G3X-2G4phiX)ph^2)r^2 + 2h ph(G4X - h G4XX ph^2) r] | gleiches Ausdruck | IDENTISCH |
| mu | Eq.(4.30): mu=2(phi'a1+2 r a4)/sqrt(fh) | mu=2*(ph*a1+2*r*a4)/sqfh | IDENTISCH (a4-Zweig separat zu pruefen) |
| P1 | Eq.(4.29): P1=(h mu/(2 f r^2 H^2)) (f r^4 H^4/(mu^2 h))' | P1=h*mu/(2 f r^2 H^2)*dY/dr, Y=f r^4 H^4/(mu^2 h) | IDENTISCH |
| Ghost-Bedingung | Eq.(4.31): K=2P1-F>0 | K=2*P1-F | IDENTISCH |

Methodik-Hinweis: pdftotext zerstueckelt die Brueche in Appendix A (1/2-Faktoren
landen isoliert in separaten Zeilen). Die endgueltige Verifikation erfolgte
VISUELL auf der gerenderten PDF-Seite 22 (pdftoppm 200dpi + Vision). Eine
fruehere Sympy-Gegenprobe mit der Text-Lesart war ein pdftotext-Artefakt,
KEIN Codefehler.

## Konsolidierter Befundstand (Pocket x~0.294-0.324, min K=-0.05586576)

1. Stencil-Robustheit: bereits auditert (5x4 ... 15x12, NEGATIVE_POCKET_STENCIL_ROBUST).
2. Routen-Robustheit: NEU (tools/recompute_kscalar_independent_route.py):
   Cubic-Spline-Ableitung und 2x-Grid-Spline weichen nur 3.8e-8 bzw. 2.4e-6
   relativ vom jet9d8-Wert ab. POCKET_ROUTE_ROBUST.
3. Formel-Treue: JETZT GEGEN DIE CANONISCHE QUELLE VERIFIZIERT (Tabelle oben).
   Die P1-Konstruktion ist eine getreue Kase-Tsujikawa-Implementierung.

## Logische Konsequenz

Die drei Erklärungsklassen "Ableitungsstencil", "Ableitungsroute" und
"Formelimplementierung" sind JETZT alle ausgeschlossen. Die verbleibenden
Hypothesen fuer die negative Tasche sind nur noch:

(H1) Die EINGANGSPRIMITIVEN (G4_background, G3, G4X, G4XX, G4phi, ... aus der
     Action-Rekonstruktion) oder die Background-Loesung selbst sind im
     x~0.30-Band inkorrekt.  [passt zum 88.46%-c2-transversal-Befund]
(H2) Der archivierte positive K_scalar-Oracle wurde mit ANDEREN Inputs oder
     einer anderen (nicht mehr rekonstruierbaren) Konvention erzeugt.
     [Oracle-Provenienzproblem, nicht Formelproblem]
(H3) Die negative Tasche ist PHYSIKALISCH: der quartische Core hat im Band
     x~0.30 einen echten Even-Parity-Ghost (K<0, Paper Eq. 4.31 verletzt).

Naechster entscheidender Schritt: (H1) pruefen, indem die Primitiven
UNABHAENGIG (symbolisch/AD) aus der eingefrorenen Action bestimmt und mit den
eingespeisten Spalten verglichen werden. Bleibt der Vergleich konsistent, ist
(H3) als ernstzunehmender Physik-Befund stehenzulassen und (H2) als
Dokumentationsluecke zu kennzeichnen.

## Status des physikalischen Claims

WEITERHIN GILT: physical_qnm_claim_allowed=false. Aber die Beweislast hat sich
verschoben: nicht mehr "vielleicht falsche P1-Formel", sondern
"Input-Primitiven-Provenienz oder echte Physik".
