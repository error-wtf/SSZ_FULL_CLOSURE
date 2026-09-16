# SSZ P5 / HSVT — Vollständiger technischer Status bis zur KRGM- und QNM-Schwelle
**Stand: 16.09.2026**

## 0. Zweck dieses Dokuments

Dieses Dokument ist ein technisches externes Gedächtnis für den aktuellen
Stand der P5-Horndeski–SVT-Rekonstruktion. Es trennt strikt:

- **CLOSED/PASS**: direkt bewiesen oder numerisch regressiert,
- **CONSTRUCTIVE PASS**: Existenz innerhalb einer kontrollierten Familie bewiesen,
- **SELECTED NUMERICAL PASS**: ein konkreter numerischer Vertreter wurde geprüft,
- **OPEN**: noch nicht vollständig numerisch abgeschlossen,
- **SUPERSEDED/REJECTED**: frühere Rechenroute, die nicht mehr als physikalische
  Evidenz verwendet werden darf.

Das Endziel bleibt:

1. ein einziger center-to-infinity unreduced same-action Vertreter;
2. vollständige profile-aware Constraint-Elimination;
3. physischer finite-l Operator `(K,R,G,M)` in einer gemeinsamen Basis;
4. gekoppelte QNM-/Resonanzrechnung mit einer numerisch stabilen Spektralmethode.

---

# 1. Gefrorener P5-Hintergrund

Die P5-Geometrie ist durch die Segmentierungsvariable `Xi` primär definiert.

Mit

    C = r_s/r = 1/x

und `r_s=1` gilt

    phi = Xi(C).

Die source-locked Suszeptibilität lautet

    chi(Xi) = dXi/dC,

wobei die P5-Closure autonom als

    dXi/dC = (1-Xi)^(3/2) P5(Xi)

geschrieben wird.

Die temporale Metrik ist

    f = (1+Xi)^(-2),

und die radiale Closure wird über das Verhältnis zur GR-Suszeptibilität

    chi_GR = 1/2 (1+Xi)^3

gebildet.

Der P5-Hintergrund besitzt:
- äußeren instabilen Light Ring bei `C=2/3`, `x=1.5`;
- inneren stabilen Light Ring bei `C≈0.706135`, `x≈1.41616064`;
- horizonlosen regulären Core.

Wichtig: der historische 1.8–2.2-r_s-Hermite-Bridge ist **nicht** die P5-Feldgleichung.

---

# 2. Exakte Regimevariable kappa und E00/E11-Basiswechsel

Da

    phi_r = -C^2 chi,

folgt identisch

    kappa = h phi_r^2
          = h C^4 chi^2
          = -2 X.

Diese Identität wurde auf dem archivierten Zentralprofil bis Floating-Point-
Genauigkeit regressiert.

Für den Metric-Backgroundblock `(f2,f2X)` ist

    J = d(E00,E11)/d(f2,f2X)
      = [[-A,  0],
         [-A, -B]],

mit

    A = r^2 f,
    B = r^2 f h phi_r^2 = A kappa.

Daher ist die einfache Zeilenpermutation

    E00 <-> E11

keine Konditionsverbesserung: sie verändert die Singularwerte nicht.

Die algebraisch natürliche Basis ist dagegen

    E0     = E00,
    EDelta = E11 - E00,

mit

    d(E0,EDelta)/d(f2,f2X) = diag(-A,-A kappa).

Der schwache linke SVD-Winkel ist exakt eine Funktion von kappa,

    theta(kappa) = -1/2 atan(2/kappa^2),

auf dem kontinuierlichen Ast mit `theta -> -45°` für `kappa -> 0`.

Ferner

    dtheta/dlnC
      = [2 kappa^2/(kappa^4+4)] dlnkappa/dlnC,

und

    dlnkappa/dlnC
      = dlnh/dlnC + 4 + 2 C chi_Xi.

Damit ist analytisch gezeigt:
die natürliche Gleichungsbasis hängt direkt von `X=-kappa/2` ab; ihre
Rotationsgeschwindigkeit hängt zusätzlich von radialer Geometrie und
Segmentierungs-Suszeptibilität ab.

---

# 3. Verbindung zum SVT-Stabilitätskontrollraum

Der gleiche Parameter

    kappa = h phi_r^2

steht im background-null Hessian-Kontrollmap:

    M_SVT =
      [[0,       0, 1],
       [0,      -1, 0],
       [kappa/2,-F, 0]],

also

    det M_SVT = kappa/2.

Damit verbindet **dieselbe** Größe:
- Kondition/Ausrichtung des inversen E00/E11-Problems,
- scalar radial gradient invariant `X`,
- Rang des einfachen transverse SVT-Stabilitätskontrollraums.

Reproduzierte Mindestwerte:
- outer handover: kappa_min ~ 0.098838,
- central SVT lobe: kappa_min ~ 0.207862,
- inner handover: kappa_min ~ 0.632753.

Daher gibt es in den Handover-Bändern keinen Verlust dieses Kontrollrangs.

Richtung exaktem Zentrum gilt `phi_r -> 0`, also `kappa -> 0`.
Dies ist ein Chartwechsel des inversen Problems. Eine künstliche Division
durch kappa macht den Chart nicht physisch regulär. Der separate Taylor-Center-
Branch bleibt notwendig.

---

# 4. Regime-/Knicktest

Mehrere dimensionslose Diagnosen wurden parallel untersucht:
- max chi,
- max |D chi|,
- max |D^2 chi|,
- max |D ln chi|,
- max |D^2 ln chi|,
- max kappa,
- max |D theta|,
- max |D^2 theta|,

mit `D=d/dlnC`.

Es gibt keinen mathematisch ausgezeichneten einzigen "Knickradius", der unter
allen Definitionen identisch ist.

Robust ist ein Regimecluster ungefähr

    C ~ 0.69 ... 0.74
    x ~ 1.35 ... 1.46.

Einige Paare liegen extrem nahe:
- stärkste logarithmische theta-Krümmung vs. chi-Krümmung: ~0.14 % in C;
- schnellste Rotation vs. max |D chi|: ~0.16 %.

Andere vernünftige Definitionen unterscheiden sich deutlich stärker.
Daher ist die gültige Aussage:

**P5 besitzt einen gemeinsamen Strong-Field-Regimecluster, in dem
Segmentierungsantwort, kappa/X und natürliche Background-Gleichungsbasis
besonders stark variieren. Eine eindeutige kausale Identifikation eines
einzigen "Knickpunktes" ist noch nicht bewiesen.**

---

# 5. P0/Middle-Bridge-Vergleich

Konventionslock:

P5 verwendet

    C = r_s/r.

Das Middle-Bridge-Paper verwendete

    u_MB = GM/(r c^2) = C/2.

Daher müssen Suszeptibilitäten vor einem Vergleich umgerechnet werden.

Nach dieser Umrechnung teilen P0 und P5:
- äußere lineare Suszeptibilitätsnormalisierung `dXi/dC = 1/2`;
- innere saturierende Antwort `dXi/dC -> 0`.

Sie besitzen jedoch unterschiedliche Transition-Dynamik.

Die alte Hermite-Bridge hat eine deutlich stärkere erzwungene
Suszeptibilitätsverstärkung als P5. Sie ist deshalb strukturelle Provenienz,
aber keine numerische Bestätigung des P5-Regimes.

---

# 6. Same-action Architektur und shared baseline

Die zentrale Hybridregel lautet

    C_HSVT = C_MH + (C_SVT - C_shared),

und muss auf der **unreduced quadratic-action Ebene** angewendet werden.

`C_shared` bezeichnet gemeinsame Einstein/Maxwell/Scalar-Anteile, die nicht
zweimal gezählt werden dürfen.

Wichtig:

    Reduce(A+B) != Reduce(A)+Reduce(B)

im Allgemeinen, weil Constraint-Elimination nicht linear in den
Action-Koeffizienten ist.

Deshalb darf ein gespeicherter genuine-SVT Split-Block nicht automatisch als
vollständige unsplit Zhang–Kase-Action in sämtliche Backgroundidentitäten
eingesetzt werden.

---

# 7. Backgroundgleichungen und E00/E11-Resolves

Unabhängige statische Gleichungen:

    E00 = 0,
    E11 = 0,
    E_phi = 0,
    E_A = 0.

Im starken SVT-Bereich:
- Metric equations bestimmen on-shell `f2,f2X`,
- vector current fixiert die elektrische Branch,
- scalar equation liefert eine nichttriviale Integrabilitätsgleichung für
  `f3X`.

Das verhindert, dass perturbativ bequeme Punktwerte ohne zugrunde liegende
Action akzeptiert werden.

## 7.1 Outer handover

Der explizit re-solved outer same-action background member liefert:

    max |E00_resolved| = 1.804112e-16
    max |E11_resolved| = 3.226586e-16
    max |J_A|          = 1.632543e-15

Die ältere einfache Standardsektorzerlegung reproduziert hier ebenfalls
annähernd die volle Gleichung:

    max |E00_H + E00_SVTg| ~ 2.16e-12
    max |E11_H + E11_SVTg| ~ 3.20e-12.

## 7.2 Inner handover — wichtige Korrektur

Der explizit re-solved inner member liefert:

    max |E00_resolved| = 8.881784e-16
    max |E11_resolved| = 4.662937e-15
    max |J_A|          = 3.072771e-15.

Damit ist der tatsächliche Background-PASS eindeutig.

Die älteren Diagnose-Spalten

    E00_H + E00_SVTg,
    E11_H + E11_SVTg

dürfen hier **nicht** als vollständige Residuen verwendet werden.
Sie erreichen über den gesamten inneren Handover sogar Werte von ungefähr

    9.90 und 3.73.

Das ist kein Scheitern des re-solved Members.
Im inneren Handover entstehen zusätzliche Partitions-/Handover-Jets und die
lower Horndeski jets werden neu gelöst. Diese Beiträge sind nicht in der
naiven Zweiteilung enthalten.

Frühere Aussage, die naive inner standard sum sei überall machine-zero:
**SUPERSEDED.**

Gültig ist:
**die tatsächlichen re-solved Background-Residuals sind machine-zero.**

---

# 8. Expliziter innerer Partitionscheck

Im inneren ausgewählten Member gilt numerisch

    f2 = S_SVT * f2_full_outer_reference

mit max Fehler

    1.67e-16,

und

    f2X = S_SVT * f2X_full_outer_reference

mit max Fehler

    5.33e-15.

Außerdem

    f2F ≈ S_SVT

bis ~3.0e-7.

Damit ist sichtbar, dass die Datei selbst den genuinen SVT-Anteil als
partitionierten Anteil eines vollständigen unsplit Referenzwertes speichert.

Dies erklärt, warum ein naiver "full-SVT-background-resolver auf stored f2"
falsche Schlussfolgerungen erzeugt.

---

# 9. Holonomizität des alten Zentralmembers

Der zentrale Action-Jet-Satz wurde direkt auf Kettenregeln getestet:

    df3/dr = f3_phi phi' + f3_X X',

sowie analog für

    f3_X,
    f4,
    f4_X,
    f4_XX,
    f2,
    f2_X.

Mit dem akzeptierten JET9D8-Differentiator schließen diese Identitäten im
Inneren typischerweise auf ~1e-6 oder besser; `f4_XX` erreicht ~1e-11.

Stencil-Variationen bestätigen die Stabilität.

Damit ist die zwischenzeitliche Hypothese

    "der alte Zentralmember ist fundamental nicht holonomisch"

**REJECTED / SUPERSEDED.**

---

# 10. Higher radial jets: JET9D8

Akzeptierter radialer derivative service:

- local polynomial jet,
- window = 9,
- degree = 8,
- auf dem nichtuniformen `r=x`-Grid.

Zentrale Regression:
- `v12 = -v6/(2h)` bis ~1e-14 absolute;
- `v13` median scaled relative ~6e-6;
- `e4` isolated median ~few 1e-7;
- `e4` end-to-end ebenfalls stark regressiert.

Daraus folgt:

    v12 = -v6/(2h)

ist die aktuell akzeptierte Appendix-A-Branch.

Die frühere F3-Auswahl

    v12 = +v6/(2h)

ist **SUPERSEDED**.

---

# 11. a5-Korrektur

Die korrekte Formel ist

    a5 = a2' - a1''
         - (A0' v4/2)'
         + A0' v5/2.

Der frühere Builder ließ `-a1''` weg und ist nicht mehr autoritativ.

Für pure Horndeski, `A0'=0`:

    a5 = a2' - a1''.

Die holonomische pure-H-Neufassung wurde gegen den smooth strong-H carrier
regressiert und ist akzeptiert.

---

# 12. Eq. A2 / K2'-Konflikt: was wirklich herauskam

Zwischenzeitlich schien die Zhang–Kase-v6'-Identität A2 gegen den
gespeicherten Zentralmember zu scheitern.

Daraufhin wurden folgende unabhängige Tests gemacht:

1. analytische P5-Hintergrundableitungen statt CSV-Splines;
2. reverse factor audit;
3. direkte Holonomizitäts-Kettenregeln;
4. direkte JET9D8-Ableitung von K2;
5. unabhängige Quotientenregel für K2=N/D;
6. interne algebraische Konsistenz von A2 und Eq. 4.33.

Ergebnis:

- JET9D8 und Quotientenregel für K2' stimmen extrem gut überein.
- Zhang–Kase A2 und Eq. 4.33 sind unter ihren vollständigen Annahmen intern
  konsistent.
- Der gespeicherte P5-Zentralblock darf im Hybridprojekt nicht blind als
  vollständige unsplit SVT-Action behandelt werden.
- Deshalb ist ein A2-Rest auf einem Split-Block nicht automatisch ein
  physikalischer Widerspruch.

Zwischenzeitlich erzeugter "Eq.131/A2 repair member":
**nicht Produktionsmember; nur Existenz-/Diagnoseprobe.**

Ebenso ist eine zwischenzeitliche Vermutung `+4h` statt `+2h` in A2
**retired**; sie entstand aus einem manuellen Ableitungsfehler, nicht aus
dem Paper.

---

# 13. Warum frühere negative low-l K-Ausgaben nicht als Ghost gelten

Die constrained even action erzeugt beim Eliminieren von `h1` und den übrigen
Hilfsfeldern Radialableitungen von Koeffizienten.

Daher ist

    freeze coefficients
    -> Fourier transform
    -> eliminate constraints

nicht äquivalent zu

    keep coefficient profiles/jets
    -> eliminate constraints
    -> form symbol.

Frühere negative Shortcut-Ausgaben sind deshalb Regressiondiagnostik,
keine akzeptierte Instabilitätsevidenz.

Zusätzlich wurde später sichtbar, dass ein Teil des A2/K2'-Vergleichs auf
einem split statt unsplit Backgroundblock stattfand.

---

# 14. Exakte common even constraint architecture

Vor der Reduktion:

    (H0,H1,H2,h1,dphi,dA0,dA1,V).

Physical basis:

    chi = (psi,dphi,V)^T,

mit

    psi = H2 + (L a4/a3) h1 + (a1/a3) dphi'.

Setze

    p = a1/a3,
    q = L a4/a3,

    B     = a2 - v2 v4/(2v1),
    Cphi  = a5 + L a6 - v2 v5/(2v1),
    CH2   = a7 + L a8 - v2 v3/(2v1),
    Ch1   = L(a9 - v2 v6/(2v1)),

    Beff  = B - a3 p' - CH2 p,
    Dh1   = Ch1 - a3 q' - CH2 q.

Dann

    h1 =
      -(a3 psi' + Beff dphi' + Cphi dphi + CH2 psi + v2 V)/Dh1.

Nach dieser Eliminierung:

    J1 = b3 dot(dphi) + b4 dot(H2) + L b5 dot(h1).

Für `(H1,dA1)`:

    [[2b1, v11],
     [v11,2v10]] [H1,dA1]^T = [R1,R2]^T,

mit

    R1 = -J1/L,
    R2 = -(2v1/L) dot(V) + (v6/2) dot(h1),

und

    DeltaV = 4 b1 v10 - v11^2.

Also

    H1   = (2v10 R1-v11 R2)/DeltaV,
    dA1  = (-v11 R1+2b1 R2)/DeltaV.

Schließlich

    dA0 =
       (v1 V)'/(L v9)
       -(v8 h1+v12 H2+v13 dphi)/(2v9).

---

# 15. Neuer vollständiger Constraint-Map-Code

Am 16.09. wurde die komplette algebraische Hilfsfeld-Elimination erstmals in
einen gemeinsamen JET9D8-Code übertragen:

    ssz_hybrid_full_constraint_maps_JET9D8.py

Er exportiert lineare Maps

    h1  = H10 y + H11 y',
    H2  = H20 y + H21 y',
    H1  = T10 dot(y) + T11 dot(y'),
    dA1 = A10 dot(y) + A11 dot(y'),
    dA0 = A00 y + A01 y',

mit

    y=(psi,dphi,V).

Auf dem 4000-Punkt-Zentralprofil wurden
`L=6,12,20,30,42,72,110,210,420,1000` getestet.

Für L=6:

    min |Dh1|     = 11.987940
    min |DeltaV|  = 3.042884
    min |2 L v9|  = 0.119951

und die beiden strukturellen höheren Ableitungs-Cancellations sind

    dphi'' cancellation = 0
    h1' cancellation ~ 1e-16.

Diese Zahlen reproduzieren die bereits im Closure-Manuskript angegebenen
Constraint-Margen.

Damit ist die **algebraische vollständige Constraint-Map-Implementierung PASS**.

Die `dA0`-Map kann große numerische Koeffizienten erreichen, besonders nahe
dem inneren Rand. Das ist bei nichtverschwindendem Pivot kein Ghost-Beweis,
aber ein wichtiges Stiffness-/Normierungsdiagnostikum für den späteren
Spektraloperator.

---

# 16. Principal stability — etablierte lokale Referenz

Für den akzeptierten zentralen genuine-SVT Principal-Vertreter liegen die
archivierten exakten Zhang–Kase-Regressionen vor:

    K_even_min ~ 10.5415 > 0
    c_r,T^2   >= 1.00507
    c_r,S^2   >= 0.903745
    c_r,V^2   >= 9.45887

und coupled angular:

    c_Omega,-^2 >= 84.912
    c_Omega,+^2 >= 643.142.

Diese Aussagen sind lokale principal-symbol Aussagen des ausgewählten
Referenzvertreters.

---

# 17. Global canonical K/R/G principal export

Ein center-to-infinity kanonischer principal-characteristic Export mit
45,617 Punkten existiert.

Er ist ein dokumentierter globaler K/R/G Principal-Vertreter.

Er darf jedoch nicht ohne Prüfung mit einem beliebigen neuen finite-l M-Block
als ein finaler same-basis physical KRGM-Operator bezeichnet werden.

Grund:
- basis/convention consistency,
- shared-baseline assembly,
- lower-order member selection,
- profile-aware product-rule terms.

---

# 18. Massmatrix M: was bereits geschlossen ist

Der direct-action High-L-Audit wurde korrigiert.

Entscheidend:

    m5_minus = a4 v13 - a6 v6.

Außerdem existieren zwei verschiedene `m1`-artige Kombinationen:
- `m1_plus` im c4-Sektor,
- `m1_minus` im dA0/v9-Schur-Sektor.

Sie dürfen nicht zusammengelegt werden.

Mit JET9D8 konvergieren insbesondere:

    M12,
    M22/L,
    M23

sehr stark gegen die Action-derived High-L-Identitäten
(median errors im Bereich ~1e-8 bei L=1e9).

Damit sind die früheren spline-basierten negativen Mass-/Angular-Ausgaben
**SUPERSEDED**.

Was noch fehlt:
ein einziger finaler finite-l M(r,L)-Block, der gemeinsam mit K,R,G aus genau
demselben globally assembled selected member erzeugt wird.

---

# 19. Aktueller KRGM-Status

## CLOSED / PASS

- P5 background geometry.
- constructive same-action existence through outer and inner overlaps.
- background-null control ranks.
- nondynamical pivot invertibility.
- selected outer background handover residuals.
- selected inner background handover residuals.
- JET9D8 higher radial jets.
- Appendix v12 minus branch.
- holonomic a5 correction.
- complete algebraic constraint maps.
- central High-L direct-action M regressions.
- global canonical principal K/R/G representation.

## CONSTRUCTIVE / NON-UNIQUE

- global finite-l KRGM operator existence.

## NOCH NICHT FINAL NUMERISCH EINGEFROREN

- ein einziger center-to-infinity 41/41 **same convention / same basis /
  same lower-order member** coefficient stream;
- daraus ein physischer finite-l `(K,R,G,M)`-Dump für ausgewählte l;
- coupled QNM spectrum.

Das ist keine neue Existenzlücke.
Es ist jetzt eine **representative-level assembly/reducer reproducibility
task**.

---

# 20. Der jetzt verbleibende Reducer-Schritt

Die neue Constraint-Map-Schicht muss in den profile-level unreduced quadratic
action kernel eingesetzt werden.

Wichtig:
- kein frozen coefficient symbol;
- alle radialen Produktregeln behalten;
- alle IBP-Terme behalten;
- JET9D8 für coefficient derivatives;
- shared-baseline pieces nur einmal;
- `v12=-v6/(2h)`;
- korrigiertes a5;
- ausgewählte `v5,c3,e3` action-realisiert und global smooth.

Danach erhält man

    S^(2) = 1/2 ∫dt dr [
       dot(chi)^T K dot(chi)
       + dot(chi)^T R chi'
       - chi'^T G chi'
       - chi^T M chi
    ].

Jeder Hilfsfeldschritt ist formal ein Schur-Komplement, aber die
profile-dependent field transformations erzeugen zusätzliche derivative terms.
Genau diese müssen im Produktionsreducer enthalten sein.

---

# 21. Acceptance gates für den finalen finite-l KRGM-Dump

Vor QNM müssen mindestens folgende Regressionen PASS sein:

1. **Constraint pivots**
   - Dh1 != 0
   - DeltaV != 0
   - 2 L v9 != 0

2. **Central exact K**
   - exakte no-ghost Kombination reproduzieren.

3. **Radial characteristics**
   - tensor/scalar/vector Reference branches reproduzieren.

4. **High-L M**
   - M11,M12,M13,M22,M23,M33 limits reproduzieren.
   - insbesondere die korrigierten M12/M23 routes.

5. **Pure-sector limits**
   - pure Horndeski/Maxwell-Horndeski limit;
   - pure genuine-SVT limit.

6. **Radial derivative convergence**
   - window/degree variation JET service.

7. **Basis consistency**
   - alle transformation derivatives mitgeführt.

8. **Shared baseline**
   - keine doppelte Einstein/Maxwell/Scalar-Zählung.

Erst nach diesen Gates wird ein finite-l KRGM-Dump als final physischer
selected member bezeichnet.

---

# 22. QNM-Status

Die asymptotische outgoing Jost-Seite ist bereits stark validiert.

Eine high-order Jost expansion wurde im weak exterior transportiert; der
Order-10 Fehler liegt ungefähr bei

    2.43e-12.

Nicht akzeptiert:
- inward shooting gedämpfter Moden entlang der reellen Radialachse.

Grund:
exponentiell wachsende contamination; Kandidaten verschoben sich beim
Outer-Radius-Test.

Produktionsmethoden:

    exterior complex scaling (ECS),
    continued fraction,
    compactified Jost / spectral determinant.

Für den gekoppelten 3-Kanal-Operator wird nach finalem KRGM ein komplexer
Resonanzdeterminant aufgebaut.

Vor einer Aussage über `Im omega` müssen geprüft werden:
- radial resolution convergence,
- complex-scaling angle convergence,
- outer/domain convergence,
- basis normalization,
- branch tracking,
- determinant/root convergence.

Es wird **kein** Stabilitätsclaim aus einem einzelnen numerischen Root gemacht.

---

# 23. Was heute ausdrücklich verworfen/korrigiert wurde

## SUPERSEDED 1
`v12=+v6/(2h)` als Appendix branch.

Aktuell akzeptiert:

    v12=-v6/(2h).

## SUPERSEDED 2
a5 ohne `-a1''`.

## SUPERSEDED 3
Spline-basierte negative angular/mass scans, die gleichzeitig bekannte High-L
Regressions verletzten.

## SUPERSEDED 4
"alter Zentralmember ist fundamental nicht holonomisch".

Direkte chain-rule tests sprechen dagegen.

## SUPERSEDED 5
A2 müsse durch eine einfache +4h statt +2h Änderung korrigiert werden.

Das war ein manueller Derivationsfehler.

## SUPERSEDED 6
Naive inner standard-sector sum sei durch den gesamten inner handover das
vollständige Backgroundresiduum.

Nein:
nur die tatsächlich re-solved residual columns sind dort vollständig.

## NICHT AKZEPTIERT
Experimenteller Eq.-A2/131-"repair member" als Produktions-QNM-member.

Er bleibt nur Diagnose-/Existenzprobe.

---

# 24. Gesamturteil des technischen Audits

Der aktuelle Stand unterstützt **keinen bekannten lokalen No-Go** der P5
same-action Konstruktion.

Die stärksten Ergebnisse sind:

1. Der Background und die same-action Handover sind konstruktiv bzw. explizit
   kontrolliert.
2. Der Regimeparameter `kappa=-2X` besitzt eine reale mathematische Doppelrolle
   im Background-Inverseproblem und im perturbativen Kontrollrang.
3. Die empfindlichen radialen Jets sind mit einem regressierten derivative
   service beherrscht.
4. Alle nondynamischen Even-Hilfsfelder sind jetzt als profile-aware lineare
   Maps maschinenimplementiert.
5. Die Massmatrix-High-L-Struktur ist action-level regressiert.
6. Die verbleibende Arbeit ist keine neue Geometrie- oder Existenzsuche,
   sondern die **einheitliche vollständige profile-level Substitution in den
   unreduced quadratic kernel**.

Die korrekte aktuelle Pipeline lautet damit:

    selected same-action action jets
        ->
    shared-baseline-safe 41-slot stream
        ->
    JET9D8 full constraint maps
        ->
    profile-level product-rule / IBP reduction
        ->
    finite-l physical K,R,G,M
        ->
    regression gates
        ->
    frozen KRGM member
        ->
    ECS / compactified-Jost coupled QNM spectrum.

---

# 25. Nächster konkrete Rechenschritt

Nicht:
- neues P5 fitten,
- outer/inner existence neu beweisen,
- A2 künstlich erzwingen,
- reduced matrices per Hand zusammenkleben.

Sondern:

**Den Termkernel `ssz_hybrid_unreduced_even_kernel.py` auf Profilebene mit den
neuen vollständigen JET9D8-Constraint-Maps substituieren und automatisiert die
bilinearen Koeffizienten K,R,G,M inklusive sämtlicher coefficient derivatives
und Integrationen durch Teile sammeln.**

Das ist der letzte große numerische Reducer-Schritt vor dem tatsächlichen
QNM-Problem.
