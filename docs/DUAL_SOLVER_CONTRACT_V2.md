# DUAL_SOLVER_CONTRACT_V2 (Phase S, nach dem T-Entscheid)

Zwei unabhängige Produktionssolver müssen denselben komplexen Pol finden:

  S1  compactified Jost / spectral determinant
  S2  exterior complex scaling (ECS)

## Akzeptanzachsen (MODE_CERTIFICATE_V2, globale Achsen)
1. Pole-Übereinstimmung: |Δω_R|, |Δγ| unter Solver-Konvergenzgrenzen
2. Radiale Auflösung: Stabil unter Gitterverfeinerung (Richardson)
3. Äußere Domain: R_max-Extrapolation konvergiert
4. ECS-Winkel: Pol ortstreu unter θ_Variation
5. Basiswechsel: gleicher Pol in zwei Diskretisierungsbasen
6. Root-Continuation: derselbe Pol unter Parameter-Tracking vom
   Bekannten (Schwarzschild-Anker) zum Member

## Pfad
Jeder Solver reportet {ω_R, ω_I, method, discretization, convergence}
als PREDICTED_MODE_CATALOG-Erweiterung; Zertifikate erst nach
beidseitiger Übereinstimmung.  Die gesunde lokale Zone (v3) bleibt
unberührt — dieser Vertrag aktiviert NUR den globalen Pfad.

## Trigger
Freigegeben erst durch den T-Entscheid (N1–N4): positive unabhängige
K_scalar-Bestätigung ODER Member-Verwerfung/-Reparatur.
