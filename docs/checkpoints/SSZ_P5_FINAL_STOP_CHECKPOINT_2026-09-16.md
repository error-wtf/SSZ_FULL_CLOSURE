# SSZ P5 — Final Stop Checkpoint, 2026-09-16

## Why this checkpoint exists
This file freezes the result of the Sep-16 finite-l closure work and explicitly stops further iterative representative hunting in this run.

## What is closed
- P5 frozen geometry and local sector witnesses.
- Constructive local same-action Horndeski/SVT closure.
- 41/41 unreduced coefficient language for selected sector representatives.
- Complete Zhang–Kase Appendix-A 41-slot emitter.
- Correct Appendix branch v12 = -v6/(2h).
- Correct holonomic a5 identity.
- Full algebraic constraint maps and nonzero elimination pivots on the certified domains.
- JET9D8 radial product-rule service.
- Profile-aware Euler-operator reducer, including the antisymmetric radial S block.
- Reducer high-L mass regression.
- Maxwell–Horndeski variable-G4(phi), G5=0 emitter regression against the archived strong-H 39/41 carrier. After the Sep-16 parser fixes, c6, d3 and e4 regress to numerical/JET accuracy.

## Important distinction: 41/41 CLOSED does not mean global finite-l same-action PASSED
The existing 41/41 closure means that complete unreduced coefficient tables exist for selected local/sector representatives.

The stronger statement required for a final coupled QNM calculation is:

> one center-to-infinity holonomic action member, with the handovers evaluated from that same member, whose finite-l reduced kinetic/gradient operator is healthy.

That stronger statement is NOT established by the 41/41 slot closure alone.

## Final explicit outer same-action attempt
A final outer member was built from:

1. the healthy variable-G4(phi), G5=0 Maxwell–Horndeski principal carrier reconstructed from the archived strong-H action data;
2. the resolved A0' != 0 outer background;
3. the reconstructed genuine-SVT Delta-C block, assembled before constraint elimination;
4. common auxiliary-V recanonicalization after assembly;
5. v12 = -v6/(2h);
6. v5 = c3 = e3 = 0 for the selected lower-order member;
7. holonomic a5 regenerated with JET9D8.

Structural checks PASS:
- H0 quadratic residual = 0 after common recanonicalization;
- v12 identity residual = 0;
- no operator derivatives above second order;
- R = 0;
- G symmetry residual is at floating-point level;
- S antisymmetry residual is ~1e-8;
- all common constraint pivots remain nonzero.

## Decisive finite-l result
The physical 3x3 kinetic matrix K is NOT positive definite for this explicit outer same-action member.

Representative results:

| L | min eig(K) | rows with negative eig(K) |
|---:|---:|---:|
| 6 | -2.116849e3 | 2113 / 2201 |
| 12 | -6.688485e2 | 2065 / 2201 |
| 20 | -2.363717e2 | 2025 / 2201 |
| 42 | -1.286435e1 | 1903 / 2201 |
| 110 | -1.884584 | 1576 / 2201 |
| 420 | -2.50559e-1 | 1061 / 2201 |
| 1000 | -1.06480e-1 | 629 / 2201 |

Where K is positive, the radial generalized eigenvalues remain positive in this audit, but that does not rescue a ghost-sign failure of K itself.

## Interpretation
This does NOT invalidate:
- the P5 geometry;
- the local SVT light-ring witness;
- the local Horndeski carrier;
- the constructive same-action existence theorem;
- the 41/41 coefficient completion;
- the reducer architecture.

It DOES mean that the particular explicit outer global same-action member selected in this Sep-16 run is not an acceptable final finite-l representative.

The failure is therefore a representative-selection failure at the stronger global finite-l level, not a failure of the previously established local closure statements.

## QNM status
A coupled QNM spectrum must NOT be claimed from this member. QNM remains blocked until an explicit global same-action member passes K>0 and radial/angular gradient gates at finite l.

## Stop rule
Do not continue ad-hoc tuning in this run. A future continuation should start from a deliberately formulated optimization/control problem for the remaining Horndeski/SVT transverse jet freedom, with K positivity imposed as a target, rather than by trying additional hand-picked representatives.
