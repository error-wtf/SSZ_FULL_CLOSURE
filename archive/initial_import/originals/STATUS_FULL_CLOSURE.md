# SSZ P5 / Horndeski + U(1)-SVT — Release Status

**Release:** 2026-09-16

## Frozen verdict

- **FULL CONSTRUCTIVE CLOSURE: PASS**
- **DIRECT GLOBAL KRGM EXPORT: OPEN IMPLEMENTATION GATE**
- **QNM numerical spectrum: do not publish until the direct regenerated KRGM gate is PASS**

The production action member is the globally patched on-shell Horndeski P5 representative plus

\[
\Delta f_2=10^{-2}Y,
\qquad
Y=\nabla_\mu\phi\nabla_\nu\phi F^{\mu\alpha}F^\nu{}_{\alpha},
\]

on the zero-vector background branch `A0prime=0`.

This deformation is exactly background-null but perturbatively nontrivial. Its global vector kinetic/radial margin is

\[
Z_A=1-2\kappa\epsilon_Y,
\quad \kappa=h\phi'^2,
\quad Z_{A,\min}=0.9872911806592305>0.
\]

The direct-global-KRGM gate is deliberately separate: a repository must regenerate all needed lower-order coefficient profiles from the covariant action and then reduce them. It must not obtain a PASS by concatenating historical rounded CSV coefficient tables.
