# SSZ P5 — Eq. (4.33) / K2-prime reverse audit
**Date:** 2026-09-16

## Purpose
This audit separates three questions that had been mixed together:

1. Is `dK2/dr` from JET9D8 numerically trustworthy?
2. Does the previously transcribed Eq.-(4.33)-style identity reproduce the derivative of the literal K2 definition?
3. If not, is the failure in the L-dependent or the L-independent part?

The accepted central P5 genuine-SVT representative and the canonical analytic P5 geometry are used.

## Algebraic decomposition

Write

    K2 = N / D,

with

    N = 4 r^2 a4/f,
    D = D0 + L Q,
    D0 = f'/f - 2/r + A0' v6/(2 a4),
    Q = 1/(r h).

The previously used identity has the form

    K2' = (A'/A) K2 - B K2^2/(8 r f^2 a4),
    A = r^2 f h,

so the B required by the literal K2 definition is exactly

    B_req = (2 f^3/r) [ D' + D (A'/A - N'/N) ].

Therefore B_req is affine in L:

    B_req = B_req,0 + L B_req,L.

For the EH+genuine-SVT branch, a4=sqrt(fh)/2.  Then

    2 a4'/a4 = f'/f + h'/h,

and the L coefficient simplifies identically to

    B_req,L = d_r[f^3/(r^2 h)].

Hence the L-dependent part of the transcribed identity is not the problem.

## Numerical results

- max |a4-sqrt(fh)/2| = 2.146833e-11.
- p95 residual of `2 a4'/a4 - f'/f - h'/h` =
  3.815002e-08.
- p95 scaled residual in the L coefficient =
  3.295075e-09.
- p95 scaled residual in the L-independent coefficient =
  1.624423e+00.

JET9D8 differentiation of K2 was also compared to differentiating the
quotient N/D term by term.  The two direct derivative routes agree tightly;
see the accompanying summary table.  Thus the low-L discrepancy is not
explained by a failure of the local derivative service.

## Interpretation

The old Eq.-(4.33) discrepancy is isolated to the L-independent term.
That is precisely the term containing the electric/SVT quantities
`A0'`, `v6`, `v10`, `alpha7` and the derivative of `A0' v6/(2a4)` when
the literal K2 definition is differentiated.

This strongly links the K2-prime mismatch to the already unresolved
Eq.-131 / v6-prime source-convention issue.  It does **not** support
altering the P5 action merely to force the current Eq.-131 transcription.

The archived central member remains holonomically consistent under direct
chain-rule tests.  Therefore the safe status is:

- direct radial differentiation: numerically validated;
- L-dependent Eq.-(4.33) structure: validated;
- L-independent source-reduced shortcut: unresolved;
- finite-L low-L K positivity cannot yet be certified by mixing the two routes.

No physical low-L ghost is claimed from the direct route until the
source-reduced identity and the unreduced same-action reduction are in
one verified convention.
