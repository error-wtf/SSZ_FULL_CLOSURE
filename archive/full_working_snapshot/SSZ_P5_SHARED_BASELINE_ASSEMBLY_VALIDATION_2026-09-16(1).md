# SSZ P5 — Shared-baseline / Background-Assembly Validation
**Date:** 2026-09-16

## Purpose

This audit distinguishes two different objects that had been conflated in earlier
intermediate checks:

1. the **actual re-solved same-action background residuals** stored as
   `metric_residual_14`, `metric_residual_15` and `JA`;
2. the older diagnostic columns `E00_H + E00_SVTg` and
   `E11_H + E11_SVTg`, which are a **naive standard-sector sum** and do not,
   in the inner handover, contain every partition-derivative / re-solved
   lower-Horndeski contribution.

Only (1) is a validity test of the selected same-action background member.

## Actual re-solved residuals

Outer handover:
- max |E00 residual| = 1.804112e-16
- max |E11 residual| = 3.226586e-16
- max |J_A|          = 1.632543e-15

Inner handover:
- max |E00 residual| = 8.881784e-16
- max |E11 residual| = 4.662937e-15
- max |J_A|          = 3.072771e-15

These are the relevant background-equation gates.

## Naive standard-sector sums

Outer:
- max |E00_H + E00_SVTg| = 2.161105e-12
- max |E11_H + E11_SVTg| = 3.203313e-12

Here the simple decomposition happens to reproduce the full residual to machine
precision.

Inner:
- max |E00_H + E00_SVTg| = 9.900747e+00
- max |E11_H + E11_SVTg| = 3.733177e+00

These large values are **not** a failure of the actual inner handover, because
the selected inner member re-solves lower Horndeski jets and contains handover
jet terms not represented by the naive two-column decomposition.  The actual
stored metric residuals remain at floating-point scale.

## Partition identities

The selected inner handover explicitly satisfies

    f2  = S_SVT * f2_full_outer_reference
    f2X = S_SVT * f2X_full_outer_reference

to

- max residual f2  = 1.665335e-16
- max residual f2X = 5.329071e-15

and `f2F ≈ S_SVT` with max residual 3.012899e-07.

The outer handover similarly has `f2F ≈ S_SVT` with max residual
9.368634e-06.

## Claim correction

An earlier statement that the *naive* `E_H + E_SVTg` total is machine-zero
through the entire inner overlap was too strong and is superseded.

The correct statement is:

**The explicitly selected inner same-action background member satisfies the
actual metric equations and vector current to machine precision; the older
two-piece standard-sum diagnostic is incomplete away from its endpoint and
must not be used as the full inner residual.**
