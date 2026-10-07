# V4 Physical Boundary Conditions — 3-DOF production branch

**Branch:** (a1, eps) = (-0.5, -0.3) — healthy 3-channel
production branch (see `V4_FULL_3DOF_BRANCH_HEALTH_V1.json`).
**Physical basis:** Psi = (psi, dphi, V)^T — three dynamical channels.
**Operator:** `V4_PHYSICAL_RESONANCE_EXPORT_V2_3DOF.npz`
(sha256 `ff2087354a78826a021a4ffe0c6d18f0…`).

## What this supersedes

All boundary-condition notes written for the old (+0.5, −0.3) branch and
its 2-channel (psi, V) reduction. That reduction is superseded (see
`V4_DOF_ARCHITECTURE_REAUDIT_V2.json`): dphi is dynamical, the physical
basis has THREE channels. Schwarzschild-style ingoing-horizon conditions
are inapplicable: the branch is horizonless with a regular interior.

## Measured asymptotic state (r → 60)

| quantity | value |
|---|---|
| f_inf | 0.24617405 |
| h_inf | 0.99945765 |
| phi_inf | 0.98746435 |

## Conventions (MANDATORY for the catalogue)

| convention | definition |
|---|---|
| time | T = sqrt(f_inf) * t (sqrt f_inf = 0.49615929) |
| tortoise | R_star = sqrt(f_inf) * r_star |
| frequency | Omega_inf = omega_t / sqrt(f_inf) = 2.01548174 * omega_t |

The catalogue reports **Omega_inf** in T-time; raw `omega_t` is stored
separately. The old factor 1.9953 belongs to the superseded +0.5 branch.

## Inner boundary condition

Regularity of (psi, dphi, V) at the inner window edge r = 0.05, matching
the smooth integrated background (f_min = 0.2482,
h_min = 0.9984, all finite).

## Outer boundary condition

Outgoing waves in every channel against the continued exterior:

    Psi_ch ~ A_ch * exp(+i omega_t r_star),   ch in {psi, dphi, V}

with the common metric-driven tortoise coordinate
r_star = ∫ dr/sqrt(f h) (tail slope 2.015726). The coupling between
channels decays with the 1/r geometry tail; the leading asymptotic
fundamental matrix is diagonal.

## Gate consumers

`COUPLED_RESONANCE_SOLVER_V1` — C-R1 (Jost/ECS agreement) through
C-R5 (matching-radius stability) consume THIS artifact. Any solver run
bound to a different BC artifact is invalid for this branch.
