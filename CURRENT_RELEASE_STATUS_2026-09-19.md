# Current release status — 2026-09-19

**Absolute Full Closure is not yet certified.**  The repository now distinguishes
three different things that must not be conflated: (1) reusable verified
infrastructure, (2) rejected concrete action members, and (3) the currently active
electric-hybrid action search.

## Software and reusable infrastructure

All **115/115 tests PASS** in bounded groups:

- 67 unit,
- 4 regression,
- 1 integration,
- 42 negative/guard,
- 1 release-auditor.

The following work remains valid and is retained independent of the final member:
P5 geometry/background, 41-slot emitters, JET9D8, conditioned Core principal
invariants, the generalized-psi descriptor, Strong-H descriptor pullback
equivalence, and the singular radial DAE pencil.

## Member decisions

### Regional electric member (2026-09-17)

**Rejected for Absolute Closure.**  Its scoped 41-slot normalization remains useful,
but the selected Central member fails finite-L kinetic positivity at
L=6,12,20,42 and its frozen c3/e3 pair is not reproduced by the same direct
action jets that reproduce the other lower-order channels.

### Zero-vector Horndeski + epsilon_Y member (2026-09-18)

**Rejected as a global production member, but retained as a useful witness.**
Its many finite-L/descriptor PASS results were real tests of those components;
they are not erased.  The stronger background test is local and decisive:

At the inner P5 light ring,

```text
r/rs = 1.4161606552723742
2 f - r f' = 0
```

so the derivative control in the necessary Maxwell-Horndeski Eq.85 vanishes.
For A0prime=0 with positive tensor_H and a4, the actual P5 geometry requires

```text
tensor_F / tensor_H = -1.6386760917...
```

and hence cannot simultaneously keep tensor_F>0 and tensor_H>0.  At the outer
light ring r/rs≈1.5 the same audit gives tensor_F/tensor_H≈1, providing a useful
control check.  The machine report is
`data/generated/absolute_attempt_2026-09-19/ZERO_VECTOR_LIGHT_RING_AUDIT.json`.

## Active path: electric Horndeski/SVT hybrid

The denominator-free Eq.85 identity also gives the first quantitative seed for
the required electric support.  If one uses the old carrier tensor target
`tensor_F=tensor_H` only as a local diagnostic target, the inner light ring needs

```text
q := A0prime^2 * v8_background_identity ≈ 1.56365148818.
```

The full pointwise seed profile is exported as
`data/generated/absolute_attempt_2026-09-19/ELECTRIC_HYBRID_EQ85_SEED.csv`.
This **is not yet a holonomic action reconstruction** and its `v8` convention must
be mapped explicitly to the selected action representation before promotion.

The active search definition is
`SSZ_P5_ELECTRIC_HYBRID_SEARCH_CANDIDATE_2026-09-19.json`.
The intended member combines a nonzero electric SVT background around the inner
light ring with general Horndeski controls for the scalar/tensor kinetic sector.
No coefficient-level transplantation is accepted as a final result.

## Remaining Absolute-Closure chain

```text
on-shell electric hybrid action
    -> common-action Direct-41 stream
    -> finite-L K / radial / angular gates
    -> action-level C-infinity Core/exterior handovers
    -> generalized-psi descriptor + DAE pencil for that same member
    -> global same-operator KRGSM
    -> coupled QNM convergence / instability search
    -> strict release PASS
```

No v1.0.0 or Absolute-Full-Closure certificate is issued at this checkpoint.

## Electric-hybrid principal feasibility checkpoint

The Central Horndeski principal-control search is now reproducible without hidden
workspace state.  The primitive basis, exact search baseline, recipe, replay and
machine audit are stored under `data/generated/absolute_attempt_2026-09-19/`.

The current recipe is **not** an Absolute-Closure member.  It is a controlled
feasibility result: the old O(1e2) low-L Central kinetic failure is reduced to
`-3.242946` (L=6), `-0.130994` (L=12) and `-0.0128774` (L=20), while L>=42 is
positive and the L=6 constraint pivots remain large (`min|Dh1|≈151`,
`min|aux det|≈442`).  This shows that the scalar kinetic obstruction is highly
controllable, but the remaining low-L edge modes should be addressed by the
joint electric/background action solve rather than by further coefficient-level
principal shaping.
