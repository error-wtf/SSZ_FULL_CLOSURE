# Electric-hybrid checkpoint — 2026-09-19

## What changed

This checkpoint does **not** roll back the verified P5/perturbation infrastructure.
It changes only the global action-member decision.

### Retained

- P5 geometry and analytic center.
- 41-slot framework and JET9D8.
- Core conditioned principal/kinetic invariants.
- Generalized-psi descriptor and its H0/h1 constraint replay.
- Strong-H descriptor pullback equivalence to the established 3-field KRGSM operator.
- Singular radial DAE pencil construction.

### Rejected as final global members

1. The 2026-09-17 regional electric member: finite-L Central kinetic failure and
   incomplete common-action lower replay.
2. The 2026-09-18 global zero-vector epsilon-Y member: a local necessary Eq.85
   condition at the inner P5 light ring excludes simultaneous positive tensor F/H
   on A0prime=0.

The second rejection is now executable rather than interpretive.  Run:

```bash
python tools/audit_zero_vector_light_ring.py
```

The two resolved light rings are

```text
inner: r/rs = 1.4161606552723742
outer: r/rs = 1.4999999937556716
```

At the inner ring the zero-vector condition requires

```text
tensor_F/tensor_H = -1.638676091723501
```

whereas at the outer ring the control result is

```text
tensor_F/tensor_H = 0.9999997699724016.
```

## Electric hybrid seed

Solving the same denominator-free Eq.85 pointwise for the electric background
combination

```text
q = A0prime^2 * v8_background_identity
```

with `tensor_F=tensor_H` used only as a diagnostic target gives

```text
q(inner light ring) = 1.5636514881801942
q(outer light ring) = 1.4907001578681684e-7.
```

The full seed profile is
`data/generated/absolute_attempt_2026-09-19/ELECTRIC_HYBRID_EQ85_SEED.csv`.
This seed is **not** a final action: its background `v8` combination must be mapped
to the chosen SVT action convention and then reconstructed holonomically.

## Current active path

The target is now one on-shell electric Horndeski/SVT hybrid that combines:

- nonzero electric support through the inner light-ring region;
- Horndeski scalar/tensor controls sufficient to remove the Central negative
  finite-L scalar kinetic direction;
- one common holonomic action for background and 41-slot perturbations;
- action-level C-infinity handovers;
- the already validated generalized-psi/DAE representation in the deep Core.

No QNM calculation is promoted until this same member has a global Direct-41/KRGSM
certificate.

## Reproduction status

Software: **109/109 tests PASS** in bounded groups.

Representation audits retained as scoped witnesses:

```bash
python tools/audit_descriptor_pullback_equivalence.py
python tools/audit_radial_descriptor_pencil.py
```

They validate the descriptor/pencil machinery; they do not promote the rejected
zero-vector witness back to a global production member.
