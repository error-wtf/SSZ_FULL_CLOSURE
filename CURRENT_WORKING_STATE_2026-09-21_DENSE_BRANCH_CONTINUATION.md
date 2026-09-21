# Current Working State — Dense Strong-Field Branch Continuation

Date: 2026-09-21

## Status

This is a checkpoint, not Absolute Full Closure.

The strong-field action continuation has produced the first finite dense-grid
principal-sector branch across the stable inner light ring.

- Domain: `0.7002175543885971 <= u <= 0.7079894973743436`
- Dense rows: `297`
- Tested multipoles: `L = 6, 12, 20, 42, 110, 420, 1000`
- Dense Even kinetic gate: `PASS`
- Dense radial hyperbolicity gate: `PASS`
- Inner light ring `u ~= 0.7061346` lies inside the passed domain.

## Finite continuation member

Final smooth Horndeski null-bundle coefficients:

```
(-0.08409032140882212,
 -0.2254425096952444,
 -0.3223223184527862,
  0.09666464883889786)
```

The member was obtained by finite predictor-corrector continuation with
re-linearization on the deformed member, not by scaling one frozen tangent.

## Dense kinetic result

Worst finite value:

```
min lambda(K) = 1.2603117420917508e-05
L             = 1000
u             = 0.7048124531132783
```

No negative dense-grid kinetic rows occur for any tested L.

## Dense radial result

Worst radial characteristic:

```
min c_r^2 = 3.647320775286289e-04
L         = 1000
u         = 0.7002175543885971
```

No negative dense-grid radial rows occur for any tested L.

## Scientific scope

This checkpoint establishes a finite action-derived continuation step through
the stable inner light ring with positive Even kinetic matrix and positive
radial characteristics over the tested finite-L set.

It does NOT yet certify:

- the validated angular gate on this member,
- continuation from `u ~= 0.708` to the Horndeski-core anchor at `u = 0.715`,
- the analytic high-L principal-symbol limit,
- one final combined absolute H+SVT emitter provenance for the complete branch,
- global Direct-41 / KRGSM,
- same-operator coupled QNM,
- Absolute Full Closure.

The historical Inner coefficient stream remains rejected as a final member and
must not be used to fill the unevaluated core-side continuation.

## Primary evidence

- `data/generated/strong_field_transition_2026-09-21/TRANSITION_DENSE_H_CONTINUATION_REPORT.json`
- `data/generated/strong_field_transition_2026-09-21/TRANSITION_DENSE_H_CONTINUATION_RADIAL.json`
- `data/generated/strong_field_transition_2026-09-21/TRANSITION_DENSE_H_CONTINUATION_STATE.npz`
- `data/generated/strong_field_transition_2026-09-21/TRANSITION_DENSE_H_CONTINUATION_TOTAL41.csv`
