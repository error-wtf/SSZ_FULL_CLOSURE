# ANTI-CIRCULAR FORWARD-PREDICTION LAW

**Status: BINDING for the entire spectroscopy chain of this repository.**
Adopted 2026-10-08, in force from commit adoption onward, applies to every
future preregistration, scan, catalogue and comparison artifact.

---

This project is FORWARD-PREDICTIVE and anti-circular.

Observational data MUST NOT be used to choose, tune, optimize or rescue:

- action parameters
- branch selection
- background solution
- physical DOF selection
- boundary conditions
- solver conventions
- scan domains
- angular mode L
- resonance acceptance thresholds
- frequency windows
- mass scaling
- candidate-selection thresholds

The complete theory-side prediction must be frozen first.

## Required order

```
ACTION
-> BRANCH
-> OPERATOR
-> HEALTH GATES
-> SOLVER CERTIFICATION
-> PREREGISTERED THEORY DOMAIN
-> RESONANCE SEARCH
-> FROZEN THEORY OUTPUT
-> SHA256 / COMMIT
-> ONLY THEN OBSERVATIONAL UNBLINDING
```

and NEVER:

```
data -> interesting peak -> move branch/scan/domain -> "prediction"
```

## Calibration vs. fitting

Numerical calibration is allowed; physical fitting is not.

Pöschl–Teller, square-well, RW/Leaver and similar controls may calibrate
NUMERICS (solver error, grid resolution, residual thresholds, ECS geometry).
An observational target may NOT calibrate PHYSICS.

NICER/GWOSC products may be fetched and prepared beforehand, but their
spectral contents must not influence theory construction or scan selection.

## Failures are results

If a theory prediction disagrees with observations:

**RECORD THE FAILURE.**

Do not retune the same prediction against the same data. Any changed
model/branch/domain becomes a NEW hypothesis version and must be tested
against an untouched holdout dataset or a future observation. No post-hoc
rescue is allowed to inherit the label "prediction".

Null results (e.g. `NO_CERTIFIED_RESONANCE_IN_EXPANDED_DOMAIN_V1`) are
forward predictions of full rank and are recorded as such.

## Mass handling

When an external dynamical mass measurement is required for frequency
scaling, the order is:

```
independent mass measurement
-> selection rule fixed in advance
-> M frozen
-> SSZ frequency scaled physically
-> ONLY THEN compare to the observation
```

Never: choose M so that the frequency fits. Never use Kerr masses or LVK
posteriors to move SSZ poles toward known ringdown locations. SSZ produces
its poles independently first; only afterwards may the blinded comparison
be opened.

## Mandatory comparison-artifact fields

Every comparison artifact MUST state:

- `theory_freeze_commit`
- `theory_freeze_timestamp`
- `observational_unblind_timestamp`
- `data_access_boundary`
- `preregistration_hash`
- `prediction_hash`

and MUST verify:

```
theory_freeze_timestamp < observational_unblind_timestamp
```

If this fails:

`ANTI_CIRCULARITY_GATE = FAIL`.

## Domain-expansion protocol (already in force)

Domain changes follow the preregistered ladder only (L=6 complete, then
L=12, L=20, L=42, then the eps=+0.3 branch from the original small box).
A larger scan may never be motivated by an observational peak location.
