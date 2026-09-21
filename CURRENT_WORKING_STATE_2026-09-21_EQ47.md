# Angular Eq47 projection and dense-action checkpoint

This is the canonical continuation state. Earlier M23-first-divergence reports
are historical and superseded by this audit. Absolute Full Closure is false.

## Derived and numerically tested

The Eq46 common K/M ratio and the pure kinetic identity imply the mixed identity
with M33 on its right-hand side. The literal printed K33 relation differs by
`(R-1) r h phi' K22 K33`. The symbolic and numerical checks pass on
`0.62 <= u < 0.70`. The reference scaled residual for this prediction is
`2.52e-14`; the maximum across the trust region is `8.09e-9`.

The primary source is https://arxiv.org/html/2404.11910v3, HTML equations 84–85
(the user's Eq4.46–4.47 labels). Classification:
`PUBLISHED_REDUCED_SHORTCUT_CONFLICT`; this is not an official erratum.

Cancellation-aware selection finds epsilon cubed at the genuine-SVT reference.
At u=0.6500150037509377 the action roots are
`(-164.97552195395, 1.00000000000, 4218.66284407651)`.
The independently projected null-sector quadratic gives
`B1=4053.68732212255`, `B2=-695976.104649254`.
Finite-L results converge at L=1e3,1e4,1e5,1e6; maximum scaled root errors are
approximately `0.380, 0.0390, 0.00391, 0.000391`.

The coupled derivation is Q=(ms-z*ks)(mg-z*kg)-(msg-z*ksg)^2. Thus
q2=ks*kg-ksg^2, B1=(ms*kg+ks*mg-2*msg*ksg)/q2 and
B2=(ms*mg-msg^2)/q2. CSV outputs separate the action determinant, projected
quadratic, printed shortcut comparator and m5-minus/m5-plus routes. The
Eq47-projected reconstruction changes the identity at the determinant level;
it does not edit a printed B1/B2 factor by hand. Any remaining discrepancy is
retained. The alpha7 kinetic identity also passes (maximum scaled error 8.42e-12).

## Dense branch and remaining limits

The 297 existing trusted points on 0.7002175544 <= u <= 0.7079894974 were
evaluated with all original radial guards retained. Finite-L angular positivity
fails for the tested set 6,12,20,42,110,420,1000. The formal dense Laurent
selection is nonuniform and remains unresolved; its extreme roots are diagnostic
values, not certified physical limits. Finite-L results are exported separately.

The known action tangent bundle fails the baseline replay guard at 1.56e-5,
above the 1e-7 criterion. Its sensitivity/repair is therefore held. Coordinate
eigenvector weights are provided but cannot identify invariant physical mode
fractions across differential basis transformations.

`EXACT_GM_GHS_REPLAY_NOT_PRESENT_IN_INPUT_CHECKPOINT` remains explicit.
No synthetic background has been substituted for that missing oracle.

Continuation 0.708 -> 0.715: **HELD**. Global QNM/stable trapping: **HELD**.
`ABSOLUTE_FULL_CLOSURE = false`.

## Reproduce

```bash
pip install -r requirements.lock
OPENBLAS_NUM_THREADS=1 python tools/run_angular_universal_pipeline.py
pytest -q tests/unit/test_eq47_projection.py tests/unit/test_angular_laurent_series.py tests/unit/test_angular_universal_provenance.py
python tools/release_manifest.py --check
```

The method pipeline may pass while the physical dense-angular gate fails.
Its machine-readable result is FULL_CLOSURE_WORKING_STATUS_2026-09-21_EQ47.json.
Full baseline testing retained one failed historical Electric-Hybrid recipe
comparison (4.20667736755e-7 versus 5e-8); no tolerance was weakened.

Test groups: 51 non-unit tests passed; the original unit batch had 81 passes
and one failure. Two additional regressions passed in the final seven-test
Angular batch: aggregate unique tests 134 PASS / 1 FAIL. Syntax checks cover
all added or modified Python files.

The comparator called `literal_shortcut` evaluates the printed B1/B2 expression
on action-extracted coefficients; it is not a complete independent replay of
every printed mass identity. The projected m5 routes also use printed subleading
mass identities. Their remaining disagreement with the primary determinant is
open. Attribution of the full discrepancy to Eq47 alone has NOT been established.
