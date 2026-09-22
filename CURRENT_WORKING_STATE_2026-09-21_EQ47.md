# Angular Eq47 projection and dense-action checkpoint

This is the canonical continuation state. Earlier M23-first-divergence reports
are historical and superseded by this audit. Absolute Full Closure is false.

## Direct epsilon^3 angular determinant on the dense branch (2026-09-22)

The exact rank-1 null-vector elimination (the project's authoritative route to
the first nonzero angular characteristic) was applied to all 297 trusted dense
points; the raw leading-power cubic diagnostic is superseded by this
factorized reduction. Route oracles: witness anchor 1.12e-15 (B1/B2 at the
canonical representative), rank-1 condition 2.34e-15, Eq. 4.46 residuals
1.39e-10, vector-root consistency 1.01e-15. The alpha7 kinetic identity holds
on the central genuine-SVT member (4.37e-13) and is classified as diagnostic
only on the hybrid dense member (max scaled deviation 0.80, healthy magnitudes,
structural not numerical).

Result: the scalar-gravity angular branch is negative on ALL 297 trusted
points; minimum cminus^2 = -1473.2079063534215 at u = 0.7063353338334584,
adjacent to the stable light ring. The vector branch stays positive
(min +0.048412752850779846 at the upper edge) and the third branch positive
(min +2514.1881725342664). All discriminants positive; no root degeneracies;
mode tracking via the exact vector-root identity, gaps 231.09/2514.09. The
independent finite-L route (L=1000, existing DENSE_FINITE_L_ANGULAR.csv) is
negative at exactly the same u points (Jaccard agreement 1.0): eikonal and
finite-L agree qualitatively; the instability is not a high-L artifact.

Classification: CANDIDATE_PHYSICAL_FALSIFICATION_GATE for the CURRENT member.
This is not yet a falsification of the theory family: the H+SVT transverse
completion is underdetermined by the background. Under the NO-FIT contract a
repair may NOT select jets by stability margins; an independent selection
principle must be registered before any repair attempt. Before any tangent
sensitivity statement the tangent baseline replay mismatch (1.5626839589798224e-05
against 1e-7) must be resolved. Evidence:
data/generated/dense_angular_direct_2026-09-22/.

## Tangent baseline resolved; sensitivity evaluated (2026-09-22, autopilot run)

The 1.5626839589798224e-05 replay mismatch was localized: the frozen TOTAL41
was emitted with an intermediate holonomic re-emission chain (e3 channel,
dr-composition in absolute_reemit) that predates the finalized sources of the
same commit; the stored STATE.npz carried no source hashes. The baseline was
re-frozen from the committed canonical sources using the HISTORICAL predictor
coefficients c (no re-optimization, no fitting); frozen originals are kept as
*.2026-09-21_frozen; provenance with full source hashes is in
REFREEZE_PROVENANCE.json. Post-refreeze replay error: 0.0 (guard 1e-7 PASS).
The dense K gate is unchanged (min_K = 1.2603117195246649e-05 at L=1000,
u=0.7048124531132783; DENSE_K_PASS).

With the guard satisfied, the tangent sensitivity was evaluated for the first
time (DENSE_ACTION_TANGENT_SENSITIVITY.csv): all four background-null axes
reach the angular margin at the worst point u=0.7063353338334584
(d(root)/d(axis) = -2429.83, +1074.15, +269.03, -193.94 at L=1000). The
unstable finite-L mode is a pure project-psi branch (psi weight 0.99997,
root -1481.48 at L=1000), agreeing with the direct eikonal value -1473.21 at
the same u. Controllability is therefore demonstrated, but per the NO-FIT
contract these gradients may not be used to choose a repair. The honest
end states are: (A) PHYSICAL_FALSIFICATION_GATE for the certified member, or
(B) declaring angular positivity a CONSTRUCTION_CONSTRAINT with a pre-registered,
stability-independent transverse selection principle -- which demotes the
angular gate from independent prediction in any future closure ledger.
This choice belongs to the theory owner, not the agent.

## Prospective member v2 (TRANSVERSE_ZERO): refuted at the K gate (2026-09-22)

Under the owner's stricter-than-B decision the failed member was frozen as a
falsification witness (commit 1be3551). A new prospective version was built
under the rule PROSPECTIVE_V2_TRANSVERSE_ZERO, frozen before any stability
evaluation (commit 587d776): c = 0, the unique parameter-free completion that
exercises none of the background-null transverse freedom (minimal operator
content). The background is identical to the failed member by the null-space
property (verified to 1e-12 on all background columns).

Forward-chain result: the member fails at the FIRST gate. K < 0 on the dense
branch at every L in {6,12,20,42,110,420,1000} (min_K = -18.52 at L=1000,
u=0.70799; -0.108 at L=6). The chain is terminated there per the
first-failing-layer discipline; a kinetic ghost has no downstream angular or
QNM semantics. Evidence: PROSPECTIVE_V2_K_GATE.json.

Structural conclusion (both constraints now mapped):

1. The background-forced action ALONE has a kinetic ghost: the historical
   K>0 repair was necessary. The minimal-completion hypothesis is refuted.
2. The K-repaired historical member is angular-unstable on all 297 trusted
   points (robust, two independent reductions, pure psi mode).
3. Therefore, within this background-null bundle, viable members require
   BOTH K>0 and angular c_Omega>0 as construction constraints. Any such
   member's angular stability is CONSTRUCTION_CONSTRAINT, never
   INDEPENDENT_PREDICTION. Absolute Full Closure in the strong sense (all
   gates independent) is unreachable within this bundle.

Breaking this circularity requires additional independent theory input: a
covariant completion principle that fixes the transverse jets from theory
(e.g. a fundamental Horndeski+SVT Lagrangian family in the coupling functions
G2..G5, F, rather than a stability-selected expansion) -- not more
optimization. This is precisely the point where continuing requires a new
model decision rather than debugging; per contract, work on the strong-field
member chain halts here pending that theory input.

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
