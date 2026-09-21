# Current Working State — Projected Strong-Field Controllability — 2026-09-21

Absolute Full Closure is **not certified**.  The active proof object is the
single action-first strong-field transition on `0.70<=u<=0.715`.

## New result

The pre-continuation local controllability question has now been answered on the
part of the right transition for which a certifiable absolute Electric-Hybrid
action seed exists.

For action coordinates `p`, background Jacobian `B=dE/dp`, and the weakest
finite-L kinetic mode gradient `g_L`, the physical tangent is `ker(B)`.  At
`u~=0.7013`, `0.7030`, the inner light ring, and `0.7080`, the projected
gradients are nonzero and a single local action-tangent direction improves the
weakest kinetic mode simultaneously for `L=6,20,42,1000`.  The same statement
remains true when the present algebraic scalar-identity residual is included as
an additional tangent row.

At the inner light ring the common-gain result is finite-difference stable:
using action perturbation steps `2e-5`, `1e-5`, and `5e-6` changes the
scalar-identity-constrained common gain only at the few-parts-in-1e6 relative
level.

## Horndeski a1 clarification

The raw positive primitive `a1` sensitivity survives projection onto the full
local quartic-Horndeski background-null tangent spanned by
`(G4_XX,G3_X,G2,G2_X)`.  Thus `a1` is a genuine on-shell tangent clue, not only
an off-shell coefficient direction.

It is **not** promoted as a one-dimensional transition repair.  A compact
`a1`-maximizing null bundle has an amplitude-converged linear response, but the
high-L sign depends on the radial localization width.  This is exactly the
reason the next predictor is a **joint Horndeski+SVT null-bundle feasibility
solve**, not an `a1` patch.

## Deliberate boundary of the result

`u=0.712` remains unevaluated in the projected-controllability audit because the
current accepted absolute Hybrid action seed stops below `u=0.71`.  The old
Inner coefficient member is rejected as a final member and is not used to fake
an absolute continuation seed.

## Next mathematical step

Construct the smallest joint H+SVT tangent bundle that simultaneously satisfies

```
B dp = 0,
d lambda_K(L) > 0,
```

for the required modes while adding radial and validated-angular inequality
sensitivities.  Use that tangent only as the predictor.  Then solve the full
nonlinear background corrector, absolutely re-emit `C41`, and re-evaluate the
physical gates before advancing in `u`.

## Primary evidence

- `data/generated/strong_field_transition_2026-09-20/TRANSITION_PROJECTED_CONTROLLABILITY.json`
- `data/generated/strong_field_transition_2026-09-20/TRANSITION_HORNDESKI_A1_PROJECTION.json`
- `data/generated/strong_field_transition_2026-09-20/TRANSITION_PROJECTED_CONTROLLABILITY_CONVERGENCE.json`
- `docs/STRONG_FIELD_TRANSITION_CONTRACT_2026-09-20.md`

The previous 2026-09-20 release and working-state files remain historical
snapshots and are not rewritten as if they had contained these later results.

## Software verification for this checkpoint

- Unit: `77/77 PASS`.
- Negative + regression + integration + release: `52/52 PASS`.
- Complete grouped matrix: `129/129 PASS`.
- No test tolerance was relaxed.
