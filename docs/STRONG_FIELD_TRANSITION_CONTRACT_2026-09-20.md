# P5 Strong-Field Transition Contract — 2026-09-20

## Purpose

This document freezes the anti-target-fitting closure architecture.  It replaces
"repair the historical Inner coefficient member" as the active final-closure
strategy.  Historical artifacts remain regression and provenance material; they
are not deleted or relabelled as physical production members.

## 1. Proof objects must not be conflated

A reproducible/normalized 41-slot stream is not by itself a physical closure
certificate.  A final physical member must be produced by one covariant action,
carry the target background on shell, pass all required local stability and
hyperbolicity gates, and be re-emitted absolutely from that same action.

The historical `CENTRAL_DIRECT_41=PASS` remains valid only in its declared
scope: authoritative coefficient import and common-convention normalization.
It is not promoted to `CENTRAL_PHYSICAL_MEMBER=PASS`.

## 2. The physical right transition is one interval

The final strong-field solve is

```
0.70 <= u <= 0.715
```

with

- left anchor: the healthy electric-hybrid bulk member at `u=0.70`;
- right anchor: the pure-Horndeski core at `u=0.715`;
- stable inner light ring: `u_LR = 0.7061345733` inside the same solve.

The historical seam `u=0.71` is **not** a physical matching surface and must not
introduce an independent action or extra matching conditions.

## 3. Action first; coefficients are outputs

The unknowns are smooth action functions/jets.  The perturbation coefficients
`a_i,b_i,c_i,d_i,e_i,v_i,alpha_i,...` are emitted quantities:

```
A_candidate -> background EOM -> C41_absolute -> reduced physical operator.
```

Historical coefficient profiles are diagnostics only.  They are never hard
interior targets in the final transition solve.

Delta/Jacobian response operators may be used for Newton steps, rank tests, and
sensitivity diagnostics, but a certificate may never be based on
`C_historical + delta C`.  Every promoted member must be independently and
absolutely re-emitted from the candidate action.

## 4. Background-preserving directions are preferred

Where possible, continuation should first use transverse/background-null action
directions.  Example:

```
Delta G2 = 1/2 q(phi) [X-X_b(phi)]^2
```

which vanishes with its first background variation while retaining a nonzero
`G2_XX` second variation.  Such directions can alter perturbation physics while
preserving an already-good P5 background.

A direction that has zero sensitivity to the failing physical mode is not an
optimization target.  The current audit establishes that the existing `G2_XX`
lift is a kinetic null direction for the right-transition obstruction.

## 5. Electric support and nondegeneracy

The electric branch must remain nontrivial through the stable inner light ring
and may smoothly approach the pure-Horndeski core only after it:

```
A0'(u_LR) != 0,
A0'(0.715) = 0.
```

The vector equation `J_A=0` with `A0'!=0` must additionally remain on a
nondegenerate constitutive root.  The diagnostic

```
chi_E = partial J_A / partial A0'
```

must therefore remain finite and bounded away from zero on the electric branch.
It is a branch-regularity gate, not a replacement for perturbative stability.

## 6. Required physical gates

After constraint elimination the physical transition must satisfy, for all
relevant propagating modes,

```
K > 0,
c_r^2 > 0,
c_Omega^2 > 0.
```

The finite training set is

```
L = 6, 12, 20, 42, 110, 420, 1000.
```

This finite set is not an eikonal proof.  Final closure additionally requires an
analytic high-L analysis of the unreduced principal symbol / canonically
normalized characteristics.

Eigenmodes must be tracked by adjacent eigenvector overlap so that crossings or
avoided crossings are not hidden by naive sorting of eigenvalues.

## 7. Angular gate ordering

The experimental angular evaluator is not a certificate until it reproduces:

1. the historical exact positive SVT angular oracle;
2. the Maxwell-Horndeski limit;
3. only then the electric-hybrid / transition candidate.

A negative number from an evaluator that fails item 1 is a reducer-regression
failure, not a physical instability claim.

## 8. Scalar background gate

The scalar equation must not be certified by a noisy standalone numerical
`J_phi'` differentiation if an exact Noether/Bianchi identity can be derived in
the same action conventions.  The final scalar gate must either use that exact
identity or a simultaneous action-level solve in which all background EOM are
satisfied by the same member.  Fixing the scalar residual while moving `E11`
off shell is explicitly non-promotable.

## 9. Continuation logic

Continuation starts from the known healthy left anchor and advances inward in
small steps.  At each step:

1. solve the background/action constraints;
2. absolutely emit the local 41-slot stream;
3. reduce the physical operator;
4. track modes;
5. enforce stability/hyperbolicity and electric nondegeneracy;
6. monitor the constraint Jacobian rank.

If an eigenvalue/characteristic reaches zero or the constraint Jacobian loses
rank, the event is treated as physical/mathematical information (possible
branch endpoint, bifurcation, or missing action direction), not something to
hide with larger regularization.

## 10. QNM and final promotion

Only after one action passes from the left anchor through the inner light ring
to the core may the regional absolute streams be assembled globally.  The
sequence is then

```
absolute regional action members
 -> center-to-infinity C41_absolute
 -> Direct Global KRGSM
 -> same-operator coupled QNM
 -> fresh replay / manifest / release archive.
```

For the horizonless P5 object the spectral boundary conditions are a regular
center and outgoing infinity, not an ingoing horizon condition.

## Current status

The transition is **OPEN**.  The current unmodified hybrid extension loses
positive even-sector kinetic definiteness before the inner light ring, while
its background remains numerically close to on shell.  The old historical
Inner member is not eligible as the right physical neighbor.  The existing
`G2_XX` lift has zero measured kinetic sensitivity to this new obstruction, so
additional action-normal directions must be identified by mode-projected
sensitivity/rank analysis rather than blind coefficient fitting.

## 11. First mode-projected direction map

The first primitive sensitivity audit at the inner light ring identifies the
compact `a1` Horndeski primitive direction as a strong positive lever on the
failing even-sector kinetic mode.  `F_tensor` is nearly null and `H_tensor` is
weakly negative at the sampled point.  This does **not** promote an `a1` patch:
the primitive direction must be realized by a joint background-preserving
Horndeski+SVT action bundle.  It only supplies the rank information needed to
choose that bundle rationally.

Two tempting one-channel shortcuts have now been falsified explicitly:

- the existing background-null `G2_XX` lift has zero measured kinetic response
  to the right-transition obstruction even at very large diagnostic amplitude;
- a standalone quartic `G4_XX` tubular compensation has a rank-deficient
  background pivot and cannot be used as a one-channel repair.

The next continuation implementation should therefore construct a jointly
compensated normal bundle containing an `a1`-effective direction and verify the
background EOM and constraint rank before using its positive kinetic
sensitivity.

## 12. Projected on-shell controllability gate

Before any finite predictor/corrector continuation step, a candidate control
bundle must pass the local action-tangent test.  For action coordinates `p` and
background equations `E`, define

```
B = dE/dp,
N = ker(B),
g_L = d lambda_min(K_L) / dp.
```

A raw primitive sensitivity is not sufficient.  The physically relevant local
question is whether `g_L` has nonzero projection onto `N`.  The current audit
implements this with fresh action-level local perturbations and absolute 41-slot
reemission.  At `u ~= 0.7013, 0.7030, u_LR, 0.7080`, one background-tangent
direction exists that improves the weakest kinetic mode simultaneously for
`L=6,20,42,1000`.  The result survives inclusion of the present algebraic
scalar-identity row and is finite-difference converged at the light ring.

The requested point `u=0.712` is deliberately not evaluated yet: the accepted
absolute Electric-Hybrid action table ends below the historical `u=0.71` seam.
Using the rejected historical Inner coefficient stream as an absolute action
seed would violate this contract.  The finite continuation must first carry the
same action across the seam.

The earlier raw Horndeski `a1` clue has also been projected onto the full local
quartic `G5=0` background-null tangent generated by
`(G4_XX,G3_X,G2,G2_X)`.  A large fraction of the `a1` gradient survives this
projection, so `a1` is not merely an off-shell artifact.  However an
`a1`-maximizing one-dimensional bundle is localization dependent: sufficiently
narrow compact bundles can improve low/intermediate `L` while worsening the
sampled high-`L` mode.  Therefore `a1` is a useful component of the normal
bundle, not a promoted one-channel repair.

Promotion rule:

```
projected local controllability PASS
    -> joint H+SVT null-bundle predictor
    -> nonlinear background corrector
    -> absolute Action -> C41 reemission
    -> K / radial / validated-angular gates
```

No finite transition PASS may be inferred directly from the tangent audit.
