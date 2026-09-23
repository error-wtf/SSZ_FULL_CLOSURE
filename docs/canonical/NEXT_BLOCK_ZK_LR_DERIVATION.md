# NEXT WORK BLOCK — ZK-LR DERIVATION (frozen 2026-09-22, commit 59050778+; STEP 0 appended per owner)

Frozen protocol (owner directive; no deviations, no interpretation freedom):

## STEP 0 — ASSERT_EQ85_V8_AT_INNER_LR (TOP PRIORITY, before everything)

Prerequisites to nail down BEFORE any derivation:
  1. Is -16.824455 the correctly mapped Eq.-85 v8 coefficient (ZK-V9 per the
     V9=MH-v8 sub-scheme mapping) at the inner stable light ring?
  2. Was it evaluated on the same final same-action member (the electric central
     stream) — not on another representative?

Output (mandatory fields):
  r_LR, u_LR, v8^Eq85, W''(u_LR), a4, sgn v8.

If v8^Eq85 ~ -16.824455 is confirmed:
  hard reject for THIS member: v8 <= 0 => Delta_LR < 0 for EVERY real A0_u
  (A0_u^2 >= 0, v8 < 0 => electric term <= 0, W'' a4 > 0) — the Eq.85 electric
  channel has the WRONG SIGN, not merely insufficient magnitude.
  The search then shifts from "more electric support" to "which SVT structure
  can deliver the correct sign at all" (STEP 1 becomes decisive).

If the mapping or member identity does NOT confirm:
  the -16.82 value is misattributed; redo the evaluation on the correctly
  mapped coefficient of the final member before any conclusion.


## STEP 0a — PROVENANCE FIELDS (mandatory in the STEP-0 report)

member_hash | member_source | evaluation_radius_r | evaluation_u |
source_coefficient_name = MH_v8 | canonical_coefficient_name = ZK_V9 |
mapping_rule: ZK_V9 <- MH_v8 (TWO SOURCES: ssz_hybrid_unreduced_even_kernel.py:52
AND Monograph chunk 38, SHA a028d69b...) | value | sign |
interpolation vs exact-row provenance.

THREE-CASE LOGIC (frozen):
  A) wrong slot / wrong member          -> STEP 0 decides A vs (B|C)
  B) correct Eq.85 coefficient, v8 < 0  -> STEP 1 decides B vs C
  C) correct coefficient, but the FULL ZK-SVT LR identity contains additional
     electric support terms
No U(1)-SVT no-go may be inferred from v8 < 0 alone (case C exists).

## STEP 1 — General ZK-LR identity, symbolically, undivided

TARGET (schematic, NOT pre-reduced to Eq.85):

  C_LR^ZK = C_geom + C_A0 + C_f2 + C_f3 + C_f~3 + C_f4 + C_f~4 = 0

built from the FULL general SVT background theory. Evaluate at the stable ring
only after assembly. If after all specializations only
C_A0 ∝ A0prime^2 * v8^Eq85 survives electric support and v8 = -16.824455 is
confirmed -> THIS member is structurally dead. If an additional term
A0prime^2 * S_SVT (or another nonvanishing genuine-SVT combination) survives ->
that operator can make LR-COMPATIBILITY possible in the first place — it does
NOT by itself enable Full Closure (Joint-DAE, Same-Member, K, c_r^2, c_Omega^2
and interfaces must still pass). Do not prejudge which branch.

CASE B PRECISE WORDING: v8^Eq85 < 0 AND the full ZK-LR identity contains no
additional positive support channel => THIS MEMBER is locally excluded at the
stable LR (state exactly WHERE and WHY) — not a statement that the theory is dead.

STEP 1 REQUIRED FINAL FORM — three-channel, dimension- and sign-clear:

  Delta_LR^full = (-W'' a4)            [geometric load]
                + (A0_u^2 * Xi_MH)     [Eq.85 channel]
                + (Sigma_SVT)          [genuine SVT]

  with F = (sqrt(h)/f^(3/2)) * Delta_LR^full.

Each channel is then separately checkable: sgn Xi_MH, sgn Sigma_SVT,
Delta_LR^full > 0. This exposes WHICH operator actually works against the
positive curvature of the stable optical potential.

## STEP 1 — General ZK-LR identity, symbolically, undivided

From the GENERAL U(1)-SVT background EOM construct a combination C_LR that
contains NO factor 1/(2f - r f') BEFORE taking the light-ring limit. Only then
set 2f - r f' = 0. Then, in this order: A0prime = 0, then Delta f2 = eps_Y Y.

SUCCESS TEST (symbolic, not numeric): factor the eps_Y contributions with at
least one of A0prime^2, F, or Y, such that
    A0prime = 0  =>  Delta C_epsY = 0.

If 0: the eps_Y zero-vector branch is STRUCTURALLY CLOSED (not merely
numerically failed). If not 0: an overlooked genuine SVT support term exists —
derive its sign/regularity. In no case may anything be "made to fit".

## STEP 2 — Joint background solve as DAE (not optimization)

From E11 = 0 algebraically: f2X = A(r) + B(r) f3X.
The linearized relation delta f2X = -(2h^2 phi' A0prime^2/(r f)) delta f3X must
RE-EMERGE AS A REGRESSION TEST (not be inserted as an assumption).
Substitute into E_phi = 0 -> first-order ODE/DAE for f3X.
Integration constant: fixed by the existing action value at a certified overlap
boundary — NOT fitted.

## STEP 3 — Same-member gate (brutal, one hash) — CORRECTED per owner

UNIVERSAL gate (case-independent):
  ONE action hash must simultaneously satisfy:
    E00 ~ 0, E11 ~ 0, J_A ~ 0, E_phi ~ 0, Delta_LR^full > 0,
  plus regularity/pivot conditions.

v8_sign and the Eq.85-only S_LR are DIAGNOSTIC FIELDS, not universal gates.
Reason: if STEP 1 finds Sigma_SVT != 0, Delta_LR^full > 0 can hold while
Xi_MH ∝ v8 < 0 (Sigma_SVT overcompensates the load) — the old gate would
falsely reject a physically valid case C.

BRANCH RULE (frozen, WITH PREMISE GUARD):
  PREMISE (mandatory): f > 0, h > 0, W'' > 0, a4 > 0 at the inner ring.
  Without these the sharpened rule is applied OUTSIDE its derivation domain.
  Sigma_SVT = 0 AND premises hold:
    necessary: v8 > 0, S_LR^Eq85 > 1, Delta_LR^full > 0
    (then W'' a4 > 0 and S_LR is well-defined with its intended meaning)
  Sigma_SVT != 0:
    Delta_LR^full > 0 only; v8, S_LR^Eq85 diagnostic
  Premises not met:
    no sharpened rule; fall back to the universal gate + explicit
    UNKNOWN/Premise-violation record.

One table, one member, one hash. No "scalar solved file" vs "metric solved file".

FORM STATUS (semantic convention, frozen): until derivation, Delta_LR^full is
merely the NAME for the complete light-ring compatibility quantity. It is NOT
yet the claim that it decomposes as -W'' a4 + A0_u^2 Xi_MH + Sigma_SVT.
STEP 1 may PROVE the three-channel form — or REPLACE it by the exact form
following from the general SVT EOM. The candidate form for
F = (sqrt(h)/f^(3/2)) Delta_LR^full must not silently become a premise.
Proving or replacing it is STEP 1's job.

## STEP 4 — Only then

K -> c_r^2 -> c_Omega^2. The 41-slot angular machine stays OFF until the common
electric background satisfies the same action.

## Reference values (frozen)

outer ring: analytically F/H = 1; numerical replay 1.000186 (carrier H(2/3)=0.9720791)
inner ring: F/H = -1.639678 (documented -1.63868); W''=+2.59314; v8(current member)=-16.824455

## Already source-confirmed

Y and all first variations of Delta f2 = eps_Y Y vanish on the F-tilde = 0
background (Monograph, chunk 3, SHA a028d69b...) — the ZK step is a
primary-source/generality check, not an open physics guess.


## ORDERING FREEZE (owner, 2026-09-22 final)

STEP 1 and STEP 2 must NOT be merged logically:
  STEP 1: derive AND FREEZE the full ZK-LR identity.
  STEP 2: joint E11 + E_phi DAE on exactly THAT identity.
Rationale: otherwise the DAE could presuppose a form of the LR balance that
STEP 1 is supposed to prove.

Threshold semantics: 0.2230195 is the EQ.85-ONLY threshold under
Sigma_SVT = 0 until STEP 1 decides. v9 > 0 at the inner ring remains the
standing proof that the Eq.85 channel itself has the correct sign.

The most interesting concrete measurement after STEP 1:
  S_LR^Eq85 = (2 u^2 f h A0_u^2 v9) / (W'' a4)
for the REAL background. Under Eq.85-only: S_LR > 1 <=> A0_u^2 > 0.2230195.
If the actual on-shell value (from J_A = 0 and the joint background EOM) lies
above: first ON-SHELL electric support of the same member.

## STEP 1 TARGET — FINAL FORMULATION (owner, 2026-09-22, session close)

  C_LR^full = C_bg |_(2f - r f' = 0)

where C_bg is formed UNDIVIDED beforehand from the COMPLETE SVT background EOM.
Then group by origin:

  C_LR^full = C_geom + C_MH-electric + C_genuine SVT

Three hard checks (in order):
  1. A0prime = 0 AND Delta f2 = eps_Y Y  =>  C_epsY =? 0
  2. MH limit  =>  C_LR^full  ->  -W'' a4 + 2 u^2 f h A0_u^2 v8
  3. Real inner ring: insert v8^Eq85 = +17.521313 with the on-shell A0_u of the
     SAME member; then S_LR^Eq85 = (2u^2 f h A0_u^2 v8)/(W'' a4) > 1?

If S_LR > 1 AND STEP 1 produces no additional negative SVT load:
first genuine ON-SHELL electric support of the same member.
If additional genuine-SVT terms survive: evaluate the complete balance.

The electric channel is no longer merely formally possible: its sign is correct
and its necessary Eq.85 threshold is analytically known. What remains is pure
mathematics: does the full ZK background theory carry this channel alone, or
does it bring a further support operator?

Protocol refinement CLOSED. No further meta-work.

## FINAL EVALUATION ORDERING (owner, session close)

Primary test: evaluate the FULL on-shell LR balance
  C_LR^full (respectively Delta_LR^full if the three-channel reduction is proven).
DIAGNOSTIC (secondary, Eq.85/MH channel only): A0_u^2 >? 0.2230195, S_LR^Eq85.
Reason: if C_genuine SVT != 0, then A0_u^2 < 0.2230195 can still be carried by a
positive genuine-SVT term — or an additional negative term can consume an
apparently sufficient Eq.85 support. 0.2230195 is the Eq.85/MH-channel threshold,
never the universal one.

No further protocol architecture needed. The next relevant output is the
actual symbolic LR identity.

## STEP 1 AMENDMENT — f2Y SYMBOLIC + SOLVE RANK (owner, final)

1. f2Y stays SYMBOLIC through the C_bg construction. Order:
   C_bg(f2Y symbolic) -> MH-/eps_Y-checks -> f2Y=0 production-branch specialization
   -> real member evaluation.
   Rationale: the combination f2F - 2 h phi'^2 f2Y is the ideal vehicle; zeroing
   f2Y at the entrance would remove the eps_Y falsification test by definition.
   (Evaluator context: the current branch REJECTS nonzero f2Y — that is a branch
   postulate, registered; the symbolic derivation must not inherit it.)

2. LOCAL SOLVE RANK is part of the derivation:
   det[ d(JA, E00, E11) / d(f2F, f2, f2X) ] != 0,
   evaluated in particular at the inner light ring. Basis order:
   JA=0 -> f2F = J[remaining jets]
   E00=0, E11-E00=0 -> f2, f2X (diagonalized basis; kappa = -2X documented)
   E_phi=0 -> d f3X / dr = S[f3X, ...]
   These results enter the UNDIVIDED Eq.85-generating combination — not the
   reverse. Then ring limit, then subtract geometry + MH-electric; only the
   remainder is Sigma_SVT.

## FULL-CLOSURE END-RUN (frozen execution plan, owner 2026-09-22)

STEP 2 RESULT: the 3.44e-4 E_phi window-max is a BOUNDARY artifact (row 0, u=0.69751,
the handover edge); interior max = 1.64e-5 (u=0.70854) — derivative level.
Classified per rule: boundary/stencil effect, documented; interior clean.

Execution order (frozen):
 1. E_phi maximum localize/classify            [DONE — boundary artifact, interior 1.64e-5]
 2. C_bg from the integrated unsplit same-action member
 3. Ring limit
 4. Sigma_SVT exact extraction
 5. Delta_LR^full on shell
 6. Joint DAE over the full relevant window
 7. Same-member gate, ONE hash (E00/E11/J_A/E_phi ~ 0 + full light-ring condition)
 8. K > 0
 9. c_r^2 > 0
10. c_Omega^2 > 0
11. Interfaces / continuation
12. QNM/trapping LAST

No new branch search. No fitting. No old slot loops. From here: close the
remaining proof chain.

## STEP 2 EXECUTION ORDER — THREE SEPARATE PROJECTIONS (owner correction, frozen)

C_bg is built FULLY SYMBOLIC: C_bg[f2Y, A0prime, f3, f~3, f4, f~4, ...].
The old order ("MH projection A0prime=0, then compare") is WRONG — it kills the
electric Eq.85 term (-2 r f h A0prime^2 v8) that must be REPRODUCED.

Test A — eps_Y-nulltest:
  A0prime = 0, f2Y stays SYMBOLIC. Require algebraically: Delta C_epsY = 0.
  This tests the background-null zero-vector eps_Y branch.

Test B — full Maxwell-Horndeski limit:
  genuine-SVT structures OFF, but A0prime stays SYMBOLIC (not set to 0).
  Require: C_bg falls back EXACTLY to the undivided Eq.85 numerator:
    (2f - r f') a4' - [ (r f'' - r f'^2/f + 2f' - 2f/r) a4
                        + f^(3/2)/(r sqrt(h)) F
                        - 2 r f h A0prime^2 v8 ] = 0.
  Only this test proves the full Eq.85 structure INCLUDING the electric
  coefficient was reproduced. Slot rule: Eq.85 'v8' = ZK-canonical V9 = MH v8
  (never the 13-slot v8 column).

Test C — real production member:
  NOW apply f2Y = 0 (branch postulate) + the actual integrated A0prime(r) + the
  re-solved action jets. Then ring limit (2f - r f' = 0), and only afterwards:
  Sigma_SVT = C_LR^full - C_LR^MH.

Rationale: A0prime=0 is the right projection for the eps_Y background-null test
but NOT for the full electric Eq.85 regression check. The two statements must
never be merged under one "MH projection" step (RAG-consistent: the old
Delta f2 = 0.01 Y member is documented background-null on the zero-vector
branch; the new electric computation examines a different branch).
