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

## STEP 3 — Same-member gate (brutal, one hash)

ONE action hash must simultaneously satisfy:
  E00 ~ 0, E11 ~ 0, J_A ~ 0, E_phi ~ 0,
and at the inner stable ring:
  v8^Eq85 > 0,  S_LR > 1,  Delta_LR > 0.

One table, one member, one hash. No "scalar solved file" vs "metric solved file".

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
