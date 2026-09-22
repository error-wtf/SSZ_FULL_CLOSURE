# CHATGPT HANDOFF — SSZ P5 Full Closure (2026-09-22)

## A. CURRENT TOP-LEVEL STATUS

```
ABSOLUTE_FULL_CLOSURE          = false
Q1_HEALTHY_K_EXISTENCE         = TRUE  (strong-H carrier 1.605e-7 @ L=1000; exterior 9.89e-7)
Q2_ELECTRIC_SAME_ACTION        = OPEN  (paused on Eq85 replay + C_shared provenance)
Q3_FULL_GATE_STATUS            = NOT_STARTED_FOR_CERTIFIED_COMMON_ACTION
EPSILON_Y_MEMBER_STATUS        = REPLAY_PENDING (neither active nor superseded-as-replayed)
HISTORICAL_SPLIT_Q2            = UNKNOWN (C_shared provenance pending; archival track)
EQ85_REPLAY_STATUS             = BLOCKED_PENDING_PRIMARY_SOURCE_RETRIEVAL
```

## B. CURRENT DECISION TREE

```
Eq85 replay (spec frozen in CANONICAL_FACTS.yaml, EQ85_PRIMARY_SOURCE_AUDIT)
   |
   +-- confirms inner-ring obstruction (required tensor_F = -1.6386 vs tensor_H = +1.0)
   |      -> EPSILON_Y member = SUPERSEDED_REPLAYED (final)
   |      -> electric Q2 continues on the partition-chain architecture
   |
   +-- refutes obstruction (convention/numerics resolved)
          -> EPSILON_Y member = ACTIVE_BLIND_TEST_CANDIDATE
          -> execute the preregistered frozen forward test (eps_Y=0.01 fixed, no repairs)
```

## C. PERMANENT MEMBER REGISTRY (summary; full: ACTION_MEMBER_REGISTRY.yaml)

| MEMBER_ID | domain | K | angular | same-action | status |
|---|---|---|---|---|---|
| C0_TRANSVERSE_ZERO | u 0.7002-0.7080 | FAIL (-18.5223 @ L1000) | diagnostic | = central data | falsification witness |
| C1_MINIMAL | = C0 stream | FAIL | diagnostic | falsified hypothesis | witness |
| C3_FIELDSPACE_MINIMAL | = C0 stream | FAIL | diagnostic | falsified hypothesis (NEW_THEORY_HYPOTHESIS_WITH_FIXED_FIELDSPACE_METRIC — never rename) | witness |
| CENTRAL_MEMBER (CENTRAL_EXACT_SVT descendant) | u 0.61-0.70998 | FAIL near u~0.708-0.71 (full finite-L; consistent w/ 2026-09-17 gate) | not replayed | the hybrid base | active construction base |
| DENSE_K_REPAIRED_MEMBER | u 0.7002-0.7080 | PASS (1.26e-5) | **FAIL_ROBUST** (eps^3: cminus^2<0 on 297/297, min -1473.2 @ u=0.706335, pure psi-mode) | falsification witness | frozen |
| STRONG_H_CARRIER | u 0.5515-0.71 | PASS 1.6e-7 @ L1000 (all L) | not replayed | OPEN (Q2) | Q1 existence witness |
| GENUINE_SVT_LOBE_FULL | u 0.61-0.715 | K_ZK_exact 10.54 = reduced combination (SUPERSEDED_NUMERICS as full witness) | 84.912/643.142 SUPERSEDED_NUMERICS | different-representative (documented u=0.61) | archival regression record |
| OUTER_H_CARRIER (exterior) | outer | PASS 9.89e-7 @ L1000 | not replayed | OPEN (Q2) | existence witness |
| INNER_H_CORE | inner | deep-core FAIL_IMPLEMENTATION (2026-09-18) | — | OPEN | implementation fail ≠ physics fail |
| EPSILON_Y_ZERO_VECTOR_MEMBER | A0'=0 global | component tests passed | — | SUPERSEDED claim (REPLAY_PENDING!) | historical |

## D. FALSIFICATION WITNESSES (permanent, never erase)

1. C0/C1/C3 minimal stream → finite-L K ghost (min K = −18.522265817831975 @ L=1000, u=0.70799)
2. K-repaired dense member → direct-action angular FAIL (cminus²<0 on 297/297; min −1473.21 @ u=0.706335; pure psi-mode weight 0.99997)
3. C3's unique minimizer is the trivial zero third covariant jet (Gram PD 1.54e9..4.30e10) — per owner rule this triviality IS its prediction
4. ε_Y/A0'=0 member → Eq.85 inner-light-ring obstruction (documented, replay pending)
5. Lobe CSV = documented DIFFERENT action-jet representative at u=0.61 (190-430% naive mismatch; NOT an action incompatibility proof)

## E. HEALTHY EXISTENCE WITNESSES (prevent "healthy K does not exist")

Strong-H carrier + exterior region (see C). These are EXISTENCE witnesses only —
they do NOT constitute certified same-action closure (Q2) nor full closure (Q3).

## F. ANGULAR METHOD STATUS

- Production angular method: direct action-level ε³ determinant via exact rank-1
  null-vector elimination (tools/audit_dense_angular_direct.py; witness anchor 1.1e-15).
- ε² leading determinant cancels structurally (rank-1 (psi,V) block) — feature, not bug.
- Historical reduced angular values (84.912 / 643.142): SUPERSEDED_NUMERICS unless
  independently replayed by the current determinant route.
- Mode tracking: vector root identity-based; no naive sorting.

## G. EQ. 4.47 / m5 STATUS

- Eq. 4.46 + 4.47: PASS in same basis/scaling on 0.62≤u<0.70; the literal printed
  residual is exactly predicted by (R−1)·r·h·φ′·K22·K33 (scaled 8.09e-9).
- m5_action_minus = a4·v13 − a6·v6 = PRIMARY_ACTION_PRODUCTION.
- m5_printed_plus = PUBLISHED_COMPARATOR_ONLY (PUBLISHED_REDUCED_SHORTCUT_CONFLICT).
- Classification remains: not an official paper erratum.

## H. KNOWLEDGE-CONTINUITY LESSON (binding)

- Repo = canonical living state. ZIP = immutable snapshot. RAG = provenance archive
  (hits need member/scope/supersession mapping — never automatic current truth).
- No "latest file wins"; no formula outside its registered scope
  (see SVT_SPLIT_SCOPE_REGISTRY.json: the S_SVT relation is OUTER_HANDOVER-scoped —
  applying it to LOBE_CENTRAL was the session's documented scope error).
- tools/context_guard.py blocks overstrong claims and scope misuse — run before
  any status promotion.
- UNKNOWN_PENDING_SOURCE_RECOVERY is preferred over any reconstruction by plausibility.

## I. NEXT EXACT STEP (ONE TASK)

Implement Kase-Tsujikawa Eq. 85 (arXiv:2301.10362) independently from the paper
text — do NOT copy src/ssz_p5/action/light_ring_electric.py (historical
implementation, mechanism structurally confirmed: at the ring (2f−rf′)·a4′=0
removes the a4′ control; geom = r·f″ − r·f′²/f + 2f′ − 2f/r) — evaluate at both
light rings on data/production/ssz_p5_horndeski_carrier_through_light_rings_to_core_2026-09-12.csv
(sha 65b32db2...), confirm/refute the documented targets
(inner: required tensor_F = −1.6386 vs tensor_H = +1.0; outer: +0.9721 consistent),
then follow the decision tree in section B.
