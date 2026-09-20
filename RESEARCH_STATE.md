# SSZ P5 — consolidated research state (2026-09-19)

This is the short truth ledger.  Historical files are preserved, but status claims
below supersede isolated older wording.

## Verified infrastructure that remains valid

- P5 regular horizonless geometry, analytic center and the two light rings.
- 41-slot perturbation language, regional emitters and shared-baseline assembly rules.
- JET9D8 radial derivative service and corrected holonomic identities.
- Conditioned Maxwell-Horndeski Core principal/kinetic invariants.
- Generalized-psi unreduced descriptor that keeps H0/h1 constraints explicit.
- Strong-H differential-operator pullback equivalence between descriptor and the
  established reduced KRGSM operator.
- Singular radial DAE quadratic pencil that does not invert the constraint kinetic block.
- Existing test-field/eikonal/Jost calculations under their historical scope.

## Rejected concrete members (not rejected infrastructure)

### 2026-09-17 regional electric member
Rejected for Absolute Closure: Central finite-L kinetic FAIL at L=6,12,20,42,
and the selected Central c3/e3 pair is not a common direct-action lower replay.

### 2026-09-18 global zero-vector epsilon_Y member
Rejected as a global on-shell member.  A denominator-free Eq.85 audit at the inner
P5 light ring r/rs=1.416160655... gives, for A0prime=0,
`tensor_F/tensor_H=-1.638676...` while tensor_H>0 implies a4>0.  Thus a stable
positive-F/positive-H zero-vector continuation through that light ring is excluded
for this geometry.  The outer light ring gives F/H≈1 and acts as a control.

The many direct finite-L and descriptor PASS results obtained with this witness are
still valid tests of the software/representations; they simply do not promote that
witness to a global on-shell action.

## Active action search

The active target is an **electric Horndeski/SVT hybrid**.  The same Eq.85 audit
provides a local quantitative electric seed: if `tensor_F=tensor_H` is used only
as a diagnostic target, the inner light ring requires
`A0prime^2*v8_background_identity≈1.56365148818`.

See:
- `SSZ_P5_ELECTRIC_HYBRID_SEARCH_CANDIDATE_2026-09-19.json`
- `data/generated/absolute_attempt_2026-09-19/ZERO_VECTOR_LIGHT_RING_AUDIT.json`
- `data/generated/absolute_attempt_2026-09-19/ELECTRIC_HYBRID_EQ85_SEED.csv`
- `CONTINUE_IMPLEMENTATION_2026-09-19.md`

The next scientific gate is not another coefficient transplant.  It is one
holonomic electric action member that simultaneously satisfies the background
identities and finite-L scalar/tensor/vector stability, after which the already
validated descriptor/DAE/QNM infrastructure can be reused.

## QNM claim boundary

No final coupled HSVT QNM spectrum is certified.  Coupled QNM remains blocked until
the same action member has a global Direct-41/KRGSM certificate.

## Principal-feasibility result after primitive-control reconstruction

The previously external Central Horndeski control basis has now been extracted to
`data/generated/absolute_attempt_2026-09-19/ELECTRIC_HYBRID_HORNDESKI_PRIMITIVE_BASIS.csv`,
and the exact 4000-row search baseline is preserved as
`ELECTRIC_HYBRID_CENTRAL_SEARCH_BASELINE.csv`.  No external pickle is required to
reproduce the current search checkpoint.

A common Appendix-A primitive-response recipe using only the effective Horndeski
controls `a1`, `F_tensor`, and `H_tensor` reduces the robust Central kinetic
failure while keeping the algebraic pivots far from zero.  On the full JET9D8 grid
for `0.62 < u < 0.70`, the current reproducible recipe gives:

- L=6:  min(K) = -3.242946...; min|Dh1| = 151.06; min|aux det| = 441.88
- L=12: min(K) = -0.130994...
- L=20: min(K) = -0.0128774...
- L=42,110,420,1000: kinetic PASS.

This is a **principal-feasibility result only**.  It is not promoted to an action
member because the remaining low-L edge modes are still negative and the added
Horndeski primitives have not yet been reconstructed together with the electric
background as one on-shell holonomic action.  Attempts to force the last low-L
modes positive by unconstrained Newton steps drove the H0/h1 constraint pivot
near zero and produced spurious O(1e8) negative modes; those branches are rejected.

The practical conclusion is that further principal-only shaping is not the next
gate.  The next missing freedom must be taken from the electric/background action
sector and solved jointly with the Horndeski controls.

## Executable repository observability

The current state is now exposed by scripts rather than requiring chat context or
manual comparison of historical markdown files.  `python tools/show_repo.py` is
the primary truth view.  Focused scripts expose the gate matrix, member matrix,
evidence hashes, full file tree and Python code map.  `REPO_EVIDENCE_REGISTRY.json`
records expected audit exits so intentional rejection witnesses are not confused
with active-search failures.

Current grouped software result after adding these guards: **115/115 PASS**.
