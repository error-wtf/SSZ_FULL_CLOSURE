# Sagnac Reference Transport — Module Spec (SAG_REF_V1, 2026-10-05)

Status: SPEC (implementation work order). This module is an INDEPENDENT
reference layer. It validates the transport METHOD, not SSZ.

## 0. Anti-circularity rationale (binding)

The registered TRUE FULL CLOSURE (43/43 gates, G00..G140) validates the SSZ
transport kernel `src/ssz_p5/postclosure/transport.py` against internal
geometric identities, Schwarzschild anchors, conservation laws and negative
controls. It does NOT yet validate the transport machinery against a fully
analytically known transport observable.

This module closes exactly that gap with the cheapest exactly-solvable
rotational transport case: the Sagnac effect.

Chain of trust being constructed:

    known physics (Sagnac, closed form)
        -> validated transport mathematics (3 independent routes agree)
            -> SSZ transport (existing kernel, unchanged)

instead of

    SSZ implementation -> SSZ validates itself.

## 1. Scope & non-touch guarantees

- NO file under `src/ssz_p5/postclosure/` is modified by this module.
- NO entry in the existing gate registry (`true_closure/graph.py` GATES)
  is redefined, renumbered or reinterpreted. `43/43 PASS` keeps its exact
  current meaning.
- This module registers its own gate group `G101..G109` (subscript "Sag")
  in its own registry file and its own verdict key. It does NOT plug into
  the G140 condition.
- Future meta-gate (separate work order, NOT part of SAG_REF_V1):
  `G150 = G_Sagnac(all) AND G130 AND G140`
  ("reference-validated transport architecture").

## 2. Physical ground truth (no SSZ content)

Flat-spacetime ring transport, inertial reasoning in the rotating frame:

    t_+ = L / (c - v)          co-rotating signal
    t_- = L / (c + v)          counter-rotating signal
    dt  = t_+ - t_- = 2 L v / (c^2 - v^2)

Equivalent area form for a circular ring of radius R (L = 2 pi R, v = Omega R):

    dt = 4 Omega A / c^2 * 1/(1 - v^2/c^2),   A = pi R^2
       -> 4 Omega A / c^2  in the limit v << c.

Phase / JIF readout at carrier frequency f = omega/2pi:

    dphi = omega * dt,          dJ = dt * f = dphi / (2 pi)

Inversion (branch-selected):

    v(dt) = (c^2 L - sqrt((c^2 L)^2 - 4 L^2 dt c^2 ... )) -- solved from the
    quadratic  dt v^2 /L + (2v) - dt c^2 /L = 0;
    physical root |v| < c, sign(dt) = sign(v).

Edge contracts (mandatory negative controls):
    v -> 0  =>  dt -> 0
    v -> -v =>  dt -> -dt
    |v| -> c^- => t_+ -> infinity (co-rotator never catches the beacon)

## 3. Module layout (new, self-contained)

    src/ssz_p5/reference_transport/
        __init__.py
        sagnac.py          # closed form (Route 1) + edge contracts
        transport_pde.py   # Route 2: first-passage PDE transport on a ring
        segment_chain.py   # Route 3: N-segment discrete chain, N -> inf
        inversion.py       # dt -> v branch-selected quadratic inversion
        gates.py           # G101..G109 registry + evaluator (own verdict)
    tests/reference_transport/
        test_sagnac_closed_form.py
        test_sagnac_pde.py
        test_sagnac_segment_chain.py
        test_sagnac_inversion.py
        test_sagnac_gates.py

Import discipline: `reference_transport/*` must NOT import
`ssz_p5.postclosure`, `ssz_p5.true_closure` or any SSZ geometry module.
Enforced by a unit test that scans the module's import graph
(F821-style forbidden-import assert, same pattern as the bridge
blindness asserts).

## 4. The three independent routes

### Route 1 — closed form (sagnac.py)
Direct evaluation of t_+, t_-, dt, dphi, dJ from the formulas in §2.
Trivially exact; serves as the oracle for Routes 2 and 3.

### Route 2 — PDE first-passage (transport_pde.py)
One-dimensional transport equation on the ring circumference x in [0,L)
(galilean kinematics in the rotating frame; c is the wave speed relative
to the medium, v the receiver speed):

    d/dt psi(t,x) + s * (c - s*v) * d/dx psi(t,x) = 0,   s in {+1,-1}

emitter/receiver at x=0, periodic boundary, signal = step or narrow
Gaussian packet. Arrival time = first time the receiver sees the packet
peak exceed a fixed threshold (declared: 0.5 of max, sub-sample by
linear interpolation). Numerics: method-of-lines + upwind flux,
CFL = 0.4, grid convergence required (G107).
Expected result: t_s = L/(c - s v), i.e. dt(PDE) -> dt(closed) as
grid -> inf.

Why this is a REAL test: the PDE route solves a drift-advection problem
and measures a first-passage time — the same mathematical operation
class the SSZ kernel performs (transport then readout), but in a case
with an exact answer.

### Route 3 — segment chain (segment_chain.py)
Discretize the ring into N equal arcs. In each arc the signal propagates
at c relative to the medium; the receiver moves v per unit time. Per
segment:

    dt_seg = (L/N) / (c - s*v)

Chain total: t_s(N) = N * dt_seg = L/(c - s v) exactly for uniform
segments (trivial), BUT the non-trivial version used here randomizes
segment lengths (Dirichlet distribution, normalized to L) — the chain
must converge to the same t_s independent of the partition, and the
segment-count convergence N -> inf must be O(1/N) at worst in the
partition-randomized case (G106, G107).

Why this is a REAL test: segment chains are precisely the SSZ
discrete-phase methodology (counted segments -> counted phase). Showing
the chain reproduces an exact known transport observable validates the
counting machinery itself.

### Inversion (inversion.py)
Given dt (from ANY route), recover v via the branch-selected quadratic.
Round-trip contract: v_rec = v_in to 1e-12 relative for |v|/c in
[0.01, 0.99]. Sign contract for negative v. Ambiguity contract: no
real solution for |dt| > L*|v|/c^2 domain edge — must raise.

## 5. Gate matrix (own registry, verdict key: `SAG_REF_V1`)

| Gate | Content | Pass criterion |
|------|---------|----------------|
| G101_Sag | closed-form t_+ | == analytic, rel err < 1e-14 |
| G102_Sag | closed-form t_- | == analytic, rel err < 1e-14 |
| G103_Sag | closed-form dt  | == analytic, rel err < 1e-14 |
| G104_Sag | direction reversal v->-v | dt antisymmetric, < 1e-15 asymmetry |
| G105_Sag | null limit v->0 | dt < 1e-15 for v < 1e-15 c |
| G106_Sag | segment convergence | |t(N)-t_inf|/t_inf < 1e-12 at N=4096; partition-independence < 1e-12 |
| G107_Sag | PDE-vs-closed-form | first-passage dt matches to grid-converged tol; observed order >= 1 |
| G108_Sag | phase readout | dphi = omega dt and dJ = f dt round-trips, rel < 1e-14 |
| G109_Sag | inversion round-trip | v(dt(v)) == v, rel < 1e-12; branch+sign contracts; no-solution raise |

Negative controls (part of the gates, not optional):
- corrupted dt (one ULP flip in route 3 chain) MUST fail G103/G106;
- inverted sign convention MUST fail G104;
- forbidden import (postclosure) MUST fail the import-discipline test.

## 6. Verdict semantics

    SAG_REF_V1_PASS  := all G101..G109_Sag PASS
    (independent of and additive to TRUE_FULL_CLOSURE_PASS)

Documentation records the composition separately:

    existing TRUE FULL CLOSURE  +  SAG_REF_V1  =  two independent pillars.

A future G150 (separate spec) may require both; it is explicitly out of
scope here and must NOT be added as a condition of G140.

## 7. Implementation order

1. `sagnac.py` + closed-form tests (G101-G105) — pure function work.
2. `segment_chain.py` + tests (G106) incl. partition randomization.
3. `inversion.py` + tests (G109).
4. `transport_pde.py` + tests (G107, G108) — the only numerically
   non-trivial route; CFL and threshold sub-sampling documented.
5. `gates.py` registry + verdict JSON
   `data/generated/reference_transport/SAG_REF_V1_VERDICT.json`
   + SHA256 freeze of the verdict.
6. Forbidden-import discipline test.
7. README section "Transport method validation (external reference)"
   — additive, no change to existing gate tables.

## 8. Interface to the later SSZ bridge (documented, not built)

The abstraction the three routes share is:

    S (signal law, here: galilean c, drift v)
    -> A (architecture: ring, emitter/receiver, first-passage readout)
    -> L (transport operator: PDE or segment chain)
    -> Psi (propagated field)
    -> dt (first-passage time)
    -> dphi/dJ (phase readout)

The existing SSZ kernel instantiates A and L with
A = A[f(r), h(r), ...] and L = geodesic/phase transport.  The later
bridge (separate work order) re-expresses the SSZ kernel inside THIS
abstraction and checks: same formalism, two geometries (flat ring and
SSZ), one validated transport mathematics.
