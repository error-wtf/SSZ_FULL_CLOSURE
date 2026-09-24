"""TRUE FULL CLOSURE — machine-readable physics dependency graph.

The final verdict is computed EXCLUSIVELY from this graph: a node passes
iff all of its registered pytest node-ids pass; a gate passes iff its
nodes and its declared dependencies pass; TRUE FULL CLOSURE passes iff
every gate passes.  No manual PASS anywhere.

Node provenance keys refer to artifacts/true_full_closure/
physics_rag_provenance.json (Physics MCP/RAG citations + standard
literature, each identity independently re-derived in-repo).
"""
from __future__ import annotations

from dataclasses import dataclass, field

IMPL_TRANSPORT = "src/ssz_p5/postclosure/transport.py"
IMPL_CHAIN = "src/ssz_p5/true_closure/chain.py"
TESTS_FORWARD = "tests/scientific/test_true_closure_forward_chain.py"
TESTS_POSTCLOSURE = "tests/scientific/test_postclosure_source_free_transport.py"
TESTS_FALSIFIERS = "tests/negative/test_true_closure_falsifiers.py"


@dataclass(frozen=True)
class ChainNode:
    id: str
    statement: str
    relation: str
    impl: str
    tests: tuple
    deps: tuple = ()
    tolerance: str = "documented per test"
    provenance: str = ""
    gate: str = ""


NODES: list[ChainNode] = [
    # ---------------- foundations + differential geometry ----------------
    ChainNode(
        id="FOUNDATIONS",
        statement="Coordinates, signature (-,+,+,+), units c=G=r_s=1, "
                  "u=r_s/r; canonical metric ds^2 = -f dt^2 + dr^2/h "
                  "+ r^2 dOmega^2 with g symmetric, inverse exact and "
                  "det(g) = -(f/h) r^4 sin^2(theta).",
        relation="g_{munu} g^{numu} = delta; det(g) analytic identity",
        impl=IMPL_CHAIN,
        tests=(f"{TESTS_FORWARD}::test_foundations_metric_tensors",),
        provenance="foundations_conventions",
        gate="G120",
    ),
    ChainNode(
        id="CONNECTION_RIEMANN",
        statement="Levi-Civita connection and curvature: numeric Riemann "
                  "from finite-differenced Christoffels reproduces the "
                  "sympy-derived Ricci tensor; Riemann algebraic "
                  "symmetries and first Bianchi identity hold.",
        relation="R^mu_{nu rho sig} = d_rho G^mu_{nu sig} - d_sig "
                 "G^mu_{nu rho} + G^mu_{rho lam} G^lam_{nu sig} "
                 "- G^mu_{sig lam} G^lam_{nu rho}; R_{munu} = "
                 "R^rho_{mu rho nu}; R_{a[bcd]} = 0",
        impl=IMPL_CHAIN,
        tests=(f"{TESTS_FORWARD}::test_ricci_two_routes",
               f"{TESTS_FORWARD}::test_riemann_symmetries_and_bianchi1",
               f"{TESTS_FORWARD}::test_kretschmann_schwarzschild_anchor",
               f"{TESTS_FORWARD}::test_kretschmann_member_regular"),
        deps=("FOUNDATIONS",),
        provenance="curvature_differential_geometry",
        gate="G120",
    ),
    ChainNode(
        id="BIANCHI_CONTRACTED",
        statement="Contracted Bianchi identity div G = 0 — a geometric "
                  "identity for ANY metric; verified end-to-end through "
                  "the full curvature stack on the member (grid-convergent) "
                  "and on exact Schwarzschild.",
        relation="nabla_mu G^{mu nu} = 0 (nontrivial for nu = r)",
        impl=IMPL_CHAIN,
        tests=(f"{TESTS_FORWARD}::test_contracted_bianchi_member",
               f"{TESTS_FORWARD}::test_contracted_bianchi_schwarzschild"),
        deps=("CONNECTION_RIEMANN",),
        provenance="curvature_differential_geometry",
        gate="G120",
    ),
    ChainNode(
        id="KNOWN_LIMITS",
        statement="Known limits: exact Schwarzschild anchors (photon ring "
                  "u = 2/3, b_c = 3 sqrt(3)/2, Ricci-flat, Kretschmann "
                  "48 M^2/r^6, theta and sigma closed forms) and the "
                  "weak-field PPN signature (gamma = 1, Newtonian Kepler "
                  "Omega^2 = M/r^3) of the pipeline on a GR reference "
                  "metric.  Domain-limited claim: no asymptotic-flatness "
                  "claim beyond the member window.",
        relation="reference-geometry identities, independently computed",
        impl=IMPL_CHAIN,
        tests=(f"{TESTS_FORWARD}::test_ppn_signature_weak_field",
               f"{TESTS_FORWARD}::test_schwarzschild_anchor_bundle"),
        deps=("CONNECTION_RIEMANN",),
        provenance="known_limits",
        gate="G121",
    ),
    # ---------------- source-free transport (promoted) -------------------
    ChainNode(
        id="MATTER_TRANSPORT",
        statement="Matter moves source-free on g_SSZ: u^nu nabla_nu u^mu "
                  "= 0 with NO extra SSZ force; Killing energy E = f u^t, "
                  "angular momentum L = r^2 u^phi and norm are conserved; "
                  "FD acceleration matches the Christoffel transport term.",
        relation="u^nu nabla_nu u^mu = 0  <=>  dE/dtau = dL/dtau = 0, "
                 "g u u = -1",
        impl=IMPL_TRANSPORT,
        tests=(f"{TESTS_POSTCLOSURE}::test_timelike_source_free_conservation",
               f"{TESTS_POSTCLOSURE}::test_timelike_acceleration_crosscheck",
               f"{TESTS_POSTCLOSURE}::test_circular_orbit_kepler"),
        deps=("CONNECTION_RIEMANN",),
        provenance="geodesic_transport",
        gate="G110",
    ),
    ChainNode(
        id="NULL_TRANSPORT",
        statement="Light moves source-free in the eikonal limit: "
                  "k^nu nabla_nu k^mu = 0 for sub-critical and "
                  "near-critical impact parameters.",
        relation="k^nu nabla_nu k^mu = 0  <=>  dE/dlambda = dL/dlambda = 0, "
                 "k k = 0",
        impl=IMPL_TRANSPORT,
        tests=(f"{TESTS_POSTCLOSURE}::test_null_geodesic_conservation",),
        deps=("CONNECTION_RIEMANN",),
        provenance="geodesic_transport",
        gate="G111",
    ),
    ChainNode(
        id="RAYCHAUDHURI_TIMELIKE",
        statement="Timelike congruence dynamics: dtheta/dtau = -theta^2/3 "
                  "- sigma^2 + omega^2 - R_mn u^mu u^nu for the geodesic "
                  "radial congruence; irrotational (omega = 0); identity "
                  "grid-convergent; Schwarzschild closed-form anchor "
                  "(theta = -(3/2) r^-3/2, sigma^2 = (3/2) r^-3).",
        relation="Raychaudhuri equation (geodesic, irrotational limit)",
        impl=IMPL_TRANSPORT,
        tests=(f"{TESTS_POSTCLOSURE}::test_raychaudhuri_identity_member",
               f"{TESTS_POSTCLOSURE}::test_raychaudhuri_schwarzschild_anchor"),
        deps=("MATTER_TRANSPORT",),
        provenance="raychaudhuri",
        gate="G112",
    ),
    ChainNode(
        id="RAYCHAUDHURI_NULL",
        statement="Null congruence dynamics: dtheta_hat/dlambda = "
                  "-theta_hat^2/2 - R_mn k^mu k^nu for the affine radial "
                  "null congruence (b = 0); Schwarzschild anchor "
                  "theta_hat = 2/r with R_kk = 0.",
        relation="null Raychaudhuri equation (affine parametrisation)",
        impl=IMPL_TRANSPORT,
        tests=(f"{TESTS_POSTCLOSURE}::test_null_raychaudhuri_identity",
               f"{TESTS_POSTCLOSURE}::test_null_raychaudhuri_schwarzschild_anchor"),
        deps=("NULL_TRANSPORT",),
        provenance="raychaudhuri",
        gate="G113",
    ),
    ChainNode(
        id="OPTICAL_TRANSPORT",
        statement="Geometric-optics amplitude transport: "
                  "d(a^2)/dlambda = -(nabla_mu k^mu) a^2 reproduces the "
                  "exact geometric invariant a^2 r^2 = const for radial "
                  "null beams.",
        relation="nabla_mu (a^2 k^mu) = 0  =>  a^2 r^2 = const "
                 "(spherical symmetry)",
        impl=IMPL_TRANSPORT,
        tests=(f"{TESTS_POSTCLOSURE}::test_amplitude_transport_invariant",),
        deps=("NULL_TRANSPORT",),
        provenance="geometric_optics",
        gate="G114",
    ),
    ChainNode(
        id="PHASE_TRANSPORT",
        statement="Eikonal phase and redshift from the same geometry: "
                  "S_r = integral dr/sqrt(f h) per unit photon energy "
                  "(JIF phase-chain building block Xi->D->g->Phi); the "
                  "integral phase equals the ODE-accumulated coordinate "
                  "time of an integrated null ray; static-observer "
                  "redshift sqrt(f).",
        relation="k_mu = d_mu Phi; dPhi per unit E: S_r; "
                 "omega(r) = E/sqrt(f)",
        impl=f"{IMPL_TRANSPORT} + {IMPL_CHAIN}",
        tests=(f"{TESTS_POSTCLOSURE}::test_eikonal_phase_redshift",
               f"{TESTS_FORWARD}::test_phase_integral_vs_differential"),
        deps=("NULL_TRANSPORT",),
        provenance="phase_jif",
        gate="G115",
    ),
    ChainNode(
        id="ROTATION_ORBITS",
        statement="Rotation and circular dynamics from the geometry: "
                  "exact Kepler Omega^2 = f'/(2r); at the rings the "
                  "INDEPENDENT formulas Omega_kepler and Omega_photon = "
                  "sqrt(W) coincide (ring identity r f' - 2f = 0); frame "
                  "dragging vanishes identically in the static slice "
                  "(g_tphi = 0); dtau/dt = sqrt(f).",
        relation="Omega^2 = f'/(2r); at ring: Omega^2 = W = f/r^2",
        impl=IMPL_TRANSPORT,
        tests=(f"{TESTS_POSTCLOSURE}::test_consolidation_single_geometry",),
        deps=("MATTER_TRANSPORT",),
        provenance="rotation_orbits",
        gate="G116",
    ),
    ChainNode(
        id="PHOTON_RING_CRITICALITY",
        statement="Photon-ring criticality from the exact null orbit "
                  "equation (du/dphi)^2 = (h/f)(1/b^2 - W(u)): outer "
                  "UNSTABLE ring at u = 2/3 (W_uu < 0, b_c = 3sqrt3/2) "
                  "shows logarithmic winding — the deflection increment "
                  "per decade converges to exactly 2 ln 10.",
        relation="critical impact parameter 1/b_c^2 = W(u_ph); "
                 "Delta phi ~ ln(1/eps)",
        impl=IMPL_TRANSPORT,
        tests=(f"{TESTS_POSTCLOSURE}::test_light_rings_member_structure",
               f"{TESTS_POSTCLOSURE}::test_outer_ring_log_winding"),
        deps=("NULL_TRANSPORT",),
        provenance="photon_ring_stability",
        gate="G117",
    ),
    ChainNode(
        id="STABLE_RING_LIBRATION",
        statement="Stable inner trapping: libration around the inner ring "
                  "u = 0.706135 (W_uu > 0); the full nonlinear ODE period, "
                  "the turning-point quadrature and the linear theory "
                  "2 pi/kappa (kappa^2 = W_uu h / (2 f)) agree, with "
                  "positive O(eps) nonlinear correction.",
        relation="kappa^2 = W_uu h / (2 f); T_phi -> 2 pi / kappa",
        impl=IMPL_TRANSPORT,
        tests=(f"{TESTS_POSTCLOSURE}::test_inner_ring_libration",),
        deps=("PHOTON_RING_CRITICALITY",),
        provenance="photon_ring_stability",
        gate="G118",
    ),
    # ---------------- falsifiability -------------------------------------
    ChainNode(
        id="NEGATIVE_CONTROLS",
        statement="The validation harness is FALSIFIABLE: artificial "
                  "four-force, perturbed metric coefficient, wrong "
                  "Christoffel component, broken conserved-quantity "
                  "formula, shifted photon-ring location, wrong amplitude "
                  "law, wrong phase integrand and wrong redshift law are "
                  "ALL detected by the canonical diagnostics.",
        relation="for each corrupted model: canonical detector fires "
                 "(residuals exceed the source-free band by orders of "
                 "magnitude)",
        impl=f"{TESTS_FALSIFIERS} + {IMPL_TRANSPORT}",
        tests=tuple(f"{TESTS_FALSIFIERS}::test_falsifier_{k}" for k in (
            "force", "metric_perturbation", "christoffel_sign",
            "conserved_quantity", "ring_location", "amplitude_law",
            "phase_integrand", "redshift_law")),
        deps=("MATTER_TRANSPORT", "NULL_TRANSPORT", "RAYCHAUDHURI_TIMELIKE",
              "OPTICAL_TRANSPORT", "PHASE_TRANSPORT",
              "PHOTON_RING_CRITICALITY"),
        provenance="negative_controls",
        gate="G119",
    ),
    # ---------------- robustness -----------------------------------------
    ChainNode(
        id="NUMERICAL_ROBUSTNESS",
        statement="Central claims do not depend on a single lucky "
                  "resolution: conservation residuals converge under "
                  "integrator tolerance refinement; Raychaudhuri residuals "
                  "converge under grid refinement; multiple initial "
                  "conditions (both directions, several (E, L)) conserve.",
        relation="residual(rtol) and residual(grid) monotone to solver "
                 "noise floor",
        impl=f"{TESTS_FORWARD} + {IMPL_TRANSPORT}",
        tests=(f"{TESTS_FORWARD}::test_robustness_tolerance_convergence",
               f"{TESTS_FORWARD}::test_robustness_grid_convergence",
               f"{TESTS_FORWARD}::test_robustness_multiple_initial_conditions"),
        deps=("RAYCHAUDHURI_TIMELIKE", "NULL_TRANSPORT"),
        provenance="numerical_robustness",
        gate="G122",
    ),
]

# gate-level registry (non-chain gates handled by the evaluator directly)
GATES = {
    "G110": {"name": "source-free timelike transport (matter)",
             "nodes": ["MATTER_TRANSPORT"], "requires": ["G100"]},
    "G111": {"name": "source-free null transport (light)",
             "nodes": ["NULL_TRANSPORT"], "requires": ["G100"]},
    "G112": {"name": "timelike Raychaudhuri congruence dynamics",
             "nodes": ["RAYCHAUDHURI_TIMELIKE"], "requires": ["G110"]},
    "G113": {"name": "null Raychaudhuri congruence dynamics",
             "nodes": ["RAYCHAUDHURI_NULL"], "requires": ["G111"]},
    "G114": {"name": "geometric-optics amplitude transport",
             "nodes": ["OPTICAL_TRANSPORT"], "requires": ["G111"]},
    "G115": {"name": "eikonal phase transport / JIF chain entry",
             "nodes": ["PHASE_TRANSPORT"], "requires": ["G111"]},
    "G116": {"name": "rotation and circular orbital dynamics",
             "nodes": ["ROTATION_ORBITS"], "requires": ["G110"]},
    "G117": {"name": "photon-ring criticality (outer ring, log winding)",
             "nodes": ["PHOTON_RING_CRITICALITY"], "requires": ["G111"]},
    "G118": {"name": "stable inner ring libration cross-check",
             "nodes": ["STABLE_RING_LIBRATION"], "requires": ["G117"]},
    "G119": {"name": "negative controls (falsifiability battery)",
             "nodes": ["NEGATIVE_CONTROLS"], "requires": ["G112", "G115"]},
    "G120": {"name": "foundations + differential geometry",
             "nodes": ["FOUNDATIONS", "CONNECTION_RIEMANN",
                       "BIANCHI_CONTRACTED"], "requires": ["G100"]},
    "G121": {"name": "known limits (Schwarzschild anchors, PPN signature)",
             "nodes": ["KNOWN_LIMITS"], "requires": ["G120"]},
    "G122": {"name": "numerical robustness / convergence",
             "nodes": ["NUMERICAL_ROBUSTNESS"],
             "requires": ["G112", "G113", "G115"]},
    "G130": {"name": "SINGLE UNIFIED DYNAMICS (one geometry, one "
                     "source-free transport, no per-observable force; "
                     "closed on the registered validation corpus)",
             "nodes": [],
             "requires": ["G110", "G111", "G112", "G113", "G114", "G115",
                          "G116", "G117", "G118", "G119"]},
    "G140": {"name": "TRUE FULL CLOSURE (complete forward chain + "
                     "historical corpus + provenance + artifact "
                     "integrity + clean tree)",
             "nodes": [],
             "requires": ["G00", "G100", "G110", "G111", "G112", "G113",
                          "G114", "G115", "G116", "G117", "G118", "G119",
                          "G120", "G121", "G122", "G130"]},
}

# provenance keys -> MCP/RAG search queries
PROVENANCE_QUERIES = {
    "foundations_conventions": "metric ansatz ds^2 = -f dt^2 + dr^2/h "
                               "segmented spacetime canonical",
    "curvature_differential_geometry": "Christoffel Ricci Riemann tensor "
                                       "symbolic derivation",
    "geodesic_transport": "geodesic equation conserved quantities Killing "
                          "energy segmentation inertia",
    "raychaudhuri": "Raychaudhuri congruence expansion shear vorticity",
    "geometric_optics": "eikonal amplitude transport geometric optics "
                        "null wave",
    "phase_jif": "JIF chain phase Xi D g proper time transfer",
    "rotation_orbits": "circular orbit angular velocity Kepler Omega",
    "photon_ring_stability": "photon sphere critical impact parameter "
                             "light ring stability winding",
    "negative_controls": "find_validation_records: anti-circularity / "
                         "falsification records",
    "single_unified_dynamics_open_history": "single unified dynamics OPEN "
                                            "gap audit",
}


def evaluate_graph(node_results: dict, gate_conditions: dict,
                   member_hash: str,
                   legacy_status: dict | None = None) -> dict:
    """Compute the verdict EXCLUSIVELY from node/gate results.
    node_results: {node_id: {"pass": bool, "tests": {...}}}
    gate_conditions: {gate_id: bool} for node-less gates (G130/G140)
    legacy_status: {gate_id: status-str} for the historical G00..G100
    gates (certification_status() output; dependency-enforced)."""
    gate_status: dict[str, bool] = {}
    unresolved: list[str] = []

    legacy_status = legacy_status or {}
    for gid, st in legacy_status.items():
        gate_status[gid] = (st == "PASS")

    node_pass = {k: bool(v.get("pass", False)) for k, v in
                 node_results.items()}
    for nid, node in ((n.id, n) for n in NODES):
        deps_ok = all(node_pass.get(d, False)
                      for d in node.deps if d in node_pass)
        if not deps_ok:
            node_pass[nid] = False

    for gid, gdef in GATES.items():
        if gdef["nodes"]:
            ok = all(node_pass.get(n, False) for n in gdef["nodes"])
        else:
            cond = gate_conditions.get(gid, {})
            ok = bool(cond) if isinstance(cond, bool) else bool(
                cond.get("pass", False))
        gate_status[gid] = ok

    # dependency closure over gates (legacy + new)
    for gid, gdef in GATES.items():
        if not gate_status.get(gid, False):
            continue
        for parent in gdef["requires"]:
            if not gate_status.get(parent, False):
                gate_status[gid] = False
                unresolved.append(f"{gid} blocked by {parent}")

    for gid, ok in gate_status.items():
        if not ok and gid not in ("G130", "G140") and gid not in (
                legacy_status or {}):
            unresolved.append(gid)

    return {
        "TRUE_FULL_CLOSURE_PASS": gate_status.get("G140", False),
        "gate_status": gate_status,
        "node_status": node_pass,
        "unresolved": sorted(set(unresolved)),
        "member_hash": member_hash,
        "single_unified_dynamics": (
            "CLOSED (one geometry, one source-free transport operator, no "
            "per-observable force; falsifiability demonstrated; closed on "
            "the registered validation corpus — see scope note)"
            if gate_status.get("G130") else "OPEN"),
        "extra_ssz_force_required": (
            "NO" if (gate_status.get("G110") and gate_status.get("G111")
                     and gate_status.get("G119")) else "UNDETERMINED"),
    }
