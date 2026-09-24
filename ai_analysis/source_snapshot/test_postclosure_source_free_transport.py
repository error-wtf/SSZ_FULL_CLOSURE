"""SSZ SOURCE-FREE TRANSPORT (post-closure research block, NON-GATING).

Dedicated verification that the closed P5 full-action geometry supports
NONTRIVIAL SOURCE-FREE TRANSPORT:  u^nu nabla_nu u^mu = 0 (matter) and
k^nu nabla_nu k^mu = 0 (light) with NO external force anywhere, while the
optical/kinematical structure (theta, sigma, Raychaudhuri, trapping,
phase) is rich.

Core principle (registered): NO new force; D_SSZ/dtau := u^nu nabla_nu
[g_SSZ].  RAG constraints are respected: no flow equation is derived from
the static profile and "single unified dynamics" stays OPEN.

All numbers come from the FROZEN electric production member (member-hash
verified at load).  No fitting anywhere; tolerances below are documented
SOLVER/INTERPOLATION accuracy bounds, not physics gates.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

from ssz_p5.postclosure.transport import (
    SSZMetric,
    amplitude_transport_radial,
    circular_timelike_orbit,
    consolidation_table,
    eikonal_phase_and_redshift,
    find_light_rings,
    forced_transport_control,
    libration_period_ode,
    load_member_metric,
    null_geodesic_transport,
    radial_congruence_scalars,
    radial_null_congruence,
    repo_root,
    ring_trapping,
    timelike_geodesic,
)

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def metric():
    return load_member_metric(ROOT)


def schwarzschild_metric() -> SSZMetric:
    """Exact Schwarzschild (r_s = 1) through the same machinery —
    analytic anchor with R_mu_nu == 0 identically."""
    u = np.linspace(0.50, 0.95, 800)
    f = 1.0 - u
    return SSZMetric(
        u=u, f_u=CubicSpline(u, f), h_u=CubicSpline(u, f.copy()),
        member_hash="ANALYTIC_SCHWARZSCHILD_RS1",
        u_min=float(u.min()), u_max=float(u.max()),
    )


# ---------------------------------------------------------------------------
# member hash / provenance
# ---------------------------------------------------------------------------


def test_member_hash_pinned(metric):
    """The post-closure block consumes the SAME frozen member as all gates."""
    man = json.loads((ROOT / "data/generated/phase2_q2/"
                      "ELECTRIC_PRODUCTION_MEMBER_CURRENT.json").read_text())
    assert metric.member_hash == man["member_hash"]


# ---------------------------------------------------------------------------
# 1. timelike source-free transport  u^nu nabla_nu u^mu = 0
# ---------------------------------------------------------------------------


def test_timelike_source_free_conservation(metric):
    """Radial + angular free fall: Killing energy, angular momentum and
    norm are conserved (solver-independent statement of the vanishing
    force term).  DOP853 rtol 1e-12; norm residual < 1e-6."""
    ck = timelike_geodesic(metric, r0=1.6, E=1.0, L=0.0, tau_span=(0, 1.5))
    assert ck.max_dE < 1e-8
    assert ck.max_dL < 1e-8
    assert ck.max_norm_res < 1e-6
    ck2 = timelike_geodesic(metric, r0=1.6, E=0.95, L=0.8, tau_span=(0, 1.5))
    assert ck2.max_dE < 1e-8
    assert ck2.max_dL < 1e-10
    assert ck2.max_norm_res < 1e-6


def test_timelike_acceleration_crosscheck(metric):
    """Finite-difference acceleration along the worldline matches the
    Christoffel transport term (relative FD tolerance 5e-3; second-order
    gradient on 4001 output points)."""
    ck = timelike_geodesic(metric, r0=1.5, E=0.9, L=0.6, inward=False,
                           tau_span=(0, 1.5))
    assert ck.max_accel_res < 5e-3


def test_negative_control_force_breaks_energy_conservation(metric):
    """NEGATIVE CONTROL: with a synthetic force F^t = 1e-2 the Killing
    energy drifts by O(F_t * tau) — proves the source-free residuals are a
    genuine discriminator (forced drift is >1e6x the geodesic residual)."""
    geodesic = timelike_geodesic(metric, r0=1.6, E=1.0, L=0.0,
                                 tau_span=(0, 1.5))
    forced = forced_transport_control(metric, r0=1.6, F_t=1e-2,
                                      tau_span=(0, 1.5))
    assert forced["max_dE"] > 1e-3
    assert forced["max_dE"] > 1e6 * max(geodesic.max_dE, 1e-12)


def test_circular_orbit_kepler(metric):
    """Circular timelike orbit: exact Kepler Omega^2 = f'/(2r) and
    dtau/dt = sqrt(f - Omega^2 r^2) both finite and positive."""
    orb = circular_timelike_orbit(metric, 1.6)
    f = float(metric.f(1.6))
    fp = float(metric.fp(1.6))
    assert orb["Omega"] ** 2 == pytest.approx(fp / (2.0 * 1.6), rel=1e-12)
    assert orb["dtau_dt"] == pytest.approx(
        np.sqrt(f - 1.6**2 * orb["Omega"] ** 2), rel=1e-12)
    assert orb["dtau_dt"] > 0


# ---------------------------------------------------------------------------
# 2. kinematical scalars + Raychaudhuri (timelike radial congruence)
# ---------------------------------------------------------------------------


def test_raychaudhuri_identity_member(metric):
    """dtheta/dtau = -theta^2/3 - sigma^2 + omega^2 - R_mn u^mu u^nu holds
    on the frozen member; two independent theta routes agree; congruence is
    irrotational (omega^2 < 1e-25); residual converges with grid
    refinement (80 -> 200 points)."""
    g80 = radial_congruence_scalars(metric, np.linspace(1.42, 1.63, 80))
    g200 = radial_congruence_scalars(metric, np.linspace(1.42, 1.63, 200))
    assert g200["theta_routes_absdiff"] < 1e-5
    assert float(np.max(np.abs(g200["omega2"]))) < 1e-25
    assert g200["max_scaled_raychaudhuri_res"] < 1e-3
    assert (g200["max_scaled_raychaudhuri_res"]
            < 0.5 * g80["max_scaled_raychaudhuri_res"])


def test_raychaudhuri_schwarzschild_anchor():
    """Analytic anchor (Ricci-flat, exact closed forms):
    theta = -(3/2) r^{-3/2}, sigma^2 = (3/2) r^{-3}, omega = 0,
    R_mn u^mu u^nu = 0 — the identity must hold to spline accuracy."""
    ms = schwarzschild_metric()
    grid = np.linspace(1.2, 1.9, 200)
    res = radial_congruence_scalars(ms, grid)
    theta_exact = -1.5 * grid ** -1.5
    sigma2_exact = 1.5 * grid ** -3.0
    assert res["max_scaled_raychaudhuri_res"] < 1e-4
    assert np.max(np.abs(res["theta"] - theta_exact)) < 1e-6
    assert np.max(np.abs(res["sigma2"] - sigma2_exact)) < 1e-6
    assert float(np.max(np.abs(res["R_uu"]))) < 1e-8


# ---------------------------------------------------------------------------
# 3. null sector  k^nu nabla_nu k^mu = 0
# ---------------------------------------------------------------------------


def test_null_geodesic_conservation(metric):
    """Null rays with impact parameter (sub-critical 2.0 and near-critical
    2.59): E, L and k.k are conserved to solver accuracy."""
    for b, span in ((2.0, 12.0), (2.59, 20.0)):
        ck = null_geodesic_transport(metric, b=b, r0=1.63, lam_span=(0, span))
        assert ck.max_dE < 1e-9
        assert ck.max_dL < 1e-9
        assert ck.max_norm_res < 1e-9


def test_null_raychaudhuri_identity(metric):
    """Radial null congruence (affine, b=0): dtheta/dlambda =
    -theta^2/2 - R_mn k^mu k^nu on the member to spline accuracy."""
    res = radial_null_congruence(metric, np.linspace(1.42, 1.63, 200))
    assert res["theta_routes_absdiff"] < 1e-8
    assert res["max_scaled_raychaudhuri_res"] < 1e-4


def test_null_raychaudhuri_schwarzschild_anchor():
    """Ricci-flat anchor: theta_hat = 2/r exactly, R_kk = 0, so the
    identity reduces to dtheta/dlambda = -theta^2/2 exactly."""
    ms = schwarzschild_metric()
    res = radial_null_congruence(ms, np.linspace(1.2, 1.9, 200))
    assert np.max(np.abs(res["theta_hat"] - 2.0 / np.linspace(1.2, 1.9, 200))) \
        < 1e-8
    assert float(np.max(np.abs(res["R_kk"]))) < 1e-8
    assert res["max_scaled_raychaudhuri_res"] < 1e-5


def test_amplitude_transport_invariant(metric):
    """Beam amplitude transport d(a^2)/dlambda = -(nabla_mu k^mu) a^2
    reproduces the exact geometric invariant a^2 r^2 = const for radial
    beams (spherical symmetry; error < 1e-9)."""
    out = amplitude_transport_radial(metric, 1.45, 1.63)
    assert out["invariant_error"] < 1e-9
    assert 0.0 < out["a2_ratio"] < 1.0


def test_eikonal_phase_redshift(metric):
    """Reduced eikonal phase integral S_r = integral dr/sqrt(fh) between
    the two light rings is finite/positive; the static-observer redshift
    between the rings equals sqrt(f_a/f_b)."""
    u_in = 0.7061345809124143
    ph = eikonal_phase_and_redshift(metric, 1.0 / u_in, 1.5)
    assert ph["reduced_phase_per_E"] > 0
    f_in = float(metric.f(1.0 / u_in))
    f_out = float(metric.f(1.5))
    assert ph["redshift_a_to_b"] == pytest.approx(np.sqrt(f_in / f_out),
                                                  rel=1e-12)


# ---------------------------------------------------------------------------
# 4. ring structure: trapping from the exact null orbit equation
# ---------------------------------------------------------------------------


def test_light_rings_member_structure(metric):
    """The frozen member carries TWO light rings: the outer unstable
    (W_uu < 0) at the Schwarzschild position u = 2/3 (b_c = 3 sqrt(3)/2,
    Omega_ph^2 = 4/27) and the inner STABLE (W_uu > 0) at
    u = 0.706135 (registered value).  Reconstruction is via spline root
    finding on the SAME member data (independent-reconstruction band
    5e-4 against the registered analytic values)."""
    rings = find_light_rings(metric)
    assert len(rings) == 2
    outer = next(r for r in rings if not r["stable"])
    inner = next(r for r in rings if r["stable"])
    assert outer["u"] == pytest.approx(2.0 / 3.0, abs=5e-4)
    assert outer["W_uu"] < 0
    assert outer["b_crit"] == pytest.approx(3.0 * np.sqrt(3.0) / 2.0,
                                            rel=5e-4)
    assert outer["Omega_ph"] ** 2 == pytest.approx(4.0 / 27.0, rel=5e-4)
    assert inner["u"] == pytest.approx(0.706135, abs=5e-4)
    assert inner["W_uu"] > 0
    # exact structural identities on BOTH rings
    for rg in rings:
        assert rg["Omega_ph"] ** 2 == pytest.approx(rg["W"], rel=1e-12)
        f_r = float(metric.f(rg["r"]))
        assert f_r / rg["r"] ** 2 == pytest.approx(rg["W"], rel=1e-6)


def test_outer_ring_log_winding(metric):
    """Unstable ring: deflection diverges logarithmically — the Delta phi
    increment per decade of turning-point distance converges to
    2 ln 10 (deviation < 0.5% at the smallest decade)."""
    rings = find_light_rings(metric)
    outer = next(r for r in rings if not r["stable"])
    trap = ring_trapping(metric, [outer])
    inc = trap["outer_unstable"]["increments"]
    target = trap["outer_unstable"]["expected_increment_per_decade"]
    assert abs(inc[-1] - target) / target < 5e-3
    # monotone convergence towards the asymptotic increment
    assert abs(inc[-1] - target) < abs(inc[0] - target)


def test_inner_ring_libration(metric):
    """Stable ring: full-ODE libration period at eps = 1e-5 agrees with
    the independent turning-point quadrature to < 1e-6 (rel) and with the
    linear theory 2 pi / kappa to < 1% (nonlinear amplitude correction
    is O(eps) and positive)."""
    rings = find_light_rings(metric)
    inner = next(r for r in rings if r["stable"])
    trap = ring_trapping(metric, [inner])
    t = trap["inner_stable"]
    assert t["libration_period_ode"] == pytest.approx(
        t["delta_phi_series"][0], rel=1e-6)
    lin = t["expected_libration_period_2pi_over_kappa"]
    assert t["libration_period_ode"] == pytest.approx(lin, rel=1e-2)
    assert t["libration_period_ode"] > lin  # positive nonlinear correction


# ---------------------------------------------------------------------------
# consolidation: one geometry -> orbit + time + phase + rotation + trapping
# ---------------------------------------------------------------------------


def test_consolidation_single_geometry(metric):
    """At BOTH rings the independent formulas Omega_kepler = sqrt(f'/(2r))
    and Omega_photon = sqrt(W) coincide to machine precision (exact ring
    identity r f' - 2f = 0); frame dragging vanishes identically in the
    static central slice; time dilation dtau/dt = sqrt(f) is consistent."""
    rings = find_light_rings(metric)
    rows = consolidation_table(metric, rings)
    assert len(rows) == 2
    for row in rows:
        assert row["kepler_vs_photon_rel_diff"] < 1e-9
        assert row["frame_dragging"] == 0.0
        assert row["dtau_dt"] == pytest.approx(
            np.sqrt(float(metric.f(row["r"]))), rel=1e-12)
    by = {row["ring"]: row for row in rows}
    assert by["inner_stable"]["W_uu"] > 0 > by["outer_unstable"]["W_uu"]
    # the stable ring binds photons: larger oscillation exponent
    assert by["inner_stable"]["stability_exponent"] \
        > by["outer_unstable"]["stability_exponent"]


def test_postclosure_artifact_registered():
    """The evidence artifact exists, pins the SAME member hash and carries
    the RAG constraint declarations."""
    path = ROOT / "data/generated/postclosure/SSZ_SOURCE_FREE_TRANSPORT.json"
    assert path.exists()
    d = json.loads(path.read_text())
    man = json.loads((ROOT / "data/generated/phase2_q2/"
                      "ELECTRIC_PRODUCTION_MEMBER_CURRENT.json").read_text())
    assert d["member_hash"] == man["member_hash"]
    for key in ("no_flow_from_static_profile", "single_unified_dynamics_open",
                "early_paper_is_intuition_archive",
                "navier_stokes_structural_only"):
        assert key in d["rag_constraints_respected"]
    assert d["policy"]["fitting"] == "NONE"
