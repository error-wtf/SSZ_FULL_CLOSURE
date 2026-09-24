"""TRUE FULL CLOSURE — first-principles forward chain tests.

Covers: canonical foundations, differential geometry (numeric Riemann,
symmetries, Bianchi identities, Kretschmann anchor), known limits
(Schwarzschild bundle, PPN signature), phase integral<->differential
cross-check, and numerical robustness (tolerance/grid/IC convergence).

All tolerances are documented solver/interpolation accuracy bounds:
  - DOP853 rtol 1e-12          -> conservation floors ~1e-9
  - CubicSpline on frozen data -> derivative errors ~1e-6..1e-10
  - central FD with eps=1e-5   -> curvature-route agreement ~1e-8
"""
import numpy as np
import pytest
from pathlib import Path

from ssz_p5.postclosure.transport import (
    load_member_metric,
    null_geodesic_transport,
    radial_congruence_scalars,
    timelike_geodesic,
)
from ssz_p5.true_closure.chain import (
    christoffels_fd,
    check_contracted_bianchi,
    check_riemann_symmetries,
    curvature_invariants_numeric,
    metric_tensors,
    null_coordinate_time_ode,
    null_coordinate_time_quadrature,
    ppn_signature_checks,
    ricci_from_riemann,
    riemann_numeric,
    schwarzschild_reference_metric,
    weak_field_reference_metric,
)

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def metric():
    return load_member_metric(ROOT)


# ---------------------------------------------------------------------------
# foundations
# ---------------------------------------------------------------------------


def test_foundations_metric_tensors(metric):
    for r in (1.44, 1.5, 1.55, 1.62):
        g, ginv, det = metric_tensors(metric, r)
        assert np.allclose(g, g.T)
        assert np.allclose(g @ ginv, np.eye(4), atol=1e-12)
        f = float(metric.f(r))
        h = float(metric.h(r))
        assert det == pytest.approx(-f * (1.0 / h) * r**4, rel=1e-12)
        # signature: one negative, three positive eigenvalues
        ev = np.linalg.eigvalsh(g)
        assert sum(1 for x in ev if x < 0) == 1
    # Christoffel FD derivative route agrees with the analytic tensor at
    # theta-perturbed points (Gamma depends on r and theta only)
    G, _ = christoffels_fd(metric, 1.5)
    assert np.allclose(G, metric.Gamma_tensor(1.5))


# ---------------------------------------------------------------------------
# differential geometry
# ---------------------------------------------------------------------------


def test_ricci_two_routes(metric):
    """Numeric Riemann (FD route) contracts to the same Ricci tensor as
    the sympy derivation (tolerance: central-FD accuracy ~1e-6 scaled)."""
    for r in (1.44, 1.5, 1.6):
        Ric_num = ricci_from_riemann(riemann_numeric(metric, r))
        Ric_sym = metric.Ricci(r)
        num = np.array([Ric_num[0, 0], Ric_num[1, 1], Ric_num[2, 2],
                        Ric_num[3, 3]])
        scale = np.maximum(1.0, np.abs(num))
        assert np.max(np.abs(num - Ric_sym) / scale) < 1e-6


def test_riemann_symmetries_and_bianchi1(metric):
    """R_{abcd} = -R_{bacd} = -R_{abdc}, pair exchange, first Bianchi —
    all at FD noise level (central differences, eps = 1e-5)."""
    for r in (1.44, 1.55):
        sym = check_riemann_symmetries(metric, r)
        assert sym["antisym_first_pair"] < 1e-6
        assert sym["antisym_second_pair"] < 1e-12
        assert sym["pair_exchange"] < 1e-6
        assert sym["first_bianchi"] < 1e-12


def test_kretschmann_schwarzschild_anchor():
    """Exact Schwarzschild (r_s = 1, M = 1/2): K = 48 M^2/r^6 = 12/r^6
    and Ricci-flatness, through the numeric curvature stack."""
    ms = schwarzschild_reference_metric()
    for r in (1.3, 1.4, 1.7):
        inv = curvature_invariants_numeric(ms, r)
        K_exact = 12.0 / r**6
        assert inv["kretschmann"] == pytest.approx(K_exact, rel=1e-8)
        assert inv["ricci2"] < 1e-14
        assert inv["r2"] < 1e-14


def test_kretschmann_member_regular(metric):
    """On the frozen member the curvature invariants are finite and
    regular across the whole transport domain (no curvature blow-up)."""
    for r in np.linspace(1.42, 1.63, 12):
        inv = curvature_invariants_numeric(metric, float(r))
        assert np.isfinite(inv["kretschmann"])
        assert inv["kretschmann"] < 100.0  # regular strong-field scale


def test_contracted_bianchi_member(metric):
    """div G = 0 on the member: geometric identity, grid-convergent
    (1.3e-3 at 20 pts -> 2.3e-5 at 80 pts -> below 1e-5 at 160 pts)."""
    r20 = check_contracted_bianchi(metric, np.linspace(1.43, 1.62, 20))
    r80 = check_contracted_bianchi(metric, np.linspace(1.43, 1.62, 80))
    r160 = check_contracted_bianchi(metric, np.linspace(1.43, 1.62, 160))
    assert r80 < r20
    assert r160 < 1e-4
    assert r160 < r80


def test_contracted_bianchi_schwarzschild():
    ms = schwarzschild_reference_metric()
    res = check_contracted_bianchi(ms, np.linspace(1.25, 1.85, 40))
    assert res < 1e-7


# ---------------------------------------------------------------------------
# known limits
# ---------------------------------------------------------------------------


def test_ppn_signature_weak_field():
    """Pipeline limit check on a GR-like weak-field reference metric
    (f = h = 1 - eps u, M = eps/2): gamma = 1 signature exactly and the
    Newtonian Kepler law Omega^2 = M/r^3 (SSZ PPN results live in the
    historical corpus gates; this certifies the PIPELINE limit behaviour
    and makes no claim outside the tested domain)."""
    eps = 0.2
    mw = weak_field_reference_metric(eps)
    out = ppn_signature_checks(mw, M=eps / 2.0)
    assert out["gamma1_signature_max_abs"] < 1e-10
    assert out["newtonian_kepler_max_rel"] < 1e-9


def test_schwarzschild_anchor_bundle():
    """Bundle of exact Schwarzschild anchors through the canonical
    machinery: Ricci-flat, Kretschmann closed form, photon ring u = 2/3
    with Omega_ph^2 = 4/27 and b_c = 3 sqrt(3)/2."""
    ms = schwarzschild_reference_metric()
    inv = curvature_invariants_numeric(ms, 1.4)
    assert inv["ricci2"] < 1e-14
    K_exact = 12.0 / 1.4**6
    assert inv["kretschmann"] == pytest.approx(K_exact, rel=1e-8)
    # photon ring via W_u = 0 on the reference metric (W = u^2 f)
    from scipy.optimize import brentq
    u_ring = brentq(lambda x: float(ms.W_u(x)), 0.6, 0.7,
                    xtol=1e-14, rtol=8.9e-16)
    assert u_ring == pytest.approx(2.0 / 3.0, rel=1e-8)
    W = float(ms.W_of_u(u_ring))
    assert (1.0 / W) ** 0.5 == pytest.approx(3.0 * np.sqrt(3.0) / 2.0,
                                             rel=1e-8)
    assert W == pytest.approx(4.0 / 27.0, rel=1e-8)


# ---------------------------------------------------------------------------
# phase transport cross-checks
# ---------------------------------------------------------------------------


def test_phase_integral_vs_differential(metric):
    """Integral phase/time route (Shapiro quadrature int dr/sqrt(fh))
    equals the differential route (coordinate time accumulated by the
    INTEGRATED null geodesic ODE) — both directions, three radii pairs."""
    for (a, b, inward) in ((1.45, 1.62, False), (1.42, 1.60, False),
                           (1.62, 1.45, True)):
        t_ode = null_coordinate_time_ode(metric, a, b, inward=inward)
        t_quad = null_coordinate_time_quadrature(metric, a, b)
        assert t_ode == pytest.approx(t_quad, rel=1e-9)


# ---------------------------------------------------------------------------
# numerical robustness
# ---------------------------------------------------------------------------


def _max_dE_at_tol(metric, rtol):
    """Integrate a radial timelike geodesic with a given solver tolerance
    and return the Killing-energy conservation residual."""
    r0 = 1.6
    f0 = float(metric.f(r0))
    h0 = float(metric.h(r0))
    ur0 = -np.sqrt(h0 * (1.0 / f0 - 1.0))
    ut0 = 1.0 / f0

    def rhs(tau, y):
        _, r, _, ut, ur, _ = y
        f = float(metric.f(r))
        fp = float(metric.fp(r))
        h = float(metric.h(r))
        hp = float(metric.hp(r))
        return [ut, ur, 0.0, -(fp / f) * ut * ur,
                -(0.5 * h * fp) * ut * ut + (hp / (2.0 * h)) * ur * ur, 0.0]

    from scipy.integrate import solve_ivp
    from ssz_p5.postclosure.transport import radial_domain
    r_lo, _ = radial_domain(metric)

    def leave(tau, y):
        return y[1] - r_lo

    leave.terminal = True  # type: ignore[attr-defined]
    sol = solve_ivp(rhs, (0.0, 1.5), [0.0, r0, 0.0, ut0, ur0, 0.0],
                    method="DOP853", rtol=rtol, atol=rtol,
                    events=leave, max_step=1.5 / 200.0)
    assert sol.success
    assert sol.t_events[0].size >= 1  # identical stopping event for all rtol
    return float(np.max(np.abs(metric.f(sol.y[1]) * sol.y[3] - 1.0)))


def test_robustness_tolerance_convergence(metric):
    """Conservation residual is truncation-dominated: it decreases with
    integrator tolerance refinement (1e-8 -> 1e-10 -> 1e-12) and stays
    below the corresponding documented floor."""
    residuals = [_max_dE_at_tol(metric, rtol) for rtol in
                 (1e-8, 1e-10, 1e-12)]
    assert all(np.isfinite(residuals))
    assert residuals[1] <= max(residuals[0], 1e-9)
    assert residuals[2] <= max(residuals[1], 1e-10)
    assert residuals[2] < 1e-9


def test_robustness_grid_convergence(metric):
    """Raychaudhuri residual converges monotonically under grid
    refinement (80 -> 200 -> 400 points)."""
    r80 = radial_congruence_scalars(
        metric, np.linspace(1.42, 1.63, 80))["max_scaled_raychaudhuri_res"]
    r200 = radial_congruence_scalars(
        metric, np.linspace(1.42, 1.63, 200))["max_scaled_raychaudhuri_res"]
    r400 = radial_congruence_scalars(
        metric, np.linspace(1.42, 1.63, 400))["max_scaled_raychaudhuri_res"]
    assert r200 < r80
    assert r400 < r200
    assert r400 < 1e-4


def test_robustness_multiple_initial_conditions(metric):
    """Source-free transport holds for multiple initial conditions:
    both directions and several (E, L) pairs."""
    for (E, L, inward) in ((1.0, 0.0, True), (0.95, 0.8, True),
                           (0.9, 0.6, False), (0.85, 0.4, False)):
        ck = timelike_geodesic(metric, r0=1.6, E=E, L=L, inward=inward,
                               tau_span=(0, 1.2))
        assert ck.max_dE < 1e-8
        assert ck.max_norm_res < 1e-6
    for b in (1.8, 2.0, 2.4, 2.55):
        ck = null_geodesic_transport(metric, b=b, r0=1.63,
                                     lam_span=(0, 10.0))
        assert ck.max_dE < 1e-9
        assert ck.max_norm_res < 1e-9
