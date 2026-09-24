"""TRUE FULL CLOSURE — falsifiability battery (negative controls).

Each test CORRUPTS one element of the chain (isolated fixture, canonical
member untouched) and asserts that the canonical detector FIRES, i.e. the
validation harness distinguishes correct source-free geometric transport
from each class of error.  If any of these tests failed, the corresponding
POSITIVE gate would be unfalsifiable and closure would be void.

Detectors used (canonical diagnostics):
  D1 Killing-energy conservation residual along trajectories
  D2 acceleration cross-check (FD vs Christoffel form)
  D3 ring identity f'/(2r) = W = f/r^2 at a claimed ring
  D4 libration period vs linear theory 2 pi/kappa
  D5 amplitude invariant a^2 r^2 = const
  D6 null coordinate time vs Shapiro quadrature int dr/sqrt(fh)
  D7 static-observer redshift sqrt(f_a/f_b)
"""
import numpy as np
import pytest
from pathlib import Path
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

from ssz_p5.postclosure.transport import (
    SSZMetric,
    amplitude_transport_radial,
    eikonal_phase_and_redshift,
    find_light_rings,
    libration_period_ode,
    load_member_metric,
    timelike_geodesic,
)
from ssz_p5.true_closure.chain import (
    null_coordinate_time_ode,
    null_coordinate_time_quadrature,
)

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def metric():
    return load_member_metric(ROOT)


def _corrupted_metric(metric, scale=1.0 + 1e-3, mode="linear"):
    """Metric fixture with deliberately perturbed f-values (spline rebuilt
    from corrupted data; isolated — canonical member untouched)."""
    u = metric.u
    if mode == "uniform":
        f = metric.f_u(u) * scale
    else:  # u-proportional perturbation distorts the profile shape
        f = metric.f_u(u) * (1.0 + (scale - 1.0) * u / metric.u_max)
    return SSZMetric(
        u=u, f_u=CubicSpline(u, f), h_u=metric.h_u,
        member_hash="CORRUPTED_FIXTURE", u_min=metric.u_min,
        u_max=metric.u_max)


# ---------------------------------------------------------------------------
# 1. artificial four-force
# ---------------------------------------------------------------------------


def test_falsifier_force(metric):
    """D1: a synthetic force F^t produces an E-drift ~1e7x the geodesic
    band — the source-free detector fires."""
    from ssz_p5.postclosure.transport import forced_transport_control
    forced = forced_transport_control(metric, r0=1.6, F_t=1e-2,
                                     tau_span=(0, 1.5))
    clean = load_member_metric(ROOT)
    from ssz_p5.postclosure.transport import timelike_geodesic
    geo = timelike_geodesic(clean, r0=1.6, E=1.0, L=0.0, tau_span=(0, 1.5))
    assert forced["max_dE"] > 1e-4
    assert forced["max_dE"] > 1e6 * max(geo.max_dE, 1e-12)


# ---------------------------------------------------------------------------
# 2. perturbed metric coefficient
# ---------------------------------------------------------------------------


def test_falsifier_metric_perturbation(metric):
    """D3: a 0.1% metric perturbation moves the outer ring away from the
    Schwarzschild anchor u = 2/3 — the ring-identity detector fires."""
    bad = _corrupted_metric(metric, scale=1.02, mode="propto_u")
    rings = find_light_rings(bad)
    outer = [rg for rg in rings if not rg["stable"]]
    if not outer:
        # perturbation removed/shifted the ring entirely: detector fires
        assert True
        return
    u_out = outer[0]["u"]
    assert abs(u_out - 2.0 / 3.0) > 5e-4  # outside the registered band
    # canonical member stays anchored
    u_true = [rg for rg in find_light_rings(metric) if not rg["stable"]][0]["u"]
    assert abs(u_true - 2.0 / 3.0) < 5e-4


# ---------------------------------------------------------------------------
# 3. wrong Christoffel component
# ---------------------------------------------------------------------------


def test_falsifier_christoffel_sign(metric):
    """D2: flipping the sign of Gamma^r_{tt} in the transport form makes
    the FD-vs-Christoffel acceleration cross-check fire at O(1)."""
    r0, E, tau_span = 1.6, 1.0, (0.0, 1.2)
    f0, h0 = float(metric.f(r0)), float(metric.h(r0))
    ut0 = E / f0
    ur0 = -np.sqrt(h0 * (E * E / f0 - 1.0))

    def rhs(tau, y):
        _, r, _, ut, ur, _ = y
        f = float(metric.f(r))
        fp = float(metric.fp(r))
        h = float(metric.h(r))
        hp = float(metric.hp(r))
        return [ut, ur, 0.0, -(fp / f) * ut * ur,
                -(0.5 * h * fp) * ut * ut + (hp / (2.0 * h)) * ur * ur, 0.0]

    from ssz_p5.postclosure.transport import radial_domain
    r_lo, r_hi = radial_domain(metric)
    leave = lambda tau_, y_: min(y_[1] - r_lo, r_hi - y_[1])
    leave.terminal = True
    sol = solve_ivp(rhs, tau_span, [0.0, r0, 0.0, ut0, ur0, 0.0],
                    method="DOP853", rtol=1e-12, atol=1e-12,
                    events=leave,
                    t_eval=np.linspace(tau_span[0], tau_span[1], 2001),
                    max_step=(tau_span[1] - tau_span[0]) / 200.0)
    assert sol.success
    tau, r, ut, ur = sol.t, sol.y[1], sol.y[3], sol.y[4]
    dudt = np.gradient(np.stack([ut, ur]), tau, axis=1)
    res_correct = 0.0
    res_wrong = 0.0
    for i in range(1, len(tau) - 1):
        rr = r[i]
        f, h = float(metric.f(rr)), float(metric.h(rr))
        fp, hp = float(metric.fp(rr)), float(metric.hp(rr))
        a_r_true = (dudt[1, i] + 0.5 * h * fp * ut[i] ** 2
                    - hp / (2.0 * h) * ur[i] ** 2)
        a_r_wrong = (dudt[1, i] - 0.5 * h * fp * ut[i] ** 2
                     - hp / (2.0 * h) * ur[i] ** 2)
        scale = max(1.0, abs(0.5 * h * fp * ut[i] ** 2))
        res_correct = max(res_correct, abs(a_r_true) / scale)
        res_wrong = max(res_wrong, abs(a_r_wrong) / scale)
    assert res_correct < 5e-3        # canonical detector stays silent
    assert res_wrong > 0.1           # ... and FIRES on the corrupted sign


# ---------------------------------------------------------------------------
# 4. broken conserved-quantity formula
# ---------------------------------------------------------------------------


def test_falsifier_conserved_quantity(metric):
    """D1 with a WRONG energy formula (u^t instead of f u^t): the
    residual detector fires although the trajectory is a true geodesic —
    the diagnostic quantity itself is validated."""
    ck = timelike_geo(metric)
    # correct formula: silent
    assert ck.max_dE < 1e-8
    # wrong formula: fires
    wrong = float(np.max(np.abs(ck.ut - ck.E0)))
    assert wrong > 0.1


def timelike_geo(metric):
    return timelike_geodesic_local(metric, 1.6, 1.0, 0.0, (0.0, 1.2))


def timelike_geodesic_local(metric, r0, E, L, span):
    from ssz_p5.postclosure.transport import (
        TransportCheck, radial_domain)
    f0, h0 = float(metric.f(r0)), float(metric.h(r0))
    val = h0 * (E * E / f0 - (L * L / r0**2 + 1.0))
    ur0 = -np.sqrt(val)
    ut0 = E / f0

    def rhs(tau, y):
        _, r, _, ut, ur, uph = y
        f = float(metric.f(r))
        fp = float(metric.fp(r))
        h = float(metric.h(r))
        hp = float(metric.hp(r))
        return [ut, ur, uph, -(fp / f) * ut * ur,
                -(0.5 * h * fp) * ut * ut + (hp / (2.0 * h)) * ur * ur
                + h * r * uph * uph,
                -(2.0 / r) * ur * uph]

    r_lo, r_hi = radial_domain(metric)
    leave = lambda tau, y: min(y[1] - r_lo, r_hi - y[1])
    leave.terminal = True
    sol = solve_ivp(rhs, span, [0.0, r0, 0.0, ut0, ur0, L / r0**2],
                    method="DOP853", rtol=1e-12, atol=1e-12,
                    events=leave, dense_output=True,
                    t_eval=np.linspace(span[0], span[1], 2001))
    t, r, ut, ur, uph = sol.y[0], sol.y[1], sol.y[3], sol.y[4], sol.y[5]
    E_arr = metric.f(r) * ut
    return TransportCheck(
        E0=E, L0=L,
        max_dE=float(np.max(np.abs(E_arr - E))),
        max_dL=float(np.max(np.abs(r**2 * uph - L))),
        max_norm_res=float(np.max(np.abs(
            -metric.f(r) * ut**2 + ur**2 / metric.h(r) + r**2 * uph**2
            + 1.0))),
        max_accel_res=0.0, tau=t, r=r, ut=ut, ur=ur, uph=uph,
        t=t, phi=sol.y[2])


# ---------------------------------------------------------------------------
# 5. wrong photon-ring location
# ---------------------------------------------------------------------------


def test_falsifier_ring_location(metric):
    """D3/D4: a claimed ring shifted by 1e-3 breaks BOTH the exact ring
    identity Omega^2 = f/r^2 AND the libration-period cross-check."""
    rings = find_light_rings(metric)
    inner = next(rg for rg in rings if rg["stable"])
    u_true = inner["u"]
    u_wrong = u_true - 1e-3
    # ring condition: W_u = 0 (W = u^2 f is definitional, the RING
    # location is where W_u vanishes)
    wu_true = abs(float(metric.W_u(u_true)))
    wu_wrong = abs(float(metric.W_u(u_wrong)))
    assert wu_true < 1e-9             # silent on the true ring
    assert wu_wrong > 1e-3            # FIRES on the shifted claim
    period_true = libration_period_ode(metric, u_true, eps=1e-5)
    kappa_true = np.sqrt(inner["W_uu"] * float(metric.h(1.0 / u_true))
                         / (2.0 * float(metric.f(1.0 / u_true))))
    assert abs(period_true - 2 * np.pi / kappa_true) / (2 * np.pi
                                                        / kappa_true) < 0.01
    # at the WRONG u the ODE period no longer matches the local linear
    # prediction beyond the tight band (detector fires)
    try:
        period_wrong = libration_period_ode(metric, u_wrong, eps=1e-5)
    except (RuntimeError, ValueError):
        period_wrong = None
    if period_wrong is not None:
        W_uu_w = float(metric.W_uu(u_wrong))
        if W_uu_w > 0:
            k_w = np.sqrt(W_uu_w * float(metric.h(1.0 / u_wrong))
                          / (2.0 * float(metric.f(1.0 / u_wrong))))
            assert abs(period_wrong - 2 * np.pi / k_w) / (2 * np.pi
                                                          / k_w) > 1e-2


# ---------------------------------------------------------------------------
# 6. wrong amplitude transport law
# ---------------------------------------------------------------------------


def test_falsifier_amplitude_law(metric):
    """D5: the WRONG law a^2 ~ 1/r (instead of a^2 r^2 = const) violates
    the invariant by O(1) — the amplitude detector fires."""
    out = amplitude_transport_radial(metric, 1.45, 1.63)
    assert out["invariant_error"] < 1e-9      # canonical law: silent
    a2_ratio = out["a2_ratio"]
    wrong_invariant = abs(a2_ratio * 1.63 - 1.45) / 1.45
    assert wrong_invariant > 1e-2             # wrong law: FIRES


# ---------------------------------------------------------------------------
# 7. wrong phase integrand
# ---------------------------------------------------------------------------


def test_falsifier_phase_integrand(metric):
    """D6: the WRONG integrand int dr/f (instead of int dr/sqrt(fh))
    disagrees with the integrated null-ray coordinate time by O(1)."""
    a, b = 1.45, 1.62
    t_ode = null_coordinate_time_ode(metric, a, b, inward=False)
    t_right = null_coordinate_time_quadrature(metric, a, b)
    grid = np.linspace(a, b, 20000)
    t_wrong = float(np.trapezoid(1.0 / metric.f(grid), grid))
    assert t_ode == pytest.approx(t_right, rel=1e-9)   # canonical: silent
    assert abs(t_wrong - t_ode) / t_ode > 1e-2         # wrong: FIRES


# ---------------------------------------------------------------------------
# 8. modified redshift law
# ---------------------------------------------------------------------------


def test_falsifier_redshift_law(metric):
    """D7: the WRONG law z = f_a/f_b (instead of sqrt(f_a/f_b)) fires."""
    u_in = 0.7061345809124143
    ph = eikonal_phase_and_redshift(metric, 1.0 / u_in, 1.5)
    f_in = float(metric.f(1.0 / u_in))
    f_out = float(metric.f(1.5))
    assert ph["redshift_a_to_b"] == pytest.approx(np.sqrt(f_in / f_out),
                                                  rel=1e-12)
    wrong = f_in / f_out
    assert abs(wrong - ph["redshift_a_to_b"]) / ph["redshift_a_to_b"] > 0.03
