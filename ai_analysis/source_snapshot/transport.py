"""Source-free covariant transport on the frozen SSZ P5 geometry.

Geometry:  ds^2 = -f(r) dt^2 + dr^2/h(r) + r^2 dOmega^2,  u = r_s/r, r_s = 1,
f and h are PCHIP interpolants of the FROZEN electric production member
(ELECTRIC_PRODUCTION_MEMBER_CURRENT.csv; hash pinned in its manifest).
No parameter is fitted anywhere: every observable below is derived from
that one geometry by standard differential geometry.

Computations (registered research block):
  1. timelike source-free transport  u^nu nabla_nu u^mu = 0
     (radial + circular; Killing-energy / angular-momentum / norm residuals
     are solver-independent statements of the vanishing force term),
  2. kinematical scalars theta, sigma_mn, omega_mn and the Raychaudhuri
     identity for the radial congruence; inner (W_uu > 0) vs outer
     (W_uu < 0) light-ring comparison,
  3. eikonal/null sector  k^nu nabla_nu k^mu = 0: photon spheres, critical
     impact parameter, circular-null Raychaudhuri balance at both rings,
     amplitude transport a^2 r^2 = const (radial beams), Shapiro delay,
     static-observer redshift and the reduced eikonal phase integral
     S_r = integral dr / sqrt(f h),
  4. consolidation table: orbit + time + phase + rotation + trapping from
     ONE geometry.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline, PchipInterpolator
from scipy.optimize import brentq

MEMBER_CSV = Path("data/generated/phase2_q2/ELECTRIC_PRODUCTION_MEMBER_CURRENT.csv")
MEMBER_MANIFEST = Path(
    "data/generated/phase2_q2/ELECTRIC_PRODUCTION_MEMBER_CURRENT.json"
)


def repo_root(start: Path | None = None) -> Path:
    p = (start or Path(__file__)).resolve()
    for cand in (p, *p.parents):
        if (cand / "MODEL_LOCK.json").exists():
            return cand
    raise FileNotFoundError("repo root not found")


# ---------------------------------------------------------------------------
# geometry model
# ---------------------------------------------------------------------------


@dataclass
class SSZMetric:
    """Cubic-spline model of the frozen member metric (static spherical
    slice).  The splines interpolate the FROZEN member values exactly
    (no fitting); cubic splines are used because transport diagnostics
    need clean second derivatives."""

    u: np.ndarray
    f_u: CubicSpline
    h_u: CubicSpline
    member_hash: str
    u_min: float
    u_max: float

    # -- coordinate helpers -------------------------------------------------
    def f(self, r):
        return self.f_u(1.0 / np.asarray(r, dtype=float))

    def h(self, r):
        return self.h_u(1.0 / np.asarray(r, dtype=float))

    def fp(self, r):
        """df/dr = -u^2 df/du (chain rule, exact for r_s = 1)."""
        u = 1.0 / np.asarray(r, dtype=float)
        return -(u**2) * self.f_u.derivative()(u)

    def hp(self, r):
        u = 1.0 / np.asarray(r, dtype=float)
        return -(u**2) * self.h_u.derivative()(u)

    def fpp(self, r):
        """d2f/dr2 = u^4 f_uu + 2 u^3 f_u."""
        u = 1.0 / np.asarray(r, dtype=float)
        return (
            u**4 * self.f_u.derivative(2)(u)
            + 2.0 * u**3 * self.f_u.derivative()(u)
        )

    def hpp(self, r):
        u = 1.0 / np.asarray(r, dtype=float)
        return (
            u**4 * self.h_u.derivative(2)(u)
            + 2.0 * u**3 * self.h_u.derivative()(u)
        )

    def W_of_u(self, u):
        """W(u) = u^2 f(u): optical potential; ring condition W_u = 0."""
        u = np.asarray(u, dtype=float)
        return u**2 * self.f_u(u)

    def W_u(self, u):
        u = np.asarray(u, dtype=float)
        return 2.0 * u * self.f_u(u) + u**2 * self.f_u.derivative()(u)

    def W_uu(self, u):
        u = np.asarray(u, dtype=float)
        return (
            2.0 * self.f_u(u)
            + 4.0 * u * self.f_u.derivative()(u)
            + u**2 * self.f_u.derivative(2)(u)
        )

    # -- Christoffel symbols -------------------------------------------------
    def Gamma_tensor(self, r, theta=np.pi / 2.0) -> np.ndarray:
        """Full symmetric Gamma[lam, mu, nu] in coordinates (t, r, th, ph)."""
        f = float(self.f(r))
        h = float(self.h(r))
        fp = float(self.fp(r))
        hp = float(self.hp(r))
        s, c = float(np.sin(theta)), float(np.cos(theta))
        G = np.zeros((4, 4, 4))
        G[0, 0, 1] = G[0, 1, 0] = fp / (2.0 * f)
        G[1, 0, 0] = 0.5 * h * fp
        G[1, 1, 1] = -hp / (2.0 * h)
        G[1, 2, 2] = -h * r
        G[1, 3, 3] = -h * r * s**2
        G[2, 1, 2] = G[2, 2, 1] = 1.0 / r
        G[3, 1, 3] = G[3, 3, 1] = 1.0 / r
        G[2, 3, 3] = -s * c
        G[3, 2, 3] = G[3, 3, 2] = c / s
        return G

    # -- Ricci tensor (symbolically derived once, lambdified) ----------------
    def Ricci(self, r) -> np.ndarray:
        """(R_tt, R_rr, R_thth, R_phph) as coordinate components."""
        args = (
            r,
            float(self.f(r)), float(self.fp(r)), float(self.fpp(r)),
            float(self.h(r)), float(self.hp(r)), float(self.hpp(r)),
        )
        return np.array(_ricci_func()(*args), dtype=float)


@lru_cache(maxsize=1)
def _ricci_func():
    """Ricci tensor of ds^2 = -f dt^2 + dr^2/h + r^2 dOmega^2, derived
    symbolically with sympy (nothing hardcoded) and lambdified to numpy."""
    coords = sp.symbols("t r th ph", real=True)
    t, r, th, ph = coords
    fr = sp.Function("f")(r)
    hr = sp.Function("h")(r)
    fp, fpp, hp, hpp = sp.symbols("fp fpp hp hpp", real=True)
    fnum, hnum = sp.symbols("fnum hnum", positive=True)

    g = sp.diag(-fr, 1 / hr, r**2, r**2 * sp.sin(th) ** 2)
    g_inv = sp.simplify(g.inv())
    n = 4
    Gam: dict[tuple[int, int, int], sp.Expr] = {}
    for a in range(n):
        for b in range(n):
            for c in range(b, n):
                expr = sp.simplify(sum(
                    g_inv[a, d] * (
                        sp.diff(g[d, c], coords[b])
                        + sp.diff(g[d, b], coords[c])
                        - sp.diff(g[b, c], coords[d])
                    )
                    for d in range(n)) / 2)
                if expr != 0:
                    Gam[(a, b, c)] = expr
                    if b != c:
                        Gam[(a, c, b)] = expr

    Ric: dict[tuple[int, int], sp.Expr] = {}
    for a in range(n):
        for b in range(a, n):
            # d_rho Gamma^rho_{ab}: each summand is differentiated w.r.t.
            # its own contracted upper-index coordinate
            expr = sum((Gam.get((m, a, b), 0) if Gam.get((m, a, b), 0) != 0
                        else sp.S.Zero).diff(coords[m])
                       for m in range(n))
            # - d_b Gamma^rho_{a rho} (trace, differentiated w.r.t. x^b)
            expr -= sp.Add(
                *[Gam.get((m, a, m), 0) for m in range(n)]
            ).diff(coords[b])
            expr += sum(
                Gam.get((m, a, b), 0) * Gam.get((s2, m, s2), 0)
                for m in range(n) for s2 in range(n))
            expr -= sum(
                Gam.get((s2, a, m), 0) * Gam.get((m, b, s2), 0)
                for m in range(n) for s2 in range(n))
            expr = sp.simplify(expr)
            if expr != 0:
                Ric[(a, b)] = expr

    subs = {
        sp.Derivative(fr, r, r): fpp, sp.Derivative(fr, r): fp,
        sp.Derivative(hr, r, r): hpp, sp.Derivative(hr, r): hp,
        fr: fnum, hr: hnum,
    }
    eq = sp.pi / 2
    comps = []
    for key in ((0, 0), (1, 1), (2, 2), (3, 3)):
        e = sp.trigsimp(sp.expand_trig(Ric.get(key, 0))).subs(subs)
        e = sp.simplify(e.subs(th, eq))
        comps.append(e)
    for e in comps:
        if th in e.free_symbols:
            raise RuntimeError("equatorial reduction left angular dependence")
    func = sp.lambdify(
        (r, fnum, fp, fpp, hnum, hp, hpp),
        comps,
        "numpy",
    )
    return func


@lru_cache(maxsize=None)
def load_member_metric(root: Path | None = None) -> SSZMetric:
    """Load the frozen electric production member (hash-verified)."""
    base = (root or repo_root()).resolve()
    csv = base / MEMBER_CSV
    man = json.loads((base / MEMBER_MANIFEST).read_text())
    digest = hashlib.sha256(csv.read_text().encode("utf-8")).hexdigest()
    if digest != man["member_hash"]:
        raise ValueError("member CSV does not match pinned member_hash")
    frame = pd.read_csv(csv)
    d = frame.sort_values("u").reset_index(drop=True)
    u = d.u.to_numpy(float)
    return SSZMetric(
        u=u,
        f_u=CubicSpline(u, d.f.to_numpy(float)),
        h_u=CubicSpline(u, d.h.to_numpy(float)),
        member_hash=man["member_hash"],
        u_min=float(u.min()),
        u_max=float(u.max()),
    )


def radial_domain(m: SSZMetric, margin: float = 1e-3) -> tuple[float, float]:
    """Open radial interval covered by the frozen member (no extrapolation)."""
    return (1.0 / m.u_max + margin, 1.0 / m.u_min - margin)


# ---------------------------------------------------------------------------
# 1. timelike source-free transport: u^nu nabla_nu u^mu = 0
# ---------------------------------------------------------------------------


@dataclass
class TransportCheck:
    """Solver-independent residuals of the source-free transport equation.

    If (and only if) u^nu nabla_nu u^mu = 0 holds, the Killing energy
    E = f u^t, the axial angular momentum L = r^2 u^phi and the norm
    g_mn u^mu u^nu are constant along the worldline."""

    E0: float
    L0: float
    max_dE: float
    max_dL: float
    max_norm_res: float
    max_accel_res: float
    tau: np.ndarray = field(default=None, repr=False)
    r: np.ndarray = field(default=None, repr=False)
    t: np.ndarray = field(default=None, repr=False)
    phi: np.ndarray = field(default=None, repr=False)
    ut: np.ndarray = field(default=None, repr=False)
    ur: np.ndarray = field(default=None, repr=False)
    uph: np.ndarray = field(default=None, repr=False)


def timelike_geodesic(
    m: SSZMetric,
    r0: float,
    E: float = 1.0,
    L: float = 0.0,
    inward: bool = True,
    tau_span: tuple = (0.0, 10.0),
) -> TransportCheck:
    """Integrate a timelike geodesic in the equatorial plane and quantify
    the source-free transport property through conserved quantities.

    Integration is restricted to the radial interval covered by the frozen
    member (terminal event on domain exit; no extrapolation anywhere)."""
    r_lo, r_hi = radial_domain(m)
    if not (r_lo < r0 < r_hi):
        raise ValueError(f"r0={r0} outside member domain {r_lo, r_hi}")
    f0 = float(m.f(r0))
    h0 = float(m.h(r0))
    val = h0 * (E * E / f0 - (L * L / r0**2 + 1.0))
    if val < 0:
        raise ValueError("timelike geodesic forbidden at r0 for given (E, L)")
    ur0 = -np.sqrt(val) if inward else np.sqrt(val)
    ut0 = E / f0

    def rhs(tau, y):
        _, r, _, ut, ur, uph = y
        f = float(m.f(r))
        fp = float(m.fp(r))
        h = float(m.h(r))
        hp = float(m.hp(r))
        return [
            ut,
            ur,
            uph,
            -(fp / f) * ut * ur,
            -(0.5 * h * fp) * ut * ut
            + (hp / (2.0 * h)) * ur * ur
            + h * r * uph * uph,
            -(2.0 / r) * ur * uph,
        ]

    def leave_domain(tau, y):
        return min(y[1] - r_lo, r_hi - y[1])

    leave_domain.terminal = True  # type: ignore[attr-defined]
    leave_domain.direction = -1  # type: ignore[attr-defined]

    t_eval = np.linspace(tau_span[0], tau_span[1], 4001)
    sol = solve_ivp(
        rhs, tau_span, [0.0, r0, 0.0, ut0, ur0, L / r0**2],
        method="DOP853", rtol=1e-12, atol=1e-12, t_eval=t_eval,
        events=leave_domain,
    )
    if not sol.success or sol.y.shape[1] < 50:
        raise RuntimeError(sol.message or "too few valid steps")
    t, r, phi, ut, ur, uph = sol.y
    tau = sol.t

    E_arr = m.f(r) * ut
    L_arr = r**2 * uph
    norm = -m.f(r) * ut**2 + ur**2 / m.h(r) + r**2 * uph**2

    # independent acceleration cross-check (finite differences vs Gamma form)
    dudtau = np.gradient(np.stack([ut, ur, uph]), tau, axis=1)
    accel_res = 0.0
    for i in range(1, len(tau) - 1):
        G = m.Gamma_tensor(r[i])
        a = dudtau[:, i] + np.array([
            2.0 * G[0, 0, 1] * ut[i] * ur[i],
            G[1, 0, 0] * ut[i] ** 2 + G[1, 1, 1] * ur[i] ** 2
            + G[1, 3, 3] * uph[i] ** 2,
            2.0 * G[3, 1, 3] * ur[i] * uph[i],
        ])
        scale = max(1.0, float(np.max(np.abs(a))))
        accel_res = max(accel_res, float(np.max(np.abs(a))) / scale)

    return TransportCheck(
        E0=float(E), L0=float(L),
        max_dE=float(np.max(np.abs(E_arr - E))),
        max_dL=float(np.max(np.abs(L_arr - L))),
        max_norm_res=float(np.max(np.abs(norm + 1.0))),
        max_accel_res=float(accel_res),
        tau=tau, r=r, t=t, phi=phi, ut=ut, ur=ur, uph=uph,
    )


def circular_timelike_orbit(m: SSZMetric, r_orb: float) -> dict:
    """Exact circular timelike orbit (u^r = 0): Kepler Omega^2 = f'/(2r)."""
    f = float(m.f(r_orb))
    fp = float(m.fp(r_orb))
    Omega = float(np.sqrt(fp / (2.0 * r_orb)))
    ut = float(1.0 / np.sqrt(f - r_orb**2 * Omega**2))
    return {
        "r": r_orb,
        "Omega": Omega,
        "ut": ut,
        "L": float(r_orb**2 * Omega * ut),
        "E": float(f * ut),
        "dtau_dt": float(np.sqrt(f - r_orb**2 * Omega**2)),
    }


# ---------------------------------------------------------------------------
# covariant-derivative helper (B_{mu nu} = nabla_nu V_mu for r-dependent flows)
# ---------------------------------------------------------------------------


def _B_matrix(m: SSZMetric, r: float, v_contra: np.ndarray,
              dv_contra_dr: np.ndarray) -> np.ndarray:
    """B[mu, nu] = nabla_nu V_mu for a stationary flow V^mu(r).

    v_contra: contravariant components at r; dv_contra_dr: the exact
    analytic radial derivative of the contravariant field (no finite
    differences, no step-size tuning)."""
    g = np.diag([
        -float(m.f(r)),
        1.0 / float(m.h(r)),
        r * r,
        (r * float(np.sin(np.pi / 2))) ** 2,
    ])
    dg_dr = np.diag([
        -float(m.fp(r)),
        -float(m.hp(r)) / float(m.h(r)) ** 2,  # d(1/h)/dr = -h'/h^2
        2.0 * r,
        2.0 * r,
    ])
    v_cov = g @ v_contra
    dv_cov_dr = dg_dr @ v_contra + g @ dv_contra_dr
    G = m.Gamma_tensor(r)
    B = np.zeros((4, 4))
    for mu in range(4):
        for nu in range(4):
            B[mu, nu] = (dv_cov_dr[mu] if nu == 1 else 0.0) \
                - G[:, mu, nu] @ v_cov
    return B


def _proj_sigma2_omega2(B: np.ndarray, g: np.ndarray,
                        ginv: np.ndarray, u_contra: np.ndarray,
                        dim: int) -> tuple[float, float, float]:
    """Projected sigma^2 and omega^2 (and theta) of a unit flow.

    q_ab = g_ab + sgn * u_a u_b with sgn = +1 for timelike (-1 norm),
    q_ab = g_ab for null (using the k/s symmetric projector instead —
    handled by the caller passing qcov/qcon directly).
    """
    u_cov = g @ u_contra
    u2 = float(u_contra @ u_cov)
    sgn = 1.0 if u2 < 0 else -1.0
    qcov = g + sgn * np.outer(u_cov, u_cov)
    qcon = ginv + sgn * np.outer(u_contra, u_contra)
    theta = float(np.einsum("ab,ab->", qcon, B))
    sym = 0.5 * (B + B.T)
    asym = 0.5 * (B - B.T)
    sig = sym - (theta / dim) * qcov
    sigma2 = float(np.einsum("ab,cd,ac,bd->", sig, sig, qcon, qcon))
    omega2 = float(np.einsum("ab,cd,ac,bd->", asym, asym, qcon, qcon))
    return sigma2, omega2, theta


# ---------------------------------------------------------------------------
# 2. kinematical scalars + Raychaudhuri (radial timelike congruence)
# ---------------------------------------------------------------------------


def radial_congruence_scalars(m: SSZMetric, r_grid, E: float = 1.0) -> dict:
    """Expansion/shear/vorticity of the stationary radial free-fall
    congruence u^mu(r) with Killing energy E (geodesic, irrotational).

    The congruence field is known in CLOSED FORM:
        u^t = E / f,   u^r = -sqrt(h (E^2/f - 1)),
    so theta is computed via the exact divergence identity
        theta = (1/sqrt(-g)) d_r (sqrt(-g) u^r)
    AND via the Christoffel route with analytic derivatives; the
    Raychaudhuri identity
        dtheta/dtau = -theta^2/3 - sigma^2 + omega^2 - R_mn u^mu u^nu
    is verified between independently computed left/right sides."""
    r_grid = np.asarray(r_grid, dtype=float)
    r_lo, r_hi = radial_domain(m)
    if not (r_lo <= r_grid.min() and r_grid.max() <= r_hi):
        raise ValueError(
            f"r_grid [{r_grid.min()}, {r_grid.max()}] outside member domain "
            f"[{r_lo}, {r_hi}]")
    f = m.f(r_grid)
    h = m.h(r_grid)
    ut = E / f
    ur = -np.sqrt(np.maximum(h * (E**2 / f - 1.0), 0.0))

    sq = np.sqrt(f / h)  # sqrt(-g) ~ r^2 sqrt(f/h)
    integ = CubicSpline(r_grid, r_grid**2 * sq * ur)
    theta_A = integ.derivative()(r_grid) / (r_grid**2 * sq)
    Gamma_sum = m.fp(r_grid) / (2 * f) - m.hp(r_grid) / (2 * h) + 2.0 / r_grid
    # analytic du^r/dr:
    #   u^r = -sqrt(h A),  A = E^2/f - 1
    #   du^r/dr = -(h' A + h A') / (2 sqrt(h A)),  A' = -E^2 f'/f^2
    A = E**2 / f - 1.0
    ur_prime = -(m.hp(r_grid) * A + h * (-E**2 * m.fp(r_grid) / f**2)) / (
        2.0 * np.sqrt(h * A))
    theta_B = ur_prime + Gamma_sum * ur

    n = len(r_grid)
    sigma2 = np.empty(n)
    omega2 = np.empty(n)
    R_uu = np.empty(n)
    for i, r in enumerate(r_grid):
        f_i, h_i = float(f[i]), float(h[i])
        fpi = float(m.fp(r))
        hpi = float(m.hp(r))
        # analytic contravariant field derivatives at r
        A_i = E**2 / f_i - 1.0
        urp = -(hpi * A_i + h_i * (-E**2 * fpi / f_i**2)) / (
            2.0 * np.sqrt(h_i * A_i))
        v = np.array([ut[i], ur[i], 0.0, 0.0])
        dv = np.array([-E * fpi / f_i**2, urp, 0.0, 0.0])
        B = _B_matrix(m, float(r), v, dv)
        g = np.diag([-f_i, 1.0 / h_i, r * r, r * r])
        ginv = np.diag([-1.0 / f_i, h_i, 1.0 / (r * r), 1.0 / (r * r)])
        s2, w2, _ = _proj_sigma2_omega2(B, g, ginv, v, dim=3)
        sigma2[i], omega2[i] = s2, w2
        Ric = m.Ricci(float(r))
        # R_mn u^mu u^nu = R_tt (u^t)^2 + R_rr (u^r)^2 (contravariant square)
        R_uu[i] = Ric[0] * ut[i] ** 2 + Ric[1] * ur[i] ** 2

    dtheta_dr = CubicSpline(r_grid, theta_A).derivative()(r_grid)
    lhs = ur * dtheta_dr
    rhs = -theta_A**2 / 3.0 - sigma2 + omega2 - R_uu
    scale = np.maximum(1.0, np.abs(rhs))
    return {
        "r": r_grid,
        "theta": theta_A,
        "theta_routeB": theta_B,
        "theta_routes_absdiff": float(np.max(np.abs(theta_A - theta_B))),
        "sigma2": sigma2,
        "omega2": omega2,
        "R_uu": R_uu,
        "raychaudhuri_lhs": lhs,
        "raychaudhuri_rhs": rhs,
        "max_scaled_raychaudhuri_res": float(np.max(np.abs(lhs - rhs) / scale)),
    }


# ---------------------------------------------------------------------------
# 3. null sector: photon spheres, balance, amplitude, phase
# ---------------------------------------------------------------------------


def find_light_rings(m: SSZMetric) -> list[dict]:
    """Roots of W_u on the member domain, classified by W_uu (stability)."""
    us = np.linspace(m.u_min + 1e-6, m.u_max - 1e-6, 4000)
    Wu = m.W_u(us)
    rings = []
    for i in np.flatnonzero(np.sign(Wu[:-1]) * np.sign(Wu[1:]) < 0):
        u_r = brentq(lambda x: float(m.W_u(x)), us[i], us[i + 1],
                     xtol=1e-14, rtol=8.9e-16)
        W = float(m.W_of_u(u_r))
        Wuu = float(m.W_uu(u_r))
        r_r = 1.0 / u_r
        f_r = float(m.f(r_r))
        h_r = float(m.h(r_r))
        rings.append({
            "u": u_r,
            "r": r_r,
            "W": W,
            "W_uu": Wuu,
            "stable": bool(Wuu > 0),
            "b_crit": float(1.0 / np.sqrt(W)),
            "Omega_ph": float(np.sqrt(W)),
            # linear exponent: kappa^2 = |W_uu| h / (2 f) at the ring
            "exponent": float(np.sqrt(abs(Wuu) * h_r / (2.0 * f_r))),
        })
    return rings


def libration_period_ode(m: SSZMetric, u_ring: float, eps: float = 1e-5,
                         n_periods: int = 4) -> float:
    """Dynamical libration period (in phi) of a near-ring null orbit with
    1/b^2 = W(u_ring) + eps, measured by integrating the FULL null geodesic
    equations and detecting successive downward crossings of u = u_ring.

    Independent dynamical cross-check of the linear-theory period
    2 pi / kappa with kappa^2 = W_uu h / (2 f) at the ring."""
    b2_inv = float(m.W_of_u(u_ring)) + eps
    b = 1.0 / np.sqrt(b2_inv)
    r0 = 1.0 / u_ring
    f0 = float(m.f(r0))
    h0 = float(m.h(r0))
    # (dr/dlam)^2 = h (1/f - b^2 u^2); u_ring is near the libration turning
    # region, so the radial velocity is small and POSITIVE (r increasing,
    # u decreasing first)
    kr0 = float(np.sqrt(max(h0 * (1.0 / f0 - u_ring**2 / b2_inv), 0.0)))
    kt0 = 1.0 / f0

    def rhs(lam, y):
        _, r, _, kt, kr, kph = y
        f = float(m.f(r))
        fp = float(m.fp(r))
        h = float(m.h(r))
        hp = float(m.hp(r))
        return [
            kt, kr, kph,
            -(fp / f) * kt * kr,
            -(0.5 * h * fp) * kt * kt + (hp / (2.0 * h)) * kr * kr
            + h * r * kph * kph,
            -(2.0 / r) * kr * kph,
        ]

    def cross_ring(lam, y):
        # u crosses u_ring downward <=> r crosses r0 = 1/u_ring upward
        return y[1] - r0

    cross_ring.direction = +1.0  # type: ignore[attr-defined]

    kappa = float(np.sqrt(abs(float(m.W_uu(u_ring))) * float(m.h(1.0 / u_ring))
                          / (2.0 * float(m.f(1.0 / u_ring)))))
    phi_per_period = 2.0 * np.pi / kappa
    dphi_dlam = b * u_ring**2
    lam_span = (0.0, 1.5 * n_periods * phi_per_period / dphi_dlam)

    sol = solve_ivp(rhs, lam_span, [0.0, r0, 0.0, kt0, kr0, b / r0**2],
                    method="DOP853", rtol=1e-12, atol=1e-12,
                    dense_output=True, events=cross_ring,
                    max_step=lam_span[1] / 20000.0)
    if not sol.success or sol.t_events[0].shape[0] < 2:
        raise RuntimeError("libration period: not enough ring crossings")
    phi_events = []
    for lam_ev in sol.t_events[0]:
        phi_events.append(float(sol.sol(lam_ev)[2]))
    diffs = np.diff(phi_events)
    return float(np.mean(diffs))


def radial_null_congruence(m: SSZMetric, r_grid) -> dict:
    """Optical scalars and Raychaudhuri identity for the OUTGOING radial
    null congruence of AFFINELY parametrised geodesics (b = 0, E = 1):

        k^mu = (1/f, sqrt(h/f), 0, 0),   k^nu nabla_nu k^mu = 0 exactly,
        theta_hat = 2 sqrt(h/f) / r,     sigma_hat = omega_hat = 0.

    Raychaudhuri  d(theta_hat)/dlambda = -theta_hat^2/2 - R_mn k^mu k^nu
    is verified between independently computed sides (LHS: theta spline
    derivative along the flow; RHS: optical scalars + symbolic Ricci)."""
    r_grid = np.asarray(r_grid, dtype=float)
    r_lo, r_hi = radial_domain(m)
    if not (r_lo <= r_grid.min() and r_grid.max() <= r_hi):
        raise ValueError("r_grid outside member domain")

    f = m.f(r_grid)
    h = m.h(r_grid)
    kt = 1.0 / f
    kr = np.sqrt(h / f)

    # route A: divergence identity with density D = r^2 sqrt(f/h)
    sq = np.sqrt(f / h)
    flux = r_grid**2 * sq * kr           # = r^2 (exact for this congruence)
    flux_sp = CubicSpline(r_grid, flux)
    theta_A = flux_sp.derivative()(r_grid) / (r_grid**2 * sq)
    # route B: closed form theta = 2 sqrt(h/f) / r
    theta_B = 2.0 * kr / r_grid

    Ric_arr = np.array([m.Ricci(float(r)) for r in r_grid])
    R_kk = Ric_arr[:, 0] * kt**2 + Ric_arr[:, 1] * kr**2

    dtheta_dr = CubicSpline(r_grid, theta_A).derivative()(r_grid)
    lhs = kr * dtheta_dr                      # d/dlambda along the flow
    rhs = -0.5 * theta_A**2 - R_kk            # sigma = omega = 0 (spherical)
    scale = np.maximum(1.0, np.abs(rhs))
    return {
        "r": r_grid,
        "theta_hat": theta_A,
        "theta_hat_routeB": theta_B,
        "theta_routes_absdiff": float(np.max(np.abs(theta_A - theta_B))),
        "sigma_hat2": np.zeros_like(r_grid),
        "omega_hat2": np.zeros_like(r_grid),
        "R_kk": R_kk,
        "raychaudhuri_lhs": lhs,
        "raychaudhuri_rhs": rhs,
        "max_scaled_raychaudhuri_res": float(np.max(np.abs(lhs - rhs) / scale)),
    }


def ring_trapping(m: SSZMetric, rings: list[dict],
                  epsilons_unstable: tuple = (1e-2, 1e-3, 1e-4, 1e-5),
                  epsilons_stable: tuple = (1e-5, 1e-6, 1e-7, 1e-8)) -> dict:
    """Source-free trapping diagnostics from the exact null orbit equation

        (du/dphi)^2 = (h/f) (1/b^2 - W(u)),

    integrated by turning-point quadrature (no ODE stiffness, no fitting):

    - UNSTABLE (outer, W_uu < 0) ring: photons with b slightly above b_crit
      wind logarithmically; the total deflection angle diverges like
      ln(1/delta) as the turning point approaches the ring.  Measured by
      the increment of Delta phi under halving delta -> ln 2.
    - STABLE (inner, W_uu > 0) ring: photons LIBRATE around the ring with
      finite angular period phi_p -> 2 pi / kappa,
      kappa^2 = W_uu h / (2 f) at the ring.

    This is the matured, force-free form of the early paper's boundary
    intuition: the ring structure is a property of the geometry alone."""
    from scipy.integrate import quad

    out = {"outer_unstable": {}, "inner_stable": {}}
    for ring in rings:
        key = "inner_stable" if ring["stable"] else "outer_unstable"
        u_r, W_r = ring["u"], ring["W"]
        delta_series = []
        epsilons = epsilons_stable if ring["stable"] else epsilons_unstable
        extras = {}
        if ring["stable"]:
            # primary dynamical measurement: full-ODE libration period at
            # the largest domain-safe energy
            extras["libration_period_ode"] = libration_period_ode(
                m, u_r, eps=epsilons_stable[0])
        for eps in epsilons:
            if ring["stable"]:
                # exact spline-consistent turning points around the minimum:
                # grow the bracket geometrically until W crosses b2_inv
                b2_inv = W_r + eps

                def _turn(sign):
                    step = max(1e-9, 1e-3 * eps / abs(W_r))
                    lo = hi = u_r
                    while float(m.W_of_u(hi + sign * step)) - b2_inv < 0.0:
                        hi = hi + sign * step
                        step *= 1.6
                        if not (m.u_min < hi < m.u_max):
                            raise ValueError("libration turning point left domain")
                    return brentq(lambda x: float(m.W_of_u(x)) - b2_inv,
                                  hi, hi + sign * step,
                                  xtol=1e-15, rtol=8.9e-16)

                u_lo = _turn(-1.0)
                u_hi = _turn(+1.0)
                amp = 0.5 * (u_hi - u_lo)
                u_mid = 0.5 * (u_lo + u_hi)

                def half_arc(psi, b2_inv=b2_inv, u_mid=u_mid, amp=amp):
                    # u = u_mid + amp cos(psi): smooth at both turning points
                    u = u_mid + amp * np.cos(psi)
                    du = amp * np.sin(psi)
                    val = (float(m.h_u(u)) / float(m.f_u(u))) * (
                        b2_inv - float(m.W_of_u(u)))
                    # relative floor kills the 0/0 at psi = 0, pi exactly
                    return du / np.sqrt(max(val, 1e-14 * (1.0 + b2_inv)))

                import warnings
                with warnings.catch_warnings():
                    # roundoff warnings near the integrable turning-point
                    # endpoints; ODE period (below) is the primary check
                    warnings.simplefilter("ignore")
                    val, _ = quad(half_arc, 0.0, np.pi, limit=500,
                                  epsabs=1e-11, epsrel=1e-11)
                delta_series.append(2.0 * val)
            else:
                # winding: turning point u_tp = u_r - eps (spline-exact),
                # integrated outward to a FIXED reference radius
                u_tp = u_r - eps
                b2_inv = float(m.W_of_u(u_tp))
                u_ref = m.u_min + 1e-6

                def arc(s, b2_inv=b2_inv, u_tp=u_tp):
                    # u = u_tp - s^2 smooths the inverse-square-root endpoint
                    u = u_tp - s * s
                    val = (float(m.h_u(u)) / float(m.f_u(u))) * (
                        b2_inv - float(m.W_of_u(u)))
                    return 2.0 * s / np.sqrt(max(val, 1e-14 * (1.0 + b2_inv)))

                from scipy.integrate import quad
                import warnings
                with warnings.catch_warnings():
                    # quad reports roundoff near the (integrable) turning
                    # point; the substituted integrand is smooth and the
                    # achieved accuracy is far below test tolerance
                    warnings.simplefilter("ignore")
                    val, _ = quad(arc, 0.0, np.sqrt(u_tp - u_ref), limit=500,
                                  epsabs=1e-10, epsrel=1e-10)
                delta_series.append(2.0 * val)
        increments = [delta_series[i + 1] - delta_series[i]
                      for i in range(len(delta_series) - 1)]
        out[key] = {
            "epsilons": [float(e) for e in epsilons],
            "delta_phi_series": [float(v) for v in delta_series],
            "increments": [float(v) for v in increments],
            # decade epsilons: log-winding adds 2 ln 10 per decade
            "expected_increment_per_decade": (
                float(2.0 * np.log(10.0)) if not ring["stable"] else None),
            "expected_libration_period_2pi_over_kappa": (
                float(2.0 * np.pi / ring["exponent"]) if ring["stable"] else None),
            "kappa_linear": ring["exponent"],
        }
        out[key].update(extras)
    return out


def null_geodesic_transport(
    m: SSZMetric, b: float, r0: float, inward: bool = True,
    lam_span: tuple = (0.0, 30.0), E: float = 1.0,
) -> TransportCheck:
    """Null geodesic with impact parameter b: k^nu nabla_nu k^mu = 0."""
    r_lo, r_hi = radial_domain(m)
    if not (r_lo < r0 < r_hi):
        raise ValueError(f"r0={r0} outside member domain {r_lo, r_hi}")
    f0 = float(m.f(r0))
    h0 = float(m.h(r0))
    val = h0 * (E * E / f0 - (b * E) ** 2 / r0**2)
    if val < 0:
        raise ValueError("null ray with this b forbidden at r0")
    kr0 = -np.sqrt(val) if inward else np.sqrt(val)
    kt0 = E / f0

    def rhs(lam, y):
        _, r, _, kt, kr, kph = y
        f = float(m.f(r))
        fp = float(m.fp(r))
        h = float(m.h(r))
        hp = float(m.hp(r))
        return [
            kt, kr, kph,
            -(fp / f) * kt * kr,
            -(0.5 * h * fp) * kt * kt + (hp / (2.0 * h)) * kr * kr
            + h * r * kph * kph,
            -(2.0 / r) * kr * kph,
        ]

    def leave_domain(lam, y):
        return min(y[1] - r_lo, r_hi - y[1])

    leave_domain.terminal = True  # type: ignore[attr-defined]
    leave_domain.direction = -1  # type: ignore[attr-defined]

    t_eval = np.linspace(lam_span[0], lam_span[1], 8001)
    sol = solve_ivp(
        rhs, lam_span, [0.0, r0, 0.0, kt0, kr0, b * E / r0**2],
        method="DOP853", rtol=1e-12, atol=1e-12, t_eval=t_eval,
        events=leave_domain,
    )
    if not sol.success or sol.y.shape[1] < 50:
        raise RuntimeError(sol.message or "too few valid steps")
    t, r, phi, kt, kr, kph = sol.y
    E_arr = m.f(r) * kt
    L_arr = r**2 * kph
    nrm = -m.f(r) * kt**2 + kr**2 / m.h(r) + r**2 * kph**2
    return TransportCheck(
        E0=float(E), L0=float(b * E),
        max_dE=float(np.max(np.abs(E_arr - E))),
        max_dL=float(np.max(np.abs(L_arr - b * E))),
        max_norm_res=float(np.max(np.abs(nrm))),
        max_accel_res=0.0,
        tau=sol.t, r=r, t=t, phi=phi, ut=kt, ur=kr, uph=kph,
    )


def amplitude_transport_radial(m: SSZMetric, r_emit: float, r_obs: float,
                               n: int = 4000) -> dict:
    """d(a^2)/dlambda = -(nabla_mu k^mu) a^2 for a radial null beam.

    For radial beams nabla_mu k^mu = (1/(r^2 sqrt(f/h)))
    d_r (r^2 sqrt(f/h) k^r) and spherical symmetry forces the exact
    geometric invariant a^2 r^2 = const (beam cross-section ~ r^2).
    The transport ODE is integrated independently in r and checked
    against that invariant."""
    r_grid = np.linspace(r_emit, r_obs, n)
    f = m.f(r_grid)
    h = m.h(r_grid)
    kr = np.sqrt(h / f)  # k^r for E = 1
    sq = np.sqrt(f / h)
    dens = r_grid**2 * sq          # D = sqrt(-g) density (equatorial)
    flux = dens * kr               # D * k^r = r^2 exactly for radial null
    flux_sp = CubicSpline(r_grid, flux)
    divk = flux_sp.derivative()(r_grid) / dens   # nabla_mu k^mu (affine)
    coef = CubicSpline(r_grid, -divk / kr)
    sol = solve_ivp(
        lambda rr, y: [float(coef(rr)) * y[0]], (r_emit, r_obs), [1.0],
        rtol=1e-12, atol=1e-14, dense_output=True,
    )
    a2_obs = float(sol.y[0, -1])
    return {
        "a2_ratio": a2_obs,
        "invariant_error": float(abs(a2_obs * r_obs**2 - r_emit**2) / r_emit**2),
    }


def eikonal_phase_and_redshift(m: SSZMetric, r_a: float, r_b: float,
                               n: int = 20000) -> dict:
    """Shapiro delay, reduced eikonal phase integral and static-observer
    redshift between two radii on a radial null ray (exact quadrature).

    The eikonal phase building block per unit photon energy is
        S_r = integral dr / sqrt(f h),
    the exact quantity the JIF phase chain (Xi -> D -> g -> Phi) consumes."""
    r_grid = np.linspace(min(r_a, r_b), max(r_a, r_b), n)
    integ = float(np.trapezoid(1.0 / np.sqrt(m.f(r_grid) * m.h(r_grid)), r_grid))
    return {
        "shapiro_dt_per_E": integ,
        "reduced_phase_per_E": integ,
        "redshift_a_to_b": float(np.sqrt(
            m.f(min(r_a, r_b)) / m.f(max(r_a, r_b)))),
    }


# ---------------------------------------------------------------------------
# 4. consolidation
# ---------------------------------------------------------------------------


def forced_transport_control(m: SSZMetric, r0: float, F_t: float = 1e-2,
                             tau_span: tuple = (0.0, 1.5)) -> dict:
    """NEGATIVE CONTROL (diagnostic only): repeat the radial timelike
    integration with a synthetic external force F^t = const.

    Purpose: prove that the conservation-based source-free diagnostics
    (Killing energy) genuinely DISTINGUISH source-free transport from
    forced transport.  dE/dtau = f * F^t for this control, so E drifts
    linearly -- in sharp contrast to the vanishing geodesic residuals."""
    f0 = float(m.f(r0))
    h0 = float(m.h(r0))
    ur0 = -np.sqrt(max(h0 * (1.0 / f0 - 1.0), 0.0))
    ut0 = 1.0 / f0

    def rhs(tau, y):
        _, r, _, ut, ur, _ = y
        f = float(m.f(r))
        fp = float(m.fp(r))
        h = float(m.h(r))
        hp = float(m.hp(r))
        return [
            ut,
            ur,
            0.0,
            -(fp / f) * ut * ur + F_t,
            -(0.5 * h * fp) * ut * ut + (hp / (2.0 * h)) * ur * ur,
            0.0,
        ]

    sol = solve_ivp(rhs, tau_span, [0.0, r0, 0.0, ut0, ur0, 0.0],
                    method="DOP853", rtol=1e-12, atol=1e-12,
                    max_step=(tau_span[1] - tau_span[0]) / 200.0)
    if not sol.success:
        raise RuntimeError(sol.message)
    t, r, _, ut, ur, _ = sol.y
    E_arr = m.f(r) * ut
    return {
        "F_t": F_t,
        "max_dE": float(np.max(np.abs(E_arr - 1.0))),
        "predicted_drift_scale": float(F_t * (tau_span[1] - tau_span[0])),
        "tau_span": list(tau_span),
    }


def libration_trajectory(m: SSZMetric, u_ring: float, eps: float = 1e-5,
                         n_periods: int = 2, n_out: int = 6000) -> dict:
    """Full trajectory of a near-ring librating photon: (lam, t, r, phi)
    arrays for visualisation of the stable-ring oscillation."""
    b2_inv = float(m.W_of_u(u_ring)) + eps
    b = 1.0 / np.sqrt(b2_inv)
    r0 = 1.0 / u_ring
    f0 = float(m.f(r0))
    h0 = float(m.h(r0))
    kr0 = float(np.sqrt(max(h0 * (1.0 / f0 - u_ring**2 / b2_inv), 0.0)))
    kt0 = 1.0 / f0

    def rhs(lam, y):
        _, r, _, kt, kr, kph = y
        f = float(m.f(r))
        fp = float(m.fp(r))
        h = float(m.h(r))
        hp = float(m.hp(r))
        return [
            kt, kr, kph,
            -(fp / f) * kt * kr,
            -(0.5 * h * fp) * kt * kt + (hp / (2.0 * h)) * kr * kr
            + h * r * kph * kph,
            -(2.0 / r) * kr * kph,
        ]

    def cross_ring(lam, y):
        return y[1] - r0

    cross_ring.direction = +1.0  # type: ignore[attr-defined]
    kappa = float(np.sqrt(abs(float(m.W_uu(u_ring)))
                          * float(m.h(1.0 / u_ring))
                          / (2.0 * float(m.f(1.0 / u_ring)))))
    phi_per_period = 2.0 * np.pi / kappa
    dphi_dlam = b * u_ring**2
    lam_span = (0.0, 1.05 * n_periods * phi_per_period / dphi_dlam)
    sol = solve_ivp(rhs, lam_span, [0.0, r0, 0.0, kt0, kr0, b / r0**2],
                    method="DOP853", rtol=1e-12, atol=1e-12,
                    dense_output=True, max_step=lam_span[1] / n_out)
    lam = np.linspace(0.0, lam_span[1], n_out)
    y = sol.sol(lam)
    return {"lam": lam, "t": y[0], "r": y[1], "phi": y[2]}


def winding_trajectory(m: SSZMetric, b: float, r0: float,
                       lam_span: tuple = (0.0, 14.0),
                       n_out: int = 8000, inward: bool = True) -> dict:
    """Full null trajectory (r, phi) for visualisation of the logarithmic
    winding near the unstable ring."""
    f0 = float(m.f(r0))
    h0 = float(m.h(r0))
    val = h0 * (1.0 / f0 - b * b / r0**2)
    if val < 0:
        raise ValueError("impact parameter forbidden at r0")
    kr0 = -float(np.sqrt(val)) if inward else float(np.sqrt(val))
    kt0 = 1.0 / f0

    def rhs(lam, y):
        _, r, _, kt, kr, kph = y
        f = float(m.f(r))
        fp = float(m.fp(r))
        h = float(m.h(r))
        hp = float(m.hp(r))
        return [
            kt, kr, kph,
            -(fp / f) * kt * kr,
            -(0.5 * h * fp) * kt * kt + (hp / (2.0 * h)) * kr * kr
            + h * r * kph * kph,
            -(2.0 / r) * kr * kph,
        ]

    r_lo, r_hi = radial_domain(m)

    def leave(lam, y):
        return min(y[1] - r_lo, r_hi - y[1])

    leave.terminal = True  # type: ignore[attr-defined]
    sol = solve_ivp(rhs, lam_span, [0.0, r0, 0.0, kt0, kr0, b / r0**2],
                    method="DOP853", rtol=1e-12, atol=1e-12,
                    dense_output=True, events=leave,
                    max_step=(lam_span[1] - lam_span[0]) / n_out)
    lam_end = float(sol.t[-1])
    lam = np.linspace(0.0, lam_end, n_out)
    y = sol.sol(lam)
    return {"lam": lam, "t": y[0], "r": y[1], "phi": y[2]}


def consolidation_table(m: SSZMetric, rings: list[dict]) -> list[dict]:
    """One geometry -> orbit + time + phase + rotation + trapping."""
    rows = []
    for ring in rings:
        r = ring["r"]
        OmK = float(np.sqrt(m.fp(r) / (2.0 * r)))
        rows.append({
            "ring": "inner_stable" if ring["stable"] else "outer_unstable",
            "u": ring["u"],
            "r": r,
            "dtau_dt": float(np.sqrt(m.f(r))),
            "Omega_kepler": OmK,
            "Omega_photon": ring["Omega_ph"],
            "kepler_vs_photon_rel_diff": float(abs(OmK - ring["Omega_ph"])
                                               / ring["Omega_ph"]),
            "b_crit": ring["b_crit"],
            "W_uu": ring["W_uu"],
            "stability_exponent": ring["exponent"],
            "frame_dragging": 0.0,  # g_tphi = 0 exactly (static central slice)
        })
    return rows
