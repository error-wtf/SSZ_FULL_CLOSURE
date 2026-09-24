"""TRUE FULL CLOSURE — canonical first-principles forward chain.

ONE geometry (the frozen electric production member), ONE source-free
transport operator, MANY physical manifestations.  Every function in this
module is either

  (a) a canonical definition imported by all downstream computations, or
  (b) an independent verification route for an identity implemented
      elsewhere (symbolic <-> numeric, analytic <-> integrated).

Nothing here is fitted; all reference values are either exact closed forms
(Schwarzschild anchors) or programmatic outputs of the frozen member.

Conventions (single source of truth):
  signature  (-, +, +, +)
  coords     (t, r, theta, phi)
  units      c = G = 1, r_s = 1, u = r_s / r
  metric     ds^2 = -f dt^2 + dr^2/h + r^2 dOmega^2
"""
from __future__ import annotations

import numpy as np
import sympy as sp
from dataclasses import dataclass
from functools import lru_cache

from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

from ..postclosure.transport import SSZMetric

COORDS = ("t", "r", "th", "ph")
SIGNATURE = "(-,+,+,+)"


@dataclass(frozen=True)
class Conventions:
    signature: str = SIGNATURE
    coordinates: tuple = COORDS
    units: str = "c = G = 1, r_s = 1, u = r_s/r"
    metric_ansatz: str = "ds^2 = -f dt^2 + dr^2/h + r^2 dOmega^2"


CONVENTIONS = Conventions()


# ---------------------------------------------------------------------------
# metric tensors, determinant (foundations)
# ---------------------------------------------------------------------------


def metric_tensors(m: SSZMetric, r: float, theta: float = np.pi / 2.0):
    """Covariant metric, inverse and determinant at (r, theta)."""
    f = float(m.f(r))
    h = float(m.h(r))
    s = float(np.sin(theta))
    g = np.diag([-f, 1.0 / h, r * r, (r * s) ** 2])
    ginv = np.diag([-1.0 / f, h, 1.0 / (r * r), 1.0 / ((r * s) ** 2)])
    det = float(-f * (1.0 / h) * (r * r) * (r * s) ** 2)
    return g, ginv, det


# ---------------------------------------------------------------------------
# numeric Riemann from finite-differenced Christoffels (independent route)
# ---------------------------------------------------------------------------


def christoffels_fd(m: SSZMetric, r: float, theta: float = np.pi / 2.0,
                    eps: float = 1e-5) -> tuple[np.ndarray, dict]:
    """Gamma[lam, mu, nu] and its coordinate derivatives dG[axis] via
    central finite differences (Gamma depends on r and theta only)."""
    G0 = m.Gamma_tensor(r, theta)
    dG_dr = (m.Gamma_tensor(r + eps, theta) - m.Gamma_tensor(r - eps, theta)) \
        / (2.0 * eps)
    dG_dth = (m.Gamma_tensor(r, theta + eps)
              - m.Gamma_tensor(r, theta - eps)) / (2.0 * eps)
    return G0, {1: dG_dr, 2: dG_dth, 0: np.zeros_like(G0), 3: np.zeros_like(G0)}


def riemann_numeric(m: SSZMetric, r: float, theta: float = np.pi / 2.0,
                    eps: float = 1e-5) -> np.ndarray:
    """R^mu_{nu rho sigma} = d_rho Gamma^mu_{nu sigma}
    - d_sigma Gamma^mu_{nu rho}
    + Gamma^mu_{rho lam} Gamma^lam_{nu sigma}
    - Gamma^mu_{sig lam} Gamma^lam_{nu rho}   (numeric route)."""
    G, dG = christoffels_fd(m, r, theta, eps)
    R = np.zeros((4, 4, 4, 4))
    for mu in range(4):
        for nu in range(4):
            for rho in range(4):
                for sig in range(4):
                    val = dG[rho][mu, nu, sig] - dG[sig][mu, nu, rho]
                    for lam in range(4):
                        val += (G[mu, rho, lam] * G[lam, nu, sig]
                                - G[mu, sig, lam] * G[lam, nu, rho])
                    R[mu, nu, rho, sig] = val
    return R


def ricci_from_riemann(R: np.ndarray) -> np.ndarray:
    """R_{mu nu} = R^rho_{mu rho nu}."""
    return np.einsum("rmrn->mn", R)


def ricci_scalar_numeric(m: SSZMetric, r: float,
                         theta: float = np.pi / 2.0) -> float:
    R = riemann_numeric(m, r, theta)
    Ric = ricci_from_riemann(R)
    _, ginv, _ = metric_tensors(m, r, theta)
    return float(np.einsum("mn,mn->", ginv, Ric))


def einstein_numeric(m: SSZMetric, r: float,
                     theta: float = np.pi / 2.0) -> np.ndarray:
    """G^{mu nu} = R^{mu nu} - 1/2 g^{mu nu} R (contravariant)."""
    R = riemann_numeric(m, r, theta)
    Ric = ricci_from_riemann(R)
    _, ginv, _ = metric_tensors(m, r, theta)
    Scal = float(np.einsum("mn,mn->", ginv, Ric))
    return np.einsum("ma,nb,ab->mn", ginv, ginv, Ric) - 0.5 * ginv * Scal


def curvature_invariants_numeric(m: SSZMetric, r: float,
                                 theta: float = np.pi / 2.0) -> dict:
    """Kretschmann K, R_mn R^mn and R^2 from the numeric Riemann tensor.
    Only the upper Riemann index is lowered with g (the other three slots
    are already covariant); contractions raise with the inverse metric."""
    R = riemann_numeric(m, r, theta)
    g, ginv, _ = metric_tensors(m, r, theta)
    R_down = np.einsum("am,mbcd->abcd", g, R)
    R_up = np.einsum("ea,fb,gc,hd,abcd->efgh", ginv, ginv, ginv, ginv, R_down)
    K = float(np.einsum("abcd,abcd->", R_down, R_up))
    Ric = ricci_from_riemann(R)
    Ric_up = np.einsum("am,bn,mn->ab", ginv, ginv, Ric)
    Ric2 = float(np.einsum("ab,ab->", Ric, Ric_up))
    Scal = float(np.einsum("mn,mn->", ginv, Ric))
    return {"kretschmann": K, "ricci2": Ric2, "r2": Scal * Scal,
            "ricci_scalar": Scal}


def check_riemann_symmetries(m: SSZMetric, r: float,
                             theta: float = np.pi / 2.0,
                             eps: float = 1e-5) -> dict:
    """Riemann algebraic symmetries and the first Bianchi identity.
    R^mu_{nu rho sigma} is a (1,3) tensor: only the UPPER index is
    lowered (g_{a mu}); the other three slots are already covariant."""
    R = riemann_numeric(m, r, theta, eps)
    g, _, _ = metric_tensors(m, r, theta)
    Rd = np.einsum("am,mbcd->abcd", g, R)
    anti_13 = float(np.max(np.abs(Rd + Rd.transpose(1, 0, 2, 3))))
    anti_24 = float(np.max(np.abs(Rd + Rd.transpose(0, 1, 3, 2))))
    pair = float(np.max(np.abs(Rd - Rd.transpose(2, 3, 0, 1))))
    bianchi1 = float(np.max(np.abs(
        Rd + Rd.transpose(0, 2, 3, 1) + Rd.transpose(0, 3, 1, 2))))
    scale = max(1.0, float(np.max(np.abs(Rd))))
    return {
        "antisym_first_pair": anti_13 / scale,
        "antisym_second_pair": anti_24 / scale,
        "pair_exchange": pair / scale,
        "first_bianchi": bianchi1 / scale,
    }


def check_contracted_bianchi(m: SSZMetric, r_grid: np.ndarray,
                             eps: float = 1e-5) -> float:
    """Contracted Bianchi identity div G = 0 (geometric identity for ANY
    metric):  nabla_mu G^{mu nu} = 0, evaluated nontrivially for nu = r.

    nabla_mu G^{mu r} = (1/sqrt(-g)) d_r (sqrt(-g) G^{rr})
                        + Gamma^r_{ab} G^{ab}.
    This is a sharp end-to-end test of the curvature stack (metric
    derivatives -> Christoffels -> Riemann -> Ricci -> Einstein)."""
    r_grid = np.asarray(r_grid, dtype=float)

    def sqrt_g(r):
        f = float(m.f(r))
        h = float(m.h(r))
        return r * r * float(np.sqrt(f / h))  # equatorial density

    def G_contrav(r):
        Gc = einstein_numeric(m, float(r))
        return Gc

    dens = np.array([sqrt_g(float(r)) for r in r_grid])
    Grr = np.array([G_contrav(float(r))[1, 1] for r in r_grid])
    flux = dens * Grr
    from scipy.interpolate import CubicSpline
    term1 = CubicSpline(r_grid, flux).derivative()(r_grid) / dens
    term2 = np.empty_like(r_grid)
    for i, r in enumerate(r_grid):
        Gc = G_contrav(float(r))
        G = m.Gamma_tensor(float(r))
        term2[i] = float(np.einsum("ab,ab->", G[1], Gc))
    div = term1 + term2
    scale = max(1.0, float(np.max(np.abs(div))))
    return float(np.max(np.abs(div)) / scale)


# ---------------------------------------------------------------------------
# known limits
# ---------------------------------------------------------------------------


@lru_cache(maxsize=1)
def schwarzschild_reference_metric() -> SSZMetric:
    """Exact Schwarzschild r_s = 1 (M = 1/2) through the canonical
    machinery — the analytic reference geometry (Ricci-flat)."""
    u = np.linspace(0.50, 0.95, 800)
    f = 1.0 - u
    return SSZMetric(
        u=u, f_u=CubicSpline(u, f), h_u=CubicSpline(u, f.copy()),
        member_hash="ANALYTIC_SCHWARZSCHILD_RS1",
        u_min=float(u.min()), u_max=float(u.max()),
    )


def weak_field_reference_metric(eps: float = 0.2) -> SSZMetric:
    """Weak-field GR-like reference f = h = 1 - eps*u (M = eps/2):
    used for the PPN-signature check of the pipeline (gamma = 1,
    Newtonian Kepler law)."""
    u = np.linspace(0.50, 0.95, 400)
    f = 1.0 - eps * u
    return SSZMetric(
        u=u, f_u=CubicSpline(u, f), h_u=CubicSpline(u, f.copy()),
        member_hash=f"ANALYTIC_WEAKFIELD_EPS{eps}",
        u_min=float(u.min()), u_max=float(u.max()),
    )


def ppn_signature_checks(m: SSZMetric, M: float) -> dict:
    """gamma = 1 signature (g_tt * g_rr = -1 pointwise, exact for f = h)
    and the Newtonian Kepler law Omega^2 = M/r^3 for the weak-field
    reference metric."""
    grid = np.linspace(float(m.u_min) + 0.01, float(m.u_max) - 0.01, 50)
    r_grid = 1.0 / grid
    prod = -float(m.f(r_grid[0])) * 0  # placeholder to keep types simple
    gamma_sig = 0.0
    for r in r_grid:
        prod = float(m.f(r)) / float(m.h(r))
        gamma_sig = max(gamma_sig, abs(prod - 1.0))
    om_mach = 0.0
    for r in r_grid:
        om2 = float(m.fp(r)) / (2.0 * r)
        om_mach = max(om_mach, abs(om2 - M / r**3) / (M / r**3))
    return {
        "gamma1_signature_max_abs": gamma_sig,
        "newtonian_kepler_max_rel": om_mach,
    }


# ---------------------------------------------------------------------------
# phase transport cross-checks (integral <-> differential)
# ---------------------------------------------------------------------------


def null_coordinate_time_ode(m: SSZMetric, r_emit: float, r_obs: float,
                             inward: bool = True) -> float:
    """Coordinate time for a radial null ray r_emit -> r_obs from the
    INTEGRATED null geodesic ODE (differential route, E = 1, b = 0)."""
    r_lo, r_hi = 1.0 / m.u_max + 1e-4, 1.0 / m.u_min - 1e-4
    if not (r_lo < min(r_emit, r_obs) and max(r_emit, r_obs) < r_hi):
        raise ValueError("radii outside member domain")
    f0 = float(m.f(r_emit))
    h0 = float(m.h(r_emit))
    kt0 = 1.0 / f0
    kr0 = float(np.sqrt(h0 * (1.0 / f0)))
    kr0 = -kr0 if inward else kr0

    def rhs(lam, y):
        _, r, kt, kr = y
        f = float(m.f(r))
        fp = float(m.fp(r))
        h = float(m.h(r))
        hp = float(m.hp(r))
        return [kt, kr, -(fp / f) * kt * kr,
                -(0.5 * h * fp) * kt * kt + (hp / (2.0 * h)) * kr * kr]

    def reach(lam, y):
        return y[1] - r_obs

    reach.terminal = True  # type: ignore[attr-defined]
    span = abs(r_obs - r_emit) * 6.0 / min(1.0, abs(kr0))
    sol = solve_ivp(rhs, (0.0, span), [0.0, r_emit, kt0, kr0],
                    method="DOP853", rtol=1e-12, atol=1e-12,
                    events=reach, dense_output=True)
    if not sol.success or not sol.t_events[0].size:
        raise RuntimeError("null ray did not reach the target sphere")
    lam_hit = float(sol.t_events[0][0])
    return float(sol.sol(lam_hit)[0])  # coordinate time t at the event


def null_coordinate_time_quadrature(m: SSZMetric, r_emit: float,
                                    r_obs: float, n: int = 20000) -> float:
    """Integral route: dt = integral dr / sqrt(f h) (Shapiro delay)."""
    lo, hi = sorted((r_emit, r_obs))
    grid = np.linspace(lo, hi, n)
    return float(np.trapezoid(
        1.0 / np.sqrt(m.f(grid) * m.h(grid)), grid))
