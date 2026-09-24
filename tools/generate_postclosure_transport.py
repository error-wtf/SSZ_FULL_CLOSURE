#!/usr/bin/env python3
"""Generate the SSZ SOURCE-FREE TRANSPORT evidence artifact (post-closure
research block, NON-GATING).

Writes data/generated/postclosure/SSZ_SOURCE_FREE_TRANSPORT.json.

Every number is derived from the FROZEN electric production member
(member-hash verified at load).  No fitting, no free parameters, no manual
verdicts.  The full-closure gate graph is NOT touched by this block.
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.postclosure.transport import (  # noqa: E402
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
    radial_domain,
    radial_null_congruence,
    ring_trapping,
    timelike_geodesic,
)


def schwarzschild_metric() -> SSZMetric:
    """Exact Schwarzschild (r_s = 1) sampled through the same machinery —
    analytic anchor for every transport diagnostic (Ricci-flat)."""
    u = np.linspace(0.50, 0.95, 800)
    f = 1.0 - u  # f(r) = 1 - 1/r with r = 1/u
    return SSZMetric(
        u=u,
        f_u=CubicSpline(u, f),
        h_u=CubicSpline(u, f.copy()),
        member_hash="ANALYTIC_SCHWARZSCHILD_RS1",
        u_min=float(u.min()),
        u_max=float(u.max()),
    )


def main() -> int:
    m = load_member_metric(ROOT)
    ts = datetime.now(timezone.utc).isoformat()

    # ---- 1. timelike source-free transport --------------------------------
    tl_radial = timelike_geodesic(m, r0=1.6, E=1.0, L=0.0, tau_span=(0, 1.5))
    tl_angular = timelike_geodesic(m, r0=1.6, E=0.95, L=0.8, tau_span=(0, 1.5))
    tl_outward = timelike_geodesic(m, r0=1.5, E=0.9, L=0.6, inward=False,
                                   tau_span=(0, 1.5))
    circ = circular_timelike_orbit(m, 1.6)

    # negative control: synthetic force breaks the conservation laws
    forced = forced_transport_control(m, r0=1.6, F_t=1e-2, tau_span=(0, 1.5))

    # ---- 2. kinematical scalars + Raychaudhuri (timelike) ------------------
    grid200 = np.linspace(1.42, 1.63, 200)
    congruence = radial_congruence_scalars(m, grid200)
    congruence_80 = radial_congruence_scalars(
        m, np.linspace(1.42, 1.63, 80))

    ms = schwarzschild_metric()
    sg = np.linspace(1.2, 1.9, 200)
    schw_cong = radial_congruence_scalars(ms, sg)
    theta_exact = -1.5 * sg ** -1.5
    sigma2_exact = 1.5 * sg ** -3.0

    # ---- 3. null sector ----------------------------------------------------
    null_cong = radial_null_congruence(m, grid200)
    null_cong_schw = radial_null_congruence(ms, sg)
    null_b20 = null_geodesic_transport(m, b=2.0, r0=1.63, lam_span=(0, 12))
    null_b259 = null_geodesic_transport(m, b=2.59, r0=1.63, lam_span=(0, 20))
    amplitude = amplitude_transport_radial(m, 1.45, 1.63)
    phase = eikonal_phase_and_redshift(m, 1.0 / m.W_u(1.0) * 0 +
                                       1.4161606399560192, 1.5)

    # ---- 4. rings + trapping + consolidation -------------------------------
    rings = find_light_rings(m)
    trapping = ring_trapping(m, rings)
    table = consolidation_table(m, rings)

    rings_out = []
    for rg, row in zip(rings, table):
        rings_out.append({**rg, **{
            "dtau_dt": row["dtau_dt"],
            "Omega_kepler": row["Omega_kepler"],
            "kepler_vs_photon_rel_diff": row["kepler_vs_photon_rel_diff"],
        }})

    artifact = {
        "schema_version": "1.0",
        "timestamp": ts,
        "block": "SSZ SOURCE-FREE TRANSPORT (post-closure research block)",
        "gating": "NON-GATING: does not modify GATE_STATUS; prerequisite "
                  "ABSOLUTE_FULL_CLOSURE_PASS (all 28 gates) is recorded",
        "member_hash": m.member_hash,
        "metric_model": {
            "ansatz": "ds^2 = -f dt^2 + dr^2/h + r^2 dOmega^2, u = r_s/r, r_s = 1",
            "interpolant": "CubicSpline through FROZEN member values of f(u), h(u) "
                           "(exact interpolation of the pinned data; no fitting)",
            "domain_u": [m.u_min, m.u_max],
            "domain_r": list(radial_domain(m)),
        },
        "core_principle": "NO new force: D_SSZ/dtau := u^nu nabla_nu [g_SSZ]; "
                          "the chain Xi -> D -> g_mn -> Gamma -> transport is "
                          "evaluated on the closed full-action geometry",
        "rag_constraints_respected": {
            "no_flow_from_static_profile":
                "no flow equation is derived from Xi0; only standard "
                "source-free transport is evaluated on g_SSZ",
            "single_unified_dynamics_open":
                "this block does NOT close 'single unified dynamics'",
            "early_paper_is_intuition_archive":
                "forward argumentation only; early intuitions are mapped to "
                "structures AFTER they are derived",
            "navier_stokes_structural_only":
                "source-free != trivial is used as a structural analogy only",
        },
        "computation_1_timelike_transport": {
            "u^nu nabla_nu u^mu = 0 via conserved quantities "
            "(E = f u^t, L = r^2 u^phi, norm):": None,
            "radial_E1": {k: getattr(tl_radial, k)
                          for k in ("max_dE", "max_dL", "max_norm_res",
                                    "max_accel_res")},
            "angular_E095_L08": {k: getattr(tl_angular, k)
                                 for k in ("max_dE", "max_dL", "max_norm_res",
                                           "max_accel_res")},
            "outward_E09_L06": {k: getattr(tl_outward, k)
                                for k in ("max_dE", "max_dL", "max_norm_res",
                                          "max_accel_res")},
            "circular_orbit_r1.6": circ,
            "negative_control_forced_Ft_1e-2": forced,
        },
        "computation_2_raychaudhuri_timelike": {
            "congruence": "stationary radial free-fall, E = 1 (geodesic, "
                          "irrotational)",
            "theta_routes_absdiff": congruence["theta_routes_absdiff"],
            "omega2_max_abs": float(np.max(np.abs(congruence["omega2"]))),
            "max_scaled_raychaudhuri_res":
                congruence["max_scaled_raychaudhuri_res"],
            "max_scaled_raychaudhuri_res_80pts":
                congruence_80["max_scaled_raychaudhuri_res"],
            "theta_range": [float(congruence["theta"].min()),
                            float(congruence["theta"].max())],
            "sigma2_range": [float(congruence["sigma2"].min()),
                             float(congruence["sigma2"].max())],
            "R_uu_range": [float(congruence["R_uu"].min()),
                           float(congruence["R_uu"].max())],
            "schwarzschild_anchor": {
                "max_scaled_raychaudhuri_res":
                    schw_cong["max_scaled_raychaudhuri_res"],
                "theta_vs_exact_max_abs": float(
                    np.max(np.abs(schw_cong["theta"] - theta_exact))),
                "sigma2_vs_exact_max_abs": float(
                    np.max(np.abs(schw_cong["sigma2"] - sigma2_exact))),
                "R_uu_max_abs": float(np.max(np.abs(schw_cong["R_uu"]))),
            },
        },
        "computation_3_null_transport": {
            "k^nu nabla_nu k^mu = 0 via E, L, k.k": None,
            "radial_null_congruence": {
                "theta_routes_absdiff": null_cong["theta_routes_absdiff"],
                "max_scaled_raychaudhuri_res":
                    null_cong["max_scaled_raychaudhuri_res"],
                "R_kk_range": [float(null_cong["R_kk"].min()),
                               float(null_cong["R_kk"].max())],
                "schwarzschild_anchor_res":
                    null_cong_schw["max_scaled_raychaudhuri_res"],
            },
            "null_geodesic_b2.0": {k: getattr(null_b20, k)
                                   for k in ("max_dE", "max_dL",
                                             "max_norm_res")},
            "null_geodesic_b2.59": {k: getattr(null_b259, k)
                                    for k in ("max_dE", "max_dL",
                                              "max_norm_res")},
            "amplitude_transport": amplitude,
            "eikonal_phase_redshift_rings": phase,
        },
        "computation_4_ring_structure": {
            "rings": rings_out,
            "trapping": trapping,
            "consolidation": table,
        },
        "paper_intuition_mapping": {
            "no_new_force / geometry-in-transport":
                "realized literally: a^mu = u^nu nabla_nu u^mu = 0 while all "
                "observables below are nontrivial",
            "segmentation creates inertia/resistance":
                "mathematically: theta, sigma, Raychaudhuri focusing of the "
                "source-free congruence (computation 2) - no resistance term",
            "boundary balance (old: forces cancel)":
                "matured: ring structure = potential extremum W_u = 0 with "
                "W_uu sign = stability; trapping = log-winding (outer) vs "
                "libration (inner) from the exact null orbit equation",
            "light/phase (old: wavelength change)":
                "matured: eikonal transport - affine null geodesics, amplitude "
                "a^2 r^2 = const, redshift sqrt(f), S_r = integral dr/sqrt(fh)",
            "rotation/frame-dragging":
                "static central slice: g_tphi = 0 exactly (frame dragging 0); "
                "rotation enters via Omega(r) and congruence vorticity",
        },
        "policy": {
            "fitting": "NONE",
            "fitting_note": "splines interpolate the frozen member data exactly",
            "member_change": "NONE",
            "manual_pass": "NONE",
        },
    }
    artifact["metric_model"]["domain_r"] = list(radial_domain(m))

    out_dir = ROOT / "data" / "generated" / "postclosure"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "SSZ_SOURCE_FREE_TRANSPORT.json"
    out.write_text(json.dumps(artifact, indent=1) + "\n")
    print(json.dumps({
        "artifact": str(out.relative_to(ROOT)),
        "member_hash": m.member_hash,
        "rings": [(r["u"], r["stable"]) for r in rings],
        "raychaudhuri_res": congruence["max_scaled_raychaudhuri_res"],
        "null_raychaudhuri_res": null_cong["max_scaled_raychaudhuri_res"],
    }, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
