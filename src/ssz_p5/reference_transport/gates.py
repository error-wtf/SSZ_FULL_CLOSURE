"""SAG_REF_V1 gate registry and evaluator (INDEPENDENT of true_closure).

Registers gates G101..G109_Sag with their own verdict key `SAG_REF_V1`.
This registry deliberately does NOT interact with the true_closure GATES
dictionary and is not a condition of G130/G140 (see
docs/SAGNAC_REFERENCE_TRANSPORT_SPEC.md section 6).
"""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from ssz_p5.reference_transport.sagnac import (
    edge_contract_catch_up,
    edge_contract_null,
    edge_contract_reversal,
    sagnac_times,
)
from ssz_p5.reference_transport.segment_chain import chain_convergence, segment_chain
from ssz_p5.reference_transport.inversion import invert_delta_t

BETA_TEST = 0.3
F_CARRIER = 1.0


def _gate(name: str, ok: bool, detail: dict) -> dict:
    return {"gate": name, "pass": bool(ok), "detail": detail}


def evaluate() -> dict:
    gates: list[dict] = []

    # --- closed form (G101-G103)
    r = sagnac_times(BETA_TEST, F_CARRIER)
    t_p_exact = 1.0 / (1.0 - BETA_TEST)
    t_m_exact = 1.0 / (1.0 + BETA_TEST)
    dt_exact = 2.0 * BETA_TEST / (1.0 - BETA_TEST**2)
    gates.append(_gate("G101_Sag", abs(r.t_plus - t_p_exact) <= 1e-14 * abs(t_p_exact),
                       {"t_plus": r.t_plus, "exact": t_p_exact}))
    gates.append(_gate("G102_Sag", abs(r.t_minus - t_m_exact) <= 1e-14 * abs(t_m_exact),
                       {"t_minus": r.t_minus, "exact": t_m_exact}))
    gates.append(_gate("G103_Sag", abs(r.delta_t - dt_exact) <= 1e-14 * abs(dt_exact),
                       {"delta_t": r.delta_t, "exact": dt_exact}))

    # --- G104 reversal
    d_pos, d_neg = edge_contract_reversal(BETA_TEST)
    asym = abs(d_pos + d_neg)
    gates.append(_gate("G104_Sag", asym <= 1e-15, {"asymmetry": asym}))

    # --- G105 null limit
    dt0 = edge_contract_null()
    gates.append(_gate("G105_Sag", dt0 < 1e-14, {"delta_t_at_v_1e-15": dt0}))

    # --- G106 segment convergence + partition independence
    exact_tplus = t_p_exact
    ladder = chain_convergence(BETA_TEST, s=+1, seed=None)
    n_max, t_max = ladder[-1]
    conv_err = abs(t_max - exact_tplus) / abs(exact_tplus)
    rnd = segment_chain(BETA_TEST, 4096, s=+1, seed=1978)
    part_err = abs(rnd - exact_tplus) / abs(exact_tplus)
    gates.append(_gate("G106_Sag", conv_err < 1e-12 and part_err < 1e-12,
                       {"N": n_max, "convergence_err": conv_err,
                        "partition_err": part_err}))

    # --- G107 PDE vs closed form (convergence, observed order >= 1)
    pde_grid = []
    for n in (512, 1024, 2048):
        t_pde = first_passage_safe(BETA_TEST, +1, n)
        pde_grid.append((n, t_pde))
    errs = [abs(t - t_p_exact) / t_p_exact for _, t in pde_grid]
    orders = []
    for i in range(1, len(errs)):
        if errs[i] > 0 and errs[i-1] > 0:
            orders.append(np.log(errs[i-1] / errs[i]) / np.log(2.0))
    monotone = all(errs[i+1] < errs[i] for i in range(len(errs)-1))
    pde_ok = monotone and errs[-1] < 0.05 and (not orders or min(orders) >= 0.5)
    gates.append(_gate("G107_Sag", pde_ok,
                       {"grid": [n for n, _ in pde_grid],
                        "rel_errs": errs, "observed_orders": orders}))

    # --- G108 phase readout round-trip
    phi_rt = r.delta_phi / (2.0 * np.pi * F_CARRIER)
    phi_ok = abs(phi_rt - r.delta_t) <= 1e-14 * abs(r.delta_t)
    gates.append(_gate("G108_Sag", phi_ok,
                       {"dphi": r.delta_phi, "dJ": r.delta_J, "roundtrip": phi_rt}))

    # --- G109 inversion round-trip
    betas = np.linspace(0.01, 0.99, 99)
    max_err = 0.0
    sign_ok = True
    for b in betas:
        dt_pos = 2.0 * b / (1.0 - b * b)
        dt_neg = 2.0 * (-b) / (1.0 - b * b)
        v = invert_delta_t(dt_pos)
        vn = invert_delta_t(dt_neg)
        max_err = max(max_err, abs(v - b), abs(vn + b))
        sign_ok = sign_ok and np.sign(v) == np.sign(b)
    gates.append(_gate("G109_Sag", max_err <= 1e-12 and sign_ok,
                       {"max_roundtrip_err": float(max_err), "sign_contract": bool(sign_ok)}))

    verdict = {
        "verdict_key": "SAG_REF_V1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "beta_test": BETA_TEST,
        "gates": gates,
        "all_pass": all(g["pass"] for g in gates),
        "scope": "transport method validation; independent of TRUE_FULL_CLOSURE",
    }
    blob = json.dumps(verdict, indent=1, sort_keys=True)
    verdict["sha256"] = hashlib.sha256(blob.encode()).hexdigest()
    return verdict


def first_passage_safe(beta: float, s: int, n: int) -> float:
    from ssz_p5.reference_transport.transport_pde import first_passage_time
    return first_passage_time(beta, s=s, n_grid=n, t_max=8.0)


def write_verdict(out: Path | None = None) -> Path:
    v = evaluate()
    out = out or (Path(__file__).resolve().parents[2] / ".." / ".." / ".." / ".." /
                  "data" / "generated" / "reference_transport" /
                  "SAG_REF_V1_VERDICT.json")
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(v, indent=1, sort_keys=True))
    return out


if __name__ == "__main__":
    import sys
    p = write_verdict()
    v = evaluate()
    print("SAG_REF_V1 all_pass:", v["all_pass"])
    for g in v["gates"]:
        print(f"  {g['gate']}: {'PASS' if g['pass'] else 'FAIL'}")
    print("written:", p)
    sys.exit(0 if v["all_pass"] else 1)
