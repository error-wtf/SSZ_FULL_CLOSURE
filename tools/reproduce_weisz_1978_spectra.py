#!/usr/bin/env python3
"""Independent numerical reproduction of Weisz et al., Phys. Rev. B 18, 3275 (1978).

Implements the finite Frenkel-Kontorova rational approximants used in the paper:
  13 lambda = 10 a, beta=P = 0.3 and 0.7
  48 lambda = 49 a, beta=P = 0.005, 0.01 and 0.02

Units: a=alpha=M=1. Then Q=2*pi*N/M and beta=V0 Q^2/(2 alpha).
For a periodic displacement y_j (x_j=j a+y_j),
  U/alpha = 1/2 sum_j (y_{j+1}-y_j)^2 + (2 beta/Q^2) sum_j cos(Q x_j)
and the equilibrium equations are
  2 y_j-y_{j-1}-y_{j+1} - (2 beta/Q) sin(Q x_j)=0.
The dimensionless dynamical matrix has eigenvalues Omega^2=M omega^2/alpha:
  H_jj = 2 - 2 beta cos(Q x_j), H_j,j+/-1=-1.
The q=0 optical residue of a normalized mode v_n is
  R_n = |<1/sqrt(M),v_n>|^2,
and sum_n R_n=1.

This directly reproduces the paper's mode-frequency / optical-weight calculation,
without using any SSZ output.  It then emits machine-readable benchmark spectra
for later cross-model spectroscopy tests.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize

ROOT=Path(__file__).resolve().parents[1]
OUTDIR=ROOT/"data/generated/spectral/weisz1978"
REPORT=OUTDIR/"WEISZ1978_FK_REPRODUCTION.json"
CSV=OUTDIR/"WEISZ1978_FK_MODES.csv"


def potential_grad(y: np.ndarray, N: int, M: int, beta: float):
    j=np.arange(M,dtype=float)
    Q=2*np.pi*N/M
    x=j+y
    dy=np.roll(y,-1)-y
    U=0.5*np.sum(dy*dy)+(2*beta/Q**2)*np.sum(np.cos(Q*x))
    grad=2*y-np.roll(y,1)-np.roll(y,-1)-(2*beta/Q)*np.sin(Q*x)
    return float(U),grad


def hessian(y: np.ndarray,N:int,M:int,beta:float):
    j=np.arange(M,dtype=float)
    Q=2*np.pi*N/M
    x=j+y
    H=np.diag(2-2*beta*np.cos(Q*x))
    for k in range(M):
        H[k,(k-1)%M]-=1
        H[k,(k+1)%M]-=1
    return (H+H.T)/2


def solve_equilibrium(N:int,M:int,beta:float):
    Q=2*np.pi*N/M
    lam=2*np.pi/Q
    starts=[]
    # Translation phase is physical relative to the substrate.  Search the whole
    # substrate period plus small FK perturbative distortions.
    for s in np.linspace(0,lam,17,endpoint=False):
        y=np.full(M,s)
        y += 0.01*np.sin(Q*(np.arange(M)+s))
        starts.append(y)
    # continuation-like zero start
    starts.append(np.zeros(M))

    sols=[]
    for y0 in starts:
        res=minimize(lambda y: potential_grad(y,N,M,beta),y0,jac=True,
                     method="L-BFGS-B",options={"ftol":1e-14,"gtol":1e-11,"maxiter":20000})
        y=np.asarray(res.x,float)
        U,g=potential_grad(y,N,M,beta)
        H=hessian(y,N,M,beta)
        ev=np.linalg.eigvalsh(H)
        sols.append((U,float(np.linalg.norm(g,np.inf)),float(ev[0]),res.success,y))
    # Prefer stable extrema, then lowest energy.
    stable=[s for s in sols if s[2]>-1e-8 and s[1]<1e-7]
    pool=stable if stable else sols
    best=min(pool,key=lambda z:z[0])
    return best, sols


def spectrum(N:int,M:int,beta:float):
    (U,ginf,minh,ok,y),allsol=solve_equilibrium(N,M,beta)
    H=hessian(y,N,M,beta)
    lam,V=np.linalg.eigh(H)
    q0=np.ones(M)/np.sqrt(M)
    R=np.abs(V.T@q0)**2
    order=np.argsort(lam)
    lam=lam[order];R=R[order];V=V[:,order]
    # Numerical noise can yield tiny negative eigenvalues on the phason mode.
    omega=np.sqrt(np.maximum(lam,0))
    idx=np.argsort(-R)
    dominant=int(idx[0])
    return {
      "N":N,"M":M,"beta":beta,"Q":2*np.pi*N/M,
      "energy":U,"grad_inf":ginf,"hessian_min":minh,"optimizer_success":bool(ok),
      "omega2":lam,"omega":omega,"residue":R,
      "dominant_index":dominant,
      "dominant_omega2":float(lam[dominant]),
      "residue_sum":float(R.sum()),
      "n_eff":float(1/np.sum(R*R)),
      "entropy":float(-np.sum(np.where(R>0,R*np.log(R),0))),
      "top_indices":[int(i) for i in idx[:10]],
      "top_weights":[float(R[i]) for i in idx[:10]],
      "equilibrium_y":y,
    }


def response_grid(spec,delta=0.1,n=2500):
    lam=spec["omega2"];R=spec["residue"]
    hi=max(0.65,float(np.nanmax(lam)*1.02))
    x=np.linspace(0,hi,n)
    # Im [1/(lambda-(x+i delta))] = delta/[(lambda-x)^2+delta^2]
    im=np.sum(R[:,None]*delta/((lam[:,None]-x[None,:])**2+delta**2),axis=0)
    return x,im


def summarize(spec):
    idx=np.argsort(-spec["residue"])
    return {
      "N_lambda_eq_M_a":f"{spec['N']} lambda = {spec['M']} a",
      "beta":spec["beta"],
      "grad_inf":spec["grad_inf"],
      "hessian_min":spec["hessian_min"],
      "residue_sum":spec["residue_sum"],
      "dominant_omega2":spec["dominant_omega2"],
      "lowest_omega2":float(np.min(spec["omega2"])),
      "top_weights":[float(spec["residue"][i]) for i in idx[:8]],
      "top_omega2":[float(spec["omega2"][i]) for i in idx[:8]],
      "effective_mode_count":spec["n_eff"],
      "spectral_entropy":spec["entropy"],
    }


def main():
    cases=[(13,10,0.3),(13,10,0.7),(48,49,0.005),(48,49,0.01),(48,49,0.02)]
    specs=[spectrum(*c) for c in cases]

    rows=[]
    for s in specs:
        dom=max(float(np.max(s["residue"])),1e-300)
        for n,(o2,w,r) in enumerate(zip(s["omega2"],s["omega"],s["residue"])):
            rows.append({
              "N":s["N"],"M":s["M"],"beta":s["beta"],"mode":n,
              "omega2":float(o2),"omega":float(w),"q0_residue":float(r),
              "residue_over_dominant":float(r/dom),
            })
    OUTDIR.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(CSV,index=False)

    small03=specs[0];small07=specs[1];hi=specs[-1]
    hi_sorted=np.argsort(hi["omega2"])
    # successive low-frequency modes and their q=0 strength, the Fig.6 quantity
    fig6=[{
      "mode":int(i),"omega2":float(hi["omega2"][i]),
      "ratio":float(hi["residue"][i]/max(np.max(hi["residue"]),1e-300))
    } for i in hi_sorted[:min(16,len(hi_sorted))]]

    checks={
      "equilibrium_converged_all":all(s["grad_inf"]<1e-7 for s in specs),
      "residue_sum_rule_all":all(abs(s["residue_sum"]-1)<1e-10 for s in specs),
      "paper_unpinned_near_zero_mode_P03":bool(np.min(small03["omega2"])<5e-3),
      "paper_locked_finite_mode_P07":bool(np.min(small07["omega2"])>1e-2),
      "few_modes_dominate_13_10_P03":bool(sum(sorted(small03["residue"],reverse=True)[:3])>0.99),
      "high_order_48_49_P02_has_strong_residue_hierarchy":bool(
          np.max(hi["residue"])/max(np.partition(hi["residue"],-10)[-10],1e-300)>1e3
      ),
    }

    # Fit exponential decay of the strongest low-frequency harmonic-like residues
    # only after sorting by omega^2 and excluding the dominant zero/near-zero mode.
    low=np.argsort(hi["omega2"])[:16]
    rr=np.array([hi["residue"][i]/np.max(hi["residue"]) for i in low])
    oo=np.array([hi["omega2"][i] for i in low])
    m=(rr>1e-12)&(rr<0.8)
    if np.sum(m)>=3:
        coef=np.polyfit(oo[m],np.log10(rr[m]),1)
        decay={"log10_ratio_per_Omega2":float(coef[0]),"intercept":float(coef[1]),"points":int(np.sum(m))}
    else:
        decay={"log10_ratio_per_Omega2":None,"intercept":None,"points":int(np.sum(m))}

    report={
      "source":"Weisz, Cardarelli, Sacco, Sokoloff, Phys. Rev. B 18, 3275 (1978)",
      "implementation":"independent finite FK equilibrium + dynamical matrix + q=0 residues",
      "cases":[summarize(s) for s in specs],
      "fig6_style_P002_low_frequency_points":fig6,
      "P002_decay_fit":decay,
      "checks":checks,
      "overall_reproduction_pass":bool(all(checks.values())),
      "guard":"This reproduces the equations/parameter cases numerically; it does not treat the printed plot pixels as exact tabulated data."
    }
    REPORT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps(report,indent=2,allow_nan=False))
    return 0 if checks["equilibrium_converged_all"] and checks["residue_sum_rule_all"] else 1

if __name__=="__main__":
    raise SystemExit(main())
