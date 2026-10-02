#!/usr/bin/env python3
"""Targeted Weisz/Sokoloff high-order near-degenerate residue test.

Phys. Rev. B 18, 3275 states that in high-order rational approximants
93 lambda = 91 a and 62 lambda = 59 a, harmonics omega_{45 Q} and omega_{20 Q}
are nearly degenerate with omega_0, yet their residues are at least ten orders
of magnitude smaller than the omega_0 residue.

We independently solve the finite Frenkel-Kontorova equilibrium/dynamical
matrix and identify the eigenmode nearest the corresponding unperturbed harmonic
frequency.  This tests the key Green-function idea: near frequency degeneracy
does not imply comparable observable residue.

Units: a=alpha=M=1, beta=V0 Q^2/(2 alpha).
"""
from __future__ import annotations
import json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/generated/spectral/weisz1978/WEISZ1978_HIGH_ORDER_NEAR_DEGENERATE_RESIDUES.json"

CASES=(
    # paper examples. Weak beta values keep the rational approximants in the
    # regime where the perturbative near-degenerate argument applies.
    (93,91,45,0.02),
    (62,59,20,0.02),
)


def Ug(y,N,M,beta):
    j=np.arange(M,dtype=float);Q=2*np.pi*N/M;x=j+y
    dy=np.roll(y,-1)-y
    U=.5*np.sum(dy*dy)+(2*beta/Q**2)*np.sum(np.cos(Q*x))
    g=2*y-np.roll(y,1)-np.roll(y,-1)-(2*beta/Q)*np.sin(Q*x)
    return float(U),g


def Hmat(y,N,M,beta):
    j=np.arange(M,dtype=float);Q=2*np.pi*N/M;x=j+y
    H=np.diag(2-2*beta*np.cos(Q*x))
    for k in range(M):
        H[k,(k-1)%M]-=1;H[k,(k+1)%M]-=1
    return (H+H.T)/2


def solve(N,M,beta):
    Q=2*np.pi*N/M
    lam=2*np.pi/Q
    starts=[]
    for s in np.linspace(0,lam,13,endpoint=False):
        y=np.full(M,s)+0.005*np.sin(Q*(np.arange(M)+s))
        starts.append(y)
    starts.append(np.zeros(M))
    best=None
    for y0 in starts:
        r=minimize(lambda y:Ug(y,N,M,beta),y0,jac=True,method="L-BFGS-B",
                   options={"ftol":1e-14,"gtol":1e-10,"maxiter":30000})
        U,g=Ug(r.x,N,M,beta);H=Hmat(r.x,N,M,beta)
        mine=float(np.linalg.eigvalsh(H)[0])
        rec=(U,float(np.linalg.norm(g,np.inf)),mine,np.asarray(r.x,float))
        if best is None or ((mine>-1e-7 and rec[1]<1e-6),( -U)) > ((best[2]>-1e-7 and best[1]<1e-6),(-best[0])):
            best=rec
    U,gi,mine,y=best
    H=Hmat(y,N,M,beta)
    ev,V=np.linalg.eigh(H)
    q=np.ones(M)/np.sqrt(M)
    R=np.abs(V.T@q)**2
    return gi,ev,R


def target_unperturbed(N,M,n):
    Q=2*np.pi*N/M
    q=((n*Q+np.pi)%(2*np.pi))-np.pi
    return 2*(1-np.cos(q))


def main():
    results=[]
    for N,M,n,beta in CASES:
        gi,ev,R=solve(N,M,beta)
        dom=int(np.argmax(R))
        target=target_unperturbed(N,M,n)
        # Exclude the dominant q=0 mode when locating the high-order near-degenerate partner.
        idx=np.argsort(np.abs(ev-target))
        partner=next(int(i) for i in idx if i!=dom)
        ratio=float(R[partner]/max(R[dom],1e-300))
        results.append({
          "approximant":f"{N} lambda = {M} a",
          "beta":beta,"harmonic_n":n,
          "unperturbed_target_omega2":float(target),
          "dominant_mode_index":dom,
          "dominant_omega2":float(ev[dom]),
          "dominant_residue":float(R[dom]),
          "partner_mode_index":partner,
          "partner_omega2":float(ev[partner]),
          "partner_residue":float(R[partner]),
          "partner_to_dominant_ratio":ratio,
          "suppression_orders":float(-np.log10(max(ratio,1e-300))),
          "grad_inf":gi,
          "residue_sum":float(np.sum(R)),
          "ten_orders_suppressed":bool(ratio<1e-10),
        })
    report={
      "source_claim":"Near-degenerate high-order harmonics can have residues at least ten orders below the main q=0 mode.",
      "results":results,
      "all_residue_sum_rules":all(abs(r["residue_sum"]-1)<1e-10 for r in results),
      "all_equilibria_converged":all(r["grad_inf"]<1e-6 for r in results),
      "all_examples_ten_orders_suppressed":all(r["ten_orders_suppressed"] for r in results),
      "guard":"Mode identification uses nearest finite-system eigenvalue to the corresponding unperturbed harmonic target; this is an independent numerical proxy for the paper's harmonic label."
    }
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps(report,indent=2,allow_nan=False))
    return 0

if __name__=="__main__": raise SystemExit(main())
