#!/usr/bin/env python3
"""Common unreduced even-parity kernel for SSZ Horndeski + U(1)-SVT.

This module encodes the shared gauge-fixed quadratic-action structure of
Kase & Tsujikawa (2023), Eqs. (4.12)-(4.14), and Zhang & Kase (2024),
Eqs. (4.12)-(4.14), *before* theory-specific constraint elimination.

Fields (common gauge):
    H0, H1, H2, h1, dphi, dA0, dA1, V
where V is the auxiliary vector variable introduced to linearize
(dA0' - dot(dA1))^2.

The local symbol freezes radial background coefficients.  This is the
appropriate eikonal/principal-symbol operation: radial derivatives of the
background coefficients are lower differential order.  The full finite-l
operator/QNM problem requires restoring them later.

Important: this module does NOT claim a final hybrid PASS.  It supplies the
correct common operator architecture so Horndeski and genuine-SVT coefficient
contributions can be combined before the constraints are eliminated.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Mapping, Tuple, List
import numpy as np

FIELDS = ("H0", "H1", "H2", "h1", "dphi", "dA0", "dA1", "V")
FI = {name: i for i, name in enumerate(FIELDS)}

A_KEYS = tuple(f"a{i}" for i in range(1,10))
B_KEYS = tuple(f"b{i}" for i in range(1,6))
C_KEYS = tuple(f"c{i}" for i in range(1,7))
D_KEYS = tuple(f"d{i}" for i in range(1,5))
E_KEYS = tuple(f"e{i}" for i in range(1,5))
V_KEYS = tuple(f"v{i}" for i in range(1,14))
CANONICAL_KEYS = A_KEYS+B_KEYS+C_KEYS+D_KEYS+E_KEYS+V_KEYS

# A term is coefficient * (dt^nt dr^nr q_i) * (dt^mt dr^mr q_j).
Term = Tuple[complex, int, Tuple[int,int], int, Tuple[int,int]]


def canonicalize_mh_vector(mh: Mapping[str,float]) -> Dict[str,float]:
    """Map the Maxwell-Horndeski vector numbering into the 13-slot ZK schema.

    Kase-Tsujikawa (2023) uses v1..v9:
      v1 flux^2; v2 H0 flux; v3 H2 flux; v4 phi' flux; v5 phi flux;
      v6 H0^2; v7 h1 A0; v8 A0^2; v9 A1^2.

    Zhang-Kase (2024) inserts a new flux-h1 coefficient at canonical v6 and
    therefore shifts the old algebraic vector slots by one:
      V7=MH v6, V8=MH v7, V9=MH v8, V10=MH v9,
      with V6=V11=V12=V13=0 in Maxwell-Horndeski.
    """
    out = {k: float(mh.get(k,0.0)) for k in A_KEYS+B_KEYS+C_KEYS+D_KEYS+E_KEYS}
    for i in range(1,14): out[f"v{i}"] = 0.0
    for i in range(1,6): out[f"v{i}"] = float(mh.get(f"v{i}",0.0))
    out["v6"] = 0.0
    out["v7"] = float(mh.get("v6",0.0))
    out["v8"] = float(mh.get("v7",0.0))
    out["v9"] = float(mh.get("v8",0.0))
    out["v10"] = float(mh.get("v9",0.0))
    return out


def blank_coefficients() -> Dict[str,float]:
    return {k: 0.0 for k in CANONICAL_KEYS}


def combine_coefficients(*pieces: Mapping[str,float]) -> Dict[str,float]:
    """Add already-baseline-separated contributions in the common schema."""
    out = blank_coefficients()
    for p in pieces:
        for k,v in p.items():
            if k in out: out[k] += float(v)
    return out


def _add(terms: List[Term], coef, f1, d1=(0,0), f2=None, d2=(0,0)):
    if f2 is None: f2=f1
    if coef != 0:
        terms.append((complex(coef), FI[f1], tuple(d1), FI[f2], tuple(d2)))


def build_terms(c: Mapping[str,float], L: float) -> List[Term]:
    """Build the frozen-coefficient common auxiliary-field Lagrangian."""
    q = {k: float(c.get(k,0.0)) for k in CANONICAL_KEYS}
    T: List[Term] = []

    # Lu, common Eq. (4.12)
    _add(T,q['a1'],'H0',(0,0),'dphi',(0,2))
    _add(T,q['a2'],'H0',(0,0),'dphi',(0,1))
    _add(T,q['a3'],'H0',(0,0),'H2',(0,1))
    _add(T,L*q['a4'],'H0',(0,0),'h1',(0,1))
    _add(T,q['a5']+L*q['a6'],'H0',(0,0),'dphi',(0,0))
    _add(T,q['a7']+L*q['a8'],'H0',(0,0),'H2',(0,0))
    _add(T,L*q['a9'],'H0',(0,0),'h1',(0,0))
    _add(T,L*q['b1'],'H1',(0,0),'H1',(0,0))
    _add(T,q['b2'],'H1',(0,0),'dphi',(1,1))
    _add(T,q['b3'],'H1',(0,0),'dphi',(1,0))
    _add(T,q['b4'],'H1',(0,0),'H2',(1,0))
    _add(T,L*q['b5'],'H1',(0,0),'h1',(1,0))
    _add(T,q['c1'],'dphi',(1,0),'H2',(1,0))
    _add(T,q['c2'],'H2',(0,0),'dphi',(0,1))
    _add(T,q['c3']+L*q['c4'],'H2',(0,0),'dphi',(0,0))
    _add(T,L*q['c5'],'H2',(0,0),'h1',(0,0))
    _add(T,q['c6'],'H2',(0,0),'H2',(0,0))
    _add(T,L*q['d1'],'h1',(1,0),'h1',(1,0))
    _add(T,L*q['d2'],'h1',(0,0),'dphi',(0,1))
    _add(T,L*q['d3'],'h1',(0,0),'dphi',(0,0))
    _add(T,L*q['d4'],'h1',(0,0),'h1',(0,0))
    _add(T,q['e1'],'dphi',(1,0),'dphi',(1,0))
    _add(T,q['e2'],'dphi',(0,1),'dphi',(0,1))
    _add(T,q['e3']+L*q['e4'],'dphi',(0,0),'dphi',(0,0))

    # LA after auxiliary V, common ZK Eq. (4.14).
    v1=q['v1']
    if abs(v1) < 1e-300:
        raise ZeroDivisionError('v1=0: auxiliary-field representation is singular')

    _add(T, 2*v1, 'V',(0,0),'dA0',(0,1))
    _add(T,-2*v1, 'V',(0,0),'dA1',(1,0))
    _add(T,-v1, 'V',(0,0),'V',(0,0))

    # S=v2 H0+v3 H2+v4 dphi'+v5 dphi+L v6 h1.
    S = [
        (q['v2'],'H0',(0,0)),
        (q['v3'],'H2',(0,0)),
        (q['v4'],'dphi',(0,1)),
        (q['v5'],'dphi',(0,0)),
        (L*q['v6'],'h1',(0,0)),
    ]
    for coef,field,der in S:
        _add(T,coef,'V',(0,0),field,der)  # + V S
    for i,(ci,fi,di) in enumerate(S):
        _add(T,-ci*ci/(4*v1),fi,di,fi,di)
        for cj,fj,dj in S[i+1:]:
            _add(T,-ci*cj/(2*v1),fi,di,fj,dj)

    _add(T,0.5*L*q['v6'],'h1',(0,0),'dA1',(1,0))
    _add(T,q['v7'],'H0',(0,0),'H0',(0,0))
    _add(T,L*q['v8'],'h1',(0,0),'dA0',(0,0))
    _add(T,L*q['v9'],'dA0',(0,0),'dA0',(0,0))
    _add(T,L*q['v10'],'dA1',(0,0),'dA1',(0,0))
    _add(T,L*q['v11'],'H1',(0,0),'dA1',(0,0))
    _add(T,L*q['v12'],'H2',(0,0),'dA0',(0,0))
    _add(T,L*q['v13'],'dphi',(0,0),'dA0',(0,0))
    return T


def symbol_matrix(c: Mapping[str,float], omega: float, kr: float, L: float) -> np.ndarray:
    """Euler-Lagrange Fourier symbol of the common frozen quadratic action.

    Convention q ~ exp(-i omega t + i kr r).  The result is Hermitian up to
    floating roundoff for real coefficients and real omega,kr.
    """
    dt = -1j*omega
    dr =  1j*kr
    P = np.zeros((len(FIELDS),len(FIELDS)),complex)
    for coef,i,a,j,b in build_terms(c,L):
        na=sum(a); nb=sum(b)
        dprod=(dt**(a[0]+b[0]))*(dr**(a[1]+b[1]))
        if i==j:
            P[i,i] += coef*((-1)**na + (-1)**nb)*dprod
        else:
            P[i,j] += coef*((-1)**na)*dprod
            P[j,i] += coef*((-1)**nb)*dprod
    return P


def h0_constraint_expected(c: Mapping[str,float], omega:float, kr:float, L:float) -> np.ndarray:
    """Expected frozen H0 row from direct variation of the auxiliary action."""
    q={k:float(c.get(k,0.0)) for k in CANONICAL_KEYS}
    dt=-1j*omega; dr=1j*kr
    row=np.zeros(len(FIELDS),complex)
    v1=q['v1']
    row[FI['dphi']] = q['a1']*dr**2 + (q['a2']-q['v2']*q['v4']/(2*v1))*dr + q['a5']+L*q['a6']-q['v2']*q['v5']/(2*v1)
    row[FI['H2']] = q['a3']*dr + q['a7']+L*q['a8']-q['v2']*q['v3']/(2*v1)
    row[FI['h1']] = L*q['a4']*dr + L*(q['a9']-q['v2']*q['v6']/(2*v1))
    row[FI['V']] = q['v2']
    # H0^2 cancels iff v7=v2^2/(4v1)
    row[FI['H0']] = 2*(q['v7']-q['v2']**2/(4*v1))
    return row


def validate_common_structure(c: Mapping[str,float], r:float|None=None, f:float|None=None,
                              tol:float=1e-10) -> Dict[str,float|bool]:
    q={k:float(c.get(k,0.0)) for k in CANONICAL_KEYS}
    out={}
    out['H0_square_identity_residual']=q['v7']-q['v2']**2/(4*q['v1']) if q['v1'] else np.nan
    if r is not None:
        out['a3_plus_r_a4']=q['a3']+r*q['a4']
    if r is not None and f is not None:
        out['b4_minus_2r_a4_over_f']=q['b4']-2*r*q['a4']/f
        out['b5_plus_a4_over_f']=q['b5']+q['a4']/f
        out['d1_minus_a4_over_2f']=q['d1']-q['a4']/(2*f)
    finite=[abs(float(v)) for v in out.values() if np.isfinite(v)]
    out['structural_identities_pass']=bool(not finite or max(finite)<tol)
    return out


def self_test() -> Dict[str,float|bool]:
    rng=np.random.default_rng(20260912)
    c={k:float(rng.normal()) for k in CANONICAL_KEYS}
    c['v1']=1.7
    c['v2']=0.8
    c['v7']=c['v2']**2/(4*c['v1'])
    r=1.51; f=.42
    c['a4']=.31; c['a3']=-r*c['a4']
    c['b4']=2*r*c['a4']/f; c['b5']=-c['a4']/f; c['d1']=c['a4']/(2*f)
    om=.73; kr=1.17; L=30.0
    P=symbol_matrix(c,om,kr,L)
    expected=h0_constraint_expected(c,om,kr,L)
    herm=np.max(np.abs(P-P.conj().T))
    h0err=np.max(np.abs(P[FI['H0']]-expected))
    ids=validate_common_structure(c,r,f)
    return {
        'hermiticity_error':float(herm),
        'H0_constraint_row_error':float(h0err),
        'H0_diagonal_abs':float(abs(P[FI['H0'],FI['H0']])),
        'structural_identities_pass':bool(ids['structural_identities_pass']),
        'pass':bool(herm<1e-12 and h0err<1e-12 and abs(P[FI['H0'],FI['H0']])<1e-12 and ids['structural_identities_pass'])
    }


if __name__=='__main__':
    print(self_test())
