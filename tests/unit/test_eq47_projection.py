import numpy as np
import pandas as pd
import sympy as sp
import json
from pathlib import Path
from ssz_p5.stability.laurent import LaurentSeries as LS
from ssz_p5.stability.angular_universal_laurent import leading_characteristic, roots_from_coeffs

def test_projected_relation_and_literal_residual_are_distinct():
    R,t,s,c,b,a=sp.symbols('R t s c b a',nonzero=True)
    k12=(t*s*c+b*a)/c
    mixed=k12*R*c-b*R*a
    assert sp.simplify(mixed-t*s*R*c)==0
    assert sp.simplify(mixed-t*s*c-(R-1)*t*s*c)==0
    assert sp.simplify(mixed-t*s*c)!=0

def test_leading_order_is_scale_invariant_and_not_hardcoded():
    d=pd.DataFrame({'x':[1.0],'f':[1.0]})
    for scale in (1e-12,1,1e12):
        k0=np.diag([0,1,0])[None]*scale
        k1=np.diag([1,0,0])[None]*scale
        k2=np.diag([0,0,1])[None]*scale
        la={'K':LS({0:k0,1:k1,2:k2}),'M':LS({-1:2*k0,0:k1,1:3*k2})}
        p,c,_=leading_characteristic(d,la)
        assert p[0]==3
        np.testing.assert_allclose(roots_from_coeffs(c)[0],[1,2,3],atol=1e-10)
        la={'K':LS({0:k0,1:k1+k2}),'M':LS({-1:2*k0,0:k1+3*k2})}
        assert leading_characteristic(d,la)[0][0]==2

def test_reference_action_and_projection_checkpoint():
    root=Path(__file__).resolve().parents[2]
    report=json.loads((root/'data/generated/angular_eq47_projection_2026-09-21/EQ47_PROJECTION_AUDIT.json').read_text())
    assert all(report['numerical_gates'].values())
    assert report['first_published_shortcut_divergence']=='EQ47_MIXED_MISSING_K_TO_M_PROJECTION_FACTOR'
    r=report['representative']
    assert r['leading_eps_power']==3
    np.testing.assert_allclose([r[f'action_root{i}'] for i in range(3)],[-164.97552195,1,4218.66284408],rtol=1e-6)
    np.testing.assert_allclose([r['action_B1'],r['action_B2']],[r['null_quadratic_B1'],r['null_quadratic_B2']],rtol=1e-8)
    assert report['finite_L_to_laurent']=='PASS'
    assert report['continuation_0p708_to_0p715']=='HELD'
    assert not report['absolute_full_closure']

def test_large_root_preserves_cubic_degree():
    c=np.polynomial.polynomial.polyfromroots([-1e16,1,2])
    roots=roots_from_coeffs([c])[0]
    assert np.isfinite(roots).all()
    np.testing.assert_allclose(roots,[-1e16,1,2],rtol=1e-6)
