import numpy as np
from ssz_p5.stability.angular_universal_laurent import (
    load_genuine_svt_witness, canonical_laurent, extract_eq83_style_coefficients,
    action_m5_mass_shortcuts, leading_characteristic, roots_from_coeffs,
)

def _scaled(a,b): return np.abs(a-b)/np.maximum(1.0,np.abs(b))

def test_action_minus_branch_matches_m22_and_cross_identity():
    d=load_genuine_svt_witness()
    # thin grid keeps the unit test light; production audit runs the full profile
    d=d.iloc[::40].reset_index(drop=True)
    la=canonical_laurent(d)
    c=extract_eq83_style_coefficients(d,la)
    mn=action_m5_mass_shortcuts(d,'minus'); pl=action_m5_mass_shortcuts(d,'plus')
    assert np.nanmax(_scaled(c['M22_0'],mn['M22_0'])) < 1e-8
    ident=(d.a6*d.v6*d.v13/(d.a4*d.v9)).to_numpy(float)
    assert np.nanmax(_scaled(mn['M22_0']-pl['M22_0'],ident)) < 1e-8

def test_raw_characteristic_keeps_vector_root_near_luminal():
    d=load_genuine_svt_witness().iloc[::40].reset_index(drop=True)
    la=canonical_laurent(d); _,coef,_=leading_characteristic(d,la); roots=roots_from_coeffs(coef)
    delta=np.abs(roots-1.0)
    valid=np.any(np.isfinite(delta),axis=1)
    dist=np.min(np.where(np.isfinite(delta),delta,np.inf),axis=1)
    assert np.median(dist[valid]) < 5e-3
