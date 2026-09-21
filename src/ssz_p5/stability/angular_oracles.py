"""Published-limit angular oracle helpers.

These functions are deliberately separated from the universal reducer.  They
consume already mapped asymptotic coefficients and therefore cannot silently
repair the common-action Laurent reduction.
"""
from __future__ import annotations
import numpy as np


def cubic_from_vector_and_coupled(c_v, b1, b2):
    """Ascending coefficients of (z-c_v)(z^2-b1*z+b2)."""
    c_v=np.asarray(c_v,float); b1=np.asarray(b1,float); b2=np.asarray(b2,float)
    return np.stack((-c_v*b2, b2+c_v*b1, -(b1+c_v), np.ones_like(c_v)),axis=-1)


def normalize_cubic(coeff, eps=1e-300):
    c=np.asarray(coeff,float)
    lead=c[...,3]
    return c/np.where(np.abs(lead)>eps,lead,np.nan)[...,None]


def maxwell_horndeski_gm_ghs_expected_cubic(n=1):
    """Exact EMD/GM-GHS angular oracle: (1-z)^3.

    Returned in ascending powers of z.  This is the published exact limiting
    polynomial; a background-specific action replay must still be supplied by
    the universal reducer audit before a production PASS is claimed.
    """
    c=np.array([1.0,-3.0,3.0,-1.0])
    return np.tile(c,(int(n),1))
