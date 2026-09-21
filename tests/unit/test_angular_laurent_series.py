import numpy as np
from ssz_p5.stability.laurent import LaurentSeries

def test_formal_reciprocal():
    x=LaurentSeries({-1:np.array([2.0]),0:np.array([3.0]),1:np.array([5.0])},pmin=-4,pmax=8)
    y=x.reciprocal()
    for eps in (1e-2,2e-3):
        assert np.allclose((x*y).evaluate(eps),1.0,rtol=1e-10,atol=1e-10)
