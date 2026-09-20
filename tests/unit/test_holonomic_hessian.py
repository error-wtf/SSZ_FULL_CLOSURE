import numpy as np
import pandas as pd

from ssz_p5.production.holonomic_hessian import (
    assemble_hessian,
    audit_existing_split_controls,
    chain_residuals,
    implied_lower_from_transverse,
    response_from_hessian,
)


def background(n=41):
    x = np.linspace(1.4, 1.3, n)
    h = np.linspace(0.39, 0.41, n)
    f = np.linspace(0.29, 0.31, n)
    ph = -np.linspace(1.25, 1.20, n)
    X = -0.5 * h * ph**2
    A = np.linspace(1.2, 0.3, n)
    return pd.DataFrame(dict(u=1/x, x=x, f=f, h=h, phiprime=ph, X=X, A0prime=A))


def principal(n=41):
    return pd.DataFrame(
        dict(
            x=np.linspace(1.4, 1.3, n),
            f2XX=np.linspace(0.2, 0.4, n),
            f2XF=np.linspace(-0.1, 0.2, n),
            f2FF=np.linspace(0.3, -0.2, n),
        )
    )


def test_transverse_block_implies_exact_background_null_hessian():
    b = background()
    p = principal()
    lower = implied_lower_from_transverse(b, p)
    residual, normalized, _ = chain_residuals(b, lower, p)
    np.testing.assert_allclose(residual, 0.0, atol=2e-12)
    assert np.max(normalized) < 2e-12


def test_first_chain_rule_cannot_be_repaired_by_principal_block():
    b = background()
    p = principal()
    lower = implied_lower_from_transverse(b, p)
    lower["f2phiphi_control"] += 1.0
    r1, _, _ = chain_residuals(b, lower, p)

    # Altering XX/XF/FF cannot enter the first row of Ht.
    p2 = p.copy()
    p2[["f2XX", "f2XF", "f2FF"]] += 123.0
    r2, _, _ = chain_residuals(b, lower, p2)
    np.testing.assert_allclose(r1[:, 0], r2[:, 0], atol=1e-12)
    assert np.max(np.abs(r1[:, 0])) > 1e-3


def test_audit_flags_nonholonomic_split_controls():
    b = background()
    p = principal()
    lower = implied_lower_from_transverse(b, p)
    lower["f2phiF_control"] += 4.0
    _, report = audit_existing_split_controls(b, lower, p)
    assert report["status"] == "FAIL"
    assert report["first_chain_rule_lower_only"] is True
    assert report["principal_controls_can_repair_first_chain_rule"] is False


def test_appendix_a_lower_response_contains_cross_and_radial_terms():
    b = background(61)
    p = principal(61)
    lower = implied_lower_from_transverse(b, p)
    H = assemble_hessian(lower, p)
    got = response_from_hessian(b, H)

    r = b.x.to_numpy(float)
    f = b.f.to_numpy(float)
    h = b.h.to_numpy(float)
    A = b.A0prime.to_numpy(float)
    ph = b.phiprime.to_numpy(float)
    phiX = H[:, 0, 1]
    phiF = H[:, 0, 2]
    phiphi = H[:, 0, 0]

    expected_c3 = (
        0.5 * r**2 * np.sqrt(f * h) * ph**2 * phiX
        - 0.5 * r**2 * A**2 * np.sqrt(h / f) * phiF
    )
    np.testing.assert_allclose(got[:, 1], expected_c3, rtol=0, atol=1e-12)

    from ssz_p5.jets.jet9d8 import profile_derivative
    expected_e3 = (
        -0.5 * profile_derivative(r, r**2 * np.sqrt(f * h) * ph * phiX, 1, 9, 8)
        -0.5 * r**2 * np.sqrt(f / h) * phiphi
    )
    np.testing.assert_allclose(got[:, 2], expected_e3, rtol=0, atol=1e-12)
