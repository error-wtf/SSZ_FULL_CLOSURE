"""Canonical exact symbolic comparison service (Section 31 policy).

One shared helper for exact expression equality.  Never invent per-test
simplification chains; never use fragile AST checks (exact Mul-tree has()).
Every comparison returns (bool, canonical_residual) so callers can record
the residual as evidence.
"""
from __future__ import annotations

import sympy as sp


def canonical_simplify(expr):
    """Deterministic canonicalization pipeline for exact symbolic identities."""
    return sp.factor(
        sp.cancel(
            sp.together(
                sp.expand(expr)
            )
        )
    )


def exact_equal(lhs, rhs):
    """Exact mathematical equality of two sympy expressions.

    Returns (is_equal, canonical_residual).  The residual is the
    canonicalized difference, suitable for evidence records.
    """
    diff = canonical_simplify(sp.sympify(lhs) - sp.sympify(rhs))
    return diff == 0, diff


def coefficient_equal(expr, sym, expected_coeff):
    """Exact coefficient check expr == expected_coeff * sym + (rest without sym).

    Robust replacement for expr.has(expected_Mul_tree)."""
    expr = sp.expand(sp.sympify(expr))
    coeff = expr.coeff(sp.sympify(sym))
    return exact_equal(coeff, expected_coeff)
