from ssz_p5.production.inner_targets import (
    JET_ORDERS,
    endpoint_constraint_count,
    endpoint_derivative_orders,
)


def test_endpoint_orders_are_generated_from_contract_only():
    assert JET_ORDERS == {"v5": 1, "c3": 1, "e3": 0}
    assert endpoint_derivative_orders("v5") == (0, 1)
    assert endpoint_derivative_orders("c3") == (0, 1)
    assert endpoint_derivative_orders("e3") == (0,)


def test_endpoint_equation_count_matches_contract():
    # two sides: 2*(2 + 2 + 1) = 10 scalar endpoint equations
    assert endpoint_constraint_count() == 10
    assert endpoint_constraint_count(("e3",)) == 2
