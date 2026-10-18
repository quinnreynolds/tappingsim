"""Unit tests for tappingsim.components.bedmodels module."""

import pytest
from tappingsim.components.bedmodels import bedmodel_ergun, bedmodel_kozenycarman


class TestBedmodelKozenyCarman:
    def test_returns_zero_a_const(self):
        """Kozeny-Carman is laminar-only so the quadratic A term is always 0."""
        a, b = bedmodel_kozenycarman(0.1, 0.2, 20, 0.02, 0.8, 0.5, 0.1, 3000)
        assert a == 0

    def test_b_const_positive(self):
        a, b = bedmodel_kozenycarman(0.1, 0.2, 20, 0.02, 0.8, 0.5, 0.1, 3000)
        assert b > 0

    def test_b_const_value(self):
        a, b = bedmodel_kozenycarman(0.1, 0.2, 20, 0.02, 0.8, 0.5, 0.1, 3000)
        assert b == pytest.approx(1740.234375, rel=1e-6)

    def test_cavity_larger_than_taphole(self):
        """When bedmindiameter > tapholediameter a different branch is taken."""
        a, b = bedmodel_kozenycarman(0.1, 0.05, 20, 0.02, 0.8, 0.5, 0.1, 3000)
        assert a == 0
        assert b > 0


class TestBedmodelErgun:
    def test_both_consts_positive(self):
        a, b = bedmodel_ergun(0.1, 0.2, 20, 0.02, 0.8, 0.5, 0.1, 3000)
        assert a > 0
        assert b > 0

    def test_a_const_value(self):
        a, b = bedmodel_ergun(0.1, 0.2, 20, 0.02, 0.8, 0.5, 0.1, 3000)
        assert a == pytest.approx(683.593065, rel=1e-5)

    def test_b_const_value(self):
        a, b = bedmodel_ergun(0.1, 0.2, 20, 0.02, 0.8, 0.5, 0.1, 3000)
        assert b == pytest.approx(1450.195313, rel=1e-5)

    def test_cavity_larger_than_taphole(self):
        a, b = bedmodel_ergun(0.1, 0.05, 20, 0.02, 0.8, 0.5, 0.1, 3000)
        assert a > 0
        assert b > 0
