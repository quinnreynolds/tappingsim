"""Unit tests for tappingsim.components.fdmodels module."""

import math

import pytest
from tappingsim.components.fdmodels import (
    fdmodel_bellos,
    fdmodel_cheng,
    fdmodel_serghides1,
    fdmodel_serghides2,
)


class TestFdmodelBellos:
    def test_positive(self):
        assert fdmodel_bellos(1.0, 7000, 0.005, 0.1, 1e-3) > 0

    def test_value(self):
        f = fdmodel_bellos(1.0, 7000, 0.005, 0.1, 1e-3)
        assert f == pytest.approx(0.037441389, rel=1e-5)

    def test_low_velocity_clamp(self):
        """Near-zero velocity should not raise and should return a finite value."""
        f = fdmodel_bellos(1e-10, 7000, 0.005, 0.1, 1e-3)
        assert math.isfinite(f)
        assert f > 0


class TestFdmodelCheng:
    def test_positive(self):
        assert fdmodel_cheng(1.0, 7000, 0.005, 0.1, 1e-3) > 0

    def test_value(self):
        f = fdmodel_cheng(1.0, 7000, 0.005, 0.1, 1e-3)
        assert f == pytest.approx(0.037501995, rel=1e-5)


class TestFdmodelSerghides1:
    def test_positive(self):
        assert fdmodel_serghides1(1.0, 7000, 0.005, 0.1, 1e-3) > 0

    def test_value(self):
        f = fdmodel_serghides1(1.0, 7000, 0.005, 0.1, 1e-3)
        assert f == pytest.approx(0.038334483, rel=1e-5)

    def test_zero_velocity_returns_laminar(self):
        assert fdmodel_serghides1(0.0, 7000, 0.005, 0.1, 1e-3) == 64

    def test_laminar_regime(self):
        """Very low Re should return the laminar value f = 64/Re."""
        f = fdmodel_serghides1(1e-4, 1000, 0.001, 0.01, 1e-5)
        assert f > 0


class TestFdmodelSerghides2:
    def test_positive(self):
        assert fdmodel_serghides2(1.0, 7000, 0.005, 0.1, 1e-3) > 0

    def test_value(self):
        f = fdmodel_serghides2(1.0, 7000, 0.005, 0.1, 1e-3)
        assert f == pytest.approx(0.038334483, rel=1e-5)

    def test_zero_velocity_returns_laminar(self):
        assert fdmodel_serghides2(0.0, 7000, 0.005, 0.1, 1e-3) == 64
