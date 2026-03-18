"""Unit tests for tappingsim.launders module."""

import pytest
from tappingsim.launders import SimpleSisoLaunder


class TestSimpleSisoLaunder:
    @pytest.fixture
    def launder(self):
        return SimpleSisoLaunder()

    def test_initial_metal_flowrate_zero(self, launder):
        assert launder.vdotmetal_out == 0

    def test_initial_slag_flowrate_zero(self, launder):
        assert launder.vdotslag_out == 0

    def test_metal_passthrough(self, launder):
        launder.calc_dt(10, 0.5, 0.3)
        assert launder.vdotmetal_out == pytest.approx(0.5)

    def test_slag_passthrough(self, launder):
        launder.calc_dt(10, 0.5, 0.3)
        assert launder.vdotslag_out == pytest.approx(0.3)

    def test_zero_inflow(self, launder):
        launder.calc_dt(10, 0.0, 0.0)
        assert launder.vdotmetal_out == 0
        assert launder.vdotslag_out == 0

    def test_dt_does_not_affect_output(self, launder):
        """Launder has no accumulation so dt should not affect output flowrates."""
        launder.calc_dt(1, 0.5, 0.3)
        assert launder.vdotmetal_out == pytest.approx(0.5)
        launder.calc_dt(1000, 0.5, 0.3)
        assert launder.vdotmetal_out == pytest.approx(0.5)
