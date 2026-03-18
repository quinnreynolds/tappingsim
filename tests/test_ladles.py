"""Unit tests for tappingsim.ladles module."""

import math
import pytest
from tappingsim.ladles import CylindricalLadle, overflowmodel_exp, overflowmodel_step


# ---------------------------------------------------------------------------
# Overflow models
# ---------------------------------------------------------------------------

class TestOverflowmodelStep:
    def test_returns_zero(self):
        assert overflowmodel_step(0.5) == 0

    def test_returns_zero_negative_delta(self):
        assert overflowmodel_step(-0.5) == 0

    def test_ignores_consts(self):
        assert overflowmodel_step(0.5, 1.0, 2.0) == 0


class TestOverflowmodelExp:
    def test_positive_delta_decays(self):
        """Entrainment should decay exponentially when interface is below outlet."""
        f = overflowmodel_exp(0.5, 0.1, 2.0)
        assert f == pytest.approx(0.036787944117144235, rel=1e-8)

    def test_negative_delta_returns_premultiplier(self):
        """When interface is above outlet, return full pre-multiplier."""
        f = overflowmodel_exp(-0.5, 0.1, 2.0)
        assert f == pytest.approx(0.1, rel=1e-8)

    def test_zero_delta_returns_premultiplier(self):
        f = overflowmodel_exp(0.0, 0.1, 2.0)
        assert f == pytest.approx(0.1, rel=1e-8)

    def test_large_positive_delta_approaches_zero(self):
        f = overflowmodel_exp(100.0, 0.1, 2.0)
        assert f == pytest.approx(0.0, abs=1e-80)


# ---------------------------------------------------------------------------
# CylindricalLadle
# ---------------------------------------------------------------------------

@pytest.fixture
def empty_ladle():
    return CylindricalLadle(
        diameter=1.8,
        depth=2.0,
        hmetal_init=0.0,
        hslag_init=0.0,
        overflowmodel=overflowmodel_step,
        overflowconsts=[],
    )


class TestCylindricalLadleInit:
    def test_initial_metal_level(self, empty_ladle):
        assert empty_ladle.hmetal == 0.0

    def test_initial_slag_level(self, empty_ladle):
        assert empty_ladle.hslag == 0.0

    def test_diameter_stored(self, empty_ladle):
        assert empty_ladle.diameter == 1.8

    def test_depth_stored(self, empty_ladle):
        assert empty_ladle.depth == 2.0


class TestCylindricalLadleFilling:
    def test_metal_level_increases(self, empty_ladle):
        empty_ladle.calc_dt(60, 0.005, 0.002)
        assert empty_ladle.hmetal > 0

    def test_slag_level_increases(self, empty_ladle):
        empty_ladle.calc_dt(60, 0.005, 0.002)
        assert empty_ladle.hslag > 0

    def test_slag_above_metal(self, empty_ladle):
        empty_ladle.calc_dt(60, 0.005, 0.002)
        assert empty_ladle.hslag >= empty_ladle.hmetal

    def test_no_outflow_when_not_full(self, empty_ladle):
        empty_ladle.calc_dt(60, 0.005, 0.002)
        assert empty_ladle.vdotmetal_out == 0
        assert empty_ladle.vdotslag_out == 0

    def test_metal_level_value(self, empty_ladle):
        empty_ladle.calc_dt(60, 0.005, 0.002)
        assert empty_ladle.hmetal == pytest.approx(0.11789255043844098, rel=1e-8)

    def test_slag_level_value(self, empty_ladle):
        empty_ladle.calc_dt(60, 0.005, 0.002)
        assert empty_ladle.hslag == pytest.approx(0.16504957061381736, rel=1e-8)


class TestCylindricalLadleOverflow:
    def test_outflow_when_full(self):
        """Ladle initialised above depth should produce outflow immediately."""
        ladle = CylindricalLadle(
            diameter=1.8,
            depth=1.0,
            hmetal_init=0.0,
            hslag_init=1.5,
            overflowmodel=overflowmodel_step,
            overflowconsts=[],
        )
        ladle.calc_dt(60, 0.0, 0.001)
        assert ladle.vdotslag_out > 0

    def test_step_model_no_metal_in_slag_overflow(self):
        """Step overflow model means no metal entrainment in slag overflow."""
        ladle = CylindricalLadle(
            diameter=1.8,
            depth=1.0,
            hmetal_init=0.5,
            hslag_init=1.5,
            overflowmodel=overflowmodel_step,
            overflowconsts=[],
        )
        ladle.calc_dt(60, 0.0, 0.001)
        assert ladle.vdotmetal_out == 0
