"""Unit tests for the tappingsim.ladles.PanLadle class."""

import functools
import math
import warnings

import pytest
from scipy.constants import g
from scipy.integrate import quad
from tappingsim.components.outlets import overflowmodel_exp, overflowmodel_step
from tappingsim.ladles import CylindricalLadle, PanLadle

K = 2 / 3 * math.sqrt(g)
RHO_METAL, RHO_SLAG = 7000.0, 3000.0
DEPTH = 0.9
BASE_L, BASE_W, TOP_L, TOP_W = 0.8, 1.0, 2.2, 1.6


def exp_model(premultiplier, decay):
    return functools.partial(overflowmodel_exp, premultiplier=premultiplier, decay=decay)


def make_pan(hmetal_init=0.0, hslag_init=0.0, overflowmodel=overflowmodel_step, **overrides):
    params = dict(
        base_length=BASE_L, base_width=BASE_W, top_length=TOP_L, top_width=TOP_W, depth=DEPTH,
        hmetal_init=hmetal_init, hslag_init=hslag_init, overflowmodel=overflowmodel,
        outlet_width=0.4, outlet_slag_Cd=0.6, outlet_metal_Cd=0.55,
        density_metal=RHO_METAL, density_slag=RHO_SLAG,
    )
    params.update(overrides)
    return PanLadle(**params)


def make_pan_quietly(**kwargs):
    """A pan which may give the warning for a base larger than the top."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return make_pan(**kwargs)


def section_area(z, bl=BASE_L, bw=BASE_W, tl=TOP_L, tw=TOP_W, depth=DEPTH):
    """The horizontal cross-sectional area at height z, with the taper continuing above the
    rim, m2."""
    return (bl + (tl - bl) * z / depth) * (bw + (tw - bw) * z / depth)


def stored_volumes(pan):
    """Volumes of (metal, slag) held in the pan, m3."""
    return pan._levels_to_volumes(pan.hmetal, pan.hslag)


@pytest.fixture
def empty_pan():
    return make_pan()


class TestPanLadleInit:
    def test_initial_levels(self):
        pan = make_pan(hmetal_init=0.1, hslag_init=0.3)
        assert pan.hmetal == 0.1
        assert pan.hslag == 0.3

    def test_geometry_stored(self, empty_pan):
        assert empty_pan.base_length == BASE_L
        assert empty_pan.base_width == BASE_W
        assert empty_pan.top_length == TOP_L
        assert empty_pan.top_width == TOP_W
        assert empty_pan.depth == DEPTH

    def test_outlet_parameters_stored(self, empty_pan):
        assert empty_pan.outlet.crest_height == DEPTH
        assert empty_pan.outlet.width == 0.4
        assert empty_pan.outlet.slag_Cd == 0.6
        assert empty_pan.outlet.metal_Cd == 0.55
        assert empty_pan.outlet.density_metal == RHO_METAL
        assert empty_pan.outlet.density_slag == RHO_SLAG
        assert empty_pan.outlet.overflowmodel is overflowmodel_step

    def test_default_max_substep(self, empty_pan):
        assert empty_pan.max_substep == 1.0

    def test_initial_flowrates_zero(self, empty_pan):
        assert empty_pan.vdotmetal_out == 0
        assert empty_pan.vdotslag_out == 0

    @pytest.mark.parametrize("name", ["base_length", "base_width", "top_length", "top_width",
                                      "depth"])
    @pytest.mark.parametrize("value", [0, -0.5])
    def test_invalid_dimensions_rejected(self, name, value):
        with pytest.raises(ValueError, match=name):
            make_pan(**{name: value})

    def test_invalid_outlet_width_rejected(self):
        with pytest.raises(ValueError):
            make_pan(outlet_width=0)

    def test_invalid_max_substep_rejected(self):
        with pytest.raises(ValueError):
            make_pan(max_substep=0)

    @pytest.mark.parametrize("missing", ["base_length", "top_width", "depth", "outlet_width",
                                         "density_slag"])
    def test_parameters_are_required(self, missing):
        params = dict(
            base_length=BASE_L, base_width=BASE_W, top_length=TOP_L, top_width=TOP_W,
            depth=DEPTH, hmetal_init=0.0, hslag_init=0.0, overflowmodel=overflowmodel_step,
            outlet_width=0.4, outlet_slag_Cd=0.6, outlet_metal_Cd=0.55,
            density_metal=RHO_METAL, density_slag=RHO_SLAG,
        )
        del params[missing]
        with pytest.raises(TypeError):
            PanLadle(**params)


class TestPanLadleBaseLargerThanTopWarning:
    def test_no_warning_for_a_pan_widening_upwards(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            make_pan()

    def test_no_warning_for_a_prism(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            make_pan(top_length=BASE_L, top_width=BASE_W)

    def test_no_warning_when_only_one_direction_is_equal(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            make_pan(top_length=BASE_L)

    @pytest.mark.parametrize("overrides", [
        dict(top_length=0.5),
        dict(top_width=0.5),
        dict(top_length=0.5, top_width=0.5),
    ])
    def test_warning_when_the_base_is_larger_in_any_direction(self, overrides):
        with pytest.warns(UserWarning, match="base of the pan is larger than the top"):
            make_pan(**overrides)

    def test_warning_points_at_the_caller(self):
        with pytest.warns(UserWarning) as record:
            make_pan(top_length=0.5)
        assert record[0].filename == __file__


class TestPanLadleGeometry:
    def test_volume_below_a_level_matches_the_integrated_cross_section(self, empty_pan):
        for h in (0.05, 0.4, DEPTH, 1.3):
            expected = quad(section_area, 0, h)[0]
            assert empty_pan._levels_to_volumes(h, h)[0] == pytest.approx(expected, rel=1e-12)

    def test_volume_to_the_rim_for_similar_rectangles(self):
        """The truncated pyramid formula, for a pan whose top is a scaled copy of its base."""
        pan = make_pan(base_length=1.0, base_width=2.0, top_length=2.0, top_width=4.0, depth=1.5)
        a0, a1 = 2.0, 8.0
        expected = 1.5 / 3 * (a0 + a1 + math.sqrt(a0 * a1))
        assert pan._levels_to_volumes(1.5, 1.5)[0] == pytest.approx(expected, rel=1e-14)

    def test_slag_volume_is_the_volume_between_the_levels(self, empty_pan):
        vmetal, vslag = empty_pan._levels_to_volumes(0.2, 0.7)
        assert vmetal == pytest.approx(quad(section_area, 0, 0.2)[0], rel=1e-12)
        assert vslag == pytest.approx(quad(section_area, 0.2, 0.7)[0], rel=1e-12)

    def test_prism_volume_is_linear(self):
        pan = make_pan(top_length=BASE_L, top_width=BASE_W)
        vmetal, vslag = pan._levels_to_volumes(0.3, 1.2)
        assert vmetal == pytest.approx(BASE_L * BASE_W * 0.3, rel=1e-14)
        assert vslag == pytest.approx(BASE_L * BASE_W * 0.9, rel=1e-14)

    @pytest.mark.parametrize("hmetal, hslag", [
        (0.0, 0.0), (1e-6, 1e-4), (0.1, 0.5), (0.5, 0.9), (0.85, 0.95), (1.2, 1.5), (3.0, 4.0),
    ])
    def test_levels_survive_conversion_to_volumes_and_back(self, empty_pan, hmetal, hslag):
        levels = empty_pan._volumes_to_levels(*empty_pan._levels_to_volumes(hmetal, hslag))
        assert levels == pytest.approx((hmetal, hslag), rel=1e-12, abs=1e-15)

    @pytest.mark.parametrize("overrides", [
        dict(),
        dict(top_length=BASE_L, top_width=BASE_W),
        dict(top_length=BASE_L),
        dict(top_length=3.0, top_width=BASE_W),
        dict(top_length=0.5, top_width=2.0),
        dict(top_length=0.5, top_width=0.4),
        dict(base_length=0.02, top_length=3, base_width=0.02, top_width=3),
    ])
    def test_level_inversion_for_any_taper(self, overrides):
        pan = make_pan_quietly(**overrides)
        for h in (1e-5, 0.01, 0.3, DEPTH, 1.1):
            level = pan._volumes_to_levels(pan._levels_to_volumes(h, h)[0], 0.0)[0]
            assert level == pytest.approx(h, rel=1e-12)

    def test_level_is_the_inverse_of_the_volume_for_a_pan_that_closes(self):
        """A pan narrowing to a point at 3 m has a limited volume."""
        pan = make_pan_quietly(base_length=1.0, base_width=1.0, top_length=0.5, top_width=0.5,
                               depth=1.0)
        assert pan._volumes_to_levels(pan._levels_to_volumes(0.9, 0.9)[0], 0.0)[0] == (
            pytest.approx(0.9, rel=1e-12))
        vmax = pan._volume(2.0)
        assert pan._level(vmax) == pytest.approx(2.0, rel=1e-6)
        with pytest.raises(RuntimeError, match="narrows to a close"):
            pan._level(1.01 * vmax)

    def test_level_of_a_slightly_negative_volume(self, empty_pan):
        """A level, and not an error, is given for rounding errors in an empty pan."""
        assert empty_pan._level(-1e-18) == pytest.approx(-1e-18 / (BASE_L * BASE_W))
        assert empty_pan._level(0.0) == 0.0

    def test_cross_section_continues_to_taper_above_the_rim(self, empty_pan):
        area_above = (empty_pan._volume(1.5 + 1e-6) - empty_pan._volume(1.5)) / 1e-6
        assert area_above == pytest.approx(section_area(1.5), rel=1e-5)
        assert section_area(1.5) > section_area(DEPTH) > BASE_L * BASE_W


class TestPanLadleFilling:
    def test_levels_increase(self, empty_pan):
        empty_pan.calc_dt(60, 0.005, 0.002)
        assert empty_pan.hmetal > 0
        assert empty_pan.hslag > empty_pan.hmetal

    def test_no_outflow_when_not_full(self, empty_pan):
        empty_pan.calc_dt(60, 0.005, 0.002)
        assert empty_pan.vdotmetal_out == 0
        assert empty_pan.vdotslag_out == 0

    def test_filling_adds_the_volumes_entering(self, empty_pan):
        empty_pan.calc_dt(30, 0.004, 0.006)
        empty_pan.calc_dt(30, 0.002, 0.001)
        vmetal, vslag = stored_volumes(empty_pan)
        assert vmetal == pytest.approx(30 * 0.006, rel=1e-12)
        assert vslag == pytest.approx(30 * 0.007, rel=1e-12)

    def test_filling_below_the_crest_is_exact(self):
        """One long step, or many short ones, give the same levels."""
        one = make_pan(hmetal_init=0.1, hslag_init=0.2)
        many = make_pan(hmetal_init=0.1, hslag_init=0.2)
        one.calc_dt(100, 0.003, 0.002)
        for _ in range(100):
            many.calc_dt(1, 0.003, 0.002)
        assert (one.hmetal, one.hslag) == pytest.approx((many.hmetal, many.hslag), rel=1e-12)

    def test_levels_rise_more_slowly_as_the_pan_widens(self, empty_pan):
        heights = []
        for _ in range(4):
            empty_pan.calc_dt(60, 0.0, 0.003)
            heights.append(empty_pan.hslag)
        rises = [b - a for a, b in zip([0] + heights[:-1], heights)]
        assert rises[0] > rises[1] > rises[2] > rises[3] > 0

    def test_emptied_pan_refills(self, empty_pan):
        empty_pan.calc_dt(60, 0.005, 0.002)
        empty_pan.hmetal, empty_pan.hslag = 0, 0
        empty_pan.calc_dt(60, 0.005, 0.002)
        assert empty_pan.hslag > empty_pan.hmetal > 0

    def test_a_step_filling_beyond_the_crest_overflows(self):
        pan = make_pan(hslag_init=0.85)
        pan.calc_dt(60, 0.0, 0.05)
        assert pan.hslag > DEPTH
        assert pan.vdotslag_out > 0


class TestPanLadleOverflow:
    def test_outflow_when_full(self):
        pan = make_pan(hslag_init=1.2)
        pan.calc_dt(60, 0.0, 0.001)
        assert pan.vdotslag_out > 0

    def test_step_model_no_metal_in_slag_overflow(self):
        pan = make_pan(hmetal_init=0.5, hslag_init=1.2)
        pan.calc_dt(60, 0.0, 0.001)
        assert pan.vdotmetal_out == 0

    def test_steady_slag_head_matches_open_channel_flow(self):
        vdotslag_in = 0.012
        expected_head = (vdotslag_in / (K * 0.6 * 0.4)) ** (2 / 3)
        pan = make_pan()
        for _ in range(500):
            pan.calc_dt(60, 0.0, vdotslag_in)
        assert pan.hslag - DEPTH == pytest.approx(expected_head, rel=1e-6)
        assert pan.vdotslag_out == pytest.approx(vdotslag_in, rel=1e-6)

    def test_metal_and_slag_steady_state(self):
        pan = make_pan(overflowmodel=exp_model(0.1, 2.0))
        for _ in range(100):
            pan.calc_dt(60, 0.004, 0.006)
        assert pan.vdotmetal_out == pytest.approx(0.004, rel=1e-6)
        assert pan.vdotslag_out == pytest.approx(0.006, rel=1e-6)

    def test_volume_balance_over_step_with_entrainment(self):
        pan = make_pan(hmetal_init=0.95, hslag_init=1.1, overflowmodel=exp_model(0.1, 2.0))
        metal0, slag0 = stored_volumes(pan)
        pan.calc_dt(60, 0.02, 0.03)
        metal1, slag1 = stored_volumes(pan)
        assert pan.vdotmetal_out > 0 and pan.vdotslag_out > 0
        assert metal1 - metal0 == pytest.approx((0.02 - pan.vdotmetal_out) * 60, abs=1e-9)
        assert slag1 - slag0 == pytest.approx((0.03 - pan.vdotslag_out) * 60, abs=1e-9)

    def test_entrainment_moves_metal_into_slag_outflow(self):
        step = make_pan(hmetal_init=0.5, hslag_init=1.0)
        exp = make_pan(hmetal_init=0.5, hslag_init=1.0, overflowmodel=exp_model(0.2, 2.0))
        step.calc_dt(1, 0.0, 0.0)
        exp.calc_dt(1, 0.0, 0.0)
        assert step.vdotmetal_out == 0
        assert exp.vdotmetal_out > 0
        assert exp.hmetal < step.hmetal

    def test_very_wide_outlet_holds_level_at_crest(self):
        pan = make_pan(outlet_width=1e4, outlet_slag_Cd=1.0, outlet_metal_Cd=1.0)
        for _ in range(120):
            pan.calc_dt(60, 0.004, 0.006)
        assert pan.hslag == pytest.approx(DEPTH, abs=1e-3)
        assert pan.vdotmetal_out == pytest.approx(0.004, rel=1e-3)
        assert pan.vdotslag_out == pytest.approx(0.006, rel=1e-3)

    def test_wider_pan_holds_more_for_the_same_head(self):
        """The same steady head is reached at the same flow in any pan, but a pan with a larger
        top takes longer to reach it."""
        small = make_pan(top_length=1.0, top_width=1.0)
        large = make_pan(top_length=2.5, top_width=2.5)
        for pan in (small, large):
            for _ in range(10):
                pan.calc_dt(60, 0.0, 0.01)
        assert small.vdotslag_out > large.vdotslag_out


class TestPanLadleAgainstReferences:
    # Levels after 25 s, from an independent solution of the volume balance with the exact
    # level-volume relationship (SciPy's DOP853 integrator at a relative tolerance of 1e-12).
    REFERENCE_LEVELS = (0.9955064630686271, 1.5655938711044397)

    @staticmethod
    def state_after(max_substep):
        pan = make_pan(hmetal_init=0.1, hslag_init=0.25, overflowmodel=exp_model(0.3, 8.0),
                       max_substep=max_substep)
        for _ in range(int(round(25 / max_substep))):
            pan.calc_dt(max_substep, 0.12, 0.2)
        return pan.hmetal, pan.hslag

    def test_converges_to_the_reference_solution(self):
        errors = [max(abs(a - b) for a, b in zip(self.state_after(dt), self.REFERENCE_LEVELS))
                  for dt in (1.0, 0.5, 0.25)]
        assert errors[0] < 1e-2
        assert errors[1] / errors[2] == pytest.approx(2, abs=0.6)
        assert errors[0] / errors[1] == pytest.approx(2, abs=0.6)
        assert errors[2] < 2.5e-3

    def test_prism_pan_is_the_same_as_a_cylinder_of_equal_area(self):
        """A prism has the same level-volume relationship as a cylinder with the same
        cross-sectional area."""
        area = 0.25 * math.pi * 1.8**2
        pan = make_pan(base_length=1.5, base_width=area / 1.5, top_length=1.5,
                       top_width=area / 1.5, depth=1.0, overflowmodel=exp_model(0.3, 8.0))
        cylinder = CylindricalLadle(
            diameter=1.8, depth=1.0, hmetal_init=0, hslag_init=0,
            overflowmodel=exp_model(0.3, 8.0), outlet_width=0.4, outlet_slag_Cd=0.6,
            outlet_metal_Cd=0.55, density_metal=RHO_METAL, density_slag=RHO_SLAG)
        for step in range(300):
            vdotmetal, vdotslag = (0.02, 0.03) if step < 200 else (0.0, 0.0)
            pan.calc_dt(10, vdotmetal, vdotslag)
            cylinder.calc_dt(10, vdotmetal, vdotslag)
            assert pan.hmetal == pytest.approx(cylinder.hmetal, abs=1e-9)
            assert pan.hslag == pytest.approx(cylinder.hslag, abs=1e-9)
            assert pan.vdotmetal_out == pytest.approx(cylinder.vdotmetal_out, abs=1e-9)
            assert pan.vdotslag_out == pytest.approx(cylinder.vdotslag_out, abs=1e-9)
        assert pan.vdotslag_out > 0

    def test_volumes_are_conserved_over_a_tapping(self):
        pan = make_pan(hmetal_init=0.1, hslag_init=0.25, overflowmodel=exp_model(0.3, 8.0))
        vmetal0, vslag0 = stored_volumes(pan)
        vmetal_out = vslag_out = 0
        for _ in range(400):
            pan.calc_dt(1, 0.12, 0.2)
            vmetal_out += pan.vdotmetal_out
            vslag_out += pan.vdotslag_out
        vmetal1, vslag1 = stored_volumes(pan)
        assert vmetal_out > 0 and vslag_out > 0
        assert vmetal1 - vmetal0 == pytest.approx(400 * 0.12 - vmetal_out, abs=1e-10)
        assert vslag1 - vslag0 == pytest.approx(400 * 0.2 - vslag_out, abs=1e-10)

    def test_a_pan_that_narrows_upwards_behaves_the_same_way(self):
        """The warning aside, an inverted pan is a valid ladle."""
        pan = make_pan_quietly(base_length=1.5, base_width=1.5, top_length=1.0, top_width=1.0)
        for _ in range(500):
            pan.calc_dt(60, 0.0, 0.012)
        expected_head = (0.012 / (K * 0.6 * 0.4)) ** (2 / 3)
        assert pan.hslag - DEPTH == pytest.approx(expected_head, rel=1e-6)
