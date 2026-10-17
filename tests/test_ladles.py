"""Unit tests for tappingsim.ladles module."""

import math
import pytest
from scipy.constants import g
from tappingsim.ladles import (
    CylindricalLadle, overflow_flowrates, overflowmodel_exp, overflowmodel_step,
)


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
# Open-channel overflow flowrates
# ---------------------------------------------------------------------------

K = 2 / 3 * math.sqrt(g)
RHO_METAL, RHO_SLAG = 7000.0, 3000.0
CREST = 1.0
AREA = math.pi * 0.25 * 1.8**2


def fraction_model(dh, fraction):
    """Entrainment model with a constant metal fraction, for testing."""
    return fraction


class TestOverflowFlowrates:
    ratio = RHO_SLAG / RHO_METAL

    def flows(self, hmetal, hslag, model=overflowmodel_step, consts=()):
        return overflow_flowrates(
            hmetal, hslag, CREST, 0.4, 0.6, 0.55, self.ratio, lambda dh: model(dh, *consts))

    def test_no_flow_below_crest(self):
        assert self.flows(0.5, 0.99) == (0, 0)

    def test_no_flow_at_crest(self):
        assert self.flows(0.5, CREST) == (0, 0)

    def test_slag_only_flow(self):
        qmetal, qslag = self.flows(0.5, 1.1)
        assert qmetal == 0
        assert qslag == pytest.approx(K * 0.6 * 0.4 * 0.1**1.5)

    def test_metal_and_slag_flow(self):
        qmetal, qslag = self.flows(1.05, 1.25)
        assert qslag == pytest.approx(K * 0.6 * 0.4 * 0.2**1.5)
        expected = K * 0.55 * 0.4 * ((self.ratio * 0.2 + 0.05) ** 1.5 - (self.ratio * 0.2) ** 1.5)
        assert qmetal == pytest.approx(expected)

    def test_metal_only_flow(self):
        """No slag layer: metal flows with its own head above the crest."""
        qmetal, qslag = self.flows(1.1, 1.1)
        assert qslag == 0
        assert qmetal == pytest.approx(K * 0.55 * 0.4 * 0.1**1.5)

    def test_entrained_metal_in_slag_only_flow(self):
        qmetal, qslag = self.flows(0.5, 1.1, overflowmodel_exp, [0.2, 2.0])
        vfrac = 0.2 * math.exp(-2.0 * 0.5)
        total = K * 0.6 * 0.4 * 0.1**1.5
        assert qmetal == pytest.approx(vfrac * total)
        assert qslag == pytest.approx((1 - vfrac) * total)

    def test_entrained_metal_adds_to_metal_flow(self):
        plain = self.flows(1.05, 1.25)
        qmetal, qslag = self.flows(1.05, 1.25, overflowmodel_exp, [0.2, 2.0])
        qslag_total = K * 0.6 * 0.4 * 0.2**1.5
        assert qmetal == pytest.approx(plain[0] + 0.2 * qslag_total)
        assert qslag == pytest.approx(0.8 * qslag_total)

    def test_no_entrainment_without_metal(self):
        """There is no metal to entrain, even if the model returns a nonzero fraction."""
        qmetal, qslag = self.flows(0.0, 1.1, fraction_model, [0.3])
        assert qmetal == 0
        assert qslag == pytest.approx(K * 0.6 * 0.4 * 0.1**1.5)

    def test_total_flow_independent_of_entrainment(self):
        """Entrainment redistributes the slag stream between phases, not its total."""
        plain = self.flows(0.5, 1.1)
        entrained = self.flows(0.5, 1.1, fraction_model, [0.3])
        assert sum(entrained) == pytest.approx(sum(plain))


# ---------------------------------------------------------------------------
# CylindricalLadle
# ---------------------------------------------------------------------------

def make_ladle(depth=CREST, hmetal_init=0.0, hslag_init=0.0, overflowmodel=overflowmodel_step,
               overflowconsts=(), **overrides):
    params = dict(
        diameter=1.8, depth=depth, hmetal_init=hmetal_init, hslag_init=hslag_init,
        overflowmodel=overflowmodel, overflowconsts=list(overflowconsts),
        outlet_width=0.4, outlet_slag_Cd=0.6, outlet_metal_Cd=0.55,
        density_metal=RHO_METAL, density_slag=RHO_SLAG,
    )
    params.update(overrides)
    return CylindricalLadle(**params)


def stored_volumes(ladle):
    """Volumes of (metal, slag) held in the ladle, m3."""
    return AREA * ladle.hmetal, AREA * (ladle.hslag - ladle.hmetal)


@pytest.fixture
def empty_ladle():
    return make_ladle(depth=2.0)


class TestCylindricalLadleInit:
    def test_initial_metal_level(self, empty_ladle):
        assert empty_ladle.hmetal == 0.0

    def test_initial_slag_level(self, empty_ladle):
        assert empty_ladle.hslag == 0.0

    def test_diameter_stored(self, empty_ladle):
        assert empty_ladle.diameter == 1.8

    def test_depth_stored(self, empty_ladle):
        assert empty_ladle.depth == 2.0

    def test_outlet_parameters_stored(self, empty_ladle):
        assert empty_ladle.outlet_width == 0.4
        assert empty_ladle.outlet_slag_Cd == 0.6
        assert empty_ladle.outlet_metal_Cd == 0.55
        assert empty_ladle.density_metal == RHO_METAL
        assert empty_ladle.density_slag == RHO_SLAG

    def test_default_max_substep(self, empty_ladle):
        assert empty_ladle.max_substep == 1.0

    def test_initial_flowrates_zero(self, empty_ladle):
        """Flowrates are defined before the first time step."""
        assert empty_ladle.vdotmetal_out == 0
        assert empty_ladle.vdotslag_out == 0

    @pytest.mark.parametrize("width", [0, -0.1])
    def test_invalid_outlet_width_rejected(self, width):
        with pytest.raises(ValueError):
            make_ladle(outlet_width=width)

    @pytest.mark.parametrize("max_substep", [0, -1.0])
    def test_invalid_max_substep_rejected(self, max_substep):
        with pytest.raises(ValueError):
            make_ladle(max_substep=max_substep)

    @pytest.mark.parametrize("missing", ["outlet_width", "outlet_slag_Cd", "outlet_metal_Cd",
                                         "density_metal", "density_slag"])
    def test_outlet_parameters_are_required(self, missing):
        params = dict(
            diameter=1.8, depth=CREST, hmetal_init=0.0, hslag_init=0.0,
            overflowmodel=overflowmodel_step, overflowconsts=[],
            outlet_width=0.4, outlet_slag_Cd=0.6, outlet_metal_Cd=0.55,
            density_metal=RHO_METAL, density_slag=RHO_SLAG,
        )
        del params[missing]
        with pytest.raises(TypeError):
            CylindricalLadle(**params)


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

    def test_filling_below_crest_is_exact(self):
        """Without overflow the levels rise by inflow volume over cross-sectional area."""
        ladle = make_ladle()
        for _ in range(5):
            ladle.calc_dt(60, 0.001, 0.0005)
        assert ladle.hmetal == pytest.approx(300 * 0.001 / AREA, rel=1e-12)
        assert ladle.hslag == pytest.approx(300 * 0.0015 / AREA, rel=1e-12)

    def test_emptied_ladle_refills(self):
        ladle = make_ladle(hmetal_init=0.5, hslag_init=1.5)
        ladle.hmetal, ladle.hslag = 0.0, 0.0
        ladle.calc_dt(60, 0.005, 0.002)
        assert ladle.hmetal > 0


class TestCylindricalLadleOverflow:
    def test_outflow_when_full(self):
        """Ladle initialised above the crest should produce outflow immediately."""
        ladle = make_ladle(hslag_init=1.5)
        ladle.calc_dt(60, 0.0, 0.001)
        assert ladle.vdotslag_out > 0

    def test_step_model_no_metal_in_slag_overflow(self):
        """Step overflow model means no metal entrainment in slag overflow."""
        ladle = make_ladle(hmetal_init=0.5, hslag_init=1.5)
        ladle.calc_dt(60, 0.0, 0.001)
        assert ladle.vdotmetal_out == 0

    def test_levels_can_exceed_depth(self):
        """A head builds up over the outlet crest."""
        ladle = make_ladle()
        for _ in range(100):
            ladle.calc_dt(60, 0.0, 0.012)
        assert ladle.hslag > CREST + 0.01

    def test_steady_slag_head_matches_open_channel_flow(self):
        vdotslag_in = 0.012
        expected_head = (vdotslag_in / (K * 0.6 * 0.4)) ** (2 / 3)
        ladle = make_ladle()
        for _ in range(500):
            ladle.calc_dt(60, 0.0, vdotslag_in)
        assert ladle.hslag - CREST == pytest.approx(expected_head, rel=1e-6)
        assert ladle.vdotslag_out == pytest.approx(vdotslag_in, rel=1e-6)

    @pytest.mark.parametrize("dt", [1.0, 10.0, 60.0])
    def test_steady_state_independent_of_time_step(self, dt):
        ladle = make_ladle()
        for _ in range(int(30000 / dt)):
            ladle.calc_dt(dt, 0.0, 0.012)
        expected_head = (0.012 / (K * 0.6 * 0.4)) ** (2 / 3)
        assert ladle.hslag - CREST == pytest.approx(expected_head, rel=1e-6)

    def test_metal_and_slag_steady_state(self):
        """With both phases flowing in, both leave at the inflow rates."""
        ladle = make_ladle(overflowmodel=overflowmodel_exp, overflowconsts=[0.1, 2.0])
        for _ in range(100):
            ladle.calc_dt(60, 0.004, 0.006)
        assert ladle.vdotmetal_out == pytest.approx(0.004, rel=1e-6)
        assert ladle.vdotslag_out == pytest.approx(0.006, rel=1e-6)

    def test_mass_balance_over_step_with_entrainment(self):
        ladle = make_ladle(hmetal_init=1.05, hslag_init=1.2, overflowmodel=overflowmodel_exp,
                           overflowconsts=[0.1, 2.0])
        metal0, slag0 = stored_volumes(ladle)
        ladle.calc_dt(60, 0.02, 0.03)
        metal1, slag1 = stored_volumes(ladle)
        assert ladle.vdotmetal_out > 0 and ladle.vdotslag_out > 0
        assert metal1 - metal0 == pytest.approx((0.02 - ladle.vdotmetal_out) * 60, abs=1e-9)
        assert slag1 - slag0 == pytest.approx((0.03 - ladle.vdotslag_out) * 60, abs=1e-9)

    def test_entrainment_moves_metal_into_slag_outflow(self):
        step = make_ladle(hmetal_init=0.5, hslag_init=1.1)
        exp = make_ladle(hmetal_init=0.5, hslag_init=1.1, overflowmodel=overflowmodel_exp,
                         overflowconsts=[0.2, 2.0])
        step.calc_dt(1, 0.0, 0.0)
        exp.calc_dt(1, 0.0, 0.0)
        assert step.vdotmetal_out == 0
        assert exp.vdotmetal_out > 0
        assert exp.hmetal < step.hmetal

    def test_long_step_equals_repeated_substeps(self):
        long_step = make_ladle(hmetal_init=1.0, hslag_init=1.1)
        repeated = make_ladle(hmetal_init=1.0, hslag_init=1.1)
        long_step.calc_dt(10.0, 0.004, 0.006)
        for _ in range(10):
            repeated.calc_dt(1.0, 0.004, 0.006)
        assert (long_step.hmetal, long_step.hslag) == pytest.approx(
            (repeated.hmetal, repeated.hslag), abs=1e-12)

    def test_first_order_convergence(self):
        """Halving the sub-step size should roughly halve the error against a fine solution."""
        def state_after(max_substep, total_time=120):
            ladle = make_ladle(hmetal_init=0.9, hslag_init=1.0, max_substep=max_substep)
            for _ in range(int(round(total_time / max_substep))):
                ladle.calc_dt(max_substep, 0.004, 0.02)
            return ladle.hmetal, ladle.hslag

        reference = state_after(0.05)
        errors = [max(abs(a - b) for a, b in zip(state_after(dt), reference))
                  for dt in (8.0, 4.0, 2.0)]
        assert errors[0] / errors[1] == pytest.approx(2, abs=0.6)
        assert errors[1] / errors[2] == pytest.approx(2, abs=0.6)

    def test_very_wide_outlet_holds_level_at_crest(self):
        """With an outlet wide enough to pass any flow at a negligible head, the level sits
        at the crest and the inflow leaves immediately."""
        ladle = make_ladle(outlet_width=1e4, outlet_slag_Cd=1.0, outlet_metal_Cd=1.0)
        for _ in range(120):
            ladle.calc_dt(60, 0.004, 0.006)
        assert ladle.hslag == pytest.approx(CREST, abs=1e-3)
        assert ladle.hmetal == pytest.approx(CREST, abs=1e-3)
        assert ladle.vdotmetal_out == pytest.approx(0.004, rel=1e-3)
        assert ladle.vdotslag_out == pytest.approx(0.006, rel=1e-3)

    def test_unphysical_step_is_halved(self):
        """A step that would give negative volumes is repeated as smaller steps."""
        ladle = make_ladle(hmetal_init=0.5, hslag_init=1.1, max_substep=60.0)
        original = ladle._step
        dts = []

        def step(dt, vdotmetal_in, vdotslag_in):
            dts.append(dt)
            flows = original(dt, vdotmetal_in, vdotslag_in)
            if dt > 20:
                ladle.hslag = -1.0
            return flows

        ladle._step = step
        ladle.calc_dt(60.0, 0.0, 0.0)
        assert 60.0 not in dts[1:]
        assert max(dts) == 60.0 and min(dts) <= 15.0
        assert ladle.hslag >= ladle.hmetal >= 0

    def test_failure_to_find_physical_step_raises(self):
        ladle = make_ladle(hmetal_init=0.5, hslag_init=1.1)
        original = ladle._step

        def step(dt, vdotmetal_in, vdotslag_in):
            original(dt, vdotmetal_in, vdotslag_in)
            ladle.hslag = float("nan")

        ladle._step = step
        with pytest.raises(RuntimeError):
            ladle._advance(1.0, 0.0, 0.0, halvings=3)
