"""Unit tests for tappingsim.launders module."""

import functools

import numpy
import pytest
from scipy.constants import g
from tappingsim.components.outlets import overflowmodel_exp
from tappingsim.launders import SimpleSisoLaunder, SkimmerLaunder


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


# ---------------------------------------------------------------------------
# SkimmerLaunder
# ---------------------------------------------------------------------------

RHO_METAL = 7000.0
RHO_SLAG = 3000.0


def make_skimmer(**overrides):
    """Build a skimmer with distinct trapezoidal troughs (outlets are wide enough
    that a steady state is reached below the trough rims)."""
    params = dict(
        slag_trough_length=2.0,
        slag_trough_upper_width=1.2,
        slag_trough_lower_width=0.8,
        slag_trough_height=0.8,
        slag_trough_outlet_height=0.5,
        slag_trough_outlet_width=0.5,
        weir_underflow_height=0.3,
        weir_underflow_width=0.3,
        metal_trough_length=1.5,
        metal_trough_upper_width=1.0,
        metal_trough_lower_width=0.6,
        metal_trough_height=0.7,
        metal_trough_outlet_height=0.4,
        metal_trough_outlet_width=0.5,
        slag_trough_outlet_slag_Cd=0.6,
        slag_trough_outlet_metal_Cd=0.55,
        metal_trough_outlet_slag_Cd=0.65,
        metal_trough_outlet_metal_Cd=0.5,
        density_metal=RHO_METAL,
        density_slag=RHO_SLAG,
    )
    params.update(overrides)
    return SkimmerLaunder(**params)


def set_levels(launder, ss_hmetal, ss_hslag, ms_hmetal, ms_hslag):
    launder.slag_trough_hmetal = ss_hmetal
    launder.slag_trough_hslag = ss_hslag
    launder.metal_trough_hmetal = ms_hmetal
    launder.metal_trough_hslag = ms_hslag


def calc_outlet_flows(launder):
    """Calculate the outlet flows for the current levels, and store them in the outflow
    attributes."""
    (
        launder.slag_trough_vdotmetal_out,
        launder.slag_trough_vdotslag_out,
        launder.metal_trough_vdotmetal_out,
        launder.metal_trough_vdotslag_out,
    ) = launder._calc_vdot_out(*levels(launder))


def total_volumes(launder):
    """Total (metal, slag) volume held in both troughs, m3."""
    vm_ss, vs_ss, vm_ms, vs_ms = launder._levels_to_volumes(
        launder.slag_trough_hmetal,
        launder.slag_trough_hslag,
        launder.metal_trough_hmetal,
        launder.metal_trough_hslag,
    )
    return vm_ss + vm_ms, vs_ss + vs_ms


@pytest.fixture
def skimmer():
    return make_skimmer()


class TestSkimmerLaunderInit:
    def test_initial_levels_zero(self, skimmer):
        assert skimmer.slag_trough_hmetal == 0
        assert skimmer.slag_trough_hslag == 0
        assert skimmer.metal_trough_hmetal == 0
        assert skimmer.metal_trough_hslag == 0

    def test_initial_outflows_zero(self, skimmer):
        assert skimmer.slag_trough_vdotmetal_out == 0
        assert skimmer.slag_trough_vdotslag_out == 0
        assert skimmer.metal_trough_vdotmetal_out == 0
        assert skimmer.metal_trough_vdotslag_out == 0

    def test_parameters_stored(self, skimmer):
        assert skimmer.slag_trough_length == 2.0
        assert skimmer.weir_underflow_height == 0.3
        assert skimmer.metal_trough_outlet.metal_Cd == 0.5
        assert skimmer.slag_trough_outlet.crest_height == 0.5
        assert skimmer.slag_trough_outlet.width == 0.5
        assert skimmer.density_slag == RHO_SLAG


class TestSkimmerLaunderGeometry:
    def test_rectangular_volumes(self):
        launder = make_skimmer(
            slag_trough_upper_width=0.8, slag_trough_lower_width=0.8,
            metal_trough_upper_width=0.6, metal_trough_lower_width=0.6,
        )
        vm_ss, vs_ss, vm_ms, vs_ms = launder._levels_to_volumes(0.1, 0.25, 0.2, 0.3)
        assert vm_ss == pytest.approx(2.0 * 0.8 * 0.1)
        assert vs_ss == pytest.approx(2.0 * 0.8 * 0.15)
        assert vm_ms == pytest.approx(1.5 * 0.6 * 0.2)
        assert vs_ms == pytest.approx(1.5 * 0.6 * 0.1)

    def test_full_trapezoidal_trough_volume(self, skimmer):
        """A trough full of metal holds length*height*mean width."""
        vm_ss, vs_ss, vm_ms, vs_ms = skimmer._levels_to_volumes(0.8, 0.8, 0.7, 0.7)
        assert vm_ss == pytest.approx(2.0 * 0.8 * (1.2 + 0.8) / 2)
        assert vs_ss == pytest.approx(0)
        assert vm_ms == pytest.approx(1.5 * 0.7 * (1.0 + 0.6) / 2)
        assert vs_ms == pytest.approx(0)

    def test_volume_between_levels_matches_integral(self, skimmer):
        """Slag layer volume equals the integral of trough width over its thickness."""
        _, vs_ss, _, _ = skimmer._levels_to_volumes(0.1, 0.4, 0.0, 0.0)
        h = numpy.linspace(0.1, 0.4, 100001)
        width = 0.8 + (1.2 - 0.8) * h / 0.8
        expected = 2.0 * numpy.trapezoid(width, h)
        assert vs_ss == pytest.approx(expected, rel=1e-8)

    @pytest.mark.parametrize("upper_width", [0.8, 1.2])
    def test_levels_volumes_round_trip(self, upper_width):
        launder = make_skimmer(
            slag_trough_upper_width=upper_width, slag_trough_lower_width=0.8,
            metal_trough_upper_width=upper_width / 2, metal_trough_lower_width=upper_width / 2,
        )
        levels = (0.12, 0.45, 0.2, 0.33)
        volumes = launder._levels_to_volumes(*levels)
        assert launder._volumes_to_levels(*volumes) == pytest.approx(levels)

    def test_empty_troughs_zero_levels(self, skimmer):
        assert skimmer._volumes_to_levels(0, 0, 0, 0) == pytest.approx((0, 0, 0, 0))


class TestSkimmerLaunderOutlets:
    def test_no_flow_below_crest(self, skimmer):
        set_levels(skimmer, 0.2, 0.45, 0.1, 0.35)
        calc_outlet_flows(skimmer)
        assert skimmer.slag_trough_vdotmetal_out == 0
        assert skimmer.slag_trough_vdotslag_out == 0
        assert skimmer.metal_trough_vdotmetal_out == 0
        assert skimmer.metal_trough_vdotslag_out == 0

    def test_slag_only_overflow(self, skimmer):
        """Slag above the crest with metal below it: broad-crested weir flow of slag."""
        set_levels(skimmer, 0.2, 0.6, 0.1, 0.55)
        calc_outlet_flows(skimmer)
        k = 2 / 3 * numpy.sqrt(g)
        assert skimmer.slag_trough_vdotslag_out == pytest.approx(k * 0.6 * 0.5 * 0.1**1.5)
        assert skimmer.slag_trough_vdotmetal_out == 0
        assert skimmer.metal_trough_vdotslag_out == pytest.approx(k * 0.65 * 0.5 * 0.15**1.5)
        assert skimmer.metal_trough_vdotmetal_out == 0

    def test_metal_and_slag_overflow(self, skimmer):
        """Metal above the crest: slag flows on its layer thickness, metal on the
        slag-equivalent head plus its own head above the crest."""
        set_levels(skimmer, 0.55, 0.7, 0.5, 0.6)
        calc_outlet_flows(skimmer)
        k = 2 / 3 * numpy.sqrt(g)
        ratio = RHO_SLAG / RHO_METAL
        # slag trough: crest 0.5, h1 = 0.15, h2 = 0.05
        assert skimmer.slag_trough_vdotslag_out == pytest.approx(k * 0.6 * 0.5 * 0.15**1.5)
        expected_metal = k * 0.55 * 0.5 * ((ratio * 0.15 + 0.05) ** 1.5 - (ratio * 0.15) ** 1.5)
        assert skimmer.slag_trough_vdotmetal_out == pytest.approx(expected_metal)
        # metal trough: crest 0.4, h1 = 0.1, h2 = 0.1
        assert skimmer.metal_trough_vdotslag_out == pytest.approx(k * 0.65 * 0.5 * 0.1**1.5)
        expected_metal = k * 0.5 * 0.5 * ((ratio * 0.1 + 0.1) ** 1.5 - (ratio * 0.1) ** 1.5)
        assert skimmer.metal_trough_vdotmetal_out == pytest.approx(expected_metal)

    def test_flow_resets_when_level_drops(self, skimmer):
        set_levels(skimmer, 0.2, 0.6, 0.1, 0.55)
        calc_outlet_flows(skimmer)
        set_levels(skimmer, 0.2, 0.45, 0.1, 0.35)
        calc_outlet_flows(skimmer)
        assert skimmer.slag_trough_vdotslag_out == 0
        assert skimmer.metal_trough_vdotslag_out == 0


class TestSkimmerLaunderCalcDt:
    def test_empty_launder_no_inflow_stays_empty(self, skimmer):
        skimmer.calc_dt(1.0, 0.0, 0.0)
        assert total_volumes(skimmer) == pytest.approx((0, 0))

    def test_inflow_fills_slag_trough_first(self, skimmer):
        skimmer.calc_dt(1.0, 0.01, 0.02)
        assert skimmer.slag_trough_hmetal > 0
        assert skimmer.slag_trough_hslag > skimmer.slag_trough_hmetal
        # material starts to pass the weir at once, but the metal trough fills much more slowly
        assert 0 <= skimmer.metal_trough_hmetal < 0.01 * skimmer.slag_trough_hmetal
        assert 0 <= skimmer.metal_trough_hslag < 0.01 * skimmer.slag_trough_hslag

    def test_volume_added_matches_inflow(self, skimmer):
        skimmer.calc_dt(2.0, 0.01, 0.02)
        vmetal, vslag = total_volumes(skimmer)
        assert vmetal == pytest.approx(0.02)
        assert vslag == pytest.approx(0.04)

    def test_mass_balance_with_transfer_and_overflow(self, skimmer):
        """Change in stored volume of each phase equals inflow minus outflow."""
        set_levels(skimmer, 0.55, 0.7, 0.3, 0.45)
        dt, vdotmetal_in, vdotslag_in = 0.5, 0.03, 0.05
        vmetal0, vslag0 = total_volumes(skimmer)
        skimmer.calc_dt(dt, vdotmetal_in, vdotslag_in)
        vmetal1, vslag1 = total_volumes(skimmer)
        metal_out = skimmer.slag_trough_vdotmetal_out + skimmer.metal_trough_vdotmetal_out
        slag_out = skimmer.slag_trough_vdotslag_out + skimmer.metal_trough_vdotslag_out
        assert metal_out > 0 and slag_out > 0
        assert vmetal1 - vmetal0 == pytest.approx((vdotmetal_in - metal_out) * dt)
        assert vslag1 - vslag0 == pytest.approx((vdotslag_in - slag_out) * dt)

    def test_higher_slag_trough_transfers_forward(self, skimmer):
        """Positive flow direction is slag trough to metal trough."""
        set_levels(skimmer, 0.25, 0.28, 0.0, 0.0)
        skimmer.calc_dt(0.1, 0.0, 0.0)
        assert skimmer.metal_trough_hmetal > 0
        assert skimmer.slag_trough_hmetal < 0.25

    def test_equal_levels_no_transfer(self, skimmer):
        """Identical phase levels in both troughs give no pressure difference, so
        nothing moves through the weir underflow."""
        set_levels(skimmer, 0.2, 0.28, 0.2, 0.28)
        before = (0.2, 0.28, 0.2, 0.28)
        skimmer.calc_dt(1.0, 0.0, 0.0)
        # equal levels, but troughs differ in cross-section: volumes must not change
        assert (
            skimmer.slag_trough_hmetal, skimmer.slag_trough_hslag,
            skimmer.metal_trough_hmetal, skimmer.metal_trough_hslag,
        ) == pytest.approx(before)

    def test_metal_trough_slag_layer_contributes_to_pressure(self, skimmer):
        """Regression test: the slag layer in the metal trough must add hydrostatic
        pressure, and so reduce the driving force across the weir. It was previously
        computed with a thickness of (hmetal - hmetal), i.e. always zero."""
        def metal_transferred(ms_hslag):
            launder = make_skimmer()
            set_levels(launder, 0.25, 0.28, 0.1, ms_hslag)
            launder.calc_dt(0.1, 0.0, 0.0)
            return launder.metal_trough_hmetal - 0.1

        assert metal_transferred(0.1) > 0
        assert metal_transferred(0.2) < metal_transferred(0.1)

    def test_steady_state_outflow_matches_inflow(self, skimmer):
        """Constant inflow eventually leaves entirely through the two overflows."""
        vdotmetal_in, vdotslag_in = 0.01, 0.02
        for _ in range(30000):
            skimmer.calc_dt(0.1, vdotmetal_in, vdotslag_in)
        metal_out = skimmer.slag_trough_vdotmetal_out + skimmer.metal_trough_vdotmetal_out
        slag_out = skimmer.slag_trough_vdotslag_out + skimmer.metal_trough_vdotslag_out
        assert metal_out == pytest.approx(vdotmetal_in, rel=1e-3)
        assert slag_out == pytest.approx(vdotslag_in, rel=1e-3)
        assert skimmer.slag_trough_hslag < skimmer.slag_trough_height
        assert skimmer.metal_trough_hslag < skimmer.metal_trough_height


def weir_pressure_difference(launder):
    """Hydrostatic pressure at the floor of the slag trough less that of the metal trough, Pa."""
    def pressure(hmetal, hslag):
        return (hmetal * RHO_METAL + (hslag - hmetal) * RHO_SLAG) * g

    return (
        pressure(launder.slag_trough_hmetal, launder.slag_trough_hslag)
        - pressure(launder.metal_trough_hmetal, launder.metal_trough_hslag)
    )


def levels(launder):
    return numpy.array(
        [
            launder.slag_trough_hmetal,
            launder.slag_trough_hslag,
            launder.metal_trough_hmetal,
            launder.metal_trough_hslag,
        ]
    )


class TestSkimmerLaunderTimeIntegration:
    def test_default_max_substep(self, skimmer):
        assert skimmer.max_substep == 1.0

    @pytest.mark.parametrize("max_substep", [0, -1.0])
    def test_invalid_max_substep_rejected(self, max_substep):
        with pytest.raises(ValueError):
            make_skimmer(max_substep=max_substep)

    def test_long_step_equals_repeated_substeps(self):
        """A step longer than max_substep is divided into equal sub-steps."""
        long_step, repeated = make_skimmer(max_substep=2.0), make_skimmer(max_substep=2.0)
        long_step.calc_dt(10.0, 0.01, 0.02)
        for _ in range(5):
            repeated.calc_dt(2.0, 0.01, 0.02)
        # (equal up to finite-difference noise in the flow sensitivities)
        assert levels(long_step) == pytest.approx(levels(repeated), abs=1e-9)

    def test_mass_balance_over_long_step_with_overflow(self, skimmer):
        """The outflows reported for a step spanning many sub-steps are the average over
        the step, so stored volume changes by inflow less outflow."""
        set_levels(skimmer, 0.55, 0.7, 0.3, 0.45)
        dt, vdotmetal_in, vdotslag_in = 60.0, 0.03, 0.05
        vmetal0, vslag0 = total_volumes(skimmer)
        skimmer.calc_dt(dt, vdotmetal_in, vdotslag_in)
        vmetal1, vslag1 = total_volumes(skimmer)
        metal_out = skimmer.slag_trough_vdotmetal_out + skimmer.metal_trough_vdotmetal_out
        slag_out = skimmer.slag_trough_vdotslag_out + skimmer.metal_trough_vdotslag_out
        assert metal_out > 0 and slag_out > 0
        assert vmetal1 - vmetal0 == pytest.approx((vdotmetal_in - metal_out) * dt)
        assert vslag1 - vslag0 == pytest.approx((vdotslag_in - slag_out) * dt)

    @pytest.mark.parametrize("dt", [0.5, 1.0, 2.0])
    def test_steady_state_reached_without_oscillation(self, dt):
        """Regression test: an explicit scheme oscillates between two states at these
        step sizes (its stability limit is about 0.2 s for this flow)."""
        launder = make_skimmer(max_substep=dt)
        vdotmetal_in, vdotslag_in = 0.01, 0.02
        for _ in range(int(3000 / dt)):
            launder.calc_dt(dt, vdotmetal_in, vdotslag_in)
        before = levels(launder)
        launder.calc_dt(dt, vdotmetal_in, vdotslag_in)
        # (an explicit scheme alternates between states about 1e-2 m apart)
        assert levels(launder) == pytest.approx(before, abs=1e-6)
        metal_out = launder.slag_trough_vdotmetal_out + launder.metal_trough_vdotmetal_out
        slag_out = launder.slag_trough_vdotslag_out + launder.metal_trough_vdotslag_out
        assert metal_out == pytest.approx(vdotmetal_in, rel=1e-3)
        assert slag_out == pytest.approx(vdotslag_in, rel=1e-3)

    def test_steady_state_independent_of_step_size(self):
        """The steady state is where the flows balance, so it does not depend on dt."""
        results = []
        for dt in (0.25, 2.0):
            launder = make_skimmer(max_substep=dt)
            for _ in range(int(3000 / dt)):
                launder.calc_dt(dt, 0.01, 0.02)
            results.append(levels(launder))
        assert results[0] == pytest.approx(results[1], abs=1e-6)

    def test_weir_pressure_difference_stays_small_after_shutoff(self):
        """Regression test: after the inflow stops the pressure difference across the weir
        should decay towards zero. An explicit scheme flips its sign on every step with an
        amplitude comparable to the operating value (about 40 Pa at this flow)."""
        launder = make_skimmer()
        for _ in range(3000):
            launder.calc_dt(1.0, 0.01, 0.02)
        assert abs(weir_pressure_difference(launder)) > 30
        for _ in range(600):
            launder.calc_dt(1.0, 0.0, 0.0)
        for _ in range(100):
            launder.calc_dt(1.0, 0.0, 0.0)
            assert abs(weir_pressure_difference(launder)) < 1.0

    def test_first_order_convergence(self):
        """Halving the step size should roughly halve the error against a fine solution."""
        def state_after(max_substep, total_time=60):
            launder = make_skimmer(max_substep=max_substep)
            for _ in range(int(round(total_time / max_substep))):
                launder.calc_dt(max_substep, 0.01, 0.02)
            return numpy.array(
                launder._levels_to_volumes(*levels(launder))
            )

        reference = state_after(0.02)
        errors = [numpy.abs(state_after(dt) - reference).max() for dt in (2.0, 1.0, 0.5)]
        assert errors[0] / errors[1] == pytest.approx(2, abs=0.5)
        assert errors[1] / errors[2] == pytest.approx(2, abs=0.5)


def constant_fraction_model(dh, fraction):
    """Entrainment model with a constant metal fraction, for testing."""
    return fraction


def exp_model(premultiplier, decay):
    return functools.partial(overflowmodel_exp, premultiplier=premultiplier, decay=decay)


class TestSkimmerLaunderEntrainment:
    k = 2 / 3 * numpy.sqrt(g)

    def test_default_is_no_entrainment(self, skimmer):
        set_levels(skimmer, 0.2, 0.6, 0.1, 0.55)
        calc_outlet_flows(skimmer)
        assert skimmer.slag_trough_vdotmetal_out == 0
        assert skimmer.metal_trough_vdotmetal_out == 0

    def test_overflow_models_are_callable_functions_of_distance(self):
        launder = make_skimmer(
            slag_trough_overflowmodel=exp_model(0.2, 2.0),
            metal_trough_overflowmodel=functools.partial(constant_fraction_model, fraction=0.3),
        )
        slag_model = launder.slag_trough_outlet.overflowmodel
        assert slag_model(0.5) == pytest.approx(0.2 * numpy.exp(-1.0))
        assert launder.metal_trough_outlet.overflowmodel(0.5) == 0.3

    def test_entrained_metal_in_slag_trough_slag_only_overflow(self):
        launder = make_skimmer(
            slag_trough_overflowmodel=exp_model(0.2, 2.0)
        )
        set_levels(launder, 0.2, 0.6, 0.1, 0.55)
        calc_outlet_flows(launder)
        vfrac = 0.2 * numpy.exp(-2.0 * (0.5 - 0.2))
        qslag = self.k * 0.6 * 0.5 * 0.1**1.5
        assert launder.slag_trough_vdotmetal_out == pytest.approx(vfrac * qslag)
        assert launder.slag_trough_vdotslag_out == pytest.approx((1 - vfrac) * qslag)

    def test_entrained_metal_in_metal_trough_slag_only_overflow(self):
        launder = make_skimmer(
            metal_trough_overflowmodel=exp_model(0.4, 3.0)
        )
        set_levels(launder, 0.2, 0.6, 0.1, 0.55)
        calc_outlet_flows(launder)
        vfrac = 0.4 * numpy.exp(-3.0 * (0.4 - 0.1))
        qslag = self.k * 0.65 * 0.5 * 0.15**1.5
        assert launder.metal_trough_vdotmetal_out == pytest.approx(vfrac * qslag)
        assert launder.metal_trough_vdotslag_out == pytest.approx((1 - vfrac) * qslag)

    def test_each_trough_uses_its_own_model(self):
        """Entrainment is only active at the outlet whose model is configured."""
        launder = make_skimmer(
            metal_trough_overflowmodel=functools.partial(constant_fraction_model, fraction=0.25),
        )
        set_levels(launder, 0.2, 0.6, 0.1, 0.55)
        calc_outlet_flows(launder)
        assert launder.slag_trough_vdotmetal_out == 0
        qslag = self.k * 0.65 * 0.5 * 0.15**1.5
        assert launder.metal_trough_vdotmetal_out == pytest.approx(0.25 * qslag)
        assert launder.metal_trough_vdotslag_out == pytest.approx(0.75 * qslag)

    def test_entrained_metal_adds_to_metal_overflow(self):
        plain, entraining = make_skimmer(), make_skimmer(
            slag_trough_overflowmodel=functools.partial(constant_fraction_model, fraction=0.1)
        )
        for launder in (plain, entraining):
            set_levels(launder, 0.55, 0.7, 0.3, 0.45)
            calc_outlet_flows(launder)
        qslag = self.k * 0.6 * 0.5 * 0.15**1.5
        assert entraining.slag_trough_vdotmetal_out == pytest.approx(
            plain.slag_trough_vdotmetal_out + 0.1 * qslag
        )
        assert entraining.slag_trough_vdotslag_out == pytest.approx(0.9 * qslag)

    def test_no_entrainment_without_metal(self):
        launder = make_skimmer(
            slag_trough_overflowmodel=functools.partial(constant_fraction_model, fraction=0.3)
        )
        set_levels(launder, 0.0, 0.6, 0.0, 0.0)
        calc_outlet_flows(launder)
        assert launder.slag_trough_vdotmetal_out == 0
        assert launder.slag_trough_vdotslag_out > 0

    def test_mass_balance_with_entrainment_over_long_step(self):
        launder = make_skimmer(
            slag_trough_overflowmodel=exp_model(0.2, 2.0),
            metal_trough_overflowmodel=exp_model(0.3, 1.0),
        )
        set_levels(launder, 0.4, 0.7, 0.3, 0.45)
        dt = 60.0
        vmetal0, vslag0 = total_volumes(launder)
        launder.calc_dt(dt, 0.03, 0.05)
        vmetal1, vslag1 = total_volumes(launder)
        metal_out = launder.slag_trough_vdotmetal_out + launder.metal_trough_vdotmetal_out
        slag_out = launder.slag_trough_vdotslag_out + launder.metal_trough_vdotslag_out
        assert launder.slag_trough_vdotslag_out > 0 and launder.metal_trough_vdotslag_out > 0
        assert vmetal1 - vmetal0 == pytest.approx((0.03 - metal_out) * dt)
        assert vslag1 - vslag0 == pytest.approx((0.05 - slag_out) * dt)

    def test_steady_state_with_entrainment_outflow_matches_inflow(self):
        launder = make_skimmer(
            slag_trough_overflowmodel=exp_model(0.2, 2.0),
            metal_trough_overflowmodel=exp_model(0.3, 1.0),
        )
        for _ in range(3000):
            launder.calc_dt(1.0, 0.01, 0.02)
        metal_out = launder.slag_trough_vdotmetal_out + launder.metal_trough_vdotmetal_out
        slag_out = launder.slag_trough_vdotslag_out + launder.metal_trough_vdotslag_out
        assert metal_out == pytest.approx(0.01, rel=1e-3)
        assert slag_out == pytest.approx(0.02, rel=1e-3)


class TestSkimmerLaunderWeirFlow:
    """Flow through the weir underflow channel (height 0.3 m, width 0.3 m in the fixture)."""

    area = 0.3 * 0.3

    def expected_flows(self, source_hmetal, source_hslag, deltaP, kl=0.0):
        """Metal and slag flowrates, with the phase proportions of the channel taken from the
        trough the flow comes from, and the sign of the pressure difference, for an entrance
        loss of kl velocity heads."""
        height = 0.3
        if source_hmetal > height:
            metal_fraction, slag_fraction = 1.0, 0.0
        elif source_hslag > height:
            metal_fraction, slag_fraction = source_hmetal / height, 1 - source_hmetal / height
        else:
            metal_fraction = source_hmetal / height
            slag_fraction = (source_hslag - source_hmetal) / height
        sign = numpy.sign(deltaP)
        return (
            metal_fraction * self.area * sign
            * numpy.sqrt(2 * abs(deltaP) / ((1 + kl) * RHO_METAL)),
            slag_fraction * self.area * sign
            * numpy.sqrt(2 * abs(deltaP) / ((1 + kl) * RHO_SLAG)),
        )

    @staticmethod
    def pressure_difference(hm_ss, hs_ss, hm_ms, hs_ms):
        def pressure(hmetal, hslag):
            return (hmetal * RHO_METAL + (hslag - hmetal) * RHO_SLAG) * g

        return pressure(hm_ss, hs_ss) - pressure(hm_ms, hs_ms)

    # (hmetal, hslag) of the trough the flow comes from, with the interface above the top of
    # the channel, below it with the slag surface above it, and with both below it
    sources = [(0.4, 0.5), (0.2, 0.5), (0.1, 0.25)]
    receivers = [(0.05, 0.06), (0.0, 0.0), (0.05, 0.05)]

    @pytest.mark.parametrize("source", sources)
    def test_forward_flow_takes_proportions_from_slag_trough(self, skimmer, source):
        receiver = self.receivers[0]
        levels = (*source, *receiver)
        deltaP = self.pressure_difference(*levels)
        assert deltaP > 0
        expected = self.expected_flows(*source, deltaP)
        assert skimmer._calc_vdot_weir(*levels) == pytest.approx(expected, rel=1e-12)

    @pytest.mark.parametrize("source", sources)
    def test_reversed_flow_takes_proportions_from_metal_trough(self, skimmer, source):
        receiver = self.receivers[0]
        levels = (*receiver, *source)
        deltaP = self.pressure_difference(*levels)
        assert deltaP < 0
        expected = self.expected_flows(*source, deltaP)
        assert skimmer._calc_vdot_weir(*levels) == pytest.approx(expected, rel=1e-12)

    def test_flow_is_antisymmetric_when_the_troughs_are_swapped(self, skimmer):
        """The two troughs are treated alike: swapping them reverses the flow."""
        for source in self.sources:
            for receiver in self.receivers:
                forward = skimmer._calc_vdot_weir(*source, *receiver)
                reverse = skimmer._calc_vdot_weir(*receiver, *source)
                assert reverse == pytest.approx((-forward[0], -forward[1]), rel=1e-12, abs=1e-15)

    def test_no_flow_when_pressures_are_equal(self, skimmer):
        assert skimmer._calc_vdot_weir(0.2, 0.28, 0.2, 0.28) == (0, 0)

    def test_reversed_flow_cannot_draw_slag_from_a_trough_without_slag(self, skimmer):
        """Regression test: with the interface above the channel and no slag layer in the
        metal trough, a reversed flow carries only metal. It used to carry slag in the
        proportions of the slag trough, which could remove more slag than the metal trough
        held."""
        levels = (0.15, 0.4, 0.45, 0.45)
        assert self.pressure_difference(*levels) < 0
        vdot_metal, vdot_slag = skimmer._calc_vdot_weir(*levels)
        assert vdot_metal < 0
        assert vdot_slag == 0


class TestSkimmerLaunderWeirEntryLoss:
    """The loss of kl velocity heads at the entrance to the weir underflow channel."""

    helper = TestSkimmerLaunderWeirFlow()

    # the interface above the top of the channel, below it, and the slag surface below it
    sources = TestSkimmerLaunderWeirFlow.sources
    receivers = TestSkimmerLaunderWeirFlow.receivers

    def test_default_is_lossless(self, skimmer):
        assert skimmer.weir_underflow_kl == 0.0

    def test_loss_coefficient_stored(self):
        assert make_skimmer(weir_underflow_kl=0.5).weir_underflow_kl == 0.5

    def test_negative_loss_coefficient_rejected(self):
        with pytest.raises(ValueError, match="weir_underflow_kl"):
            make_skimmer(weir_underflow_kl=-0.1)

    def test_explicit_zero_is_the_same_as_the_default(self):
        explicit, default = make_skimmer(weir_underflow_kl=0), make_skimmer()
        levels = (0.4, 0.5, 0.05, 0.06)
        assert explicit._calc_vdot_weir(*levels) == default._calc_vdot_weir(*levels)

    @pytest.mark.parametrize("kl", [0.25, 1.0, 4.0])
    @pytest.mark.parametrize("source", sources)
    def test_flows_match_the_bernoulli_equation_with_the_loss(self, source, kl):
        skimmer = make_skimmer(weir_underflow_kl=kl)
        levels = (*source, *self.receivers[0])
        deltaP = self.helper.pressure_difference(*levels)
        expected = self.helper.expected_flows(*source, deltaP, kl)
        assert skimmer._calc_vdot_weir(*levels) == pytest.approx(expected, rel=1e-12)

    @pytest.mark.parametrize("kl", [0.25, 1.0, 4.0])
    def test_flows_are_reduced_by_the_root_of_one_plus_kl(self, kl):
        lossy, lossless = make_skimmer(weir_underflow_kl=kl), make_skimmer()
        for source in self.sources:
            for receiver in self.receivers:
                lossy_flows = lossy._calc_vdot_weir(*source, *receiver)
                lossless_flows = lossless._calc_vdot_weir(*source, *receiver)
                assert lossy_flows == pytest.approx(
                    tuple(q / numpy.sqrt(1 + kl) for q in lossless_flows), rel=1e-12, abs=1e-15)

    def test_loss_applies_to_reversed_flow(self):
        skimmer = make_skimmer(weir_underflow_kl=1.0)
        for source in self.sources:
            for receiver in self.receivers:
                forward = skimmer._calc_vdot_weir(*source, *receiver)
                reverse = skimmer._calc_vdot_weir(*receiver, *source)
                assert reverse == pytest.approx((-forward[0], -forward[1]), rel=1e-12, abs=1e-15)

    def test_no_flow_when_pressures_are_equal(self):
        skimmer = make_skimmer(weir_underflow_kl=2.0)
        assert skimmer._calc_vdot_weir(0.2, 0.28, 0.2, 0.28) == (0, 0)

    def test_loss_holds_more_metal_in_the_slag_trough(self):
        """The same supply passes through a less permeable channel at a larger level
        difference, with nothing lost over the outlets."""
        lossless, lossy = make_skimmer(), make_skimmer(weir_underflow_kl=4.0)
        for skimmer in (lossless, lossy):
            for _ in range(60):
                skimmer.calc_dt(1.0, 0.003, 0.004)
            assert total_volumes(skimmer) == pytest.approx((60 * 0.003, 60 * 0.004), rel=1e-9)
        assert lossy.slag_trough_hmetal > lossless.slag_trough_hmetal
        assert lossy.metal_trough_hmetal < lossless.metal_trough_hmetal


class TestSkimmerLaunderPhysicalVolumes:
    """The phase volumes stay positive when the weir flow reverses and for large sub-steps."""

    @staticmethod
    def smallest_volume(launder):
        return min(launder._levels_to_volumes(*levels(launder)))

    def test_reversed_flow_conserves_volume_and_stays_positive(self, skimmer):
        """A much fuller metal trough drains back through the weir."""
        set_levels(skimmer, 0.05, 0.08, 0.45, 0.46)
        assert weir_pressure_difference(skimmer) < 0
        dt = 60.0
        vmetal0, vslag0 = total_volumes(skimmer)
        skimmer.calc_dt(dt, 0.0, 0.0)
        vmetal1, vslag1 = total_volumes(skimmer)
        metal_out = skimmer.slag_trough_vdotmetal_out + skimmer.metal_trough_vdotmetal_out
        slag_out = skimmer.slag_trough_vdotslag_out + skimmer.metal_trough_vdotslag_out
        assert vmetal1 - vmetal0 == pytest.approx(-metal_out * dt, abs=1e-9)
        assert vslag1 - vslag0 == pytest.approx(-slag_out * dt, abs=1e-9)
        assert self.smallest_volume(skimmer) >= 0

    def test_draining_after_tap_shut_off_stays_positive(self):
        """Regression test: with entrainment at the slag trough outlet, draining a launder
        after the inflow stops used to give a negative slag volume, as the reversed weir flow
        drew slag from a trough holding less than it removed."""
        launder = make_skimmer(
            slag_trough_overflowmodel=functools.partial(
                overflowmodel_exp, premultiplier=0.4, decay=3.0
            )
        )
        for _ in range(600):
            launder.calc_dt(1.0, 0.03, 0.06)
        for _ in range(2000):
            launder.calc_dt(1.0, 0.0, 0.0)
            assert self.smallest_volume(launder) >= -1e-12

    def test_random_inflows_and_steps_stay_positive(self):
        rng = numpy.random.default_rng(11)
        launder = make_skimmer(
            slag_trough_overflowmodel=exp_model(0.2, 2.0),
            metal_trough_overflowmodel=exp_model(0.3, 1.0),
        )
        for _ in range(300):
            dt = rng.uniform(0.05, 60)
            launder.calc_dt(dt, rng.uniform(0, 0.05), rng.uniform(0, 0.08))
            assert self.smallest_volume(launder) >= -1e-12

    def test_large_substeps_stay_positive(self):
        """Sub-steps far beyond the accurate range are repeated as smaller steps instead of
        overshooting to negative volumes."""
        rng = numpy.random.default_rng(3)
        launder = make_skimmer(max_substep=5.0)
        for _ in range(200):
            launder.calc_dt(20.0, rng.uniform(0, 0.05), rng.uniform(0, 0.08))
            assert self.smallest_volume(launder) >= -1e-12
