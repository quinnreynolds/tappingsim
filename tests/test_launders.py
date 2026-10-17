"""Unit tests for tappingsim.launders module."""

import numpy
import pytest
from scipy.constants import g
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
        assert skimmer.metal_trough_outlet_metal_Cd == 0.5
        assert skimmer.density_slag == RHO_SLAG

    def test_outlet_premult(self, skimmer):
        assert skimmer.outlet_premult == pytest.approx(2 / 3 * numpy.sqrt(g))


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
        skimmer._calc_vdot_out(1.0)
        assert skimmer.slag_trough_vdotmetal_out == 0
        assert skimmer.slag_trough_vdotslag_out == 0
        assert skimmer.metal_trough_vdotmetal_out == 0
        assert skimmer.metal_trough_vdotslag_out == 0

    def test_slag_only_overflow(self, skimmer):
        """Slag above the crest with metal below it: broad-crested weir flow of slag."""
        set_levels(skimmer, 0.2, 0.6, 0.1, 0.55)
        skimmer._calc_vdot_out(1.0)
        k = 2 / 3 * numpy.sqrt(g)
        assert skimmer.slag_trough_vdotslag_out == pytest.approx(k * 0.6 * 0.5 * 0.1**1.5)
        assert skimmer.slag_trough_vdotmetal_out == 0
        assert skimmer.metal_trough_vdotslag_out == pytest.approx(k * 0.65 * 0.5 * 0.15**1.5)
        assert skimmer.metal_trough_vdotmetal_out == 0

    def test_metal_and_slag_overflow(self, skimmer):
        """Metal above the crest: slag flows on its layer thickness, metal on the
        slag-equivalent head plus its own head above the crest."""
        set_levels(skimmer, 0.55, 0.7, 0.5, 0.6)
        skimmer._calc_vdot_out(1.0)
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
        skimmer._calc_vdot_out(1.0)
        set_levels(skimmer, 0.2, 0.45, 0.1, 0.35)
        skimmer._calc_vdot_out(1.0)
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
        assert levels(long_step) == pytest.approx(levels(repeated), abs=1e-12)

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
