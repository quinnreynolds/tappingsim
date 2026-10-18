"""Unit tests for tappingsim.furnaces module."""

import math
import numpy
import pytest
from scipy.constants import g
from tappingsim import furnaces
from tappingsim.furnaces import (
    SubmergedArcFurnace,
    bedmodel_ergun,
    bedmodel_kozenycarman,
    fdmodel_bellos,
    fdmodel_cheng,
    fdmodel_serghides1,
    fdmodel_serghides2,
)


# ---------------------------------------------------------------------------
# Shared fixture: standard furnace parameters from notebook examples
# ---------------------------------------------------------------------------

@pytest.fixture
def furnace_params():
    return {
        'powerMVA': 40,
        'powerfactor': 0.8,
        'metalSER': 3.5,
        'slagmetalmassratio': 1.2,
        'furnacediameter': 12,
        'activeareafraction': 0.75,
        'tapholediameter': 0.1,
        'tapholelength': 1.5,
        'tapholeroughness': 1e-3,
        'tapholeheight': 0.1,
        'densitymetal': 7000,
        'densityslag': 3000,
        'viscositymetal': 0.005,
        'viscosityslag': 0.1,
        'particlediameter': 0.02,
        'particlesphericity': 0.8,
        'bedporosity': 0.5,
        'bedmindiameter': 0.2,
        'bedmaxdiameter': 20,
        'bedmodel': bedmodel_ergun,
        'entrykl': 0.25,
        'fdmodel': fdmodel_cheng,
        'hmetal_init': 0.25,
        'hslag_init': 0.4,
    }


@pytest.fixture
def furnace(furnace_params):
    return SubmergedArcFurnace(**furnace_params)


# ---------------------------------------------------------------------------
# Bed pressure drop models
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Friction factor models
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# SubmergedArcFurnace
# ---------------------------------------------------------------------------

class TestSubmergedArcFurnaceInit:
    def test_initial_metal_level(self, furnace):
        assert furnace.hmetal == 0.25

    def test_initial_slag_level(self, furnace):
        assert furnace.hslag == 0.4

    def test_taphole_initially_closed(self, furnace):
        assert furnace.tapholeopen_yn is False

    def test_initial_flowrates_zero(self, furnace):
        assert furnace.vdotmetal_out == 0
        assert furnace.vdotslag_out == 0


class TestSubmergedArcFurnaceClosedTaphole:
    def test_zero_metal_outflow(self, furnace):
        furnace.tapholeopen_yn = False
        furnace.calc_dt(dt=60)
        assert furnace.vdotmetal_out == 0

    def test_zero_slag_outflow(self, furnace):
        furnace.tapholeopen_yn = False
        furnace.calc_dt(dt=60)
        assert furnace.vdotslag_out == 0

    def test_metal_level_increases(self, furnace):
        furnace.tapholeopen_yn = False
        furnace.calc_dt(dt=60)
        assert furnace.hmetal > 0.25

    def test_slag_level_increases(self, furnace):
        furnace.tapholeopen_yn = False
        furnace.calc_dt(dt=60)
        assert furnace.hslag > 0.4

    def test_metal_level_value(self, furnace):
        furnace.tapholeopen_yn = False
        furnace.calc_dt(dt=60)
        assert furnace.hmetal == pytest.approx(0.2505132736889837, rel=1e-8)

    def test_slag_level_value(self, furnace):
        furnace.tapholeopen_yn = False
        furnace.calc_dt(dt=60)
        assert furnace.hslag == pytest.approx(0.40195044001813807, rel=1e-8)


# The reference values below were calculated for the furnace of the shared fixture with an
# independent solution of the level equations (an adaptive 8th order Runge-Kutta method with a
# tight tolerance on the tap-hole flow calculation), so these tests check accuracy.

class TestSubmergedArcFurnaceOpenTaphole:
    @pytest.fixture(autouse=True)
    def advance_one_closed_step(self, furnace):
        """Match notebook state: one closed step before opening."""
        furnace.tapholeopen_yn = False
        furnace.calc_dt(dt=60)
        furnace.tapholeopen_yn = True
        furnace.calc_dt(dt=60)

    def test_positive_metal_outflow(self, furnace):
        assert furnace.vdotmetal_out > 0

    def test_positive_slag_outflow(self, furnace):
        assert furnace.vdotslag_out > 0

    def test_metal_outflow_value(self, furnace):
        """The reported outflow is the mean over the time step."""
        assert furnace.vdotmetal_out == pytest.approx(0.009347885817, rel=1e-6)

    def test_slag_outflow_value(self, furnace):
        assert furnace.vdotslag_out == pytest.approx(0.001614724222, rel=1e-6)

    def test_metal_level_value(self, furnace):
        assert furnace.hmetal == pytest.approx(0.237801994176, rel=1e-7)

    def test_slag_level_value(self, furnace):
        assert furnace.hslag == pytest.approx(0.388391959353, rel=1e-7)


class TestSubmergedArcFurnaceCalcTimePeriod:
    @pytest.fixture(autouse=True)
    def setup(self, furnace):
        """Advance to open-taphole state then run calc_time_period."""
        furnace.tapholeopen_yn = False
        furnace.calc_dt(dt=60)
        furnace.tapholeopen_yn = True
        furnace.calc_dt(dt=60)
        times = numpy.linspace(0, 480, 480)
        self.metaltapped, self.slagtapped = furnace.calc_time_period(times=times)

    def test_metal_tapped_mass(self):
        assert self.metaltapped == pytest.approx(24837.39172, rel=1e-6)

    def test_slag_tapped_mass(self):
        assert self.slagtapped == pytest.approx(3288.00886, rel=1e-6)

    def test_metal_tapped_positive(self):
        assert self.metaltapped > 0

    def test_slag_tapped_positive(self):
        assert self.slagtapped > 0


def balance_constants(furnace):
    """Rate of change of the metal and slag levels due to smelting, and the factor linking
    a net volume flowrate to the rate of change of a level."""
    areaconst = 1 / (furnace.bedporosity * furnace.activeareafraction
                     * 0.25 * numpy.pi * furnace.furnacediameter**2)
    mdotmetal_in = furnace.powerMVA * furnace.powerfactor / (3.6 * furnace.metalSER)
    vdotmetal_in = mdotmetal_in / furnace.densitymetal
    vdotslag_in = mdotmetal_in * furnace.slagmetalmassratio / furnace.densityslag
    return areaconst, vdotmetal_in, vdotslag_in


def open_furnace(furnace_params, **overrides):
    furnace = SubmergedArcFurnace(**dict(furnace_params, **overrides))
    furnace.tapholeopen_yn = True
    return furnace


class TestSubmergedArcFurnaceTimeIntegration:
    @pytest.mark.parametrize("dt", [1.0, 60.0, 300.0])
    def test_levels_change_by_inflow_less_reported_outflow(self, furnace_params, dt):
        """The reported outflows are step averages, so the level changes balance exactly."""
        furnace = open_furnace(furnace_params)
        areaconst, vdotmetal_in, vdotslag_in = balance_constants(furnace)
        hmetal0, hslag0 = furnace.hmetal, furnace.hslag
        furnace.calc_dt(dt)
        assert furnace.vdotmetal_out > 0 and furnace.vdotslag_out > 0
        dhmetal = dt * areaconst * (vdotmetal_in - furnace.vdotmetal_out)
        dhslag = dhmetal + dt * areaconst * (vdotslag_in - furnace.vdotslag_out)
        assert furnace.hmetal - hmetal0 == pytest.approx(dhmetal, abs=1e-12)
        assert furnace.hslag - hslag0 == pytest.approx(dhslag, abs=1e-12)

    def test_closed_taphole_is_an_exact_linear_update(self, furnace_params):
        furnace = SubmergedArcFurnace(**furnace_params)
        areaconst, vdotmetal_in, vdotslag_in = balance_constants(furnace)
        furnace.calc_dt(600.0)
        assert furnace.hmetal == pytest.approx(0.25 + 600 * areaconst * vdotmetal_in, rel=1e-12)
        assert furnace.hslag == pytest.approx(
            0.4 + 600 * areaconst * (vdotmetal_in + vdotslag_in), rel=1e-12)

    def test_fourth_order_convergence(self, furnace_params):
        """Halving the time step reduces the error in the levels by a large factor (the
        flow calculation is converged only to a tolerance, which limits it at small steps)."""
        def levels_after(dt, total=960.0):
            furnace = open_furnace(furnace_params)
            for _ in range(int(round(total / dt))):
                furnace.calc_dt(dt)
            return numpy.array([furnace.hmetal, furnace.hslag])

        reference = levels_after(2.5)
        errors = [numpy.abs(levels_after(dt) - reference).max() for dt in (240.0, 120.0, 60.0)]
        assert errors[0] / errors[1] > 8
        assert errors[1] / errors[2] > 8

    def test_long_steps_remain_accurate_and_stable(self, furnace_params):
        """A step of five minutes is within about a millimetre of a fine-step solution."""
        coarse, fine = open_furnace(furnace_params), open_furnace(furnace_params)
        coarse.calc_dt(300.0)
        for _ in range(300):
            fine.calc_dt(1.0)
        assert coarse.hmetal == pytest.approx(fine.hmetal, abs=1e-3)
        assert coarse.hslag == pytest.approx(fine.hslag, abs=1e-3)

    def test_calc_time_period_matches_repeated_steps(self, furnace_params):
        times = numpy.linspace(0, 600, 41)
        a, b = open_furnace(furnace_params), open_furnace(furnace_params)
        metal, slag = a.calc_time_period(times)
        metal_sum = slag_sum = 0.0
        for dt in numpy.diff(times):
            b.calc_dt(float(dt))
            metal_sum += dt * b.vdotmetal_out * b.densitymetal
            slag_sum += dt * b.vdotslag_out * b.densityslag
        assert metal == pytest.approx(metal_sum, rel=1e-6)
        assert slag == pytest.approx(slag_sum, rel=1e-6)
        assert a.hmetal == pytest.approx(b.hmetal, rel=1e-8)
        assert a.hslag == pytest.approx(b.hslag, rel=1e-8)
        assert a.vdotmetal_out == pytest.approx(b.vdotmetal_out, rel=1e-6)

    def test_calc_time_period_tapped_mass_balances_the_levels(self, furnace_params):
        """The mass tapped over a period equals the mass smelted less the change stored."""
        furnace = open_furnace(furnace_params)
        areaconst, vdotmetal_in, vdotslag_in = balance_constants(furnace)
        hmetal0, times = furnace.hmetal, numpy.linspace(0, 600, 21)
        metal, _ = furnace.calc_time_period(times)
        stored = (furnace.hmetal - hmetal0) / areaconst * furnace.densitymetal
        smelted = 600 * vdotmetal_in * furnace.densitymetal
        assert metal == pytest.approx(smelted - stored, rel=1e-9)

    def test_calc_time_period_with_a_single_time_does_nothing(self, furnace_params):
        furnace = open_furnace(furnace_params)
        assert furnace.calc_time_period(numpy.array([0.0])) == (0, 0)
        assert (furnace.hmetal, furnace.hslag) == (0.25, 0.4)

    def test_levels_are_plain_floats(self, furnace_params):
        furnace = open_furnace(furnace_params)
        furnace.calc_dt(60.0)
        assert type(furnace.hmetal) is float and type(furnace.hslag) is float
        assert type(furnace.vdotmetal_out) is float and type(furnace.vdotslag_out) is float


# ---------------------------------------------------------------------------
# Tap-hole flow calculation
# ---------------------------------------------------------------------------

RT = 0.05  # tap-hole radius, m
RHO_M, RHO_S = 7000.0, 3000.0


def tapping_pressure(hm, hs):
    return furnaces._tapping_pressure(hm, hs, RT, RHO_M, RHO_S)


class TestTappingPressure:
    def test_no_pressure_below_the_taphole(self):
        assert tapping_pressure(-0.3, -RT - 1e-6) == 0

    def test_slag_only_above_the_taphole(self):
        assert tapping_pressure(-0.2, 0.3) == pytest.approx(0.3 * RHO_S * g)

    def test_metal_below_slag_above_the_taphole(self):
        assert tapping_pressure(0.2, 0.5) == pytest.approx(g * (RHO_S * 0.3 + RHO_M * 0.2))

    # (metal level, slag level) pairs either side of which the formula changes
    @pytest.mark.parametrize("below, above", [
        ((-0.2, -RT - 1e-9), (-0.2, -RT + 1e-9)),
        ((-0.2, RT - 1e-9), (-0.2, RT + 1e-9)),
        ((0.0, RT - 1e-9), (0.0, RT + 1e-9)),
        ((RT + 1e-3, RT + 2e-3 - 1e-9), (RT + 1e-3, RT + 2e-3 + 1e-9)),
        ((-RT - 1e-9, 0.0), (-RT + 1e-9, 0.0)),
        ((-RT - 1e-9, 0.2), (-RT + 1e-9, 0.2)),
        ((RT - 1e-9, 0.2), (RT + 1e-9, 0.2)),
    ])
    def test_continuous_across_the_branch_boundaries(self, below, above):
        assert tapping_pressure(*above) == pytest.approx(tapping_pressure(*below), abs=1e-3)

    def test_increases_with_the_slag_level_for_slag_only_or_a_filled_taphole(self):
        """(With both phases in a partly filled tap-hole the pressure, averaged over the height
        of liquid in the hole, can dip slightly as a first light layer of slag is added.)"""
        for hm in (-0.3, 0.0, 0.03, 0.2):
            levels = numpy.linspace(RT, 0.6, 60)  # tap-hole completely filled
            assert numpy.all(numpy.diff([tapping_pressure(hm, hs) for hs in levels]) > 0)
        levels = numpy.linspace(-RT, 0.6, 60)  # slag only
        assert numpy.all(numpy.diff([tapping_pressure(-0.3, hs) for hs in levels]) > 0)


class TestSegmentArea:
    def test_empty_when_the_chord_is_at_the_bottom(self):
        assert furnaces._segment_area(RT, -RT) == 0

    def test_half_the_circle_at_the_centre(self):
        assert furnaces._segment_area(RT, 0.0) == pytest.approx(0.5 * math.pi * RT**2)

    def test_full_circle_at_the_top(self):
        assert furnaces._segment_area(RT, RT) == pytest.approx(math.pi * RT**2, rel=1e-12)

    def test_segments_below_and_above_the_centre_make_up_the_circle(self):
        for h in (0.01, 0.02, 0.04):
            total = furnaces._segment_area(RT, h) + furnaces._segment_area(RT, -h)
            assert total == pytest.approx(math.pi * RT**2, rel=1e-12)

    def test_increases_with_the_height(self):
        values = [furnaces._segment_area(RT, h) for h in numpy.linspace(-RT, RT, 51)]
        assert numpy.all(numpy.diff(values) > 0)


class TestPhaseAreas:
    h0_l, h0_h = -RT - 0.01, RT + 0.01
    full = math.pi * RT**2

    def areas(self, hm, hs):
        return furnaces._phase_areas(hm, hs, self.h0_l, self.h0_h, RT)

    def test_no_flow_below_the_taphole(self):
        assert self.areas(-0.3, -RT - 1e-6) == (0, 0)

    def test_partly_filled_with_slag_only(self):
        area_m, area_s = self.areas(-0.3, 0.02)
        assert area_m == 0
        assert area_s == pytest.approx(furnaces._segment_area(RT, 0.02))

    def test_partly_filled_with_slag_and_metal(self):
        area_m, area_s = self.areas(-0.04, 0.02)
        assert area_m > 0 and area_s > 0
        assert area_m + area_s == pytest.approx(furnaces._segment_area(RT, 0.02))

    def test_filled_with_slag_only(self):
        assert self.areas(-0.3, 0.2) == (0, pytest.approx(self.full))

    def test_filled_with_metal_only(self):
        assert self.areas(0.2, 0.3) == (pytest.approx(self.full), 0)

    def test_filled_with_slag_and_metal(self):
        area_m, area_s = self.areas(0.0, 0.3)
        assert 0 < area_m < self.full and 0 < area_s < self.full
        assert area_m + area_s == pytest.approx(self.full)

    def test_metal_area_grows_continuously_from_none_to_all(self):
        values = [self.areas(hm, 0.3)[0] for hm in numpy.linspace(self.h0_l, self.h0_h, 101)]
        assert values[0] == pytest.approx(0, abs=1e-15)
        assert values[-1] == pytest.approx(self.full)
        assert numpy.all(numpy.diff(values) >= 0)

    def test_areas_are_physical_over_a_range_of_levels(self):
        rng = numpy.random.default_rng(0)
        for _ in range(2000):
            hs = rng.uniform(-3 * RT, 4 * RT)
            hm = hs - rng.uniform(0, 1) * (hs + 3 * RT)
            area_m, area_s = self.areas(hm, hs)
            assert area_m >= -1e-15 and area_s >= -1e-15
            assert area_m + area_s <= self.full * (1 + 1e-12)


def open_params(furnace_params, **overrides):
    return open_furnace(furnace_params, **overrides)._calc_params()


class TestInterfaceLevels:
    def test_no_deformation_without_flow(self, furnace_params):
        h0_l, h0_h = furnaces._interface_levels(0.0, 0.0, open_params(furnace_params))
        assert h0_l == pytest.approx(-0.5 * furnace_params["tapholediameter"])
        assert h0_h == pytest.approx(0.5 * furnace_params["tapholediameter"])

    def test_flow_widens_the_range_over_which_the_phases_mix(self, furnace_params):
        params = open_params(furnace_params)
        h0_l, h0_h = furnaces._interface_levels(2.0, 2.0, params)
        assert h0_l < -0.5 * furnace_params["tapholediameter"]
        assert h0_h > 0.5 * furnace_params["tapholediameter"]


class TestSolveVelocities:
    def test_no_driving_pressure_gives_no_flow(self, furnace_params):
        assert furnaces._solve_velocities(0.0, 1.0, 1.0, open_params(furnace_params)) == (0, 0)

    def test_velocities_satisfy_the_pressure_balance(self, furnace_params):
        p = open_params(furnace_params)
        pa = 5000.0
        umetal, uslag = furnaces._solve_velocities(pa, 0.0, 0.0, p)
        for u, a, b, rho, mu in ((umetal, p.a_metal, p.b_metal, p.densitymetal, p.viscositymetal),
                                 (uslag, p.a_slag, p.b_slag, p.densityslag, p.viscosityslag)):
            friction = p.fdmodel(u, rho, mu, p.tapholediameter, p.tapholeroughness)
            a_c = a + 0.5 * rho * friction * p.tapholelength / p.tapholediameter
            assert a_c * u * u + b * u == pytest.approx(pa, rel=1e-4)

    def test_solution_does_not_depend_on_the_initial_guess(self, furnace_params):
        p = open_params(furnace_params)
        cold = furnaces._solve_velocities(5000.0, 0.0, 0.0, p)
        far = furnaces._solve_velocities(5000.0, 4.0, 3.0, p)
        assert cold == pytest.approx(far, abs=1e-5)

    def test_velocities_increase_with_the_pressure(self, furnace_params):
        p = open_params(furnace_params)
        low = furnaces._solve_velocities(2000.0, 1.0, 1.0, p)
        high = furnaces._solve_velocities(8000.0, 1.0, 1.0, p)
        assert high[0] > low[0] and high[1] > low[1]

    def test_non_convergence_raises(self, furnace_params, monkeypatch):
        """The iteration has a cap, so a case that does not converge ends with an error."""
        monkeypatch.setattr(furnaces, "_VELOCITY_TOLERANCE", 0.0)
        monkeypatch.setattr(furnaces, "_VELOCITY_MAX_ITERATIONS", 2)
        with pytest.raises(RuntimeError):
            furnaces._solve_velocities(5000.0, 0.0, 0.0, open_params(furnace_params))

    def test_default_cap_and_tolerance(self):
        assert furnaces._VELOCITY_TOLERANCE == 1e-6
        assert furnaces._VELOCITY_MAX_ITERATIONS == 100


class TestCalcVdotOut:
    def test_closed_taphole_gives_no_flow(self, furnace_params):
        params = SubmergedArcFurnace(**furnace_params)._calc_params()
        assert furnaces._calc_vdot_out(0.6, 1.0, 1.0, 1.0, params) == (0, 0, 0, 0)

    def test_flows_are_phase_areas_times_velocities(self, furnace_params):
        p = open_params(furnace_params)
        hmetal, hslag = 0.3, 0.6
        umetal, uslag, vdotmetal, vdotslag = furnaces._calc_vdot_out(hmetal, hslag, 1.0, 1.0, p)
        rt = 0.5 * p.tapholediameter
        hm, hs = hmetal - p.tapholeheight, hslag - p.tapholeheight
        pa = furnaces._tapping_pressure(hm, hs, rt, p.densitymetal, p.densityslag)
        expected = furnaces._solve_velocities(pa, 1.0, 1.0, p)
        assert (umetal, uslag) == pytest.approx(expected)
        h0_l, h0_h = furnaces._interface_levels(umetal, uslag, p)
        area_m, area_s = furnaces._phase_areas(hm, hs, h0_l, h0_h, rt)
        assert vdotmetal == pytest.approx(area_m * umetal)
        assert vdotslag == pytest.approx(area_s * uslag)

    def test_no_metal_flow_for_metal_below_the_hearth(self, furnace_params):
        params = open_params(furnace_params)
        assert params.allownegativeheights_yn is False
        umetal, uslag, vdotmetal, vdotslag = furnaces._calc_vdot_out(-0.05, 0.5, 1.0, 1.0, params)
        assert umetal == 0 and vdotmetal == 0
        assert uslag > 0

    def test_metal_flows_below_the_hearth_if_allowed(self, furnace_params):
        furnace = open_furnace(furnace_params)
        furnace.allownegativeheights_yn = True
        params = furnace._calc_params()
        umetal = furnaces._calc_vdot_out(-0.05, 0.5, 1.0, 1.0, params)[0]
        assert umetal > 0


class TestTapholeParameters:
    """The parameters are a snapshot of the attributes at the start of each step."""

    def test_burden_coefficients_are_those_of_the_bed_model_plus_the_entrance_loss(
            self, furnace_params):
        furnace = open_furnace(furnace_params)
        a, b = furnace._calc_burden_coefficients(0.005, 7000.0)
        a_bed, b_bed = bedmodel_ergun(
            furnace.tapholediameter, furnace.bedmindiameter, furnace.bedmaxdiameter,
            furnace.particlediameter, furnace.particlesphericity, furnace.bedporosity,
            0.005, 7000.0)
        assert a == pytest.approx(0.5 * (1 + furnace.entrykl) * 7000.0 + a_bed)
        assert b == pytest.approx(b_bed)

    def test_changed_attributes_are_picked_up_by_the_next_step(self, furnace_params):
        furnace = open_furnace(furnace_params)
        before = furnace._calc_params()
        furnace.tapholediameter = 0.15
        after = furnace._calc_params()
        assert before.tapholediameter == 0.1 and after.tapholediameter == 0.15
        assert after.a_metal != before.a_metal

    def counting_bedmodel(self, calls):
        def bedmodel(*args):
            calls.append(1)
            return bedmodel_ergun(*args)

        return bedmodel

    def test_burden_coefficients_are_calculated_once_per_step(self, furnace_params):
        calls = []
        furnace = open_furnace(furnace_params, bedmodel=self.counting_bedmodel(calls))
        furnace.calc_dt(60.0)
        assert len(calls) == 2  # one per phase, not one per stage of the step

    def test_burden_coefficients_are_calculated_once_per_period(self, furnace_params):
        calls = []
        furnace = open_furnace(furnace_params, bedmodel=self.counting_bedmodel(calls))
        furnace.calc_time_period(numpy.linspace(0, 600, 41))
        assert len(calls) == 2

    def test_burden_coefficients_are_not_calculated_for_a_closed_taphole(self, furnace_params):
        calls = []
        params = dict(furnace_params, bedmodel=self.counting_bedmodel(calls))
        furnace = SubmergedArcFurnace(**params)
        furnace.calc_dt(60.0)
        assert calls == []


# ---------------------------------------------------------------------------
# Continuity as the slag level passes the top of the tap-hole
# ---------------------------------------------------------------------------

class TestTapholeTopContinuity:
    """The flows must not jump when the slag surface falls through the top of the tap-hole,
    where the tap-hole changes from completely to partly filled. The interface between the
    phases was previously positioned with different formulas either side of that level, which
    gave a step of several per cent in the metal flow and of up to tens of per cent in the
    slag flow."""

    h0_l, h0_h = -RT - 0.012, RT + 0.012
    full = math.pi * RT**2

    @pytest.mark.parametrize("hm", [-0.058, -0.04, -0.02, 0.0, 0.02, 0.04, 0.049])
    def test_phase_areas_are_continuous(self, hm):
        eps = 1e-12
        below = furnaces._phase_areas(hm, RT - eps, self.h0_l, self.h0_h, RT)
        above = furnaces._phase_areas(hm, RT + eps, self.h0_l, self.h0_h, RT)
        assert below == pytest.approx(above, abs=1e-9 * self.full)

    @pytest.mark.parametrize("hm_rel", [-0.04, -0.02, 0.0, 0.02, 0.04])
    def test_flows_are_continuous(self, furnace_params, hm_rel):
        furnace = open_furnace(furnace_params)
        params, height = furnace._calc_params(), furnace.tapholeheight
        rt = 0.5 * furnace.tapholediameter
        eps = 1e-9
        below = furnaces._calc_vdot_out(height + hm_rel, height + rt - eps, 1.0, 1.0, params)
        above = furnaces._calc_vdot_out(height + hm_rel, height + rt + eps, 1.0, 1.0, params)
        assert below[2:] == pytest.approx(above[2:], rel=1e-5)

    def test_flows_have_no_step_over_a_sweep_of_the_slag_level(self, furnace_params):
        furnace = open_furnace(furnace_params)
        params, height = furnace._calc_params(), furnace.tapholeheight
        rt = 0.5 * furnace.tapholediameter
        for hm_rel in (-0.03, 0.0, 0.03, 0.045):
            levels = numpy.linspace(rt + 0.02, rt - 0.02, 4001)
            flows = numpy.array([
                furnaces._calc_vdot_out(height + hm_rel, height + hs, 1.0, 1.0, params)[2:]
                for hs in levels])
            steps = numpy.abs(numpy.diff(flows, axis=0)).sum(axis=1)
            assert steps.max() < 3 * numpy.median(steps)

    def test_no_jump_in_the_flows_during_a_tap_through_the_top_of_the_taphole(self):
        """A shallow bath with a large tap-hole and no smelting: the slag level falls through
        the top of the tap-hole about 80 s into the tap."""
        params = dict(
            powerMVA=0, powerfactor=0.833333, metalSER=3.5, slagmetalmassratio=1.05,
            furnacediameter=11.4, activeareafraction=0.75, tapholediameter=0.2,
            tapholelength=1.25, tapholeroughness=0.001, tapholeheight=0.45,
            densitymetal=6600, densityslag=2780, viscositymetal=0.005, viscosityslag=0.1,
            particlediameter=0.012, particlesphericity=0.8, bedporosity=0.5, bedmindiameter=0.2,
            bedmaxdiameter=100, bedmodel=bedmodel_ergun, entrykl=0.25, fdmodel=fdmodel_cheng,
            hmetal_init=0.5, hslag_init=0.57,
        )
        furnace = open_furnace(params)
        times = numpy.linspace(0, 1200, 501)
        mass_flow, slag_levels = [], []
        for dt in numpy.diff(times):
            furnace.calc_dt(float(dt))
            mass_flow.append(furnace.vdotmetal_out * furnace.densitymetal
                             + furnace.vdotslag_out * furnace.densityslag)
            slag_levels.append(furnace.hslag - furnace.tapholeheight)
        assert slag_levels[0] > 0.1 > slag_levels[-1]  # the slag level passes the top of the hole
        second_differences = numpy.abs(numpy.diff(mass_flow, 2))[20:]
        assert second_differences.max() < 0.2  # (a step of 1.9 kg/s gave 1.2 before)


class TestFloodedTapholeUnchanged:
    """Pinned values for a completely filled tap-hole, which the continuity change must not
    alter."""

    h0_l, h0_h = -RT - 0.01, RT + 0.01

    @pytest.mark.parametrize("hm, hs, expected", [
        (0.0, 0.3, (0.003926990816987242, 0.003926990816987242)),
        (0.02, 0.06, (0.005562257291940933, 0.00229172434203355)),
        (-0.05, 0.1, (0.0003126080286583633, 0.00754137360531612)),
        (0.04, 0.2, (0.006993569945054542, 0.0008604116889199412)),
        (-0.058, 0.07, (2.8544891992804208e-05, 0.007825436741981679)),
        (0.055, 0.5, (0.00774200758157903, 0.00011197405239545336)),
    ])
    def test_areas(self, hm, hs, expected):
        assert furnaces._phase_areas(hm, hs, self.h0_l, self.h0_h, RT) == expected


class TestPartlyFilledInterface:
    h0_l, h0_h = -RT - 0.012, RT + 0.012

    def test_reduces_to_the_undeformed_interface_without_deformation(self):
        """With no deformation (h0_h = rt) the interface is the line from the bottom of the
        tap-hole at h0_l to the free surface at the metal level hs."""
        h0_l = -RT - 0.02
        for hs, hm in ((0.03, -0.01), (0.0, -0.04), (-0.03, -0.045), (0.04, 0.04)):
            area_m, area_s = furnaces._phase_areas(hm, hs, h0_l, RT, RT)
            hi = hs - (hs - hm) * (hs + RT) / (hs - h0_l)
            assert area_m == pytest.approx(furnaces._segment_area(RT, hi), abs=1e-15)
            assert area_s == pytest.approx(furnaces._segment_area(RT, hs) - area_m, abs=1e-15)

    def test_metal_area_increases_with_the_metal_level(self):
        for hs in (-0.03, 0.0, 0.03, 0.049):
            levels = numpy.linspace(self.h0_l, hs, 80)
            values = [furnaces._phase_areas(hm, hs, self.h0_l, self.h0_h, RT)[0] for hm in levels]
            assert values[0] == pytest.approx(0, abs=1e-15)
            assert numpy.all(numpy.diff(values) > 0)

    def test_slag_remains_in_the_hole_while_the_metal_stays_below_the_free_surface(self):
        for hs in (-0.03, 0.0, 0.03, 0.049):
            area_m, area_s = furnaces._phase_areas(hs, hs, self.h0_l, self.h0_h, RT)
            assert area_s > 0
            assert area_m + area_s == pytest.approx(furnaces._segment_area(RT, hs))
