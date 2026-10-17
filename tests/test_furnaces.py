"""Unit tests for tappingsim.furnaces module."""

import math
import numpy
import pytest
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
