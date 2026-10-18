"""Unit tests for the tappingsim.furnaces.DCArcFurnace class."""

import math

import numpy
import pytest
from tappingsim.components.fdmodels import fdmodel_bellos
from tappingsim.furnaces import DCArcFurnace

# The tap-hole at 0.2 m (the "metal" tap-hole) and the one at 0.5 m (the "slag" tap-hole)
# both have a radius of 0.05 m. The initial levels put the first tap-hole in the range where
# both phases enter it, and the second one completely filled with slag only.
HMETAL_INIT, HSLAG_INIT = 0.22, 0.8


@pytest.fixture
def dc_params():
    return {
        "powerMW": 50,
        "metalSER": 3.5,
        "slagmetalmassratio": 1.1,
        "furnacediameter": 10,
        "activeareafraction": 0.9,
        "tapholediameters": [0.1, 0.1],
        "tapholelengths": [1, 1],
        "tapholeroughnesses": [0.001, 0.002],
        "tapholeheights": [0.2, 0.5],
        "entrykls": [0.25, 0.25],
        "fdmodels": [fdmodel_bellos, fdmodel_bellos],
        "densitymetal": 6500,
        "densityslag": 3000,
        "viscositymetal": 0.005,
        "viscosityslag": 0.1,
        "hmetal_init": HMETAL_INIT,
        "hslag_init": HSLAG_INIT,
    }


def make_furnace(dc_params, states=(False, False)):
    furnace = DCArcFurnace(**dc_params)
    furnace.tapholesopen_yn = list(states)
    return furnace


def balance_constants(furnace):
    """Factor linking a net volume flowrate to the rate of change of a level, and the volume
    flowrates of metal and slag into the furnace."""
    areaconst = 1 / (furnace.activeareafraction * 0.25 * math.pi * furnace.furnacediameter**2)
    mdotmetal_in = furnace.powerMW / (3.6 * furnace.metalSER)
    vdotmetal_in = mdotmetal_in / furnace.densitymetal
    vdotslag_in = mdotmetal_in * furnace.slagmetalmassratio / furnace.densityslag
    return areaconst, vdotmetal_in, vdotslag_in


class TestDCArcFurnaceInit:
    def test_initial_levels(self, dc_params):
        furnace = DCArcFurnace(**dc_params)
        assert furnace.hmetal == HMETAL_INIT
        assert furnace.hslag == HSLAG_INIT

    def test_tapholes_initially_closed(self, dc_params):
        furnace = DCArcFurnace(**dc_params)
        assert furnace.tapholesopen_yn == [False, False]

    def test_number_of_tapholes(self, dc_params):
        assert DCArcFurnace(**dc_params).ntapholes == 2

    def test_initial_velocities_zero(self, dc_params):
        furnace = DCArcFurnace(**dc_params)
        assert furnace.umetals == [0, 0]
        assert furnace.uslags == [0, 0]


class TestDCArcFurnaceClosedTapholes:
    def test_no_outflow(self, dc_params):
        furnace = make_furnace(dc_params)
        furnace.calc_dt(60.0)
        assert furnace.vdotmetal_outs == [0, 0]
        assert furnace.vdotslag_outs == [0, 0]

    def test_levels_rise_in_an_exact_linear_update(self, dc_params):
        furnace = make_furnace(dc_params)
        areaconst, vdotmetal_in, vdotslag_in = balance_constants(furnace)
        furnace.calc_dt(600.0)
        assert furnace.hmetal == pytest.approx(HMETAL_INIT + 600 * areaconst * vdotmetal_in,
                                               rel=1e-12)
        assert furnace.hslag == pytest.approx(
            HSLAG_INIT + 600 * areaconst * (vdotmetal_in + vdotslag_in), rel=1e-12)


# The reference values below were calculated for the furnace of the fixture with an independent
# solution of the level equations (an adaptive 8th order Runge-Kutta method with a tight
# tolerance), so these tests check accuracy. The tap-holes are flooded in all of them.

class TestDCArcFurnaceOpenTapholes:
    # open states, then the levels after one 60 s step, and the metal flows from the two
    # tap-holes followed by the slag flows from the two tap-holes, averaged over the step
    cases = [
        ([True, True], (0.213323062785, 0.774655505043),
         (8.476581949099e-03, 0.0, 8.970556934958e-03, 1.447666786876e-02)),
        ([True, False], (0.213292989996, 0.786865718104),
         (8.512010619361e-03, 0.0, 9.026977764149e-03, 0.0)),
        ([False, True], (0.220518209013, 0.789306861066),
         (0.0, 0.0, 0.0, 1.466308458538e-02)),
    ]

    @pytest.mark.parametrize("states, levels, flows", cases)
    def test_levels_after_a_step(self, dc_params, states, levels, flows):
        furnace = make_furnace(dc_params, states)
        furnace.calc_dt(60.0)
        assert furnace.hmetal == pytest.approx(levels[0], rel=1e-7)
        assert furnace.hslag == pytest.approx(levels[1], rel=1e-7)

    @pytest.mark.parametrize("states, levels, flows", cases)
    def test_mean_flows_over_a_step(self, dc_params, states, levels, flows):
        furnace = make_furnace(dc_params, states)
        furnace.calc_dt(60.0)
        assert furnace.vdotmetal_outs + furnace.vdotslag_outs == pytest.approx(
            list(flows), rel=1e-6, abs=1e-15)

    def test_the_high_taphole_carries_no_metal(self, dc_params):
        furnace = make_furnace(dc_params, [True, True])
        furnace.calc_dt(60.0)
        assert furnace.vdotmetal_outs[1] == 0
        assert furnace.vdotslag_outs[1] > 0

    def test_a_closed_taphole_carries_nothing(self, dc_params):
        furnace = make_furnace(dc_params, [False, True])
        furnace.calc_dt(60.0)
        assert furnace.vdotmetal_outs[0] == 0 and furnace.vdotslag_outs[0] == 0

    def test_tapped_masses_over_a_period(self, dc_params):
        furnace = make_furnace(dc_params, [True, True])
        metal, slag = furnace.calc_time_period(numpy.linspace(0, 600, 41))
        assert metal == pytest.approx(23668.422920, rel=1e-6)
        assert slag == pytest.approx(39123.716723, rel=1e-6)
        assert furnace.hmetal == pytest.approx(0.173668291798, rel=1e-7)
        assert furnace.hslag == pytest.approx(0.581523150055, rel=1e-7)


class TestDCArcFurnaceTimeIntegration:
    @pytest.mark.parametrize("states", [[True, False], [False, True], [True, True]])
    @pytest.mark.parametrize("dt", [1.0, 60.0, 300.0])
    def test_levels_change_by_inflow_less_reported_outflow(self, dc_params, states, dt):
        """The reported outflows are step averages, so the level changes balance exactly."""
        furnace = make_furnace(dc_params, states)
        areaconst, vdotmetal_in, vdotslag_in = balance_constants(furnace)
        furnace.calc_dt(dt)
        metal_out, slag_out = sum(furnace.vdotmetal_outs), sum(furnace.vdotslag_outs)
        dhmetal = dt * areaconst * (vdotmetal_in - metal_out)
        dhslag = dhmetal + dt * areaconst * (vdotslag_in - slag_out)
        assert furnace.hmetal - HMETAL_INIT == pytest.approx(dhmetal, abs=1e-12)
        assert furnace.hslag - HSLAG_INIT == pytest.approx(dhslag, abs=1e-12)

    def test_fourth_order_convergence(self, dc_params):
        """Halving the time step reduces the error in the levels by a large factor. The period
        is one in which the slag level stays above the top of the slag tap-hole, so that the flow
        calculation is smooth."""
        def levels_after(dt, total=600.0):
            furnace = make_furnace(dc_params, [True, True])
            for _ in range(int(round(total / dt))):
                furnace.calc_dt(dt)
            assert furnace.hslag > 0.55
            return numpy.array([furnace.hmetal, furnace.hslag])

        reference = levels_after(1.25)
        errors = [numpy.abs(levels_after(dt) - reference).max() for dt in (300.0, 150.0, 75.0)]
        assert errors[0] / errors[1] > 8
        assert errors[1] / errors[2] > 8

    def test_long_steps_remain_accurate(self, dc_params):
        coarse, fine = make_furnace(dc_params, [True, True]), make_furnace(dc_params, [True, True])
        coarse.calc_dt(300.0)
        for _ in range(300):
            fine.calc_dt(1.0)
        assert coarse.hmetal == pytest.approx(fine.hmetal, abs=1e-5)
        assert coarse.hslag == pytest.approx(fine.hslag, abs=1e-5)

    def test_calc_time_period_matches_repeated_steps(self, dc_params):
        times = numpy.linspace(0, 600, 41)
        a, b = make_furnace(dc_params, [True, True]), make_furnace(dc_params, [True, True])
        metal, slag = a.calc_time_period(times)
        metal_sum = slag_sum = 0.0
        for dt in numpy.diff(times):
            b.calc_dt(float(dt))
            metal_sum += dt * sum(b.vdotmetal_outs) * b.densitymetal
            slag_sum += dt * sum(b.vdotslag_outs) * b.densityslag
        assert metal == pytest.approx(metal_sum, rel=1e-6)
        assert slag == pytest.approx(slag_sum, rel=1e-6)
        assert a.hmetal == pytest.approx(b.hmetal, rel=1e-8)
        assert a.hslag == pytest.approx(b.hslag, rel=1e-8)
        assert a.vdotmetal_outs == pytest.approx(b.vdotmetal_outs, rel=1e-6, abs=1e-15)
        assert a.umetals == pytest.approx(b.umetals, rel=1e-6, abs=1e-12)

    def test_tapped_masses_balance_the_levels(self, dc_params):
        """The mass tapped over a period equals the mass smelted less the change stored."""
        furnace = make_furnace(dc_params, [True, True])
        areaconst, vdotmetal_in, vdotslag_in = balance_constants(furnace)
        metal, slag = furnace.calc_time_period(numpy.linspace(0, 600, 21))
        stored_metal = (furnace.hmetal - HMETAL_INIT) / areaconst * furnace.densitymetal
        smelted_metal = 600 * vdotmetal_in * furnace.densitymetal
        assert metal == pytest.approx(smelted_metal - stored_metal, rel=1e-9)
        stored_slag_top = (furnace.hslag - HSLAG_INIT) / areaconst
        stored_slag = ((stored_slag_top - (furnace.hmetal - HMETAL_INIT) / areaconst)
                       * furnace.densityslag)
        smelted_slag = 600 * vdotslag_in * furnace.densityslag
        assert slag == pytest.approx(smelted_slag - stored_slag, rel=1e-9)

    def test_calc_time_period_returns_the_masses(self, dc_params):
        """Regression test: the masses used to be accumulated with the wrong entries of the
        parameter tuple, which made this raise a TypeError."""
        furnace = make_furnace(dc_params, [True, False])
        metal, slag = furnace.calc_time_period(numpy.linspace(0, 300, 11))
        assert type(metal) is float and type(slag) is float
        assert metal > 0 and slag > 0

    def test_calc_time_period_with_a_single_time_does_nothing(self, dc_params):
        furnace = make_furnace(dc_params, [True, True])
        assert furnace.calc_time_period(numpy.array([0.0])) == (0, 0)
        assert (furnace.hmetal, furnace.hslag) == (HMETAL_INIT, HSLAG_INIT)

    def test_opening_and_closing_tapholes_between_steps_takes_effect(self, dc_params):
        furnace = make_furnace(dc_params, [True, False])
        furnace.calc_dt(60.0)
        assert furnace.vdotmetal_outs[0] > 0 and furnace.vdotslag_outs[1] == 0
        furnace.tapholesopen_yn = [False, True]
        furnace.calc_dt(60.0)
        assert furnace.vdotmetal_outs[0] == 0 and furnace.vdotslag_outs[1] > 0

    def test_velocities_are_those_at_the_start_of_the_step(self, dc_params):
        furnace = make_furnace(dc_params, [True, False])
        furnace.calc_dt(60.0)
        assert furnace.umetals[0] > 0 and furnace.uslags[0] > 0
        assert furnace.umetals[1] == 0 and furnace.uslags[1] == 0

    def test_results_are_lists_of_plain_floats(self, dc_params):
        furnace = make_furnace(dc_params, [True, True])
        furnace.calc_dt(60.0)
        assert type(furnace.hmetal) is float and type(furnace.hslag) is float
        for values in (furnace.vdotmetal_outs, furnace.vdotslag_outs):
            assert len(values) == 2
            assert all(type(value) is float for value in values)

    def test_slag_tapped_through_both_tapholes_is_the_total(self, dc_params):
        """A hole whose top is below the slag level passes slag as well as metal."""
        furnace = make_furnace(dc_params, [True, True])
        furnace.calc_dt(60.0)
        assert furnace.vdotslag_outs[0] > 0 and furnace.vdotslag_outs[1] > 0
