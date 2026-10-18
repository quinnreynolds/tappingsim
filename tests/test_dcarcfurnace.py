"""Unit tests for the tappingsim.furnaces.DCArcFurnace class."""

import math

import numpy
import pytest
from tappingsim import furnaces
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


class TestDCTapholeParameters:
    """The parameters of the flow calculation are a snapshot for each tap-hole."""

    def test_one_snapshot_for_each_taphole(self, dc_params):
        params = make_furnace(dc_params, [True, False])._calc_params()
        assert len(params) == 2
        assert all(isinstance(p, furnaces._TapholeParams) for p in params)

    def test_values_follow_the_taphole_lists(self, dc_params):
        params = make_furnace(dc_params, [True, False])._calc_params()
        assert [p.tapholediameter for p in params] == [0.1, 0.1]
        assert [p.tapholeheight for p in params] == [0.2, 0.5]
        assert [p.tapholeroughness for p in params] == [0.001, 0.002]
        assert [p.tapholeopen_yn for p in params] == [True, False]
        assert all(p.fdmodel is fdmodel_bellos for p in params)
        assert all(p.densitymetal == 6500 and p.densityslag == 3000 for p in params)
        assert all(p.viscositymetal == 0.005 and p.viscosityslag == 0.1 for p in params)

    def test_open_bath_has_unit_porosity_and_only_the_entrance_loss(self, dc_params):
        dc_params["entrykls"] = [0.25, 0.5]
        for p, entrykl in zip(make_furnace(dc_params)._calc_params(), (0.25, 0.5)):
            assert p.bedporosity == 1.0
            assert p.b_metal == 0 and p.b_slag == 0
            assert p.a_metal == pytest.approx(0.5 * (1 + entrykl) * 6500)
            assert p.a_slag == pytest.approx(0.5 * (1 + entrykl) * 3000)

    def test_snapshot_is_not_changed_by_later_changes(self, dc_params):
        furnace = make_furnace(dc_params, [True, False])
        before = furnace._calc_params()
        furnace.tapholesopen_yn = [False, True]
        assert [p.tapholeopen_yn for p in before] == [True, False]
        assert [p.tapholeopen_yn for p in furnace._calc_params()] == [False, True]

    def test_changed_attributes_are_used_by_the_next_step(self, dc_params):
        furnace = make_furnace(dc_params, [True, False])
        furnace.calc_dt(60.0)
        before = furnace.vdotmetal_outs[0]
        furnace.entrykls[0] = 2.0
        furnace.calc_dt(60.0)
        assert furnace.vdotmetal_outs[0] < before


class TestDCFlowCalculation:
    # (metal level, slag level) and the metal velocities, slag velocities, metal flows and slag
    # flows from the two tap-holes, both open, starting the iteration from 0.5 m/s. Both
    # tap-holes are completely filled. These pinned values were calculated before the velocity
    # iteration was replaced by the shared one, which must not change them.
    pinned = [
        ((0.22, 0.8), ([1.8698590665899828, 1.250389293831332],
                       [2.783269981948702, 1.883992190130447],
                       [0.008800820584131949, 0.0],
                       [0.008759801604620154, 0.014796840059835894])),
        ((0.3, 0.9), ([2.1310514148759436, 1.4436482550725358],
                      [3.1784985213907144, 2.1763368354537493],
                      [0.014673886894439788, 0.0],
                      [0.0030775233919942547, 0.017092909534995895])),
        ((0.15, 0.7), ([1.6679050553370796, 1.0211702932754265],
                       [2.4771516603102013, 1.536619928839739],
                       [0.0031208350372109246, 0.0],
                       [0.01482047860575409, 0.012068584699506488])),
        ((0.5, 1.0), ([2.528496562421828, 1.637274918771857],
                      [3.7781406344327033, 2.46807159540749],
                      [0.019858765562828652, 0.006539708417610767],
                      [0.0, 0.009526059230394906])),
        ((0.1, 0.65), ([1.5824255288413915, 0.884545121786582],
                       [2.3474815617636215, 1.33121537918555],
                       [0.00020861847025351073, 0.0],
                       [0.018127597729595234, 0.010455341138987688])),
        ((0.25, 1.2), ([2.425701726016871, 1.9094396315591884],
                       [3.6232622711237297, 2.8761816268312312],
                       [0.01278489122305204, 0.0],
                       [0.0093602875699371, 0.02258947767310734])),
    ]

    @pytest.mark.parametrize("levels, expected", pinned)
    def test_flows_are_unchanged(self, dc_params, levels, expected):
        furnace = make_furnace(dc_params, [True, True])
        result = furnace._calc_vdot_out(*levels, [0.5, 0.5], [0.5, 0.5], furnace._calc_params())
        assert tuple(result) == tuple(expected)

    def test_each_taphole_is_calculated_independently(self, dc_params):
        """A tap-hole in a furnace with two gives the flows of a furnace with only that one."""
        both = make_furnace(dc_params, [True, True])
        umetals, uslags, vdotmetals, vdotslags = both._calc_vdot_out(
            0.3, 0.9, [0.5, 0.5], [0.5, 0.5], both._calc_params())
        for i in (0, 1):
            single_params = dict(dc_params)
            for key in ("tapholediameters", "tapholelengths", "tapholeroughnesses",
                        "tapholeheights", "entrykls", "fdmodels"):
                single_params[key] = [dc_params[key][i]]
            single = make_furnace(single_params, [True])
            expected = single._calc_vdot_out(0.3, 0.9, [0.5], [0.5], single._calc_params())
            assert (umetals[i], uslags[i], vdotmetals[i], vdotslags[i]) == tuple(
                values[0] for values in expected)

    # (metal level, slag level) in the regimes of the first tap-hole (at 0.2 m, radius 0.05 m)
    # and the second one (at 0.5 m): completely filled with both phases, completely filled with
    # slag only, partly filled with slag only, partly filled with metal in it, and with the
    # slag surface below the tap-hole.
    regimes = [(0.22, 0.8), (0.3, 0.9), (0.1, 0.52), (0.21, 0.23), (0.19, 0.24), (0.05, 0.1)]

    @pytest.mark.parametrize("levels", regimes)
    def test_each_taphole_matches_the_shared_calculation(self, dc_params, levels):
        """The flows of each tap-hole are those of the calculation that the submerged-arc
        furnace uses for its tap-hole, with the snapshot of the parameters of that tap-hole."""
        furnace = make_furnace(dc_params, [True, True])
        params = furnace._calc_params()
        umetals, uslags, vdotmetals, vdotslags = furnace._calc_vdot_out(
            *levels, [0.5, 0.5], [0.5, 0.5], params)
        for i in (0, 1):
            expected = furnaces._calc_vdot_out(*levels, 0.5, 0.5, params[i])
            assert (umetals[i], uslags[i], vdotmetals[i], vdotslags[i]) == expected

    def test_no_metal_below_the_hearth_unless_allowed(self, dc_params):
        dc_params["tapholeheights"] = [0.0, 0.5]
        furnace = make_furnace(dc_params, [True, False])
        umetals, _, vdotmetals, _ = furnace._calc_vdot_out(
            -0.02, 0.3, [0.5, 0.5], [0.5, 0.5], furnace._calc_params())
        assert umetals[0] == 0 and vdotmetals[0] == 0
        furnace.allownegativeheights_yn = True
        umetals, _, vdotmetals, _ = furnace._calc_vdot_out(
            -0.02, 0.3, [0.5, 0.5], [0.5, 0.5], furnace._calc_params())
        assert umetals[0] > 0 and vdotmetals[0] > 0

    def test_closed_taphole_has_no_velocity_or_flow(self, dc_params):
        furnace = make_furnace(dc_params, [False, True])
        umetals, uslags, vdotmetals, vdotslags = furnace._calc_vdot_out(
            0.3, 0.9, [0.5, 0.5], [0.5, 0.5], furnace._calc_params())
        assert (umetals[0], uslags[0], vdotmetals[0], vdotslags[0]) == (0, 0, 0, 0)
        assert uslags[1] > 0


class TestDCVelocitySolver:
    def test_the_shared_solver_is_used_for_each_open_taphole(self, dc_params, monkeypatch):
        calls = []
        original = furnaces._solve_velocities

        def recording(pa, umetal0, uslag0, params):
            calls.append(params.tapholeheight)
            return original(pa, umetal0, uslag0, params)

        monkeypatch.setattr(furnaces, "_solve_velocities", recording)
        furnace = make_furnace(dc_params, [False, True])
        furnace.calc_dt(60.0)
        assert calls == [0.5] * 4  # one for each of the four stages of the step

    def test_non_convergence_raises(self, dc_params, monkeypatch):
        """The iteration has a cap, so a case that does not converge ends with an error
        instead of looping for ever."""
        monkeypatch.setattr(furnaces, "_VELOCITY_TOLERANCE", 0.0)
        monkeypatch.setattr(furnaces, "_VELOCITY_MAX_ITERATIONS", 2)
        with pytest.raises(RuntimeError):
            make_furnace(dc_params, [True, False]).calc_dt(60.0)

    def test_closed_tapholes_do_not_run_the_solver(self, dc_params, monkeypatch):
        monkeypatch.setattr(furnaces, "_VELOCITY_TOLERANCE", 0.0)
        monkeypatch.setattr(furnaces, "_VELOCITY_MAX_ITERATIONS", 1)
        furnace = make_furnace(dc_params, [False, False])
        furnace.calc_dt(60.0)
        assert furnace.vdotmetal_outs == [0, 0]


class TestDCPartlyFilledTapholes:
    # open states, (metal level, slag level), and the metal velocities, slag velocities, metal
    # flows and slag flows from the two tap-holes. The slag surface is below the top of an open
    # tap-hole and there is no metal in it. These pinned values were calculated before the
    # interface and area calculations were replaced by the shared ones, which must not change
    # them.
    pinned = [
        ([False, True], (0.3, 0.52), ([0, 0.42818215700404955], [0, 0.6583848055272],
                                      [0, 0.0], [0, 0.003866231823569957])),
        ([False, True], (0.3, 0.505), ([0, 0.3797742935948482], [0, 0.576468625139181],
                                       [0, 0.0], [0, 0.0025515401960795815])),
        ([True, False], (0.05, 0.22), ([0.4434306520893389, 0], [0.6584046781471846, 0],
                                       [0.0, 0], [0.003866348521517959, 0])),
        ([True, False], (0.0, 0.2), ([0, 0], [0.5464208098727793, 0],
                                     [0, 0], [0.002145789502581136, 0])),
    ]

    @pytest.mark.parametrize("states, levels, expected", pinned)
    def test_flows_without_metal_in_the_taphole_are_unchanged(
            self, dc_params, states, levels, expected):
        furnace = make_furnace(dc_params, states)
        result = furnace._calc_vdot_out(*levels, [0.5, 0.5], [0.5, 0.5], furnace._calc_params())
        assert tuple(result) == tuple(expected)

    def test_the_shared_interface_and_area_calculations_are_used(self, dc_params, monkeypatch):
        calls = []
        for name in ("_interface_levels", "_phase_areas"):
            original = getattr(furnaces, name)

            def recording(*args, original=original, name=name):
                calls.append(name)
                return original(*args)

            monkeypatch.setattr(furnaces, name, recording)
        furnace = make_furnace(dc_params, [True, True])
        furnace.calc_dt(60.0)
        assert calls.count("_interface_levels") == 8  # two tap-holes, four stages
        assert calls.count("_phase_areas") == 8

    def test_closed_tapholes_do_not_use_them(self, dc_params, monkeypatch):
        def fail(*args):
            raise AssertionError("called for a closed tap-hole")

        monkeypatch.setattr(furnaces, "_interface_levels", fail)
        monkeypatch.setattr(furnaces, "_phase_areas", fail)
        make_furnace(dc_params, [False, False]).calc_dt(60.0)


class TestDCTapholeTopContinuity:
    """The flows must not jump when the slag surface falls through the top of a tap-hole, where
    the tap-hole changes from completely to partly filled. The interface between the phases was
    previously positioned with different formulas either side of that level, which gave a step in
    both flows when there was metal in the tap-hole. The deformation of the interface is large
    for an open bath, so the step was large (about 10% of the metal flow here)."""

    @pytest.mark.parametrize("hole", [0, 1])
    @pytest.mark.parametrize("hm_rel", [-0.045, -0.02, 0.0, 0.02, 0.04])
    def test_flows_are_continuous(self, dc_params, hole, hm_rel):
        furnace = make_furnace(dc_params, [hole == 0, hole == 1])
        params, height, rt = furnace._calc_params(), dc_params["tapholeheights"][hole], 0.05
        eps = 1e-9
        below = furnace._calc_vdot_out(height + hm_rel, height + rt - eps, [1.5, 1.5], [1.5, 1.5],
                                       params)
        above = furnace._calc_vdot_out(height + hm_rel, height + rt + eps, [1.5, 1.5], [1.5, 1.5],
                                       params)
        assert below[2][hole] == pytest.approx(above[2][hole], rel=1e-5)
        assert below[3][hole] == pytest.approx(above[3][hole], rel=1e-5)

    def test_flows_have_no_step_over_a_sweep_of_the_slag_level(self, dc_params):
        furnace = make_furnace(dc_params, [True, False])
        params, rt = furnace._calc_params(), 0.05
        for hm_rel in (-0.03, 0.0, 0.03, 0.045):
            levels = numpy.linspace(rt + 0.02, rt - 0.02, 4001)
            flows = numpy.array([
                furnace._calc_vdot_out(0.2 + hm_rel, 0.2 + hs, [1.5, 1.5], [1.5, 1.5], params)[2:]
                for hs in levels], dtype=object)
            totals = numpy.array([[m[0], sl[0]] for m, sl in flows], dtype=float)
            steps = numpy.abs(numpy.diff(totals, axis=0)).sum(axis=1)
            assert steps.max() < 3 * numpy.median(steps)

    def test_no_jump_while_the_slag_level_drains_through_the_top_of_the_taphole(self, dc_params):
        """No smelting, metal in the low tap-hole, and its top (0.25 m) crossed by the slag
        surface after about 225 s."""
        dc_params.update(powerMW=0, hmetal_init=0.21, hslag_init=0.27)
        furnace = make_furnace(dc_params, [True, False])
        mass_flow, slag_levels = [], []
        for _ in range(300):
            furnace.calc_dt(1.0)
            mass_flow.append(furnace.vdotmetal_outs[0] * furnace.densitymetal
                             + furnace.vdotslag_outs[0] * furnace.densityslag)
            slag_levels.append(furnace.hslag)
        assert slag_levels[0] > 0.25 > slag_levels[-1]
        assert numpy.abs(numpy.diff(mass_flow, 2)).max() < 0.05  # (0.41 kg/s with the step)
