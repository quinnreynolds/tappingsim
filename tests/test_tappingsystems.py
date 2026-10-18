"""Unit tests for tappingsim.tappingsystems module."""

import numpy
import pytest
from tappingsim.components.bedmodels import bedmodel_ergun
from tappingsim.components.fdmodels import fdmodel_cheng
from tappingsim.furnaces import DCArcFurnace, SubmergedArcFurnace
from tappingsim.ladles import CylindricalLadle
from tappingsim.components.outlets import overflowmodel_step
from tappingsim.launders import SimpleSisoLaunder
from tappingsim.tappingsystems import DCAF, DCAFWithLadles, SAF, SAFWithLadles


# ---------------------------------------------------------------------------
# Shared fixtures
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
        'tapholediameter': 0.075,
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
        'bedmindiameter': 0.15,
        'bedmaxdiameter': 20,
        'bedmodel': bedmodel_ergun,
        'entrykl': 0.25,
        'fdmodel': fdmodel_cheng,
        'hmetal_init': 0.25,
        'hslag_init': 0.4,
    }


@pytest.fixture
def saf_system(furnace_params):
    furnace = SubmergedArcFurnace(**furnace_params)
    return SAF(furnace=furnace)


@pytest.fixture
def saf_with_ladles_system(furnace_params):
    furnace = SubmergedArcFurnace(**furnace_params)
    launder = SimpleSisoLaunder()
    ladles = [
        CylindricalLadle(
            diameter=1.8, depth=2.0, hmetal_init=0, hslag_init=0,
            overflowmodel=overflowmodel_step,
            outlet_width=0.3, outlet_slag_Cd=0.6, outlet_metal_Cd=0.6,
            density_metal=7000, density_slag=3000,
        )
        for _ in range(3)
    ]
    return SAFWithLadles(furnace=furnace, launder=launder, ladles=ladles)


# ---------------------------------------------------------------------------
# SAF
# ---------------------------------------------------------------------------

class TestSAFInit:
    def test_totalisers_start_at_zero(self, saf_system):
        assert saf_system.timetotaliser == 0
        assert saf_system.powertotaliserkWh == 0
        assert saf_system.metalmasstotaliser == 0
        assert saf_system.slagmasstotaliser == 0


class TestSAFTapholeControl:
    def test_open_taphole(self, saf_system):
        saf_system.open_taphole()
        assert saf_system.furnace.tapholeopen_yn is True

    def test_close_taphole(self, saf_system):
        saf_system.open_taphole()
        saf_system.close_taphole()
        assert saf_system.furnace.tapholeopen_yn is False


class TestSAFTotalisers:
    def test_time_totaliser_increments(self, saf_system):
        saf_system.calc_dt(60)
        assert saf_system.timetotaliser == pytest.approx(60)

    def test_power_totaliser_increments(self, saf_system):
        saf_system.calc_dt(60)
        assert saf_system.powertotaliserkWh > 0

    def test_mass_totalisers_zero_when_closed(self, saf_system):
        saf_system.close_taphole()
        saf_system.calc_dt(60)
        assert saf_system.metalmasstotaliser == 0
        assert saf_system.slagmasstotaliser == 0

    def test_mass_totalisers_positive_when_open(self, saf_system):
        saf_system.open_taphole()
        for _ in range(5):
            saf_system.calc_dt(60)
        assert saf_system.metalmasstotaliser > 0
        assert saf_system.slagmasstotaliser > 0

    def test_reset_time_totaliser(self, saf_system):
        saf_system.calc_dt(60)
        saf_system.reset_time_totaliser()
        assert saf_system.timetotaliser == 0

    def test_reset_power_totaliser(self, saf_system):
        saf_system.calc_dt(60)
        saf_system.reset_power_totaliser()
        assert saf_system.powertotaliserkWh == 0

    def test_reset_mass_totaliser(self, saf_system):
        saf_system.open_taphole()
        for _ in range(5):
            saf_system.calc_dt(60)
        saf_system.reset_mass_totaliser()
        assert saf_system.metalmasstotaliser == 0
        assert saf_system.slagmasstotaliser == 0

    def test_reset_all_totalisers(self, saf_system):
        saf_system.open_taphole()
        for _ in range(5):
            saf_system.calc_dt(60)
        saf_system.reset_all_totalisers()
        assert saf_system.timetotaliser == 0
        assert saf_system.powertotaliserkWh == 0
        assert saf_system.metalmasstotaliser == 0
        assert saf_system.slagmasstotaliser == 0


# ---------------------------------------------------------------------------
# SAFWithLadles
# ---------------------------------------------------------------------------

class TestSAFWithLadlesInit:
    def test_has_three_ladles(self, saf_with_ladles_system):
        assert len(saf_with_ladles_system.ladles) == 3

    def test_totalisers_start_at_zero(self, saf_with_ladles_system):
        assert saf_with_ladles_system.timetotaliser == 0
        assert saf_with_ladles_system.metalmasstotaliser == 0
        assert saf_with_ladles_system.slagmasstotaliser == 0


class TestSAFWithLadlesEmptyLadles:
    def test_empty_ladles_resets_levels(self, saf_with_ladles_system):
        saf_with_ladles_system.open_taphole()
        for _ in range(10):
            saf_with_ladles_system.calc_dt(60)
        saf_with_ladles_system.empty_ladles()
        for ldl in saf_with_ladles_system.ladles:
            assert ldl.hmetal == 0
            assert ldl.hslag == 0


class TestSAFWithLadlesLaunderFlow:
    def test_ladles_fill_when_taphole_open(self, saf_with_ladles_system):
        """Material should flow from furnace through launder into first ladle."""
        saf_with_ladles_system.open_taphole()
        for _ in range(20):
            saf_with_ladles_system.calc_dt(60)
        assert saf_with_ladles_system.ladles[0].hslag > 0

    def test_ladles_do_not_fill_when_taphole_closed(self, saf_with_ladles_system):
        saf_with_ladles_system.close_taphole()
        for _ in range(10):
            saf_with_ladles_system.calc_dt(60)
        for ldl in saf_with_ladles_system.ladles:
            assert ldl.hslag == 0


class TestSAFWithLadlesLadleMasses:
    def test_returns_correct_list_lengths(self, saf_with_ladles_system):
        metal_masses, slag_masses = saf_with_ladles_system.ladle_masses()
        assert len(metal_masses) == 3
        assert len(slag_masses) == 3

    def test_zero_masses_when_empty(self, saf_with_ladles_system):
        metal_masses, slag_masses = saf_with_ladles_system.ladle_masses()
        assert all(m == 0 for m in metal_masses)
        assert all(m == 0 for m in slag_masses)

    def test_positive_masses_after_tapping(self, saf_with_ladles_system):
        saf_with_ladles_system.open_taphole()
        for _ in range(20):
            saf_with_ladles_system.calc_dt(60)
        metal_masses, slag_masses = saf_with_ladles_system.ladle_masses()
        assert metal_masses[0] > 0
        assert slag_masses[0] > 0


class TestSAFWithLadlesMultiTapCycle:
    def test_ten_tap_cycles_produce_material(self, saf_with_ladles_system, furnace_params):
        """Integration test: 10 tap cycles should produce metal and slag."""
        timeopen = numpy.linspace(0, 2400, 300)
        timeclosed = numpy.linspace(0, 10800, 100)
        total_metal, total_slag = 0, 0

        for _ in range(10):
            saf_with_ladles_system.reset_all_totalisers()
            saf_with_ladles_system.open_taphole()
            for dt in timeopen[1:] - timeopen[:-1]:
                saf_with_ladles_system.calc_dt(dt=dt)
            saf_with_ladles_system.close_taphole()
            total_metal += saf_with_ladles_system.metalmasstotaliser
            total_slag += saf_with_ladles_system.slagmasstotaliser
            saf_with_ladles_system.empty_ladles()
            for dt in timeclosed[1:] - timeclosed[:-1]:
                saf_with_ladles_system.calc_dt(dt=dt)

        assert total_metal > 0
        assert total_slag > 0
        # Slag-to-metal ratio should be in the right ballpark
        assert 0.5 < total_slag / total_metal < 2.0


class TestSAFWithOpenChannelLadles:
    def test_ladle_chain_conserves_volume(self, furnace_params):
        """Everything entering the first ladle is either stored in the ladles or has
        left over the outlet of the last ladle."""
        furnace = SubmergedArcFurnace(**furnace_params)
        ladles = [
            CylindricalLadle(
                diameter=1.8, depth=0.5, hmetal_init=0, hslag_init=0,
                overflowmodel=overflowmodel_step,
                outlet_width=0.3, outlet_slag_Cd=0.6, outlet_metal_Cd=0.6,
                density_metal=7000, density_slag=3000,
            )
            for _ in range(2)
        ]
        system = SAFWithLadles(furnace=furnace, launder=SimpleSisoLaunder(), ladles=ladles)
        system.open_taphole()
        area = numpy.pi * 0.25 * 1.8**2
        metal_in = slag_in = metal_lost = slag_lost = 0
        for _ in range(120):
            system.calc_dt(60)
            metal_in += 60 * furnace.vdotmetal_out
            slag_in += 60 * furnace.vdotslag_out
            metal_lost += 60 * ladles[-1].vdotmetal_out
            slag_lost += 60 * ladles[-1].vdotslag_out
        assert ladles[-1].vdotslag_out + ladles[-1].vdotmetal_out > 0
        metal_stored = sum(area * ldl.hmetal for ldl in ladles)
        slag_stored = sum(area * (ldl.hslag - ldl.hmetal) for ldl in ladles)
        assert metal_stored + metal_lost == pytest.approx(metal_in, rel=1e-9)
        assert slag_stored + slag_lost == pytest.approx(slag_in, rel=1e-9)


# ---------------------------------------------------------------------------
# DCAF and DCAFWithLadles
# ---------------------------------------------------------------------------

DC_AREA = 0.9 * 0.25 * numpy.pi * 10**2


@pytest.fixture
def dc_furnace_params():
    return {
        'powerMW': 50,
        'metalSER': 3.5,
        'slagmetalmassratio': 1.1,
        'furnacediameter': 10,
        'activeareafraction': 0.9,
        'tapholediameters': [0.1, 0.1],
        'tapholelengths': [1, 1],
        'tapholeroughnesses': [0.001, 0.002],
        'tapholeheights': [0.2, 0.5],
        'entrykls': [0.25, 0.25],
        'fdmodels': [fdmodel_cheng, fdmodel_cheng],
        'densitymetal': 6500,
        'densityslag': 3000,
        'viscositymetal': 0.005,
        'viscosityslag': 0.1,
        'hmetal_init': 0.22,
        'hslag_init': 0.8,
    }


def make_ladle(depth=0.5, density_metal=6500, density_slag=3000):
    return CylindricalLadle(
        diameter=1.8, depth=depth, hmetal_init=0, hslag_init=0,
        overflowmodel=overflowmodel_step,
        outlet_width=0.3, outlet_slag_Cd=0.6, outlet_metal_Cd=0.6,
        density_metal=density_metal, density_slag=density_slag,
    )


@pytest.fixture
def dcaf_system(dc_furnace_params):
    return DCAF(DCArcFurnace(**dc_furnace_params))


@pytest.fixture
def dcaf_with_ladles_system(dc_furnace_params):
    furnace = DCArcFurnace(**dc_furnace_params)
    launders = [SimpleSisoLaunder(), SimpleSisoLaunder()]
    ladles = [[make_ladle(), make_ladle()], [make_ladle()]]
    return DCAFWithLadles(furnace, launders, ladles)


class TestDCAFInit:
    def test_totalisers_start_at_zero(self, dcaf_system):
        assert dcaf_system.timetotaliser == 0
        assert dcaf_system.powertotaliserkWh == 0
        assert dcaf_system.metalmasstotalisers == [0, 0]
        assert dcaf_system.slagmasstotalisers == [0, 0]
        assert dcaf_system.metalmasstotaliser == 0
        assert dcaf_system.slagmasstotaliser == 0

    def test_holds_furnace(self, dc_furnace_params):
        furnace = DCArcFurnace(**dc_furnace_params)
        assert DCAF(furnace).furnace is furnace


class TestDCAFTapholeControl:
    def test_open_and_close_by_index(self, dcaf_system):
        dcaf_system.open_taphole(1)
        assert dcaf_system.furnace.tapholesopen_yn == [False, True]
        dcaf_system.open_taphole(0)
        assert dcaf_system.furnace.tapholesopen_yn == [True, True]
        dcaf_system.close_taphole(1)
        assert dcaf_system.furnace.tapholesopen_yn == [True, False]

    @pytest.mark.parametrize('index', [2, -1, 1.0, '0', None, True])
    def test_bad_index_rejected(self, dcaf_system, index):
        with pytest.raises(IndexError, match='tap-hole index'):
            dcaf_system.open_taphole(index)
        with pytest.raises(IndexError, match='tap-hole index'):
            dcaf_system.close_taphole(index)
        assert dcaf_system.furnace.tapholesopen_yn == [False, False]

    def test_numpy_integer_index_accepted(self, dcaf_system):
        dcaf_system.open_taphole(numpy.int64(1))
        assert dcaf_system.furnace.tapholesopen_yn == [False, True]


class TestDCAFTotalisers:
    def test_time_and_power_totalisers(self, dcaf_system):
        dcaf_system.calc_dt(60)
        dcaf_system.calc_dt(30)
        assert dcaf_system.timetotaliser == pytest.approx(90)
        assert dcaf_system.powertotaliserkWh == pytest.approx(50 * 90 / 3.6)

    def test_mass_totalisers_zero_when_closed(self, dcaf_system):
        dcaf_system.calc_dt(60)
        assert dcaf_system.metalmasstotalisers == [0, 0]
        assert dcaf_system.slagmasstotalisers == [0, 0]

    def test_only_open_tapholes_accumulate(self, dcaf_system):
        dcaf_system.open_taphole(1)
        for _ in range(5):
            dcaf_system.calc_dt(60)
        assert dcaf_system.metalmasstotalisers[0] == 0
        assert dcaf_system.slagmasstotalisers[0] == 0
        assert dcaf_system.slagmasstotalisers[1] > 0

    def test_tapholes_carry_mixtures(self, dcaf_system):
        """The lower tap-hole sits in the metal and the slag layers, so it passes
        both phases, while the upper one is above the metal."""
        dcaf_system.open_taphole(0)
        dcaf_system.open_taphole(1)
        for _ in range(5):
            dcaf_system.calc_dt(60)
        assert dcaf_system.metalmasstotalisers[0] > 0
        assert dcaf_system.slagmasstotalisers[0] > 0
        assert dcaf_system.metalmasstotalisers[1] == 0
        assert dcaf_system.slagmasstotalisers[1] > 0

    def test_totals_are_sums_over_tapholes(self, dcaf_system):
        dcaf_system.open_taphole(0)
        dcaf_system.open_taphole(1)
        for _ in range(5):
            dcaf_system.calc_dt(60)
        assert dcaf_system.metalmasstotaliser == sum(dcaf_system.metalmasstotalisers)
        assert dcaf_system.slagmasstotaliser == sum(dcaf_system.slagmasstotalisers)

    def test_totalisers_match_furnace_inventory(self, dcaf_system):
        """The mass tapped is the mass supplied less the change in the furnace
        contents, as the totalisers use the step-averaged flows."""
        furnace = dcaf_system.furnace
        hm0, hs0 = furnace.hmetal, furnace.hslag
        dcaf_system.open_taphole(0)
        dcaf_system.open_taphole(1)
        total = 600
        for _ in range(10):
            dcaf_system.calc_dt(total / 10)
        mdotmetal_in = furnace.powerMW / (3.6 * furnace.metalSER)
        metal_in = total * mdotmetal_in
        slag_in = metal_in * furnace.slagmetalmassratio
        metal_held = DC_AREA * (furnace.hmetal - hm0) * furnace.densitymetal
        slag_held = (DC_AREA * ((furnace.hslag - furnace.hmetal) - (hs0 - hm0))
                     * furnace.densityslag)
        assert dcaf_system.metalmasstotaliser == pytest.approx(metal_in - metal_held, rel=1e-9)
        assert dcaf_system.slagmasstotaliser == pytest.approx(slag_in - slag_held, rel=1e-9)

    def test_resets(self, dcaf_system):
        dcaf_system.open_taphole(0)
        dcaf_system.open_taphole(1)
        for _ in range(5):
            dcaf_system.calc_dt(60)
        dcaf_system.reset_time_totaliser()
        assert dcaf_system.timetotaliser == 0
        assert dcaf_system.powertotaliserkWh > 0
        dcaf_system.reset_power_totaliser()
        assert dcaf_system.powertotaliserkWh == 0
        assert dcaf_system.metalmasstotaliser > 0
        dcaf_system.reset_mass_totaliser()
        assert dcaf_system.metalmasstotalisers == [0, 0]
        assert dcaf_system.slagmasstotalisers == [0, 0]
        dcaf_system.calc_dt(60)
        dcaf_system.reset_all_totalisers()
        assert dcaf_system.timetotaliser == 0
        assert dcaf_system.powertotaliserkWh == 0
        assert dcaf_system.metalmasstotaliser == 0

    def test_reset_keeps_a_list_per_taphole(self, dcaf_system):
        dcaf_system.reset_mass_totaliser()
        dcaf_system.open_taphole(1)
        dcaf_system.calc_dt(60)
        assert dcaf_system.slagmasstotalisers[1] > 0

    def test_furnace_state_matches_the_furnace_alone(self, dc_furnace_params):
        system = DCAF(DCArcFurnace(**dc_furnace_params))
        alone = DCArcFurnace(**dc_furnace_params)
        system.open_taphole(0)
        alone.tapholesopen_yn = [True, False]
        for _ in range(10):
            system.calc_dt(30)
            alone.calc_dt(30)
        assert (system.furnace.hmetal, system.furnace.hslag) == (alone.hmetal, alone.hslag)


class TestDCAFWithLadlesInit:
    def test_holds_units(self, dc_furnace_params):
        furnace = DCArcFurnace(**dc_furnace_params)
        launders = [SimpleSisoLaunder(), SimpleSisoLaunder()]
        ladles = [[make_ladle(), make_ladle()], [make_ladle()]]
        system = DCAFWithLadles(furnace, launders, ladles)
        assert system.furnace is furnace
        assert system.launders == launders
        assert system.ladles == ladles
        assert system.metalmasstotalisers == [0, 0]

    def test_none_launder_and_missing_chains_allowed(self, dc_furnace_params):
        furnace = DCArcFurnace(**dc_furnace_params)
        system = DCAFWithLadles(furnace, [None, SimpleSisoLaunder()], [None, []])
        assert system.ladles == [[], []]

    def test_wrong_lengths_rejected(self, dc_furnace_params):
        furnace = DCArcFurnace(**dc_furnace_params)
        with pytest.raises(ValueError, match='launders and ladles'):
            DCAFWithLadles(furnace, [SimpleSisoLaunder()], [[], []])
        with pytest.raises(ValueError, match='launders and ladles'):
            DCAFWithLadles(furnace, [None, None], [[]])

    def test_ladles_without_launder_rejected(self, dc_furnace_params):
        furnace = DCArcFurnace(**dc_furnace_params)
        with pytest.raises(ValueError, match='tap-hole 1 has ladles but no launder'):
            DCAFWithLadles(furnace, [None, None], [[], [make_ladle()]])

    def test_shared_units_rejected(self, dc_furnace_params):
        furnace = DCArcFurnace(**dc_furnace_params)
        ladle, launder = make_ladle(), SimpleSisoLaunder()
        with pytest.raises(ValueError, match='only one tap-hole'):
            DCAFWithLadles(furnace, [SimpleSisoLaunder(), SimpleSisoLaunder()],
                           [[ladle], [ladle]])
        with pytest.raises(ValueError, match='only one tap-hole'):
            DCAFWithLadles(furnace, [launder, launder], [[], []])


class TestDCAFWithLadlesFlow:
    def test_chains_fill_only_when_their_taphole_is_open(self, dcaf_with_ladles_system):
        system = dcaf_with_ladles_system
        system.open_taphole(1)
        for _ in range(10):
            system.calc_dt(60)
        assert all(ldl.hslag == 0 for ldl in system.ladles[0])
        assert system.launders[0].vdotslag_out == 0
        assert system.ladles[1][0].hslag > 0

    def test_slag_taphole_chain_gets_no_metal(self, dcaf_with_ladles_system):
        system = dcaf_with_ladles_system
        system.open_taphole(1)
        for _ in range(10):
            system.calc_dt(60)
        assert system.ladles[1][0].hmetal == 0

    def test_lower_taphole_chain_gets_both_phases(self, dcaf_with_ladles_system):
        system = dcaf_with_ladles_system
        system.open_taphole(0)
        for _ in range(10):
            system.calc_dt(60)
        ladle = system.ladles[0][0]
        assert ladle.hmetal > 0
        assert ladle.hslag > ladle.hmetal

    def test_closed_taphole_chain_does_not_move(self, dcaf_with_ladles_system):
        system = dcaf_with_ladles_system
        system.open_taphole(0)
        for _ in range(5):
            system.calc_dt(60)
        system.close_taphole(0)
        before = [(ldl.hmetal, ldl.hslag) for ldl in system.ladles[0]]
        for _ in range(5):
            system.calc_dt(60)
        assert [(ldl.hmetal, ldl.hslag) for ldl in system.ladles[0]] == before

    def test_chain_without_launder_only_totalises(self, dc_furnace_params):
        furnace = DCArcFurnace(**dc_furnace_params)
        system = DCAFWithLadles(furnace, [None, None], [[], []])
        system.open_taphole(0)
        system.calc_dt(60)
        assert system.metalmasstotalisers[0] > 0

    def test_furnace_and_totalisers_match_dcaf(self, dc_furnace_params, dcaf_with_ladles_system):
        plain = DCAF(DCArcFurnace(**dc_furnace_params))
        system = dcaf_with_ladles_system
        for s in (plain, system):
            s.open_taphole(0)
            s.open_taphole(1)
        for _ in range(10):
            plain.calc_dt(30)
            system.calc_dt(30)
        assert (system.furnace.hmetal, system.furnace.hslag) == (
            plain.furnace.hmetal, plain.furnace.hslag)
        assert system.metalmasstotalisers == plain.metalmasstotalisers
        assert system.slagmasstotalisers == plain.slagmasstotalisers

    def test_ladle_chains_conserve_volume(self, dc_furnace_params):
        """For each tap-hole, everything entering the first ladle is either stored in
        its ladles or has left over the outlet of the last one."""
        furnace = DCArcFurnace(**dc_furnace_params)
        ladles = [[make_ladle(), make_ladle()], [make_ladle(), make_ladle()]]
        system = DCAFWithLadles(furnace, [SimpleSisoLaunder(), SimpleSisoLaunder()], ladles)
        system.open_taphole(0)
        system.open_taphole(1)
        area = numpy.pi * 0.25 * 1.8**2
        metal_in, slag_in = [0, 0], [0, 0]
        metal_lost, slag_lost = [0, 0], [0, 0]
        for _ in range(120):
            system.calc_dt(60)
            for i in range(2):
                metal_in[i] += 60 * furnace.vdotmetal_outs[i]
                slag_in[i] += 60 * furnace.vdotslag_outs[i]
                metal_lost[i] += 60 * ladles[i][-1].vdotmetal_out
                slag_lost[i] += 60 * ladles[i][-1].vdotslag_out
        for i in range(2):
            assert ladles[i][-1].vdotslag_out + ladles[i][-1].vdotmetal_out > 0
            metal_stored = sum(area * ldl.hmetal for ldl in ladles[i])
            slag_stored = sum(area * (ldl.hslag - ldl.hmetal) for ldl in ladles[i])
            assert metal_stored + metal_lost[i] == pytest.approx(metal_in[i], rel=1e-9)
            assert slag_stored + slag_lost[i] == pytest.approx(slag_in[i], rel=1e-9)


class TestDCAFWithLadlesLadleMasses:
    def test_zero_when_empty(self, dcaf_with_ladles_system):
        metalmasses, slagmasses = dcaf_with_ladles_system.ladle_masses()
        assert metalmasses == [[0, 0], [0]]
        assert slagmasses == [[0, 0], [0]]

    def test_masses_follow_the_ladle_levels(self, dcaf_with_ladles_system):
        system = dcaf_with_ladles_system
        ladle = system.ladles[0][0]
        ladle.hmetal, ladle.hslag = 0.1, 0.3
        area = numpy.pi * 0.25 * 1.8**2
        metalmasses, slagmasses = system.ladle_masses()
        assert metalmasses[0][0] == pytest.approx(area * 0.1 * 6500)
        assert slagmasses[0][0] == pytest.approx(area * 0.2 * 3000)

    def test_uses_the_ladle_densities(self, dc_furnace_params):
        furnace = DCArcFurnace(**dc_furnace_params)
        ladle = make_ladle(density_metal=7000, density_slag=2800)
        system = DCAFWithLadles(furnace, [SimpleSisoLaunder(), None], [[ladle], []])
        ladle.hmetal, ladle.hslag = 0.1, 0.3
        area = numpy.pi * 0.25 * 1.8**2
        metalmasses, slagmasses = system.ladle_masses()
        assert metalmasses == [[pytest.approx(area * 0.1 * 7000)], []]
        assert slagmasses == [[pytest.approx(area * 0.2 * 2800)], []]

    def test_positive_after_tapping(self, dcaf_with_ladles_system):
        system = dcaf_with_ladles_system
        system.open_taphole(0)
        system.open_taphole(1)
        for _ in range(10):
            system.calc_dt(60)
        metalmasses, slagmasses = system.ladle_masses()
        assert metalmasses[0][0] > 0
        assert slagmasses[1][0] > 0

    def test_empty_ladles(self, dcaf_with_ladles_system):
        system = dcaf_with_ladles_system
        system.open_taphole(0)
        system.open_taphole(1)
        for _ in range(10):
            system.calc_dt(60)
        system.empty_ladles()
        assert all(ldl.hmetal == 0 and ldl.hslag == 0
                   for chain in system.ladles for ldl in chain)
