"""Unit tests for tappingsim.tappingsystems module."""

import numpy
import pytest
from tappingsim.furnaces import SubmergedArcFurnace, bedmodel_ergun, fdmodel_cheng
from tappingsim.ladles import CylindricalLadle, overflowmodel_step
from tappingsim.launders import SimpleSisoLaunder
from tappingsim.tappingsystems import SAF, SAFWithLadles


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
            overflowmodel=overflowmodel_step, overflowconsts=[],
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
