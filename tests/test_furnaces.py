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
        assert furnace.vdotmetal_out == pytest.approx(0.00957043121552626, rel=1e-6)

    def test_slag_outflow_value(self, furnace):
        assert furnace.vdotslag_out == pytest.approx(0.0015383077470994438, rel=1e-6)

    def test_metal_level_value(self, furnace):
        assert furnace.hmetal == pytest.approx(0.23748715684043864, rel=1e-8)

    def test_slag_level_value(self, furnace):
        assert furnace.hslag == pytest.approx(0.38818522921477017, rel=1e-8)


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
        assert self.metaltapped == pytest.approx(24812.07235, rel=1e-5)

    def test_slag_tapped_mass(self):
        assert self.slagtapped == pytest.approx(3294.08917, rel=1e-5)

    def test_metal_tapped_positive(self):
        assert self.metaltapped > 0

    def test_slag_tapped_positive(self):
        assert self.slagtapped > 0
