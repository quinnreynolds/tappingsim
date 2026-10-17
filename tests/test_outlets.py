"""Unit tests for tappingsim.outlets module."""

import functools
import math

import pytest
from scipy.constants import g
from tappingsim.outlets import OUTLET_PREMULT, Outlet, overflowmodel_exp, overflowmodel_step

K = 2 / 3 * math.sqrt(g)
RHO_METAL, RHO_SLAG = 7000.0, 3000.0
CREST = 1.0


def fraction_model(dh, fraction):
    """Entrainment model with a constant metal fraction, for testing."""
    return fraction


# ---------------------------------------------------------------------------
# Overflow models
# ---------------------------------------------------------------------------

class TestOutletPremult:
    def test_value(self):
        assert OUTLET_PREMULT == pytest.approx(K)


class TestOverflowmodelStep:
    def test_returns_zero(self):
        assert overflowmodel_step(0.5) == 0

    def test_returns_zero_negative_delta(self):
        assert overflowmodel_step(-0.5) == 0


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

    def test_configured_with_partial(self):
        """A model with parameters is a callable of the interface distance alone."""
        model = functools.partial(overflowmodel_exp, premultiplier=0.1, decay=2.0)
        assert model(0.5) == pytest.approx(0.036787944117144235, rel=1e-8)
        assert model(-0.5) == pytest.approx(0.1)


# ---------------------------------------------------------------------------
# Outlet
# ---------------------------------------------------------------------------

def make_outlet(model=overflowmodel_step, crest_height=CREST):
    return Outlet(crest_height, 0.4, 0.6, 0.55, RHO_METAL, RHO_SLAG, model)


class TestOutlet:
    ratio = RHO_SLAG / RHO_METAL

    def test_parameters_stored(self):
        outlet = make_outlet()
        assert outlet.crest_height == CREST
        assert outlet.width == 0.4
        assert outlet.slag_Cd == 0.6
        assert outlet.metal_Cd == 0.55
        assert outlet.density_metal == RHO_METAL
        assert outlet.density_slag == RHO_SLAG
        assert outlet.overflowmodel is overflowmodel_step

    def test_default_model_is_no_entrainment(self):
        outlet = Outlet(CREST, 0.4, 0.6, 0.55, RHO_METAL, RHO_SLAG)
        assert outlet.overflowmodel is overflowmodel_step
        assert outlet.flowrates(0.5, 1.1)[0] == 0

    def test_no_flow_below_crest(self):
        assert make_outlet().flowrates(0.5, 0.99) == (0, 0)

    def test_no_flow_at_crest(self):
        assert make_outlet().flowrates(0.5, CREST) == (0, 0)

    def test_slag_only_flow(self):
        qmetal, qslag = make_outlet().flowrates(0.5, 1.1)
        assert qmetal == 0
        assert qslag == pytest.approx(K * 0.6 * 0.4 * 0.1**1.5)

    def test_metal_and_slag_flow(self):
        qmetal, qslag = make_outlet().flowrates(1.05, 1.25)
        assert qslag == pytest.approx(K * 0.6 * 0.4 * 0.2**1.5)
        expected = K * 0.55 * 0.4 * ((self.ratio * 0.2 + 0.05) ** 1.5 - (self.ratio * 0.2) ** 1.5)
        assert qmetal == pytest.approx(expected)

    def test_metal_only_flow(self):
        """No slag layer: metal flows with its own head above the crest."""
        qmetal, qslag = make_outlet().flowrates(1.1, 1.1)
        assert qslag == 0
        assert qmetal == pytest.approx(K * 0.55 * 0.4 * 0.1**1.5)

    def test_flow_depends_only_on_levels_relative_to_crest(self):
        low, high = make_outlet(crest_height=1.0), make_outlet(crest_height=3.0)
        assert high.flowrates(3.05, 3.25) == pytest.approx(low.flowrates(1.05, 1.25))

    def test_entrained_metal_in_slag_only_flow(self):
        model = functools.partial(overflowmodel_exp, premultiplier=0.2, decay=2.0)
        qmetal, qslag = make_outlet(model).flowrates(0.5, 1.1)
        vfrac = 0.2 * math.exp(-2.0 * 0.5)
        total = K * 0.6 * 0.4 * 0.1**1.5
        assert qmetal == pytest.approx(vfrac * total)
        assert qslag == pytest.approx((1 - vfrac) * total)

    def test_entrained_metal_adds_to_metal_flow(self):
        model = functools.partial(overflowmodel_exp, premultiplier=0.2, decay=2.0)
        plain = make_outlet().flowrates(1.05, 1.25)
        qmetal, qslag = make_outlet(model).flowrates(1.05, 1.25)
        qslag_total = K * 0.6 * 0.4 * 0.2**1.5
        assert qmetal == pytest.approx(plain[0] + 0.2 * qslag_total)
        assert qslag == pytest.approx(0.8 * qslag_total)

    def test_model_receives_distance_of_interface_below_crest(self):
        distances = []

        def model(distance):
            distances.append(distance)
            return 0.0

        make_outlet(model).flowrates(0.7, 1.1)
        assert distances == [pytest.approx(0.3)]

    def test_no_entrainment_without_metal(self):
        """There is no metal to entrain, even if the model returns a nonzero fraction."""
        model = functools.partial(fraction_model, fraction=0.3)
        qmetal, qslag = make_outlet(model).flowrates(0.0, 1.1)
        assert qmetal == 0
        assert qslag == pytest.approx(K * 0.6 * 0.4 * 0.1**1.5)

    def test_total_flow_independent_of_entrainment(self):
        """Entrainment redistributes the slag stream between phases, not its total."""
        model = functools.partial(fraction_model, fraction=0.3)
        plain = make_outlet().flowrates(0.5, 1.1)
        entrained = make_outlet(model).flowrates(0.5, 1.1)
        assert sum(entrained) == pytest.approx(sum(plain))

    def test_changing_a_parameter_changes_the_flow(self):
        outlet = make_outlet()
        narrow = outlet.flowrates(0.5, 1.1)[1]
        outlet.width = 0.8
        assert outlet.flowrates(0.5, 1.1)[1] == pytest.approx(2 * narrow)
