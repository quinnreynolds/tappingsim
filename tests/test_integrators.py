"""Unit tests for tappingsim.integrators module."""

import numpy
import pytest
from tappingsim import integrators
from tappingsim.integrators import integrate, linearly_implicit_step, rk4_step

# One volume drained by one flow
F1 = numpy.array([[-1.0]])
# Two volumes in series: the flow out of the first fills the second, which drains
F2 = numpy.array([[-1.0, 0.0], [1.0, -1.0]])


def linear_outflow(k):
    return lambda V: numpy.array([k * V[0]])


def series_flows(V):
    """Nonlinear flows: a square root law out of the first volume, a 3/2 power law out of
    the second."""
    return numpy.array([0.3 * numpy.sqrt(V[0]), 0.2 * V[1] ** 1.5])


class TestLinearlyImplicitStep:
    def test_linear_decay_is_backward_euler(self):
        V, q = linearly_implicit_step(linear_outflow(0.5), [2.0], [0.0], F1, 0.4)
        assert V[0] == pytest.approx(2.0 / (1 + 0.5 * 0.4), rel=1e-6)

    def test_effective_flows_close_the_volume_balance(self):
        volumes, inflows, dt = numpy.array([1.0, 0.5]), numpy.array([0.1, 0.0]), 0.7
        new_volumes, q = linearly_implicit_step(series_flows, volumes, inflows, F2, dt)
        assert new_volumes - volumes == pytest.approx(dt * (inflows + F2 @ q), abs=1e-12)

    def test_stable_for_very_large_steps(self):
        """An explicit step would be wildly unstable for k*dt = 1e7."""
        V, q = linearly_implicit_step(linear_outflow(1e3), [2.0], [0.0], F1, 1e4)
        assert numpy.isfinite(V[0])
        assert 0 <= V[0] < 1e-6

    def test_zero_step_leaves_volumes_unchanged(self):
        V, q = linearly_implicit_step(linear_outflow(0.5), [2.0], [0.1], F1, 0.0)
        assert V[0] == 2.0
        assert q[0] == pytest.approx(1.0)

    def test_steady_state_is_a_fixed_point(self):
        """Volumes at which the flows balance the inflows do not change, for any step."""
        V, q = linearly_implicit_step(linear_outflow(0.5), [0.4], [0.2], F1, 50.0)
        assert V[0] == pytest.approx(0.4, rel=1e-6)
        assert q[0] == pytest.approx(0.2, rel=1e-6)

    def test_inputs_are_not_modified(self):
        volumes, inflows = numpy.array([1.0, 0.5]), numpy.array([0.1, 0.0])
        linearly_implicit_step(series_flows, volumes, inflows, F2, 0.7)
        assert volumes.tolist() == [1.0, 0.5]
        assert inflows.tolist() == [0.1, 0.0]

    def test_first_order_in_time(self):
        """Halving the step halves the error against the exact exponential decay."""
        def error(n):
            V = numpy.array([1.0])
            for _ in range(n):
                V, _ = linearly_implicit_step(linear_outflow(1.0), V, [0.0], F1, 1.0 / n)
            return abs(V[0] - numpy.exp(-1.0))

        assert error(10) / error(20) == pytest.approx(2, abs=0.2)
        assert error(20) / error(40) == pytest.approx(2, abs=0.2)


class TestIntegrate:
    def test_substeps_equal_repeated_steps(self):
        flows, inflows = series_flows, numpy.array([0.1, 0.0])
        volumes = numpy.array([1.0, 0.5])
        V, q = integrate(flows, volumes, inflows, F2, 10.0, max_substep=2.0)
        expected = volumes
        for _ in range(5):
            expected, _ = linearly_implicit_step(flows, expected, inflows, F2, 2.0)
        assert V == pytest.approx(expected, abs=1e-14)

    def test_step_shorter_than_max_substep_is_a_single_step(self):
        volumes, inflows = numpy.array([1.0, 0.5]), numpy.array([0.1, 0.0])
        V, q = integrate(series_flows, volumes, inflows, F2, 0.5, max_substep=1.0)
        expected, expected_q = linearly_implicit_step(series_flows, volumes, inflows, F2, 0.5)
        assert V.tolist() == expected.tolist()
        assert q.tolist() == expected_q.tolist()

    def test_default_max_substep_is_one_second(self):
        volumes, inflows = numpy.array([1.0, 0.5]), numpy.array([0.1, 0.0])
        V, _ = integrate(series_flows, volumes, inflows, F2, 3.0)
        expected = volumes
        for _ in range(3):
            expected, _ = linearly_implicit_step(series_flows, expected, inflows, F2, 1.0)
        assert V == pytest.approx(expected, abs=1e-14)

    def test_averaged_flows_close_the_volume_balance(self):
        volumes, inflows, dt = numpy.array([1.0, 0.5]), numpy.array([0.1, 0.0]), 10.0
        V, q = integrate(series_flows, volumes, inflows, F2, dt, max_substep=2.0)
        assert V - volumes == pytest.approx(dt * (inflows + F2 @ q), abs=1e-9)

    def test_input_volumes_are_not_modified(self):
        volumes = numpy.array([1.0, 0.5])
        integrate(series_flows, volumes, numpy.array([0.1, 0.0]), F2, 10.0)
        assert volumes.tolist() == [1.0, 0.5]


class TestIntegrateUnphysicalSteps:
    @pytest.fixture
    def patched(self, monkeypatch):
        """Make steps longer than 1.5 s overshoot to negative volumes, recording the steps."""
        original = integrators.linearly_implicit_step
        steps = []

        def step(flows, volumes, inflows, flows_to_rates, dt, eps=1e-7):
            steps.append(dt)
            new_volumes, q = original(flows, volumes, inflows, flows_to_rates, dt, eps)
            if dt > 1.5:
                new_volumes = new_volumes - 10.0
            return new_volumes, q

        monkeypatch.setattr(integrators, "linearly_implicit_step", step)
        return original, steps

    def test_unphysical_step_is_repeated_as_half_steps(self, patched):
        original, steps = patched
        volumes, inflows = numpy.array([1.0]), numpy.array([0.1])
        V, q = integrate(linear_outflow(0.3), volumes, inflows, F1, 8.0, max_substep=8.0)
        assert steps[0] == 8.0
        assert min(steps) == 1.0
        expected = volumes
        for _ in range(8):
            expected, _ = original(linear_outflow(0.3), expected, inflows, F1, 1.0)
        assert V == pytest.approx(expected, abs=1e-14)
        assert V[0] > 0

    def test_failure_to_find_a_physical_step_raises(self, monkeypatch):
        def step(flows, volumes, inflows, flows_to_rates, dt, eps=1e-7):
            return numpy.array([numpy.nan]), numpy.array([0.0])

        monkeypatch.setattr(integrators, "linearly_implicit_step", step)
        with pytest.raises(RuntimeError):
            integrate(linear_outflow(0.3), [1.0], [0.1], F1, 1.0, halvings=3)

    def test_rounding_sized_negative_volumes_are_set_to_zero(self, monkeypatch):
        def step(flows, volumes, inflows, flows_to_rates, dt, eps=1e-7):
            return numpy.array([-1e-15, 2.0]), numpy.array([0.0])

        monkeypatch.setattr(integrators, "linearly_implicit_step", step)
        V, q = integrate(linear_outflow(0.3), [1.0, 1.0], [0.1, 0.0], F2, 1.0)
        assert V.tolist() == [0.0, 2.0]

    def test_positivity_not_checked_when_not_requested(self, monkeypatch):
        steps = []

        def step(flows, volumes, inflows, flows_to_rates, dt, eps=1e-7):
            steps.append(dt)
            return numpy.array([-1e-15, -0.5]), numpy.array([0.0])

        monkeypatch.setattr(integrators, "linearly_implicit_step", step)
        V, q = integrate(linear_outflow(0.3), [1.0, 1.0], [0.1, 0.0], F2, 1.0, positive=False)
        assert steps == [1.0]
        assert V.tolist() == [-1e-15, -0.5]

    def test_persistent_failure_raises_after_bounded_work(self, monkeypatch):
        steps = []

        def step(flows, volumes, inflows, flows_to_rates, dt, eps=1e-7):
            steps.append(dt)
            return numpy.array([-1.0]), numpy.array([0.0])

        monkeypatch.setattr(integrators, "linearly_implicit_step", step)
        with pytest.raises(RuntimeError):
            integrate(linear_outflow(0.3), [1.0], [0.1], F1, 1.0)
        assert len(steps) <= 100


class TestRk4Step:
    def test_effective_flows_close_the_balance(self):
        volumes, inflows, dt = numpy.array([1.0, 0.5]), numpy.array([0.1, 0.0]), 0.7
        new_volumes, q = rk4_step(series_flows, volumes, inflows, F2, dt)
        assert new_volumes - volumes == pytest.approx(dt * (inflows + F2 @ q), abs=1e-14)

    def test_exact_for_a_constant_rate(self):
        """With no dependence on the state the step is an exact linear update."""
        V, q = rk4_step(lambda x: numpy.array([0.2]), [1.0], [0.5], F1, 10.0)
        assert V[0] == pytest.approx(1.0 + 10.0 * (0.5 - 0.2))
        assert q[0] == pytest.approx(0.2)

    def test_fourth_order_in_time(self):
        """Halving the step reduces the error against the exact exponential decay by a
        factor of about 16."""
        def error(n):
            x = numpy.array([1.0])
            for _ in range(n):
                x, _ = rk4_step(linear_outflow(1.0), x, [0.0], F1, 1.0 / n)
            return abs(x[0] - numpy.exp(-1.0))

        assert error(8) / error(16) == pytest.approx(16, rel=0.1)
        assert error(16) / error(32) == pytest.approx(16, rel=0.05)

    def test_far_more_accurate_than_the_implicit_step(self):
        x_rk, x_im = numpy.array([1.0]), numpy.array([1.0])
        for _ in range(10):
            x_rk, _ = rk4_step(linear_outflow(1.0), x_rk, [0.0], F1, 0.1)
            x_im, _ = linearly_implicit_step(linear_outflow(1.0), x_im, [0.0], F1, 0.1)
        assert abs(x_rk[0] - numpy.exp(-1.0)) < 1e-3 * abs(x_im[0] - numpy.exp(-1.0))

    def test_unstable_for_steps_beyond_its_stability_limit(self):
        """Explicit schemes have a stability limit: for dx/dt = -k*x, k*dt must be below
        about 2.78. This is why rk4_step is for non-stiff flows only."""
        x = numpy.array([1.0])
        for _ in range(20):
            x, _ = rk4_step(linear_outflow(1.0), x, [0.0], F1, 4.0)
        assert abs(x[0]) > 1e3

    def test_stages_evaluate_the_flows_four_times(self):
        calls = []

        def counting(x):
            calls.append(x.copy())
            return numpy.array([0.2 * x[0]])

        rk4_step(counting, [1.0], [0.0], F1, 1.0)
        assert len(calls) == 4
        assert calls[0].tolist() == [1.0]

    def test_inputs_are_not_modified(self):
        volumes, inflows = numpy.array([1.0, 0.5]), numpy.array([0.1, 0.0])
        rk4_step(series_flows, volumes, inflows, F2, 0.7)
        assert volumes.tolist() == [1.0, 0.5]
        assert inflows.tolist() == [0.1, 0.0]

    def test_zero_step_leaves_the_state_unchanged(self):
        V, q = rk4_step(linear_outflow(0.5), [2.0], [0.1], F1, 0.0)
        assert V[0] == 2.0
        assert q[0] == pytest.approx(1.0)
