"""Functions for the time integration of the state of tapping units.

The state x of a unit, such as the phase volumes in a ladle or launder trough
or the bath levels in a furnace, obeys the balance::

    dx/dt = inflows + F q(x)

where q(x) are the flowrates, which depend on the state, and F is a constant
matrix of coefficients linking the flows to the rates of change of the state.

The flows over overflow outlets and through underflow channels are strongly
nonlinear, with a stiffness which can be very large, so explicit integration
is only stable for small time steps. For these the functions here provide a
linearly-implicit (Rosenbrock-Euler) scheme, which has no stability limit for
the linearised problem. For smooth flows with slow dynamics, such as the
drainage of a furnace through its tap-hole, they provide an explicit
fourth-order Runge-Kutta scheme.

"""

import numpy


def linearly_implicit_step(flows, volumes, inflows, flows_to_rates, dt,
                           eps=1e-7):
    """Advance the phase volumes over a single linearly-implicit step.

    With the flows q(V), the step is
    dV = dt*(I - dt*J)^-1 * (inflows + F*q), where F is the matrix of mass
    balance coefficients and J = F*dq/dV is the Jacobian of the volume rates,
    with dq/dV estimated by forward differences. The scheme is first order
    in time.

    Parameters
    ----------
    flows : function
        Function of the phase volumes, as a numpy array, returning the
        flowrates as a numpy array, m3/s.
    volumes : array_like of float
        The phase volumes at the start of the step, m3.
    inflows : array_like of float
        The constant rate of addition to each phase volume over the step,
        m3/s.
    flows_to_rates : array_like of float
        The matrix of mass balance coefficients, with a row for each phase
        volume and a column for each flow.
    dt : float
        Length of time step, in s.
    eps : float, optional
        The volume perturbation used to estimate the sensitivity of the
        flowrates to the volumes, m3.

    Returns
    -------
    volumes : numpy.ndarray
        The phase volumes at the end of the step, m3.
    flows : numpy.ndarray
        The effective flowrates over the step, m3/s, which are the flowrates
        that satisfy the volume balance exactly, so that the change in volumes
        is dt*(inflows + F*flows).

    """
    volumes = numpy.asarray(volumes, dtype=float)
    q = flows(volumes)
    dq_dV = numpy.empty((len(q), len(volumes)))
    for j in range(len(volumes)):
        perturbed = volumes.copy()
        perturbed[j] += eps
        dq_dV[:, j] = (flows(perturbed) - q) / eps
    dV_dt = flows_to_rates @ q + inflows
    jacobian = flows_to_rates @ dq_dV
    V_delta = dt * numpy.linalg.solve(
        numpy.eye(len(volumes)) - dt * jacobian, dV_dt)
    return volumes + V_delta, q + dq_dV @ V_delta


def _advance(flows, volumes, inflows, flows_to_rates, dt, halvings, positive,
             eps, tol=1e-12):
    """Advance over a time step, repeating the step as two half-steps if the
    result is not physical.

    A very large step relative to the stiffness of the flows can overshoot to
    negative phase volumes. If positive is True, the step is then repeated as
    two half-steps, to any depth up to the given number of halvings.
    Rounding-sized negative volumes are set to zero.

    """
    new_volumes, q = linearly_implicit_step(
        flows, volumes, inflows, flows_to_rates, dt, eps)
    if not positive:
        return new_volumes, q
    if (numpy.all(numpy.isfinite(new_volumes))
            and numpy.all(new_volumes > -tol)):
        return numpy.maximum(new_volumes, 0), q
    if halvings == 0:
        raise RuntimeError("could not advance the phase volumes over the "
                           "time step with positive volumes")
    half_volumes, q_first = _advance(
        flows, volumes, inflows, flows_to_rates, 0.5 * dt, halvings - 1,
        positive, eps)
    new_volumes, q_second = _advance(
        flows, half_volumes, inflows, flows_to_rates, 0.5 * dt, halvings - 1,
        positive, eps)
    return new_volumes, 0.5 * (q_first + q_second)


def integrate(flows, volumes, inflows, flows_to_rates, dt, max_substep=1.0,
              positive=True, halvings=12, eps=1e-7):
    """Advance the phase volumes over a time step, which is divided into
    equal sub-steps no longer than max_substep.

    Each sub-step is a linearly-implicit step (see linearly_implicit_step). If
    positive is True and a sub-step would produce negative volumes, which can
    happen for very large steps relative to the stiffness of the flows, it is
    repeated as half-steps.

    Parameters
    ----------
    flows : function
        Function of the phase volumes, as a numpy array, returning the
        flowrates as a numpy array, m3/s.
    volumes : array_like of float
        The phase volumes at the start of the time step, m3.
    inflows : array_like of float
        The constant rate of addition to each phase volume over the time
        step, m3/s.
    flows_to_rates : array_like of float
        The matrix of mass balance coefficients, with a row for each phase
        volume and a column for each flow.
    dt : float
        Length of time step, in s.
    max_substep : float, optional
        The longest sub-step, in s. The default is 1 s.
    positive : bool, optional
        Whether to require the phase volumes to stay positive. If so, a
        sub-step that gives negative volumes is repeated as half-steps, and
        rounding-sized negative volumes are set to zero. Otherwise the
        volumes are not checked. The default is True.
    halvings : int, optional
        The number of times a sub-step can be halved before giving up. The
        default of 12 limits the work for a persistently unphysical state to
        a few thousand steps per sub-step.
    eps : float, optional
        The volume perturbation used to estimate the sensitivity of the
        flowrates to the volumes, m3.

    Returns
    -------
    volumes : numpy.ndarray
        The phase volumes at the end of the time step, m3.
    flows : numpy.ndarray
        The effective flowrates averaged over the time step, m3/s, so that the
        change in volumes is dt*(inflows + F*flows).

    Raises
    ------
    RuntimeError
        If positive is True and a step with positive volumes cannot be found.

    """
    n_substeps = max(1, int(numpy.ceil(dt / max_substep)))
    volumes = numpy.asarray(volumes, dtype=float)
    total_flows = 0
    for _ in range(n_substeps):
        volumes, q = _advance(flows, volumes, inflows, flows_to_rates,
                              dt / n_substeps, halvings, positive, eps)
        total_flows = total_flows + q
    return volumes, total_flows / n_substeps


def rk4_step(flows, state, inflows, flows_to_rates, dt):
    """Advance the state over a single explicit fourth-order Runge-Kutta step.

    This is the classical Runge-Kutta scheme for dx/dt = inflows + F*q(x),
    where F is the matrix of mass balance coefficients. It is fourth order in
    time, but only stable for time steps that are short relative to the
    fastest time scale of the flows, so it should not be used for stiff flows
    (see linearly_implicit_step for those).

    Parameters
    ----------
    flows : function
        Function of the state, as a numpy array, returning the flowrates as a
        numpy array, m3/s.
    state : array_like of float
        The state at the start of the step.
    inflows : array_like of float
        The constant rate of addition to each component of the state over the
        step.
    flows_to_rates : array_like of float
        The matrix of mass balance coefficients, with a row for each component
        of the state and a column for each flow.
    dt : float
        Length of time step, in s.

    Returns
    -------
    state : numpy.ndarray
        The state at the end of the step.
    flows : numpy.ndarray
        The effective flowrates over the step, m3/s, which are the weighted
        average of the flowrates at the four stages of the step, so that the
        change in the state is dt*(inflows + F*flows).

    """
    state = numpy.asarray(state, dtype=float)
    q1 = flows(state)
    k1 = flows_to_rates @ q1 + inflows
    q2 = flows(state + 0.5 * dt * k1)
    k2 = flows_to_rates @ q2 + inflows
    q3 = flows(state + 0.5 * dt * k2)
    k3 = flows_to_rates @ q3 + inflows
    q4 = flows(state + dt * k3)
    k4 = flows_to_rates @ q4 + inflows
    return (state + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4),
            (q1 + 2 * q2 + 2 * q3 + q4) / 6)
