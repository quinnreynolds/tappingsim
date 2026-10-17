"""Classes to create ladle objects.

Molten products tapped from metallurgical furnaces are often stored in ladles,
large refractory-lined containers which can then be transported on to secondary
processing steps or product separation. Multiple ladles may be connected in a
strand during tapping, to prevent spillage and to assist with separation of
product and waste phases.

"""

import math

import numpy
from scipy.constants import pi

from tappingsim.outlets import Outlet


class CylindricalLadle():
    """Tapping ladle class. This version is a cylindrical ladle with an
    open-channel overflow outlet and a choice of empirical overflow models.

    Metal and slag overflow the ladle as open-channel flow over a
    broad-crested weir outlet with its crest at the ladle depth (see
    tappingsim.outlets.Outlet). The levels rise above the depth until the
    flowrate
    over the outlet matches the inflow. The metal entrained in the slag
    flowing over the outlet is described by the chosen overflow model.

    Parameters
    ----------
    diameter : float
        The internal diameter of the ladle, m.
    depth : float
        The internal depth of the ladle, m. This is the height of the overflow
        outlet crest.
    hmetal_init : float
        The initial level of metal in the ladle, m.
    hslag_init : float
        The initial level of slag in the ladle, m.
    overflowmodel : function
        Function to be used to model the entrainment of metal in the slag
        overflow, taking the distance of the slag-metal interface below the
        crest, in m, and returning the metal volume fraction in the slag
        stream (see tappingsim.outlets).
    outlet_width : float
        The width of the overflow outlet, m.
    outlet_slag_Cd : float
        Discharge coefficient for the slag phase over the outlet,
        dimensionless.
    outlet_metal_Cd : float
        Discharge coefficient for the metal phase over the outlet,
        dimensionless.
    density_metal : float
        The density of the metal phase, kg/m3.
    density_slag : float
        The density of the slag phase, kg/m3.
    max_substep : float, optional
        The longest internal time step used by calc_dt, in s. Longer time
        steps are divided into equal sub-steps no longer than this. The
        default is 1 s.

    Attributes
    ----------
    diameter : float
        The internal diameter of the ladle, m.
    depth : float
        The internal depth of the ladle, m. This is the height of the
        overflow outlet crest.
    hmetal : float
        The current level of metal in the ladle, m.
    hslag : float
        The current level of slag in the ladle, m.
    outlet : tappingsim.outlets.Outlet
        The overflow outlet, with its width, discharge coefficients, phase
        densities and overflow model.
    max_substep : float
        The longest internal time step used by calc_dt, in s.
    vdotmetal_out : float
        The current outlet volume flowrate of metal from the unit, m3/s.
    vdotslag_out : float
        The current outlet volume flowrate of slag from the unit, m3/s.

    Note
    ----
    The state is integrated in time with a linearly-implicit (Rosenbrock-Euler)
    scheme, as used by the skimmer launder. It is first order in time and
    stable for any time step, as the flowrate over the outlet becomes stiffer
    for large heads. The outlet flowrates reported after calc_dt are the
    effective flowrates averaged over the time step, so that the change in the
    volumes held in the ladle is exactly the volume added by the inflows less
    the volume leaving over the outlet. If a step would produce negative phase
    volumes, which can happen for very wide outlets relative to the time step,
    it is repeated as smaller steps. The ladle walls are assumed to be high
    enough that the contents do not spill elsewhere.

    """

    # Volume perturbation used to estimate the sensitivity of the flows to
    # the state, m3.
    _JACOBIAN_EPS = 1e-7

    def __init__(self, diameter, depth, hmetal_init, hslag_init,
                 overflowmodel, outlet_width, outlet_slag_Cd, outlet_metal_Cd,
                 density_metal, density_slag, max_substep=1.0):
        if outlet_width <= 0:
            raise ValueError("outlet_width must be positive")
        if max_substep <= 0:
            raise ValueError("max_substep must be positive")
        self.diameter = diameter
        self.hmetal = hmetal_init
        self.hslag = hslag_init
        self.outlet = Outlet(depth, outlet_width, outlet_slag_Cd,
                             outlet_metal_Cd, density_metal, density_slag,
                             overflowmodel)
        self.max_substep = max_substep
        self.vdotmetal_out = 0
        self.vdotslag_out = 0

    @property
    def depth(self):
        """The internal depth of the ladle, m, which is the height of the
        overflow outlet crest."""
        return self.outlet.crest_height

    def _calc_flows(self, hmetal, hslag):
        """Calculate the open-channel flowrates of metal and slag over the
        outlet for given levels, without changing the state.

        Parameters
        ----------
        hmetal : float
            Level of metal in the ladle, m.
        hslag : float
            Level of the top of the slag layer in the ladle, m.

        Returns
        -------
        tuple of float
            Volumetric flowrates of metal and slag over the outlet, m3/s.

        """
        return self.outlet.flowrates(hmetal, hslag)

    def _step(self, dt, vdotmetal_in, vdotslag_in):
        """Advance the state over a single linearly-implicit (Rosenbrock-Euler)
        step.

        With the phase volumes V, the outflows q(V) and the inflows, the step
        is dV = dt*(I + dt*dq/dV)^-1 * (inflow - q), with dq/dV estimated by
        forward differences. The flowrates over the step that exactly satisfy
        the volume balance are q + dq/dV*dV.

        Parameters
        ----------
        dt : float
            Length of time step, in s.
        vdotmetal_in : float
            Flowrate of metal into the ladle, m3/s.
        vdotslag_in : float
            Flowrate of slag into the ladle, m3/s.

        Returns
        -------
        numpy.ndarray
            Effective flowrates of metal and slag over the outlet, m3/s.

        """
        xarea = 0.25 * pi * self.diameter**2

        def flows(V):
            hmetal = V[0] / xarea
            return numpy.array(self._calc_flows(hmetal, hmetal + V[1] / xarea))

        V = numpy.array([xarea * self.hmetal,
                         xarea * (self.hslag - self.hmetal)])
        q = flows(V)
        dq_dV = numpy.empty((2, 2))
        for j in range(2):
            V_perturbed = V.copy()
            V_perturbed[j] += self._JACOBIAN_EPS
            dq_dV[:, j] = (flows(V_perturbed) - q) / self._JACOBIAN_EPS
        inflow = numpy.array([vdotmetal_in, vdotslag_in])
        V_delta = dt * numpy.linalg.solve(numpy.eye(2) + dt * dq_dV,
                                          inflow - q)
        V_new = V + V_delta
        self.hmetal = V_new[0] / xarea
        self.hslag = self.hmetal + V_new[1] / xarea
        return q + dq_dV @ V_delta

    def _advance(self, dt, vdotmetal_in, vdotslag_in, halvings=30):
        """Advance the state over a time step, repeating the step as two
        half-steps if the result is not physical.

        A very large step relative to the stiffness of the flow over the
        outlet can overshoot to negative phase volumes. The step is then
        repeated as two half-steps, to any depth up to the given number of
        halvings. Rounding-sized negative volumes are set to zero.

        Parameters
        ----------
        dt : float
            Length of time step, in s.
        vdotmetal_in : float
            Flowrate of metal into the ladle, m3/s.
        vdotslag_in : float
            Flowrate of slag into the ladle, m3/s.
        halvings : int
            The number of times the step can be halved before giving up.

        Returns
        -------
        numpy.ndarray
            Effective flowrates of metal and slag over the outlet, m3/s.

        """
        hmetal_0, hslag_0 = self.hmetal, self.hslag
        flows = self._step(dt, vdotmetal_in, vdotslag_in)
        tol = 1e-12
        if (numpy.isfinite(self.hslag) and self.hmetal > -tol
                and self.hslag - self.hmetal > -tol):
            self.hmetal = max(self.hmetal, 0)
            self.hslag = max(self.hslag, self.hmetal)
            return flows
        if halvings == 0:
            raise RuntimeError("could not advance the ladle state over the "
                               "time step with positive volumes")
        self.hmetal, self.hslag = hmetal_0, hslag_0
        return 0.5 * (
            self._advance(0.5 * dt, vdotmetal_in, vdotslag_in, halvings - 1)
            + self._advance(0.5 * dt, vdotmetal_in, vdotslag_in, halvings - 1))

    def calc_dt(self, dt, vdotmetal_in, vdotslag_in):
        """
        Integrate object state forward over a single time step, and update
        output flowrates.

        The time step is divided into equal sub-steps no longer than
        max_substep, each integrated with a linearly-implicit scheme. The
        inflow rates are held constant over the time step. The outlet
        flowrates are set to the average of the flowrates over the sub-steps.
        A ladle that stays below the crest over the whole time step is filled
        exactly in a single step.

        Parameters
        ----------
        dt : float
            Length of time step, in s.
        vdotmetal_in : float
            Flowrate of metal into the ladle, m3/s.
        vdotslag_in : float
            Flowrate of slag into the ladle, m3/s.

        Returns
        -------
        None.

        """
        xarea = 0.25 * pi * self.diameter**2
        dhslag_noflow = dt * (vdotmetal_in + vdotslag_in) / xarea
        if (self.hslag <= self.depth
                and self.hslag + dhslag_noflow <= self.depth):
            # the levels only rise and stay below the crest: there is no
            # outflow over the time step, and the filling is exact
            self.hmetal += dt * vdotmetal_in / xarea
            self.hslag += dhslag_noflow
            self.vdotmetal_out = 0
            self.vdotslag_out = 0
            return
        n_substeps = max(1, int(math.ceil(dt / self.max_substep)))
        flows = sum(
            self._advance(dt / n_substeps, vdotmetal_in, vdotslag_in)
            for _ in range(n_substeps)) / n_substeps
        self.vdotmetal_out, self.vdotslag_out = flows
