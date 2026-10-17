"""Classes to create ladle objects.

Molten products tapped from metallurgical furnaces are often stored in ladles,
large refractory-lined containers which can then be transported on to secondary
processing steps or product separation. Multiple ladles may be connected in a
strand during tapping, to prevent spillage and to assist with separation of
product and waste phases.

This module also includes support functions for describing phase entrainment
effects in the ladle outflow, and for calculating open-channel flow over an
overflow outlet. These are also used by the launder classes.

"""

import math

import numpy
from scipy.constants import g, pi

# Pre-multiplier in the expression for open-channel flow over a
# broad-crested weir.
OUTLET_PREMULT = 2 / 3 * numpy.sqrt(g)


def overflowmodel_step(interfacedeltah, *consts):
    """Function to calculate metal volume fraction entrained in slag outflow
    stream as a function of the location of the slag-metal interface.

    Parameters
    ----------
    interfacedeltah : float
        Not used in this overflow model.
    consts : list of float
        Not used in this overflow model.

    Returns
    -------
    float
        Metal volume fraction in slag outlet stream.

    Note
    ----
    This model describes perfectly separated phase model with no metal
    entrainment into the slag phase.

    Fraction = 0

    """
    return 0


def overflowmodel_exp(interfacedeltah, *consts):
    r"""Function to calculate metal volume fraction entrained in slag outflow
    stream as a function of the location of the slag-metal interface.

    Parameters
    ----------
    interfacedeltah : float
        The distance of the slag-metal interface from the outlet, height of
        the outlet crest minus interface position, in m (for a ladle, the
        ladle depth minus the interface position).
    consts : list of float
        The exponential pre-multiplier and exponential factor, dimensionless
        and 1/m respectively.

    Returns
    -------
    float
        Metal volume fraction in slag outlet stream.

    Note
    ----
    This model describes an exponential decay of metal entrainment as a
    function of the position of the interface relative to the outlet. [1]_

    Fraction = consts[0]*exp(-consts[1]*interfacedeltah) if interfacedeltah > 0
             = consts[0]                                 if interfacedeltah < 0

    References
    ----------
    .. [1] Q.G. Reynolds, J.E. Olsen, J.D. Steenkamp, Variability in Ferroalloy
       Furnace Tapping - Insights from Modelling. Proceedings of the 16th
       International Ferro-Alloys Congress (INFACON XVI) 2021, Available at
       SSRN: doi:10.2139/ssrn.3926222 or https://ssrn.com/abstract=3926222.

    """
    if interfacedeltah > 0:
        return consts[0] * math.exp(-consts[1] * interfacedeltah)
    else:
        return consts[0]


def overflow_flowrates(hmetal, hslag, crest_height, width, slag_Cd, metal_Cd,
                       density_ratio, overflowmodel):
    """Function to calculate the open-channel flowrates of metal and slag over
    an overflow outlet.

    Parameters
    ----------
    hmetal : float
        The level of metal, measured from the same datum as the crest, m.
    hslag : float
        The level of the top of the slag layer, measured from the same datum
        as the crest, m.
    crest_height : float
        The height of the outlet crest, m.
    width : float
        The width of the outlet, m.
    slag_Cd : float
        Discharge coefficient for the slag phase, dimensionless.
    metal_Cd : float
        Discharge coefficient for the metal phase, dimensionless.
    density_ratio : float
        The density of the slag phase divided by that of the metal phase.
    overflowmodel : function
        Function of the distance of the slag-metal interface below the crest,
        in m, returning the volume fraction of metal entrained in the slag
        flowing over the crest.

    Returns
    -------
    tuple of float
        Volumetric flowrates of metal and slag over the outlet, m3/s.

    Note
    ----
    There is no outflow while the top of the slag layer is below the crest.
    Otherwise each phase flows over the crest as over a broad-crested weir,
    with flowrate proportional to Cd*width*head**1.5. If the metal level is
    below the crest, only slag flows, with a head equal to the height of the
    slag surface above the crest. If the metal level is above the crest, the
    slag flows with a head equal to the thickness of the slag layer, and the
    metal flows with the head of an equivalent column of metal, which is the
    slag layer thickness scaled by the density ratio plus the height of the
    metal level above the crest.

    The volume fraction of metal entrained in the slag flow is calculated
    using the overflow model. It is zero when there is no metal present.

    """
    if hslag < crest_height:
        return 0, 0
    # TODO: need to correct slag phase density if entrained metal is present?
    vfrac = overflowmodel(crest_height - hmetal) if hmetal > 0 else 0
    if hmetal < crest_height:
        # slag is above outlet, metal is below outlet
        h1 = hslag - crest_height
        qslag = OUTLET_PREMULT * slag_Cd * width * h1**1.5
        return vfrac * qslag, (1 - vfrac) * qslag
    # slag and metal are both above outlet
    h1 = hslag - hmetal
    h2 = hmetal - crest_height
    qslag = OUTLET_PREMULT * slag_Cd * width * h1**1.5
    qmetal = (
        OUTLET_PREMULT
        * metal_Cd
        * width
        * ((density_ratio * h1 + h2) ** 1.5 - (density_ratio * h1) ** 1.5)
    )
    return qmetal + vfrac * qslag, (1 - vfrac) * qslag


class CylindricalLadle():
    """Tapping ladle class. This version is a cylindrical ladle with a
    choice of empirical overflow models.

    By default, the overflow is instantaneous: any volume above the ladle
    depth is removed within the time step, so the levels never exceed the
    depth at the end of a step. Alternatively, if the outlet_width parameter
    is provided, the overflow is calculated as open-channel flow over a
    broad-crested weir outlet with its crest at the ladle depth (see
    overflow_flowrates). The levels then rise above the depth until the
    flowrate over the outlet matches the inflow.

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
        Function to be used to model the phase overflow behaviour.
    overflowconsts : list of float
        Empirical parameters for the chosen overflow model.
    outlet_width : float, optional
        The width of the overflow outlet, m. If provided, the open-channel
        overflow calculation is used, and outlet_slag_Cd, outlet_metal_Cd,
        density_metal and density_slag must also be provided.
    outlet_slag_Cd : float, optional
        Discharge coefficient for the slag phase over the outlet,
        dimensionless.
    outlet_metal_Cd : float, optional
        Discharge coefficient for the metal phase over the outlet,
        dimensionless.
    density_metal : float, optional
        The density of the metal phase, kg/m3.
    density_slag : float, optional
        The density of the slag phase, kg/m3.
    max_substep : float, optional
        For the open-channel overflow calculation, the longest internal time
        step used by calc_dt, in s. Longer time steps are divided into equal
        sub-steps no longer than this. The default is 1 s.

    Attributes
    ----------
    diameter : float
        The internal diameter of the ladle, m.
    depth : float
        The internal depth of the ladle, m.
    hmetal : float
        The current level of metal in the ladle, m.
    hslag : float
        The current level of slag in the ladle, m.
    overflowmodel : function
        Reference to internal lambda function created for overflow modelling.
    outlet_width : float or None
        The width of the overflow outlet, m, or None if the instantaneous
        overflow calculation is used.
    outlet_slag_Cd : float or None
        Discharge coefficient for the slag phase over the outlet.
    outlet_metal_Cd : float or None
        Discharge coefficient for the metal phase over the outlet.
    density_metal : float or None
        The density of the metal phase, kg/m3.
    density_slag : float or None
        The density of the slag phase, kg/m3.
    max_substep : float
        The longest internal time step used by calc_dt for the open-channel
        overflow calculation, in s.
    vdotmetal_out : float
        The current outlet volume flowrate of metal from the unit, m3/s.
    vdotslag_out : float
        The current outlet volume flowrate of slag from the unit, m3/s.

    Note
    ----
    The open-channel overflow calculation integrates the state with a
    linearly-implicit (Rosenbrock-Euler) scheme, as used by the skimmer
    launder. It is first order in time and stable for any time step, as the
    flowrate over the outlet becomes stiffer for large heads. The outlet
    flowrates reported after calc_dt are the effective flowrates averaged over
    the time step, so that the change in the volumes held in the ladle is
    exactly the volume added by the inflows less the volume leaving over the
    outlet. If a step would produce negative phase volumes, which can happen
    for very wide outlets relative to the time step, it is repeated as
    smaller steps. The ladle walls are assumed to be high enough that the
    contents do not spill elsewhere.

    """

    # Volume perturbation used to estimate the sensitivity of the flows to
    # the state, m3.
    _JACOBIAN_EPS = 1e-7

    def __init__(self, diameter, depth, hmetal_init, hslag_init,
                 overflowmodel, overflowconsts, outlet_width=None,
                 outlet_slag_Cd=None, outlet_metal_Cd=None,
                 density_metal=None, density_slag=None, max_substep=1.0):
        open_channel_args = (outlet_slag_Cd, outlet_metal_Cd, density_metal,
                             density_slag)
        if outlet_width is None:
            if any(arg is not None for arg in open_channel_args):
                raise ValueError("outlet_width is required for the "
                                 "open-channel overflow calculation")
        else:
            if any(arg is None for arg in open_channel_args):
                raise ValueError("outlet_slag_Cd, outlet_metal_Cd, "
                                 "density_metal and density_slag are "
                                 "required with outlet_width")
            if outlet_width <= 0:
                raise ValueError("outlet_width must be positive")
        if max_substep <= 0:
            raise ValueError("max_substep must be positive")
        self.diameter = diameter
        self.depth = depth
        self.hmetal = hmetal_init
        self.hslag = hslag_init
        self.overflowmodel = lambda dh: overflowmodel(dh, *overflowconsts)
        self.outlet_width = outlet_width
        self.outlet_slag_Cd = outlet_slag_Cd
        self.outlet_metal_Cd = outlet_metal_Cd
        self.density_metal = density_metal
        self.density_slag = density_slag
        self.max_substep = max_substep
        self.vdotmetal_out = 0
        self.vdotslag_out = 0

    def _calc_vdot_out(self, dt):
        """Calculate volumetric outflow rates of metal and slag from the ladle
        for the current state, and update the vdotmetal_out and vdotslag_out
        attributes.

        Overflow occurs when the total bath level exceeds the ladle depth.
        Metal entrainment in the slag outflow is calculated using the
        configured overflow model. If the open-channel overflow calculation
        is not used, all of the volume above the ladle depth is removed over
        the time step, otherwise the flowrates are calculated by
        overflow_flowrates.

        Parameters
        ----------
        dt : float
            Length of time step, in s.

        Returns
        -------
        None.

        """
        if self.outlet_width is not None:
            self.vdotmetal_out, self.vdotslag_out = self._calc_flows(
                self.hmetal, self.hslag)
            return
        xarea = 0.25 * pi * self.diameter**2
        if self.hslag < self.depth:
            self.vdotmetal_out = 0
            self.vdotslag_out = 0
        else:
            vfrac = self.overflowmodel(self.depth - self.hmetal)
            if self.hmetal < self.depth:
                vslagout = (self.hslag - self.depth) * xarea
                self.vdotmetal_out = vslagout * vfrac / dt
                self.vdotslag_out = vslagout * (1 - vfrac) / dt
            else:
                vslagout = (self.hslag - self.hmetal) * xarea
                self.vdotmetal_out = ((self.hmetal - self.depth) * xarea
                                      + vslagout * vfrac) / dt
                self.vdotslag_out = vslagout * (1 - vfrac) / dt

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
        return overflow_flowrates(
            hmetal, hslag, self.depth, self.outlet_width, self.outlet_slag_Cd,
            self.outlet_metal_Cd, self.density_slag / self.density_metal,
            self.overflowmodel)

    def _step_open_channel(self, dt, vdotmetal_in, vdotslag_in):
        """Advance the state over a single linearly-implicit (Rosenbrock-Euler)
        step with open-channel overflow.

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

    def _advance_open_channel(self, dt, vdotmetal_in, vdotslag_in,
                              halvings=30):
        """Advance the state over a time step with open-channel overflow,
        repeating the step as two half-steps if the result is not physical.

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
        flows = self._step_open_channel(dt, vdotmetal_in, vdotslag_in)
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
            self._advance_open_channel(0.5 * dt, vdotmetal_in, vdotslag_in,
                                       halvings - 1)
            + self._advance_open_channel(0.5 * dt, vdotmetal_in, vdotslag_in,
                                         halvings - 1))

    def calc_dt(self, dt, vdotmetal_in, vdotslag_in):
        """
        Integrate object state forward over a single time step, and update
        output flowrates.

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
        if self.outlet_width is not None:
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
                self._advance_open_channel(dt / n_substeps, vdotmetal_in,
                                           vdotslag_in)
                for _ in range(n_substeps)) / n_substeps
            self.vdotmetal_out, self.vdotslag_out = flows
            return
        dhmetal = dt * vdotmetal_in / xarea
        dhslag = dt * vdotslag_in / xarea
        self.hmetal += dhmetal
        self.hslag += dhmetal + dhslag
        self._calc_vdot_out(dt)
        dhmetal = -dt * self.vdotmetal_out / xarea
        dhslag = -dt * self.vdotslag_out / xarea
        self.hmetal += dhmetal
        self.hslag += dhmetal + dhslag
