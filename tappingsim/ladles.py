"""Classes to create ladle objects.

Molten products tapped from metallurgical furnaces are often stored in ladles,
large refractory-lined containers which can then be transported on to secondary
processing steps or product separation. Multiple ladles may be connected in a
strand during tapping, to prevent spillage and to assist with separation of
product and waste phases.

"""

import numpy
from scipy.constants import pi

from tappingsim.components.integrators import integrate
from tappingsim.components.outlets import Outlet


class _OverflowLadle():
    """Base class for ladles with an open-channel overflow outlet.

    A ladle class derived from this one describes the geometry of the vessel
    only, with the methods that convert between the levels of the phases and
    the volumes held (_levels_to_volumes and _volumes_to_levels) and the one
    that fills the vessel when nothing is overflowing (_fill). The overflow
    outlet, the time integration and the reported flowrates are common to all
    of the ladles, and are described in the derived classes.

    """

    # Mass balance coefficients linking the flows over the outlet to the rates
    # of change of the metal and slag volumes.
    _FLOWS_TO_VOLUME_RATES = -numpy.eye(2)

    def __init__(self, depth, hmetal_init, hslag_init, overflowmodel,
                 outlet_width, outlet_slag_Cd, outlet_metal_Cd, density_metal,
                 density_slag, max_substep):
        if outlet_width <= 0:
            raise ValueError("outlet_width must be positive")
        if max_substep <= 0:
            raise ValueError("max_substep must be positive")
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
        overflow outlet crest. This is the depth given when the ladle was
        created, and it cannot be changed."""
        return self.outlet.crest_height

    def _levels_to_volumes(self, hmetal, hslag):
        """Volumes of metal and slag held in the ladle, m3, for given levels
        of the top of the metal and of the slag, m."""
        raise NotImplementedError

    def _volumes_to_levels(self, vmetal, vslag):
        """Levels of the top of the metal and of the slag, m, for given
        volumes of metal and slag held in the ladle, m3."""
        raise NotImplementedError

    def _fill(self, dt, vdotmetal_in, vdotslag_in):
        """Levels of the top of the metal and of the slag, m, after the
        inflows have filled the ladle for a time step, with no outflow."""
        raise NotImplementedError

    def _calc_flows(self, V):
        """Calculate the open-channel flowrates of metal and slag over the
        outlet for given phase volumes.

        Parameters
        ----------
        V : array_like of float
            Volumes of metal and slag in the ladle, m3.

        Returns
        -------
        numpy.ndarray
            Volumetric flowrates of metal and slag over the outlet, m3/s.

        """
        hmetal, hslag = self._volumes_to_levels(V[0], V[1])
        return numpy.array(self.outlet.flowrates(hmetal, hslag))

    def calc_dt(self, dt, vdotmetal_in, vdotslag_in):
        """
        Integrate object state forward over a single time step, and update
        output flowrates.

        The time step is divided into equal sub-steps no longer than
        max_substep, each integrated with a linearly-implicit scheme (see
        tappingsim.components.integrators.integrate). The inflow rates are
        held constant over the time step. The outlet flowrates are set to the
        average of the flowrates over the sub-steps. A ladle that stays below
        the crest over the whole time step is filled exactly in a single step.

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
        if self.hslag <= self.depth:
            hmetal, hslag = self._fill(dt, vdotmetal_in, vdotslag_in)
            if hslag <= self.depth:
                # the levels only rise and stay below the crest: there is no
                # outflow over the time step, and the filling is exact
                self.hmetal, self.hslag = hmetal, hslag
                self.vdotmetal_out = 0
                self.vdotslag_out = 0
                return
        V, flows = integrate(
            self._calc_flows,
            numpy.array(self._levels_to_volumes(self.hmetal, self.hslag)),
            numpy.array([vdotmetal_in, vdotslag_in]),
            self._FLOWS_TO_VOLUME_RATES,
            dt,
            self.max_substep,
        )
        self.hmetal, self.hslag = self._volumes_to_levels(V[0], V[1])
        self.vdotmetal_out, self.vdotslag_out = flows


class CylindricalLadle(_OverflowLadle):
    """Tapping ladle class. This version is a cylindrical ladle with an
    open-channel overflow outlet and a choice of empirical overflow models.

    Metal and slag overflow the ladle as open-channel flow over a
    broad-crested weir outlet with its crest at the ladle depth (see
    tappingsim.components.outlets.Outlet). The levels rise above the depth
    until the flowrate over the outlet matches the inflow. The metal entrained
    in the slag flowing over the outlet is described by the chosen overflow
    model.

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
        stream (see tappingsim.components.outlets).
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
    hmetal : float
        The current level of metal in the ladle, m.
    hslag : float
        The current level of slag in the ladle, m.
    outlet : tappingsim.components.outlets.Outlet
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
    scheme (see tappingsim.components.integrators), as used by the skimmer
    launder. It is first order in time and stable for any time step, as the
    flowrate over the outlet becomes stiffer for large heads. The outlet
    flowrates reported after calc_dt are the effective flowrates averaged over
    the time step, so that the change in the volumes held in the ladle is
    exactly the volume added by the inflows less the volume leaving over the
    outlet. If a step would produce negative phase volumes, which can happen
    for very wide outlets relative to the time step, it is repeated as smaller
    steps. The ladle walls are assumed to be high enough that the contents do
    not spill elsewhere.

    """

    def __init__(self, diameter, depth, hmetal_init, hslag_init,
                 overflowmodel, outlet_width, outlet_slag_Cd, outlet_metal_Cd,
                 density_metal, density_slag, max_substep=1.0):
        super().__init__(depth, hmetal_init, hslag_init, overflowmodel,
                         outlet_width, outlet_slag_Cd, outlet_metal_Cd,
                         density_metal, density_slag, max_substep)
        self.diameter = diameter

    def _levels_to_volumes(self, hmetal, hslag):
        xarea = 0.25 * pi * self.diameter**2
        return xarea * hmetal, xarea * (hslag - hmetal)

    def _volumes_to_levels(self, vmetal, vslag):
        xarea = 0.25 * pi * self.diameter**2
        hmetal = vmetal / xarea
        return hmetal, hmetal + vslag / xarea

    def _fill(self, dt, vdotmetal_in, vdotslag_in):
        xarea = 0.25 * pi * self.diameter**2
        return (self.hmetal + dt * vdotmetal_in / xarea,
                self.hslag + dt * (vdotmetal_in + vdotslag_in) / xarea)

