"""Classes to create ladle objects.

Molten products tapped from metallurgical furnaces are often stored in ladles,
large refractory-lined containers which can then be transported on to secondary
processing steps or product separation. Multiple ladles may be connected in a
strand during tapping, to prevent spillage and to assist with separation of
product and waste phases.

"""

import math
import warnings

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
        inflows have filled the ladle for a time step, with no outflow. If the
        ladle would then be filled above the crest, only the slag level is used
        by the caller, and it need only be above the crest."""
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


class PanLadle(_OverflowLadle):
    """Tapping ladle class. This version is a pan with the shape of a truncated
    pyramid with a rectangular base and top, which is usually wider than the
    base, with an open-channel overflow outlet and a choice of empirical
    overflow models.

    Metal and slag overflow the ladle as open-channel flow over a
    broad-crested weir outlet with its crest at the ladle depth (see
    tappingsim.components.outlets.Outlet), which is usually a slot cut into the
    rim of the pan. The levels rise above the depth until the flowrate over
    the outlet matches the inflow. The metal entrained in the slag flowing
    over the outlet is described by the chosen overflow model.

    Parameters
    ----------
    base_length : float
        The internal length of the base of the pan, m.
    base_width : float
        The internal width of the base of the pan, m.
    top_length : float
        The internal length of the pan at the top, where the outlet crest is,
        m.
    top_width : float
        The internal width of the pan at the top, where the outlet crest is,
        m.
    depth : float
        The internal depth of the pan, m. This is the height of the overflow
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
    base_length : float
        The internal length of the base of the pan, m.
    base_width : float
        The internal width of the base of the pan, m.
    top_length : float
        The internal length of the pan at the top, m.
    top_width : float
        The internal width of the pan at the top, m.
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
    The length and the width of the pan change linearly with the height above
    the base, so the horizontal cross-sectional area is a quadratic function of
    the height and the volume below a level is a cubic one. The taper
    continues above the depth, where the levels rise while the pan overflows.
    A warning is given when the base is larger than the top in either
    direction, which would be unusual for a pan. If the pan narrows so much
    that it closes up above the rim, a level cannot be found for contents
    larger than the volume of the closed pan, and a RuntimeError is raised.

    The state is integrated in time as for the CylindricalLadle: with a
    linearly-implicit (Rosenbrock-Euler) scheme (see
    tappingsim.components.integrators), first order in time and stable for any
    time step. The outlet flowrates reported after calc_dt are the effective
    flowrates averaged over the time step, so that the change in the volumes
    held in the ladle is exactly the volume added by the inflows less the
    volume leaving over the outlet. If a step would produce negative phase
    volumes, it is repeated as smaller steps.

    """

    def __init__(self, base_length, base_width, top_length, top_width, depth,
                 hmetal_init, hslag_init, overflowmodel, outlet_width,
                 outlet_slag_Cd, outlet_metal_Cd, density_metal, density_slag,
                 max_substep=1.0):
        for name, value in (("base_length", base_length),
                            ("base_width", base_width),
                            ("top_length", top_length),
                            ("top_width", top_width),
                            ("depth", depth)):
            if value <= 0:
                raise ValueError(f"{name} must be positive")
        super().__init__(depth, hmetal_init, hslag_init, overflowmodel,
                         outlet_width, outlet_slag_Cd, outlet_metal_Cd,
                         density_metal, density_slag, max_substep)
        if base_length > top_length or base_width > top_width:
            warnings.warn("the base of the pan is larger than the top in at "
                          "least one direction", stacklevel=2)
        self.base_length = base_length
        self.base_width = base_width
        self.top_length = top_length
        self.top_width = top_width
        # The length and width are base + rate*h, so that the volume below
        # the height h is h*(c1 + h*(c2 + h*c3)).
        lrate = (top_length - base_length) / depth
        wrate = (top_width - base_width) / depth
        self._c1 = base_length * base_width
        self._c2 = 0.5 * (base_length * wrate + base_width * lrate)
        self._c3 = lrate * wrate / 3
        # the height at which the pan would close up, if it narrows with height
        self._hclose = min((-base / rate for base, rate in
                            ((base_length, lrate), (base_width, wrate))
                            if rate < 0), default=math.inf)
        self._vcrest = self._volume(depth)

    def _volume(self, h):
        """Volume of the pan below a given height above the base, m3."""
        return h * (self._c1 + h * (self._c2 + h * self._c3))

    def _area(self, h):
        """Horizontal cross-sectional area of the pan at a given height above
        the base, m2."""
        return self._c1 + h * (2 * self._c2 + 3 * h * self._c3)

    def _level(self, volume):
        """Height above the base of the top of the contents of the pan, m, for
        a given volume, m3.

        The volume is an increasing function of the height, which is inverted
        with Newton's method, within an interval that is narrowed in every
        iteration so that the method cannot fail to converge.

        """
        if volume <= 0:
            # an empty pan, or negative volumes from rounding errors
            return volume / self._c1
        # when the pan narrows with height, or at least does not widen, the
        # height is no larger than that of a prism with the area of the base
        hhigh = volume / self._c1
        hlow = 0.0
        if hhigh >= self._hclose:
            hhigh = self._hclose
            if volume > self._volume(hhigh):
                raise RuntimeError("the contents of the pan are larger than "
                                   "the volume of the pan, which narrows to "
                                   "a close above the rim")
        while self._volume(hhigh) < volume:
            hlow, hhigh = hhigh, min(2 * hhigh, self._hclose)
            if hhigh == hlow:
                raise RuntimeError("the contents of the pan are larger than "
                                   "the volume of the pan, which narrows to "
                                   "a close above the rim")
        h = hhigh
        for _ in range(100):
            excess = self._volume(h) - volume
            if excess == 0:
                return h
            if excess > 0:
                hhigh = h
            else:
                hlow = h
            area = self._area(h)
            if area > 0:
                step = excess / area
                if abs(step) <= 1e-14 * h:
                    return h - step
                hnext = h - step
            else:
                hnext = hlow
            if not hlow < hnext < hhigh:
                hnext = 0.5 * (hlow + hhigh)
                if hhigh - hlow <= 1e-14 * hhigh:
                    return hnext
            h = hnext
        raise RuntimeError("could not find the level of the contents of the "
                           "pan")

    def _levels_to_volumes(self, hmetal, hslag):
        vmetal = self._volume(hmetal)
        return vmetal, self._volume(hslag) - vmetal

    def _volumes_to_levels(self, vmetal, vslag):
        return self._level(vmetal), self._level(vmetal + vslag)

    def _fill(self, dt, vdotmetal_in, vdotslag_in):
        vmetal = self._volume(self.hmetal) + dt * vdotmetal_in
        vtotal = (self._volume(self.hslag)
                  + dt * (vdotmetal_in + vdotslag_in))
        if vtotal > self._vcrest:
            return self.hmetal, math.inf
        return self._level(vmetal), self._level(vtotal)
