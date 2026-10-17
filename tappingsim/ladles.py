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
    """Tapping ladle class. This version is a cylindrical ladle with an
    choice of empirical overflow models.

    Parameters
    ----------
    diameter : float
        The internal diameter of the ladle, m.
    depth : float
        The internal depth of the ladle, m.
    hmetal_init : float
        The initial level of metal in the ladle, m.
    hslag_init : float
        The initial level of slag in the ladle, m.
    overflowmodel : function
        Function to be used to model the phase overflow behaviour.
    overflowconsts : list of float
        Empirical parameters for the chosen overflow model.

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
    vdotmetal_out : float
        The current outlet volume flowrate of metal from the unit, m3/s.
    vdotslag_out : float
        The current outlet volume flowrate of slag from the unit, m3/s.

    """

    def __init__(self, diameter, depth, hmetal_init, hslag_init,
                 overflowmodel, overflowconsts):
        self.diameter = diameter
        self.depth = depth
        self.hmetal = hmetal_init
        self.hslag = hslag_init
        self.overflowmodel = lambda dh: overflowmodel(dh, *overflowconsts)

    def _calc_vdot_out(self, dt):
        """Calculate volumetric outflow rates of metal and slag from the ladle
        for the current state, and update the vdotmetal_out and vdotslag_out
        attributes.

        Overflow occurs when the total bath level exceeds the ladle depth.
        Metal entrainment in the slag outflow is calculated using the
        configured overflow model.

        Parameters
        ----------
        dt : float
            Length of time step, in s.

        Returns
        -------
        None.

        """
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
        dhmetal = dt * vdotmetal_in / xarea
        dhslag = dt * vdotslag_in / xarea
        self.hmetal += dhmetal
        self.hslag += dhmetal + dhslag
        self._calc_vdot_out(dt)
        dhmetal = -dt * self.vdotmetal_out / xarea
        dhslag = -dt * self.vdotslag_out / xarea
        self.hmetal += dhmetal
        self.hslag += dhmetal + dhslag
