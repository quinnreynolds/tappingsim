"""Classes and functions to describe overflow outlets.

Molten products leave ladles, launder troughs and similar vessels over
overflow outlets, which are open-channel gravity flows over a weir crest. The
Outlet class calculates the flowrates of metal and slag over such an outlet.

This module also includes support functions for describing phase entrainment
effects in the outflow. An overflow model is any callable which takes the
distance of the slag-metal interface below the outlet crest, in m, and returns
the volume fraction of metal entrained in the slag flowing over the crest.
Models with empirical parameters can be configured with functools.partial, for
example::

    functools.partial(overflowmodel_exp, premultiplier=0.1, decay=2.0)

"""

import math
from dataclasses import dataclass
from typing import Callable

import numpy
from scipy.constants import g

# Pre-multiplier in the expression for open-channel flow over a
# broad-crested weir.
OUTLET_PREMULT = 2 / 3 * numpy.sqrt(g)


def overflowmodel_step(interfacedeltah):
    """Function to calculate metal volume fraction entrained in slag outflow
    stream as a function of the location of the slag-metal interface.

    Parameters
    ----------
    interfacedeltah : float
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


def overflowmodel_exp(interfacedeltah, premultiplier, decay):
    r"""Function to calculate metal volume fraction entrained in slag outflow
    stream as a function of the location of the slag-metal interface.

    Parameters
    ----------
    interfacedeltah : float
        The distance of the slag-metal interface from the outlet, height of
        the outlet crest minus interface position, in m (for a ladle, the
        ladle depth minus the interface position).
    premultiplier : float
        The exponential pre-multiplier, dimensionless.
    decay : float
        The exponential decay factor, 1/m.

    Returns
    -------
    float
        Metal volume fraction in slag outlet stream.

    Note
    ----
    This model describes an exponential decay of metal entrainment as a
    function of the position of the interface relative to the outlet. [1]_

    Fraction = premultiplier*exp(-decay*interfacedeltah) if interfacedeltah > 0
             = premultiplier                             if interfacedeltah < 0

    References
    ----------
    .. [1] Q.G. Reynolds, J.E. Olsen, J.D. Steenkamp, Variability in Ferroalloy
       Furnace Tapping - Insights from Modelling. Proceedings of the 16th
       International Ferro-Alloys Congress (INFACON XVI) 2021, Available at
       SSRN: doi:10.2139/ssrn.3926222 or https://ssrn.com/abstract=3926222.

    """
    if interfacedeltah > 0:
        return premultiplier * math.exp(-decay * interfacedeltah)
    else:
        return premultiplier


@dataclass
class Outlet:
    """Open-channel overflow outlet over a broad-crested weir.

    Parameters
    ----------
    crest_height : float
        The height of the outlet crest, m.
    width : float
        The width of the outlet, m.
    slag_Cd : float
        Discharge coefficient for the slag phase, dimensionless.
    metal_Cd : float
        Discharge coefficient for the metal phase, dimensionless.
    density_metal : float
        The density of the metal phase, kg/m3.
    density_slag : float
        The density of the slag phase, kg/m3.
    overflowmodel : function, optional
        Function of the distance of the slag-metal interface below the crest,
        in m, returning the volume fraction of metal entrained in the slag
        flowing over the crest. The default is no entrainment.

    Note
    ----
    The parameters are stored as attributes of the same name.

    """

    crest_height: float
    width: float
    slag_Cd: float
    metal_Cd: float
    density_metal: float
    density_slag: float
    overflowmodel: Callable[[float], float] = overflowmodel_step

    def flowrates(self, hmetal, hslag):
        """Calculate the open-channel flowrates of metal and slag over the
        outlet.

        Parameters
        ----------
        hmetal : float
            The level of metal, measured from the same datum as the crest, m.
        hslag : float
            The level of the top of the slag layer, measured from the same
            datum as the crest, m.

        Returns
        -------
        tuple of float
            Volumetric flowrates of metal and slag over the outlet, m3/s.

        Note
        ----
        There is no outflow while the top of the slag layer is below the
        crest. Otherwise each phase flows over the crest as over a
        broad-crested weir, with flowrate proportional to Cd*width*head**1.5.
        If the metal level is below the crest, only slag flows, with a head
        equal to the height of the slag surface above the crest. If the metal
        level is above the crest, the slag flows with a head equal to the
        thickness of the slag layer, and the metal flows with the head of an
        equivalent column of metal, which is the slag layer thickness scaled
        by the density ratio plus the height of the metal level above the
        crest.

        The volume fraction of metal entrained in the slag flow is calculated
        using the overflow model. It is zero when there is no metal present.

        """
        if hslag < self.crest_height:
            return 0, 0
        # TODO: need to correct slag phase density if entrained metal is
        # present?
        vfrac = (self.overflowmodel(self.crest_height - hmetal)
                 if hmetal > 0 else 0)
        if hmetal < self.crest_height:
            # slag is above outlet, metal is below outlet
            h1 = hslag - self.crest_height
            qslag = OUTLET_PREMULT * self.slag_Cd * self.width * h1**1.5
            return vfrac * qslag, (1 - vfrac) * qslag
        # slag and metal are both above outlet
        density_ratio = self.density_slag / self.density_metal
        h1 = hslag - hmetal
        h2 = hmetal - self.crest_height
        qslag = OUTLET_PREMULT * self.slag_Cd * self.width * h1**1.5
        qmetal = (
            OUTLET_PREMULT
            * self.metal_Cd
            * self.width
            * ((density_ratio * h1 + h2) ** 1.5 - (density_ratio * h1) ** 1.5)
        )
        return qmetal + vfrac * qslag, (1 - vfrac) * qslag
