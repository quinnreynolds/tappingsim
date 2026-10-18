"""Models of the friction factor of the flow in a tap-hole channel.

The furnace models are configured with one of these functions to describe the
pressure drop due to friction in the tap-hole channel. Each takes the velocity,
the density and viscosity of the fluid, and the diameter and roughness of the
channel, and returns the Darcy friction factor.

"""

import math


def fdmodel_bellos(velocity, density, viscosity, diameter, roughness):
    """Function to calculate the friction factor of fluid flow in the tap-hole
    channel.

    Parameters
    ----------
    velocity : float
        Velocity of fluid flow in the direction of the channel length, m/s.
    density : float
        Density of the fluid, kg/m3.
    viscosity : float
        Viscosity of the fluid, Pa.s.
    diameter : float
        Diameter of tap-hole channel, m.
    roughness : float
        The roughness of the tap-hole channel surface.

    Returns
    -------
    float
        The Darcy friction factor for the tap-hole channel.

    Note
    ----
    dP/L = f*(rho*u^2)/(2*D), where u is the fluid velocity. The Bellos et al
    correlation is valid in all flow regimes including transitional flow.
    [Bellos2018]_ [Bellos2020]_

    References
    ----------
    .. [Bellos2018] V. Bellos, I. Nalbantis, G. Tsakiris. Friction Modeling of
       Flood Flow Simulations. Journal of Hydraulic Engineering, 144(12):
       04018073, 2018. doi:10.1061/(asce)hy.1943-7900.0001540.
    .. [Bellos2020] V. Bellos, I. Nalbantis, G. Tsakiris. Erratum for Friction
       Modeling of Flood Flow Simulations. Journal of Hydraulic Engineering,
       146(10): 08220005, 2020. doi:10.1061/(ASCE)HY.1943-7900.0001802.
    """
    d_over_e = diameter / roughness
    NRe = diameter * velocity * density / viscosity
    if NRe < 1:
        NRe = 1
    a = 1 / (1 + (0.000368732 * NRe) ** 8.4)
    b = 1 / (1 + (0.006666667 * NRe / d_over_e) ** 1.8)
    return ((0.015625 * NRe) ** (-a)
            * (0.75 * math.log(0.186219739 * NRe)) ** (2 * (a - 1) * b)
            * (0.88 * math.log(3.41 * d_over_e)) ** (2 * (a - 1) * (1 - b)))


def fdmodel_cheng(velocity, density, viscosity, diameter, roughness):
    """Function to calculate the friction factor of fluid flow in the tap-hole
    channel.

    Parameters
    ----------
    velocity : float
        Velocity of fluid flow in the direction of the channel length, m/s.
    density : float
        Density of the fluid, kg/m3.
    viscosity : float
        Viscosity of the fluid, Pa.s.
    diameter : float
        Diameter of tap-hole channel, m.
    roughness : float
        The roughness of the tap-hole channel surface.

    Returns
    -------
    float
        The Darcy friction factor for the tap-hole channel.

    Note
    ----
    dP/L = f*(rho*u^2)/(2*D), where u is the fluid velocity. The Cheng
    correlation is valid in all flow regimes including transitional flow.
    [Cheng2008]_

    References
    ----------
    .. [Cheng2008] N.-S. Cheng. Formulas for Friction Factor in Transitional
       Regimes. Journal of Hydraulic Engineering. 134(9): 1357-1362, 2008.
    """
    d_over_e = diameter / roughness
    NRe = diameter * velocity * density / viscosity
    if NRe < 1:
        NRe = 1
    a = 1 / (1 + (0.000367647 * NRe) ** 9)
    b = 1 / (1 + (0.00625 * NRe / d_over_e) ** 2)
    return ((0.015625 * NRe) ** (-a)
            * (1.8 * math.log10(0.147058824 * NRe)) ** (2 * (a - 1) * b)
            * (2.0 * math.log10(3.7 * d_over_e)) ** (2 * (a - 1) * (1 - b)))


def fdmodel_serghides1(velocity, density, viscosity, diameter, roughness):
    """Function to calculate the friction factor of fluid flow in the tap-hole
    channel.

    Parameters
    ----------
    velocity : float
        Velocity of fluid flow in the direction of the channel length, m/s.
    density : float
        Density of the fluid, kg/m3.
    viscosity : float
        Viscosity of the fluid, Pa.s.
    diameter : float
        Diameter of tap-hole channel, m.
    roughness : float
        The roughness of the tap-hole channel surface.

    Returns
    -------
    float
        The Darcy friction factor for the tap-hole channel.

    Note
    ----
    dP/L = f*(rho*u^2)/(2*D), where u is the fluid velocity. The Serghides
    correlation is valid for laminar and turbulent flow regimes. [Serghides1984]_

    References
    ----------
    .. [Serghides1984] T.K. Serghides. Estimate friction factor accurately.
       Chemical Engineering Journal. 91(5): 63-64, 1984.
    """
    if velocity < 1e-6:
        return 64
    e_over_d = roughness / diameter
    invNRe = viscosity / (diameter * velocity * density)
    if invNRe > 1:
        return 64
    elif invNRe > 0.000333333:
        return 64 * invNRe
    else:
        eod = 0.27027027 * e_over_d
        psi1 = -2 * math.log10(eod + 12 * invNRe)
        psi2 = -2 * math.log10(eod + 2.51 * psi1 * invNRe)
        psi3 = -2 * math.log10(eod + 2.51 * psi2 * invNRe)
        invsqrtf = psi1 - (psi2 - psi1) * (psi2 - psi1) / \
            (psi3 - 2 * psi2 + psi1)
        return 1 / (invsqrtf * invsqrtf)


def fdmodel_serghides2(velocity, density, viscosity, diameter, roughness):
    """Function to calculate the friction factor of fluid flow in the tap-hole
    channel.

    Parameters
    ----------
    velocity : float
        Velocity of fluid flow in the direction of the channel length, m/s.
    density : float
        Density of the fluid, kg/m3.
    viscosity : float
        Viscosity of the fluid, Pa.s.
    diameter : float
        Diameter of tap-hole channel, m.
    roughness : float
        The roughness of the tap-hole channel surface.

    Returns
    -------
    float
        The Darcy friction factor for the tap-hole channel.

    Note
    ----
    dP/L = f*(rho*u^2)/(2*D), where u is the fluid velocity. The Serghides
    correlation is valid laminar and turbulent flow regimes. [Serghides1984b]_

    References
    ----------
    .. [Serghides1984b] T.K. Serghides. Estimate friction factor accurately.
       Chemical Engineering Journal. 91(5): 63-64, 1984.
    """
    if velocity < 1e-6:
        return 64
    e_over_d = roughness / diameter
    invNRe = viscosity / (diameter * velocity * density)
    if invNRe > 1:
        return 64
    elif invNRe > 0.000333333:
        return 64 * invNRe
    else:
        eod = 0.27027027 * e_over_d
        psi1 = -2 * math.log10(eod + 12 * invNRe)
        psi2 = -2 * math.log10(eod + 2.51 * psi1 * invNRe)
        invsqrtf = 4.781 - (psi1 - 4.781) * (psi1 - 4.781) / \
            (psi2 - 2 * psi1 + 4.781)
        return 1 / (invsqrtf * invsqrtf)
