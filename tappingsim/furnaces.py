"""Classes to create furnace objects.

Metallurgical furnaces typically produce one or more immiscible molten phase,
and do so continuously provided the furnace is supplied with power and feed
material. The molten products are removed from the vessel by tapping through
one or more tap-holes.

This module also includes support functions for describing flow drag effects in
packed beds and tap-hole channels.
"""

import math
from collections import namedtuple

import numpy
from scipy.constants import g, pi

from tappingsim.integrators import rk4_step


def bedmodel_kozenycarman(tapholediameter, bedmindiameter, bedmaxdiameter,
                          bedparticlediameter, bedparticlesphericity,
                          bedporosity, viscosity, density):
    """Function to calculate pressure drop through packed burden as flow enters
    the tap-hole.

    Parameters
    ----------
    tapholediameter : float
        Diameter of tap-hole channel, m.
    bedmindiameter : float
        Diameter of cavity in burden in front of tap-hole, m (set to zero
        for no cavity).
    bedmaxdiameter : float
        Extent of burden from tap-hole, m (set to a high value for entire
        furnace).
    bedparticlediameter : float
        Diameter of the constituent particles of the burden layer, m.
    bedparticlesphericity : float
        Sphericity of the constituent particles of the burden layer.
    bedporosity : float
        Porosity of the burden layer.
    viscosity : float
        Viscosity of the fluid, Pa.s.
    density : float
        Density of the fluid, kg/m3.

    Returns
    -------
    tuple of float
        A and B constants in pressure drop expression.

    Note
    ----
    dP = A*u^2 + B*u, where u is the fluid velocity. The Kozeny-Carman
    correlation is valid for laminar flow conditions only. [Kozeny1927]_
    [Carman1937]_

    References
    ----------
    .. [Kozeny1927] J. Kozeny. Ueber kapillare Leitung des Wassers im Boden.
       Sitzungsber Akad. Wiss., Wien, 136(2a): 271-306, 1927.
    .. [Carman1937] P.C. Carman. Fluid flow through granular beds. Transactions,
       Institution of Chemical Engineers, London, 15: 150-166, 1937.

    """
    rt, rmin = 0.5 * tapholediameter, 0.5 * bedmindiameter
    eff_d = bedparticlediameter * bedparticlesphericity
    rrmax = rt / (0.5 * bedmaxdiameter)
    if rmin > rt:
        rrmin = rt / rmin
        bconst = (180 * viscosity * rt * (1 - bedporosity) * (1 - bedporosity)
                  / (2 * eff_d * eff_d * bedporosity * bedporosity * bedporosity)
                  * (rrmin - rrmax))
    else:
        rrmin = rmin / rt
        bconst = (180 * viscosity * rt * (1 - bedporosity) * (1 - bedporosity)
                  / (2 * eff_d * eff_d * bedporosity * bedporosity * bedporosity)
                  * (2 - rrmin - rrmax))
    return 0, bconst


def bedmodel_ergun(tapholediameter, bedmindiameter, bedmaxdiameter,
                   bedparticlediameter, bedparticlesphericity, bedporosity,
                   viscosity, density):
    """Function to calculate pressure drop through packed burden as flow enters
    the tap-hole.

    Parameters
    ----------
    tapholediameter : float
        Diameter of tap-hole channel, m.
    bedmindiameter : float
        Diameter of cavity in burden in front of tap-hole, m (set to zero
        for no cavity).
    bedmaxdiameter : float
        Extent of burden from tap-hole, m (set to a high value for entire
        furnace).
    bedparticlediameter : float
        Diameter of the constituent particles of the burden layer, m.
    bedparticlesphericity : float
        Sphericity of the constituent particles of the burden layer.
    bedporosity : float
        Porosity of the burden layer.
    viscosity : float
        Viscosity of the fluid, Pa.s.
    density : float
        Density of the fluid, kg/m3.

    Returns
    -------
    tuple of float
        A and B constants in pressure drop expression.

    Note
    ----
    dP = A*u^2 + B*u, where u is the fluid velocity. The Ergun correlation is
    valid for both laminar and turbulent flow regimes. [Ergun1952]_

    References
    ----------
    .. [Ergun1952] S. Ergun. Fluid flow through packed columns. Chem. Eng.
       Prog. 48: 89-94, 1952.
    """
    rt, rmin = 0.5 * tapholediameter, 0.5 * bedmindiameter
    eff_d = bedparticlediameter * bedparticlesphericity
    rrmax = rt / (0.5 * bedmaxdiameter)
    cmult = rt * (1 - bedporosity) / (eff_d * eff_d *
                                      bedporosity * bedporosity * bedporosity)
    amult = 0.145833333 * density * eff_d * cmult
    bmult = 75 * viscosity * (1 - bedporosity) * cmult
    if rmin > rt:
        rrmin = rt / rmin
        aconst = amult * (rrmin * rrmin * rrmin - rrmax * rrmax * rrmax)
        bconst = bmult * (rrmin - rrmax)
    else:
        rrmin = rmin / rt
        aconst = amult * (4 - 3 * rrmin * rrmin *
                          rrmin - rrmax * rrmax * rrmax)
        bconst = bmult * (2 - rrmin - rrmax)
    return aconst, bconst


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


# Convergence tolerance of the iteration for the tap-hole velocities, as the
# sum of the changes in the metal and slag velocities between iterations, m/s,
# and the number of iterations after which the iteration is abandoned.
_VELOCITY_TOLERANCE = 1e-6
_VELOCITY_MAX_ITERATIONS = 100

# Parameters of the tap-hole flow calculation at the start of a time step. The
# burden coefficients (a_ and b_) are those of the pressure drop through the
# burden and the tap-hole entrance, for each phase.
_TapholeParams = namedtuple('_TapholeParams', [
    'tapholediameter', 'tapholelength', 'tapholeheight', 'tapholeroughness',
    'bedporosity', 'densitymetal', 'densityslag', 'viscositymetal',
    'viscosityslag', 'fdmodel', 'tapholeopen_yn', 'allownegativeheights_yn',
    'a_metal', 'b_metal', 'a_slag', 'b_slag'])


def _tapping_pressure(hm, hs, rt, densitymetal, densityslag):
    """Calculate the pressure driving the flow through the tap-hole.

    The pressure is applied to both phases, based on undeformed interfaces. It
    is the integrated pressure field over the height of the tap-hole, averaged
    with the total height in the tap-hole.

    Parameters
    ----------
    hm : float
        Metal level relative to the tap-hole centreline, m.
    hs : float
        Slag level relative to the tap-hole centreline, m.
    rt : float
        Radius of the tap-hole, m.
    densitymetal : float
        Density of the molten metal phase, kg/m3.
    densityslag : float
        Density of the molten slag phase, kg/m3.

    Returns
    -------
    float
        The driving pressure, Pa.

    """
    if hs < -rt:
        # liquid level below taphole - no flow
        return 0
    if hs < rt:
        # liquid level between taphole limits
        if hm < -rt:
            # slag only
            return 0.5 * g * (hs + rt) * densityslag
        # slag and metal
        return 0.5 * g * ((densitymetal * (hm + rt) * (hm + rt)
                           + densityslag * (hs * (hs + 2 * rt)
                                            - hm * (hm + 2 * rt))) / (hs + rt))
    # liquid level above taphole
    if hm < -rt:
        # slag only
        return hs * densityslag * g
    if hm < rt:
        # slag and metal
        return 0.25 * g * ((densitymetal * (hm + rt) * (hm + rt)
                            + densityslag * (rt * (4 * hs - rt)
                                             - hm * (hm + 2 * rt))) / rt)
    # metal only
    return g * (densityslag * (hs - hm) + densitymetal * hm)


def _solve_velocities(pa, umetal0, uslag0, params):
    """Solve for the phase velocities in the tap-hole.

    The velocity of each phase satisfies a quadratic pressure balance, with
    the friction factor of the tap-hole channel depending on the velocity, so
    the solution is found by iterating on the friction factors.

    Parameters
    ----------
    pa : float
        The pressure driving the flow, Pa.
    umetal0 : float
        Initial guess for metal velocity in tap-hole, m/s.
    uslag0 : float
        Initial guess for slag velocity in tap-hole, m/s.
    params : _TapholeParams
        Parameters of the tap-hole flow calculation.

    Returns
    -------
    umetal : float
        Velocity of metal through the tap-hole, m/s.
    uslag : float
        Velocity of slag through the tap-hole, m/s.

    Raises
    ------
    RuntimeError
        If the iteration has not converged after the maximum number of
        iterations.

    """
    fdmodel, diameter, roughness = (params.fdmodel, params.tapholediameter,
                                    params.tapholeroughness)
    length = params.tapholelength
    rho_m, rho_s = params.densitymetal, params.densityslag
    mu_m, mu_s = params.viscositymetal, params.viscosityslag
    a_m, b_m, a_s, b_s = (params.a_metal, params.b_metal, params.a_slag,
                          params.b_slag)
    umetal, uslag = umetal0, uslag0
    for _ in range(_VELOCITY_MAX_ITERATIONS):
        fdm = fdmodel(umetal, rho_m, mu_m, diameter, roughness)
        fds = fdmodel(uslag, rho_s, mu_s, diameter, roughness)
        a_mc = a_m + (0.5 * rho_m * fdm * length / diameter)
        a_sc = a_s + (0.5 * rho_s * fds * length / diameter)
        nvm = (-b_m + math.sqrt(b_m * b_m + 4 * pa * a_mc)) / (2 * a_mc)
        nvs = (-b_s + math.sqrt(b_s * b_s + 4 * pa * a_sc)) / (2 * a_sc)
        residual = abs(nvm - umetal) + abs(nvs - uslag)
        umetal, uslag = nvm, nvs
        if not residual > _VELOCITY_TOLERANCE:
            return umetal, uslag
    raise RuntimeError("the tap-hole velocity iteration did not converge in "
                       f"{_VELOCITY_MAX_ITERATIONS} iterations")


def _interface_levels(umetal, uslag, params):
    """Calculate the deformation of the slag-metal interface by the flow.

    Parameters
    ----------
    umetal : float
        Velocity of metal through the tap-hole, m/s.
    uslag : float
        Velocity of slag through the tap-hole, m/s.
    params : _TapholeParams
        Parameters of the tap-hole flow calculation.

    Returns
    -------
    h0_l : float
        The level, relative to the tap-hole centreline, below which only slag
        enters the tap-hole, m.
    h0_h : float
        The level, relative to the tap-hole centreline, above which only metal
        enters the tap-hole, m.

    """
    p = params
    rt = 0.5 * p.tapholediameter
    pcorr = p.bedporosity ** 2
    h0_l = (-rt - p.densityslag * uslag * uslag
            / (8 * g * pcorr * (p.densitymetal - p.densityslag)))
    h0_h = (rt + p.densitymetal * umetal * umetal
            / (8 * g * pcorr * (p.densitymetal - p.densityslag)))
    return h0_l, h0_h


def _segment_area(rt, h):
    """Calculate the area of a circular segment of the tap-hole.

    Parameters
    ----------
    rt : float
        Radius of the tap-hole, m.
    h : float
        The height of the chord relative to the centre, m, between -rt and rt.
        The segment is the part of the circle below the chord.

    Returns
    -------
    float
        The area of the segment, m2.

    """
    theta = 2 * math.acos(-h / rt)
    return 0.5 * rt * rt * (theta - math.sin(theta))


def _phase_areas(hm, hs, h0_l, h0_h, rt):
    """Calculate the areas of the tap-hole cross-section occupied by each
    phase.

    Parameters
    ----------
    hm : float
        Metal level relative to the tap-hole centreline, m.
    hs : float
        Slag level relative to the tap-hole centreline, m.
    h0_l : float
        The level below which only slag enters the tap-hole, m.
    h0_h : float
        The level above which only metal enters the tap-hole, m.
    rt : float
        Radius of the tap-hole, m.

    Returns
    -------
    area_m : float
        Area of the tap-hole occupied by metal, m2.
    area_s : float
        Area of the tap-hole occupied by slag, m2.

    Note
    ----
    The chord height of the slag-metal interface in the tap-hole is a linear
    function of the metal level, from -rt (only slag enters the tap-hole) at
    h0_l to the top of the liquid in the tap-hole (only metal enters) at the
    metal level h0_h for a completely filled tap-hole. For a partly filled
    tap-hole that level is hs + (h0_h - rt)*(hs + rt)/(2*rt), which reduces to
    hs when there is no deformation of the interface (h0_h = rt), and to h0_h
    when the tap-hole is full, so that the areas, and with them the flows, are
    continuous as the slag level passes the top of the tap-hole.

    """
    if hs < -rt:
        # slag below taphole - no flow
        return 0, 0
    if hs < rt:
        # taphole partially filled
        if hm < h0_l:
            # slag only
            return 0, _segment_area(rt, hs)
        # slag and metal. The interface rises linearly with the metal level
        # from the bottom of the tap-hole at h0_l to the free surface at the
        # metal level h_top. This moves from hs for a nearly empty tap-hole,
        # where there is no deformation, to h0_h for a full one, so that the
        # interface joins the one below for a completely filled tap-hole.
        h_top = hs + (h0_h - rt) * (hs + rt) / (2 * rt)
        hi = -rt + (hs + rt) * (hm - h0_l) / (h_top - h0_l)
        area_m = _segment_area(rt, hi)
        return area_m, _segment_area(rt, hs) - area_m
    # taphole completely filled
    if hm < h0_l:
        # slag only
        return 0, pi * rt * rt
    if hm > h0_h:
        # metal only
        return pi * rt * rt, 0
    # slag and metal
    hi = rt * (1 - 2 * (h0_h - hm) / (h0_h - h0_l))
    area_m = _segment_area(rt, hi)
    return area_m, pi * rt * rt - area_m


def _calc_vdot_out(hmetal, hslag, umetal0, uslag0, params):
    """Calculate volumetric outflow rates of metal and slag through the
    tap-hole for the given furnace state.

    Solves for phase velocities accounting for burden pressure drop, tap-hole
    friction, and interface deformation effects. Returns zero flow when the
    tap-hole is closed.

    Parameters
    ----------
    hmetal : float
        Current metal bath level above hearth, m.
    hslag : float
        Current slag bath level above hearth, m.
    umetal0 : float
        Initial guess for metal velocity in tap-hole, m/s.
    uslag0 : float
        Initial guess for slag velocity in tap-hole, m/s.
    params : _TapholeParams
        Parameters of the tap-hole flow calculation.

    Returns
    -------
    umetal : float
        Velocity of metal through the tap-hole, m/s.
    uslag : float
        Velocity of slag through the tap-hole, m/s.
    vdotmetal_out : float
        Volumetric outflow rate of metal, m3/s.
    vdotslag_out : float
        Volumetric outflow rate of slag, m3/s.

    """
    if not params.tapholeopen_yn:
        return 0, 0, 0, 0
    rt = 0.5 * params.tapholediameter
    hs = hslag - params.tapholeheight
    hm = hmetal - params.tapholeheight
    pa = _tapping_pressure(hm, hs, rt, params.densitymetal, params.densityslag)
    umetal, uslag = _solve_velocities(pa, umetal0, uslag0, params)
    h0_l, h0_h = _interface_levels(umetal, uslag, params)
    area_m, area_s = _phase_areas(hm, hs, h0_l, h0_h, rt)
    if not params.allownegativeheights_yn:
        if hslag <= 0:
            uslag = 0
        if hmetal <= 0:
            umetal = 0
    return umetal, uslag, area_m * umetal, area_s * uslag


class SubmergedArcFurnace():
    """Submerged-arc furnace class. This version assumes a slag and metal
    smelting process, with a porous burden layer and a single tap-hole for
    both phases.

    Parameters
    ----------
    powerMVA : float
        The furnace power level, MVA.
    powerfactor : float
        The furnace power factor for conversion between VA and W.
    metalSER : float
        The specific energy requirement of the smelting process, MWh/ton
        metal produced.
    slagmetalmassratio : float
        The ratio of slag to metal produced during smelting.
    furnacediameter : float
        The furnace vessel inner diameter, m.
    activeareafraction : float
        The fraction of the cross-sectional area of the furnace vessel
        occupied by the molten bath.
    tapholediameter : float
        The diameter of the tap-hole channel, m.
    tapholelength : float
        The length of the tap-hole channel, m.
    tapholeroughness : float
        The roughness of the tap-hole channel surface.
    tapholeheight : float
        The position of the tap-hole centerline relative to the furnace
        hearth level, m.
    densitymetal : float
        Density of the molten metal phase, kg/m3.
    densityslag : float
        Density of the molten slag phase, kg/m3.
    viscositymetal : float
        Viscosity of the molten metal phase, Pa.s.
    viscosityslag : float
        Viscosity of the molten slag phase, Pa.s.
    particlediameter : float
        Diameter of the constituent particles of the burden layer, m.
    particlesphericity : float
        Sphericity of the constituent particles of the burden layer.
    bedporosity : float
        Porosity of the burden layer.
    bedmindiameter : float
        Diameter of cavity in burden in front of tap-hole, m (set to zero
        for no cavity).
    bedmaxdiameter : float
        Extent of burden from tap-hole, m (set to a high value for entire
        furnace).
    bedmodel : callable
        Function to be used to model the pressure drop through the burden
        section.
    entrykl : float
        Pressure loss coefficient to account for tap-hole entrance effects.
    fdmodel : callable
        Function to be used to model the pressure drop through the tap-hole
        channel.
    hmetal_init : float
        The initial level of metal in the furnace, m.
    hslag_init : float
        The initial level of slag in the furnace, m.

    Attributes
    ----------
    powerMVA : float
        The furnace power level, MVA.
    powerfactor : float
        The furnace power factor for conversion between VA and W.
    metalSER : float
        The specific energy requirement of the smelting process, MWh/ton
        metal produced.
    slagmetalmassratio : float
        The ratio of slag to metal produced during smelting.
    furnacediameter : float
        The furnace vessel inner diameter, m.
    activeareafraction : float
        The fraction of the cross-sectional area of the furnace vessel
        occupied by the molten bath.
    tapholediameter : float
        The diameter of the tap-hole channel, m.
    tapholelength : float
        The length of the tap-hole channel, m.
    tapholeroughness : float
        The roughness of the tap-hole channel surface.
    tapholeheight : float
        The position of the tap-hole centerline relative to the furnace
        hearth level, m.
    densitymetal : float
        Density of the molten metal phase, kg/m3.
    densityslag : float
        Density of the molten slag phase, kg/m3.
    viscositymetal : float
        Viscosity of the molten metal phase, Pa.s.
    viscosityslag : float
        Viscosity of the molten slag phase, Pa.s.
    particlediameter : float
        Diameter of the constituent particles of the burden layer, m.
    particlesphericity : float
        Sphericity of the constituent particles of the burden layer.
    bedporosity : float
        Porosity of the burden layer.
    bedmindiameter : float
        Diameter of cavity in burden in front of tap-hole, m (set to zero
        for no cavity).
    bedmaxdiameter : float
        Extent of burden from tap-hole, m (set to a high value for entire
        furnace).
    bedmodel : callable
        Function to be used to model the pressure drop through the burden
        section.
    entrykl : float
        Pressure loss coefficient to account for tap-hole entrance effects.
    fdmodel : callable
        Function to be used to model the pressure drop through the tap-hole
        channel.
    hmetal : float
        The current level of metal in the furnace, m.
    hslag : float
        The current level of slag in the furnace, m.
    allownegativeheights_yn : boolean
        Whether or not to allow the slag or metal height to pass below the
        level of the hearth. Default is False.
    tapholeopen_yn : boolean
        Indicate whether furnace tap-hole is currently open (True) or
        closed (False).
    vdotmetal_out : float
        The outlet volume flowrate of metal from the unit, averaged over the
        last time step, m3/s.
    vdotslag_out : float
        The outlet volume flowrate of slag from the unit, averaged over the
        last time step, m3/s.
    umetal : float
        The velocity of metal through the tap-hole at the start of the last
        time step, m/s
    uslag : float
        The velocity of slag through the tap-hole at the start of the last
        time step, m/s

    Note
    ----
    This model is based on the formulation by Olsen & Reynolds. [Olsen2020]_

    The metal and slag levels are integrated in time with the classical
    fourth-order Runge-Kutta scheme (see tappingsim.integrators.rk4_step). The
    drainage of the furnace is slow and smooth, so the scheme is stable for
    time steps of many minutes, and its error in the levels and the tapped
    masses is orders of magnitude lower than that of a first-order scheme with
    the same step. The outlet flowrates reported after a time step are the
    effective flowrates averaged over the step, so that the change in the
    amounts of metal and slag in the furnace is exactly the amount added less
    the amount tapped.

    References
    ----------
    .. [Olsen2020] J.E. Olsen, Q.G. Reynolds. Mathematical Modeling of Furnace
        Drainage While Tapping Slag and Metal Through a Single Tap-Hole.
        Metallurgical and Materials Transactions B 51(4): 1750-1759, 2020.
        doi:10.1007/s11663-020-01873-1.

    """

    def __init__(self, powerMVA, powerfactor, metalSER, slagmetalmassratio,
                 furnacediameter, activeareafraction, tapholediameter,
                 tapholelength, tapholeroughness, tapholeheight, densitymetal,
                 densityslag, viscositymetal, viscosityslag, particlediameter,
                 particlesphericity, bedporosity, bedmindiameter,
                 bedmaxdiameter, bedmodel, entrykl, fdmodel, hmetal_init,
                 hslag_init):
        self.furnacediameter = furnacediameter
        self.activeareafraction = activeareafraction
        self.tapholediameter = tapholediameter
        self.tapholelength = tapholelength
        self.tapholeroughness = tapholeroughness
        self.tapholeheight = tapholeheight
        self.densitymetal = densitymetal
        self.densityslag = densityslag
        self.viscositymetal = viscositymetal
        self.viscosityslag = viscosityslag
        self.particlediameter = particlediameter
        self.particlesphericity = particlesphericity
        self.bedporosity = bedporosity
        self.bedmindiameter = bedmindiameter
        self.bedmaxdiameter = bedmaxdiameter
        self.bedmodel = bedmodel
        self.entrykl = entrykl
        self.fdmodel = fdmodel
        self.powerMVA = powerMVA
        self.powerfactor = powerfactor
        self.metalSER = metalSER
        self.slagmetalmassratio = slagmetalmassratio
        self.hmetal = hmetal_init
        self.hslag = hslag_init
        self.allownegativeheights_yn = False
        self.tapholeopen_yn = False
        self.umetal = 0
        self.uslag = 0
        self.vdotmetal_out = 0
        self.vdotslag_out = 0

    def _calc_burden_coefficients(self, viscosity, density):
        """Calculate the coefficients of the pressure drop through the burden
        and the tap-hole entrance for one phase.

        The pressure drop is dP = A*u^2 + B*u, where u is the fluid velocity.

        Parameters
        ----------
        viscosity : float
            Viscosity of the phase, Pa.s.
        density : float
            Density of the phase, kg/m3.

        Returns
        -------
        a : float
            The coefficient of the quadratic term, including the entrance loss.
        b : float
            The coefficient of the linear term.

        """
        ab, bb = self.bedmodel(
            self.tapholediameter, self.bedmindiameter, self.bedmaxdiameter,
            self.particlediameter, self.particlesphericity, self.bedporosity,
            viscosity, density)
        return 0.5 * (1 + self.entrykl) * density + ab, 0 + bb

    def _calc_params(self):
        """Collect the parameters of the tap-hole flow calculation.

        The parameters are a snapshot of the current attributes, so they can be
        changed between time steps. The burden coefficients do not depend on
        the furnace state, so they are calculated once for a time step instead
        of for each evaluation of the flow.

        Returns
        -------
        _TapholeParams
            The parameters of _calc_vdot_out.

        """
        a_m = b_m = a_s = b_s = 0
        if self.tapholeopen_yn:
            a_m, b_m = self._calc_burden_coefficients(self.viscositymetal,
                                                      self.densitymetal)
            a_s, b_s = self._calc_burden_coefficients(self.viscosityslag,
                                                      self.densityslag)
        return _TapholeParams(
            self.tapholediameter, self.tapholelength, self.tapholeheight,
            self.tapholeroughness, self.bedporosity, self.densitymetal,
            self.densityslag, self.viscositymetal, self.viscosityslag,
            self.fdmodel, self.tapholeopen_yn, self.allownegativeheights_yn,
            a_m, b_m, a_s, b_s)

    def _calc_rates(self):
        """Calculate the coefficients linking the flows to the rates of change
        of the metal and slag levels.

        The metal level changes at the net inflow of metal over the area of the
        bath, and the level of the top of the slag layer changes at the net
        inflow of metal and slag.

        Returns
        -------
        inflows : numpy.ndarray
            The rates of change of the metal level and the slag level due to
            the smelting of feed, m/s.
        flows_to_rates : numpy.ndarray
            The matrix of rates of change of the levels per unit of tapped
            metal and slag flowrates, s/m2.

        """
        areaconst = 1 / (self.bedporosity * self.activeareafraction
                         * 0.25 * pi * self.furnacediameter**2)
        mdotmetal_in = (self.powerMVA * self.powerfactor /
                        (3.6 * self.metalSER))
        vdotmetal_in = mdotmetal_in / self.densitymetal
        vdotslag_in = mdotmetal_in * self.slagmetalmassratio / self.densityslag
        inflows = areaconst * numpy.array([vdotmetal_in,
                                           vdotmetal_in + vdotslag_in])
        flows_to_rates = -areaconst * numpy.array([[1.0, 0.0], [1.0, 1.0]])
        return inflows, flows_to_rates

    def _advance(self, dt, levels, rates, params, velocities):
        """Advance the metal and slag levels over a single Runge-Kutta step.

        Parameters
        ----------
        dt : float
            Length of time step, in s.
        levels : tuple of float
            The metal level and the slag level at the start of the step, m.
        rates : tuple of numpy.ndarray
            The inflows and flows-to-rates matrix from _calc_rates.
        params : _TapholeParams
            The parameters from _calc_params.
        velocities : list of float
            The metal and slag velocities in the tap-hole, m/s, used as initial
            guesses for the flow calculation and updated in place with the
            velocities from each evaluation.

        Returns
        -------
        levels : numpy.ndarray
            The metal level and the slag level at the end of the step, m.
        flows : numpy.ndarray
            The metal and slag outflow rates averaged over the step, m3/s.
        start_velocities : list of float
            The metal and slag velocities in the tap-hole at the start of the
            step, m/s.

        """
        start_velocities = []

        def flows(levels):
            umetal, uslag, vdotmetal_out, vdotslag_out = _calc_vdot_out(
                levels[0], levels[1], velocities[0], velocities[1], params)
            velocities[0], velocities[1] = umetal, uslag
            if not start_velocities:
                start_velocities.extend((umetal, uslag))
            return numpy.array([vdotmetal_out, vdotslag_out])

        levels, vdots = rk4_step(flows, levels, *rates, dt)
        return levels, vdots, start_velocities

    def calc_dt(self, dt):
        """
        Integrate object state forward over a single time step.

        The levels are advanced with a fourth-order Runge-Kutta step, and the
        outlet flowrates are set to the average of the flowrates over the
        step, so that the change in the amounts of metal and slag in the
        furnace is exactly the amount added less the amount tapped.

        Parameters
        ----------
        dt : float
            Length of time step, in s.

        Returns
        -------
        None.

        """
        levels, vdots, (umetal, uslag) = self._advance(
            dt, (self.hmetal, self.hslag), self._calc_rates(),
            self._calc_params(), [self.umetal, self.uslag])
        self.hmetal, self.hslag = float(levels[0]), float(levels[1])
        self.umetal, self.uslag = umetal, uslag
        self.vdotmetal_out = float(vdots[0])
        self.vdotslag_out = float(vdots[1])

    def calc_time_period(self, times):
        """
        Integrate object state forward over a series of time steps. The
        model parameters are assumed to remain fixed during this period. This
        is more performant than individual time step calculations, for large
        parameter sweeps and similar applications.

        Each step is a fourth-order Runge-Kutta step, as for calc_dt, and the
        tapped masses are accumulated from the flowrates averaged over each
        step.

        Parameters
        ----------
        times : array
            A 1D numpy array of time steps at which to perform each calculation,
            in s.

        Returns
        -------
        metalmassout : float
            The mass of metal removed from the furnace during the time period,
            kg
        slagmassout : float
            The mass of slag removed from the furnace during the time period, kg

        """
        dts = [float(dt) for dt in times[1:] - times[:-1]]
        params = self._calc_params()
        rates = self._calc_rates()
        metalmassout, slagmassout = 0, 0
        if not dts:
            return metalmassout, slagmassout

        levels, velocities = (self.hmetal, self.hslag), [1, 1]
        for dt in dts:
            levels, vdots, start_velocities = self._advance(
                dt, levels, rates, params, velocities)
            metalmassout += dt * vdots[0] * self.densitymetal
            slagmassout += dt * vdots[1] * self.densityslag

        self.hmetal, self.hslag = float(levels[0]), float(levels[1])
        self.umetal, self.uslag = start_velocities
        self.vdotmetal_out = float(vdots[0])
        self.vdotslag_out = float(vdots[1])

        return float(metalmassout), float(slagmassout)
