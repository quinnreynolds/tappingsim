"""Classes to create furnace objects.

Metallurgical furnaces typically produce one or more immiscible molten phase, 
and do so continuously provided the furnace is supplied with power and feed 
material. The molten products are removed from the vessel by tapping through 
one or more tap-holes.

This module also includes support functions for describing flow drag effects in 
packed beds and tap-hole channels.
"""

import math
from scipy.constants import g, pi

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
    correlation is valid for laminar flow conditions only. [1]_ [2]_
    
    References
    ----------
    .. [1] J. Kozeny. Ueber kapillare Leitung des Wassers im Boden. Sitzungsber 
       Akad. Wiss., Wien, 136(2a): 271-306, 1927.
    .. [2] P.C. Carman. Fluid flow through granular beds. Transactions, 
       Institution of Chemical Engineers, London, 15: 150-166, 1937.

    """
    rt, rmin = 0.5*tapholediameter, 0.5*bedmindiameter
    eff_d = bedparticlediameter * bedparticlesphericity
    rrmax = rt / (0.5*bedmaxdiameter)
    if rmin > rt:
        rrmin = rt / rmin
        bconst = (180 * viscosity * rt * (1-bedporosity)*(1-bedporosity) 
                  / (2*eff_d*eff_d * bedporosity*bedporosity*bedporosity) 
                  * (rrmin - rrmax))
    else:
        rrmin = rmin / rt
        bconst = (180 * viscosity * rt * (1-bedporosity)*(1-bedporosity)
                  / (2*eff_d*eff_d * bedporosity*bedporosity*bedporosity) 
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
    valid for both laminar and turbulent flow regimes. [1]_
    
    References
    ----------
    .. [1] S. Ergun. Fluid flow through packed columns. Chem. Eng. Prog. 48: 
       89-94, 1952.
    """
    rt, rmin = 0.5*tapholediameter, 0.5*bedmindiameter
    eff_d = bedparticlediameter * bedparticlesphericity
    rrmax = rt / ( 0.5*bedmaxdiameter)
    cmult = rt * (1-bedporosity) / (eff_d*eff_d * 
                                    bedporosity*bedporosity*bedporosity)
    amult = 0.145833333 * density * eff_d * cmult
    bmult = 75 * viscosity * (1-bedporosity) * cmult
    if rmin > rt:
        rrmin = rt / rmin
        aconst = amult * (rrmin*rrmin*rrmin - rrmax*rrmax*rrmax)
        bconst = bmult * (rrmin - rrmax)
    else:
        rrmin = rmin / rt
        aconst = amult * (4 - 3*rrmin*rrmin*rrmin - rrmax*rrmax*rrmax)
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
    correlation is valid in all flow regimes including transitional flow. [1]_ 
    [2]_
    
    References
    ----------
    .. [1] V. Bellos, I. Nalbantis, G. Tsakiris. Friction Modeling of Flood 
       Flow Simulations. Journal of Hydraulic Engineering, 144(12): 04018073, 
       2018. doi:10.1061/(asce)hy.1943-7900.0001540.
    .. [2] V. Bellos, I. Nalbantis, G. Tsakiris. Erratum for Friction Modeling 
       of Flood Flow Simulations. Journal of Hydraulic Engineering, 146(10): 
       08220005, 2020. doi:10.1061/(ASCE)HY.1943-7900.0001802.
    """
    d_over_e = diameter / roughness
    NRe = diameter * velocity * density / viscosity
    if NRe < 1:
        NRe = 1
    a = 1 / (1 + (0.000368732*NRe) ** 8.4)
    b = 1 / (1 + (0.006666667*NRe/d_over_e)** 1.8)
    return ((0.015625*NRe) ** (-a)
            * (0.75 * math.log(0.186219739*NRe)) ** (2*(a-1)*b) 
            * (0.88 * math.log(3.41*d_over_e)) ** (2*(a-1)*(1-b)))

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
    correlation is valid in all flow regimes including transitional flow. [1]_
    
    References
    ----------
    .. [1] N.-S. Cheng. Formulas for Friction Factor in Transitional Regimes. 
       Journal of Hydraulic Engineering. 134(9): 1357-1362, 2008.
    """
    d_over_e = diameter / roughness
    NRe = diameter * velocity * density / viscosity
    if NRe < 1:
        NRe = 1
    a = 1 / (1 + (0.000367647*NRe) ** 9)
    b = 1 / (1 + (0.00625*NRe/d_over_e) ** 2)
    return ((0.015625*NRe) ** (-a) 
            * (1.8 * math.log10(0.147058824*NRe)) ** (2*(a-1)*b) 
            * (2.0 * math.log10(3.7*d_over_e)) ** (2*(a-1)*(1-b)))

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
    correlation is valid for laminar and turbulent flow regimes. [1]_
    
    References
    ----------
    .. [1] T.K. Serghides. Estimate friction factor accurately. Chemical 
       Engineering Journal. 91(5): 63-64, 1984.
    """
    if velocity < 1e-6:
        return 64
    e_over_d = roughness / diameter
    invNRe = viscosity / (diameter * velocity * density)
    if invNRe > 1:
        return 64
    elif invNRe > 0.000333333:
        return 64*invNRe
    else:
        eod = 0.27027027*e_over_d
        psi1 = -2 * math.log10(eod + 12*invNRe)
        psi2 = -2 * math.log10(eod + 2.51*psi1*invNRe)
        psi3 = -2 * math.log10(eod + 2.51*psi2*invNRe)
        invsqrtf = psi1 - (psi2-psi1)*(psi2-psi1)/(psi3-2*psi2+psi1)
        return 1/(invsqrtf*invsqrtf)

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
    correlation is valid laminar and turbulent flow regimes. [1]_
    
    References
    ----------
    .. [1] T.K. Serghides. Estimate friction factor accurately. Chemical 
       Engineering Journal. 91(5): 63-64, 1984.
    """
    if velocity < 1e-6:
        return 64
    e_over_d = roughness / diameter
    invNRe = viscosity / (diameter * velocity * density)
    if invNRe > 1:
        return 64
    elif invNRe > 0.000333333:
        return 64*invNRe
    else:
        eod = 0.27027027*e_over_d
        psi1 = -2 * math.log10(eod + 12*invNRe)
        psi2 = -2 * math.log10(eod + 2.51*psi1*invNRe)
        invsqrtf = 4.781 - (psi1-4.781)*(psi1-4.781)/(psi2-2*psi1+4.781)
        return 1/(invsqrtf*invsqrtf)


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
    bedmodel : float
        Function to be used to model the pressure drop through the burden
        section.
    entrykl : float
        Pressure loss coefficient to account for tap-hole entrance effects.
    fdmodel : float
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
    bedmodel : float
        Function to be used to model the pressure drop through the burden
        section.
    entrykl : float
        Pressure loss coefficient to account for tap-hole entrance effects.
    fdmodel : float
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
        The current outlet volume flowrate of metal from the unit, m3/s.
    vdotslag_out : float
        The current outlet volume flowrate of slag from the unit, m3/s.
    umetal : float
        The current velocity of metal through the tap-hole, m/s
    uslag : float
        The current velocity of slag through the tap-hole, m/s
        
    Note
    ----
    This model is based on the formulation by Olsen & Reynolds. [1]_
    
    References
    ----------
    .. [1] J.E. Olsen, Q.G. Reynolds. Mathematical Modeling of Furnace 
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
    
    def _calc_vdot_out(self, hmetal, hslag, umetal0, uslag0, 
                       tapholediameter, tapholelength, tapholeheight, 
                       tapholeroughness, entrykl, bedmindiameter, 
                       bedmaxdiameter, bedporosity, particlediameter, 
                       particlesphericity, densitymetal, densityslag, 
                       viscositymetal, viscosityslag, tapholeopen_yn, 
                       allownegativeheights_yn, bedmodel, fdmodel):
        umetal, uslag, vdotmetal_out, vdotslag_out = 0, 0, 0, 0
        if tapholeopen_yn:
            rt = 0.5 * tapholediameter
            a_m = 0.5 * (1 + entrykl) * densitymetal
            a_s = 0.5 * (1 + entrykl) * densityslag
            b_m, b_s = 0, 0
            
            # Bed pressure drop coefficient calculation
            ab, bb = bedmodel(tapholediameter, bedmindiameter, bedmaxdiameter, 
                              particlediameter, particlesphericity, bedporosity,
                              viscositymetal, densitymetal)
            a_m += ab
            b_m += bb
            ab, bb = bedmodel(tapholediameter, bedmindiameter, bedmaxdiameter, 
                              particlediameter, particlesphericity, bedporosity,
                              viscosityslag, densityslag)
            a_s += ab
            b_s += bb
                
            # Tapping pressure calculation (applied to both phases, based on 
            # undeformed interfaces). Integrated pressure field in z, averaged 
            # with total height in taphole.
            hs = hslag - tapholeheight
            hm = hmetal - tapholeheight
            if hs < -rt:
                # liquid level below taphole - no flow
                pa = 0
            elif hs < rt:
                # liquid level between taphole limits
                if hm < -rt:
                    # slag only
                    pa = 0.5*g*(hs+rt)*densityslag
                else:
                    # slag and metal
                    pa = 0.5*g*((densitymetal*(hm+rt)*(hm+rt)
                                 + densityslag*(hs*(hs+2*rt) 
                                                - hm*(hm+2*rt))) / (hs+rt))
            else:
                # liquid level above taphole
                if hm < -rt:
                    # slag only
                    pa = hs * densityslag * g
                elif hm < rt:
                    # slag and metal
                    pa = 0.25*g*((densitymetal*(hm+rt)*(hm+rt)
                                  + densityslag*(rt*(4*hs-rt) 
                                                 - hm*(hm+2*rt))) / rt)
                else:
                    # metal only
                    pa = g * (densityslag*(hs-hm) + densitymetal*hm)
            
            # calculate phase velocities
            converged, umetal, uslag = 1, umetal0, uslag0
            while converged > 1e-6:
                fdm = fdmodel(umetal, densitymetal, viscositymetal, 
                              tapholediameter, tapholeroughness)
                fds = fdmodel(uslag, densityslag, viscosityslag, 
                              tapholediameter, tapholeroughness)
                a_mc = a_m + (0.5 * densitymetal * fdm 
                              * tapholelength/tapholediameter)
                a_sc = a_s + (0.5 * densityslag * fds
                              * tapholelength/tapholediameter)
                nvm = (-b_m+math.sqrt(b_m*b_m + 4*pa*a_mc)) / (2*a_mc)
                nvs = (-b_s+math.sqrt(b_s*b_s + 4*pa*a_sc)) / (2*a_sc)
                converged = abs(nvm-umetal) + abs(nvs-uslag)
                umetal, uslag = nvm, nvs
            
            # interface deformations
            pcorr = self.bedporosity ** 2
            h0_l = (-rt - densityslag*uslag*uslag 
                    / (8*g*pcorr*(densitymetal-densityslag)))
            h0_h = (rt + densitymetal*umetal*umetal 
                    / (8*g*pcorr*(densitymetal-densityslag)))
            
            if hs < -rt:
                # slag below taphole - no flow
                area_m = 0
                area_s = 0
            elif hs < rt:
                # taphole partially filled
                if hm < h0_l:
                    # slag only
                    theta = 2 * math.acos(-hs/rt)
                    area_m = 0
                    area_s = 0.5 * rt*rt * (theta - math.sin(theta))
                else:
                    # slag and metal
                    hi = hs - (hs-hm) * (hs+rt) / (hs-h0_l)
                    theta = 2*math.acos(-hi/rt)
                    area_m = 0.5 * rt*rt * (theta - math.sin(theta))
                    theta = 2*math.acos(-hs/rt)
                    area_s = 0.5 * rt*rt * (theta - math.sin(theta)) - area_m
            else:
                # taphole completely filled
                if hm < h0_l:
                    # slag only
                    area_m = 0
                    area_s = pi * rt*rt
                elif hm > h0_h:
                    # metal only
                    area_m = pi * rt*rt
                    area_s = 0
                else:
                    # slag and metal
                    hi = rt * (1 - 2 * (h0_h-hm) / (h0_h-h0_l))
                    theta = 2 * math.acos(-hi/rt)
                    area_m = 0.5 * rt*rt * (theta - math.sin(theta))
                    area_s = pi * rt*rt - area_m                
            if not allownegativeheights_yn:
                if hslag <= 0:
                    uslag = 0
                if hmetal <= 0:
                    umetal = 0
            vdotmetal_out = area_m * umetal
            vdotslag_out = area_s * uslag
            
        return umetal, uslag, vdotmetal_out, vdotslag_out

    def calc_dt(self, dt):
        """
        Integrate object state forward over a single time step.

        Parameters
        ----------
        dt : float
            Length of time step, in s.

        Returns
        -------
        None.
        
        """
        areaconst = 1 / (self.bedporosity * self.activeareafraction 
                         * 0.25*pi*self.furnacediameter**2)
        mdotmetal_in = (self.powerMVA * self.powerfactor / (3.6*self.metalSER))
        vdotmetal_in = mdotmetal_in / self.densitymetal
        vdotslag_in = mdotmetal_in*self.slagmetalmassratio / self.densityslag
        
        umetal, uslag, vdotmetal_out, vdotslag_out = self._calc_vdot_out(
            self.hmetal, self.hslag, self.umetal, self.uslag, 
            self.tapholediameter, self.tapholelength, self.tapholeheight, 
            self.tapholeroughness, self.entrykl, self.bedmindiameter, 
            self.bedmaxdiameter, self.bedporosity, self.particlediameter, 
            self.particlesphericity, self.densitymetal, self.densityslag, 
            self.viscositymetal, self.viscosityslag, self.tapholeopen_yn, 
            self.allownegativeheights_yn, self.bedmodel, self.fdmodel)
        
        dhmetal = dt * areaconst * (vdotmetal_in-vdotmetal_out)
        dhslag = dt * areaconst * (vdotslag_in-vdotslag_out)
        self.hmetal += dhmetal
        self.hslag += dhmetal + dhslag
        self.umetal, self.uslag = umetal, uslag
        self.vdotmetal_out, self.vdotslag_out = vdotmetal_out, vdotslag_out
        
    def calc_time_period(self, times):
        """
        Integrate object state forward over a series of time steps. The 
        model parameters are assumed to remain fixed during this period. This 
        is more performant than individual time step calculations, for large 
        parameter sweeps and similar applications.

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
        dts = [float(dt) for dt in times[1:]-times[:-1]]
        params = (self.tapholediameter, self.tapholelength, self.tapholeheight,
                  self.tapholeroughness, self.entrykl, self.bedmindiameter, 
                  self.bedmaxdiameter, self.bedporosity, self.particlediameter,
                  self.particlesphericity, self.densitymetal, self.densityslag,
                  self.viscositymetal, self.viscosityslag, self.tapholeopen_yn,
                  self.allownegativeheights_yn, self.bedmodel, self.fdmodel)
        areaconst = 1 / (self.bedporosity * self.activeareafraction 
                         * 0.25*pi*self.furnacediameter**2)
        mdotmetal_in = (self.powerMVA * self.powerfactor / (3.6*self.metalSER))
        vdotmetal_in = mdotmetal_in / self.densitymetal
        vdotslag_in = mdotmetal_in*self.slagmetalmassratio / self.densityslag
        
        metalmassout, slagmassout = 0, 0
        hmetal, hslag, umetal, uslag = self.hmetal, self.hslag, 1, 1
        for dt in dts:
            umetal, uslag, vdotmetal_out, vdotslag_out = self._calc_vdot_out(
                hmetal, hslag, umetal, uslag, *params)
            dhmetal = dt * areaconst * (vdotmetal_in-vdotmetal_out)
            dhslag = dt * areaconst * (vdotslag_in-vdotslag_out)
            hmetal += dhmetal
            hslag += dhmetal + dhslag
            metalmassout += dt * (vdotmetal_out * params[10])
            slagmassout += dt * (vdotslag_out * params[11])
        
        self.hmetal, self.hslag = hmetal, hslag
        self.umetal, self.uslag = umetal, uslag
        self.vdotmetal_out, self.vdotslag_out = vdotmetal_out, vdotslag_out
        
        return metalmassout, slagmassout
        
        
