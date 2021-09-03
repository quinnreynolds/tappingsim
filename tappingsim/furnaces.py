import math
from scipy.constants import g, pi

_POWER_TIME_FACTOR = 1000/3600

def bedmodel_carmenkozeny(tapholediameter, bedmindiameter, bedmaxdiameter,
                          bedparticlediameter, bedparticlesphericity, 
                          bedporosity, viscosity, density):
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
    d_over_e = diameter / roughness
    NRe = diameter * velocity * density / viscosity
    if NRe < 1:
        NRe = 1
    a = 1 / (1 + (0.000368732*NRe) ** 8.4)
    b = 1 / (1 + (NRe/(150*d_over_e)) ** 1.8)
    return ((64 / NRe) ** a
            * 0.75 * math.log(0.186219739*NRe) ** (2*(a-1)*b) 
            * 0.88 * math.log(3.41*d_over_e) ** (2*(a-1)*(1-b)))

def fdmodel_serghides1(velocity, density, viscosity, diameter, roughness):
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

def fdmodel_eck(velocity, density, viscosity, diameter, roughness):
    e_over_d = roughness / diameter
    invNRe = viscosity / (diameter * velocity * density)
    if invNRe > 1:
        return 64
    elif invNRe > 0.000333333:
        return 64*invNRe
    else:
        invsqrtf = -2 * math.log10(0.269179004*e_over_d + 15*invNRe)
        return 1/(invsqrtf*invsqrtf)


class SubmergedArcFurnace():
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
            if umetal < 1e-6:
                umetal = 1e-6
            if uslag < 1e-6:
                uslag = 1e-6
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
            h0_l = (-rt - densityslag*uslag*uslag 
                    / (8*g*(densitymetal-densityslag)))
            h0_h = (rt + densitymetal*umetal*umetal 
                    / (8*g*(densitymetal-densityslag)))
            
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
        #activearea = (self.activeareafraction 
        #              * 0.25*pi*self.furnacediameter**2)
        areaconst = 1 / (self.bedporosity * self.activeareafraction 
                         * 0.25*pi*self.furnacediameter**2)
        mdotmetal_in = (_POWER_TIME_FACTOR * self.powerMVA * self.powerfactor 
                        / self.metalSER)
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
        dts = [float(dt) for dt in times[1:]-times[:-1]]
        params = (self.tapholediameter, self.tapholelength, self.tapholeheight,
                  self.tapholeroughness, self.entrykl, self.bedmindiameter, 
                  self.bedmaxdiameter, self.bedporosity, self.particlediameter,
                  self.particlesphericity, self.densitymetal, self.densityslag,
                  self.viscositymetal, self.viscosityslag, self.tapholeopen_yn,
                  self.allownegativeheights_yn, self.bedmodel, self.fdmodel)
        areaconst = 1 / (self.bedporosity * self.activeareafraction 
                         * 0.25*pi*self.furnacediameter**2)
        mdotmetal_in = (_POWER_TIME_FACTOR * self.powerMVA * self.powerfactor 
                        / self.metalSER)
        vdotmetal_in = mdotmetal_in / self.densitymetal
        vdotslag_in = mdotmetal_in*self.slagmetalmassratio / self.densityslag
        
        hmetal, hslag, umetal, uslag = self.hmetal, self.hslag, 1, 1
        for dt in dts:
            umetal, uslag, vdotmetal_out, vdotslag_out = self._calc_vdot_out(
                hmetal, hslag, umetal, uslag, *params)
            dhmetal = dt * areaconst * (vdotmetal_in-vdotmetal_out)
            dhslag = dt * areaconst * (vdotslag_in-vdotslag_out)
            hmetal += dhmetal
            hslag += dhmetal + dhslag
        
        self.hmetal, self.hslag = hmetal, hslag
        self.umetal, self.uslag = umetal, uslag
        self.vdotmetal_out, self.vdotslag_out = vdotmetal_out, vdotslag_out
        
        
        
