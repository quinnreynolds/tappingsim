import math
#import numpy
from scipy.constants import g, pi

POWER_TIME_FACTOR = 1000/3600

def bedmodel_carmenkozeny(tapholediameter, bedmindiameter, bedmaxdiameter,
                          bedparticlediameter, bedparticlesphericity, 
                          bedporosity, viscosity, density):
    rt, rmin, rmax = 0.5*tapholediameter, 0.5*bedmindiameter, 0.5*bedmaxdiameter
    eff_d = bedparticlediameter * bedparticlesphericity
    rrmax = rt / rmax
    if rmin > rt:
        rrmin = rt / rmin
        bconst = (180 * viscosity * rt * (1-bedporosity)**2 / 
                  (2*eff_d**2 * bedporosity**3) * (rrmin - rrmax))
    else:
        rrmin = rmin / rt
        bconst = (180 * viscosity * rt * (1-bedporosity)**2 / 
               (2*eff_d**2 * bedporosity**3) * (2 - rrmin - rrmax))
    return 0, bconst

def bedmodel_ergun(tapholediameter, bedmindiameter, bedmaxdiameter,
                   bedparticlediameter, bedparticlesphericity, bedporosity,
                   viscosity, density):
    rt, rmin = 0.5*tapholediameter, 0.5*bedmindiameter
    eff_d = bedparticlediameter * bedparticlesphericity
    rrmax = rt / ( 0.5*bedmaxdiameter)
    if rmin > rt:
        rrmin = rt / rmin
        aconst = (1.75 * density * rt * (1-bedporosity) / 
                  (12*eff_d*bedporosity**3) * (rrmin**3 - rrmax**3))
        bconst = (150*viscosity*rt*(1-bedporosity)**2 / 
                  (2*eff_d**2 * bedporosity**3) * (rrmin - rrmax))
    else:
        rrmin = rmin / rt
        aconst = (1.75 * density * rt * (1-bedporosity) / 
                  (12*eff_d*bedporosity**3) * (4 - 3*rrmin**3 - rrmax**3))
        bconst = (150 * viscosity * rt * (1-bedporosity)**2 / 
                  (2*eff_d**2 * bedporosity**3) * (2 - rrmin - rrmax))
    return aconst, bconst
    
def fdmodel_cheng(velocity, density, viscosity, diameter, roughness):
    d_over_e = diameter / roughness
    Re = max(diameter * velocity * density / viscosity, 1)
    a = 1 / (1 + (Re/2712) ** 8.4)
    b = 1 / (1 + (Re/(150*d_over_e)) ** 1.8)
    return ((64 / Re) ** a 
            * (0.75 * math.log(Re/5.37)) ** (2*(a-1)*b) 
            * (0.88 * math.log(3.41*d_over_e)) ** (2*(a-1)*(1-b)))


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
    
    def _calc_vdot_out(self, hmetal, hslag, tapholediameter, tapholelength, 
                       tapholeheight, tapholeroughness, entrykl, 
                       bedmindiameter, bedmaxdiameter, bedporosity, 
                       particlediameter, particlesphericity, densitymetal, 
                       densityslag, viscositymetal, viscosityslag, 
                       tapholeopen_yn, allownegativeheights_yn, bedmodel, 
                       fdmodel):
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
                    pa = 0.5*g*((densitymetal*(hm+rt)**2 
                                 + densityslag*(hs*(hs+2*rt) 
                                                - hm*(hm+2*rt))) / (hs+rt))
            else:
                # liquid level above taphole
                if hm < -rt:
                    # slag only
                    pa = hs * densityslag * g
                elif hm < rt:
                    # slag and metal
                    pa = 0.25*g*((densitymetal*(hm+rt)**2 
                                  + densityslag*(rt*(4*hs-rt) 
                                                 - hm*(hm+2*rt))) / rt)
                else:
                    # metal only
                    pa = g * (densityslag*(hs-hm) + densitymetal*hm)
            
            # calculate phase velocities
            converged, umetal, uslag = 1, 1, 1
            while converged > 1e-6:
                fdm = fdmodel(umetal, densitymetal, viscositymetal, 
                              tapholediameter, tapholeroughness)
                fds = fdmodel(uslag, densityslag, viscosityslag, 
                              tapholediameter, tapholeroughness)
                a_mc = a_m + (0.5 * densitymetal * fdm 
                              * tapholelength/tapholediameter)
                a_sc = a_s + (0.5 * densityslag * fds
                              * tapholelength/tapholediameter)
                nvm = (-b_m+math.sqrt(b_m**2 + 4*pa*a_mc)) / (2*a_mc)
                nvs = (-b_s+math.sqrt(b_s**2 + 4*pa*a_sc)) / (2*a_sc)
                converged = abs(nvm-umetal) + abs(nvs-uslag)
                umetal, uslag = nvm, nvs
                
            # interface deformations
            h0_l = -rt - densityslag*uslag**2/(8*g*(densitymetal-densityslag))
            h0_h = rt + densitymetal*umetal**2/(8*g*(densitymetal-densityslag))
            
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
                    area_s = 0.5 * rt**2 * (theta - math.sin(theta))
                else:
                    # slag and metal
                    hi = hs - (hs-hm) * (hs+rt) / (hs-h0_l)
                    theta = 2*math.acos(-hi/rt)
                    area_m = 0.5 * rt**2 * (theta - math.sin(theta))
                    theta = 2*math.acos(-hs/rt)
                    area_s = 0.5 * rt**2 * (theta - math.sin(theta)) - area_m
            else:
                # taphole completely filled
                if hm < h0_l:
                    # slag only
                    area_m = 0
                    area_s = pi * rt**2
                elif hm > h0_h:
                    # metal only
                    area_m = pi * rt**2
                    area_s = 0
                else:
                    # slag and metal
                    hi = rt * (1 - 2 * (h0_h-hm) / (h0_h-h0_l))
                    theta = 2 * math.acos(-hi/rt)
                    area_m = 0.5 * rt**2 * (theta - math.sin(theta))
                    area_s = pi * rt**2 - area_m                
            vdotmetal_out = area_m * umetal
            vdotslag_out = area_s * uslag
            
            if not allownegativeheights_yn:
                if hslag <= 0:
                    vdotslag_out = 0
                if hmetal <= 0:
                    vdotmetal_out = 0, 0
        else:
            vdotmetal_out, vdotslag_out = 0, 0
        
        return vdotmetal_out, vdotslag_out

    def calc_dt(self, dt):
        activearea = (self.activeareafraction 
                      * 0.25*pi*self.furnacediameter**2)
        mdotmetal_in = (POWER_TIME_FACTOR * self.powerMVA * self.powerfactor 
                        / self.metalSER)
        vdotmetal_in = mdotmetal_in / self.densitymetal
        vdotslag_in = mdotmetal_in*self.slagmetalmassratio / self.densityslag
        
        vdotmetal_out, vdotslag_out = self._calc_vdot_out(self.hmetal, 
            self.hslag, self.tapholediameter, self.tapholelength, 
            self.tapholeheight, self.tapholeroughness, self.entrykl, 
            self.bedmindiameter, self.bedmaxdiameter, self.bedporosity, 
            self.particlediameter, self.particlesphericity, self.densitymetal, 
            self.densityslag, self.viscositymetal, self.viscosityslag, 
            self.tapholeopen_yn, self.allownegativeheights_yn, 
            self.bedmodel, self.fdmodel)
        
        dhmetal = dt*(vdotmetal_in-vdotmetal_out)/(activearea*self.bedporosity)
        dhslag = dt*(vdotslag_in-vdotslag_out)/(activearea*self.bedporosity)
        self.hmetal += dhmetal
        self.hslag += dhmetal + dhslag
        self.vdotmetal_out = vdotmetal_out
        self.vdotslag_out = vdotslag_out
        
    def calc_time_period(self, times):
        # TODO
        pass
        
