import numpy
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
    rt, rmin, rmax = 0.5*tapholediameter, 0.5*bedmindiameter, 0.5*bedmaxdiameter
    eff_d = bedparticlediameter * bedparticlesphericity
    rrmax = rt / rmax
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
            * (0.75 * numpy.log(Re/5.37)) ** (2*(a-1)*b) 
            * (0.88 * numpy.log(3.41*d_over_e)) ** (2*(a-1)*(1-b)))


class SubmergedArcFurnace():
    def __init__(self, powerMVA, powerfactor, metalSER, slagmetalmassratio, 
                 activearea, tapholediameter, tapholelength, tapholeroughness,
                 tapholeheight, densitymetal, densityslag, viscositymetal, 
                 viscosityslag, particlediameter, particlesphericity, 
                 bedporosity, bedmindiameter, bedmaxdiameter, bedmodel, 
                 entrykl, fdmodel, hmetal_init, hslag_init):
        self.activearea = activearea
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
    
    def calc_vdot_out(self):
        """Return outlet flowrates as a function of the current furnace state. 
        Extended model with semi-empirical interface deformation near tap-hole 
        entry.
        """
        if self.tapholeopen_yn:
            rt = 0.5 * self.tapholediameter
            a_m = 0.5 * (1 + self.entrykl) * self.densitymetal
            a_s = 0.5 * (1 + self.entrykl) * self.densityslag
            b_m, b_s = 0, 0
            
            # Bed pressure drop coefficient calculation
            ab, bb = self.bedmodel(self.tapholediameter, self.bedmindiameter, 
                                   self.bedmaxdiameter, self.particlediameter, 
                                   self.particlesphericity, self.bedporosity,
                                   self.viscositymetal, self.densitymetal)
            a_m += ab
            b_m += bb
            ab, bb = self.bedmodel(self.tapholediameter, self.bedmindiameter, 
                                   self.bedmaxdiameter, self.particlediameter, 
                                   self.particlesphericity, self.bedporosity,
                                   self.viscosityslag, self.densityslag)
            a_s += ab
            b_s += bb
                
            # Tapping pressure calculation (applied to both phases, based on 
            # undeformed interfaces). Integrated pressure field in z, averaged 
            # with total height in taphole.
            hs = self.hslag-self.tapholeheight
            hm = self.hmetal-self.tapholeheight
            if hs < -rt:
                # liquid level below taphole - no flow
                pa = 0
            elif hs < rt:
                # liquid level between taphole limits
                if hm < -rt:
                    # slag only
                    pa = 0.5*g*(hs+rt)*self.densityslag
                else:
                    # slag and metal
                    pa = 0.5*g*((self.densitymetal*(hm+rt)**2 
                                 + self.densityslag*(hs*(hs+2*rt) 
                                                     - hm*(hm+2*rt))) / (hs+rt))
            else:
                # liquid level above taphole
                if hm < -rt:
                    # slag only
                    pa = hs * self.densityslag * g
                elif hm < rt:
                    # slag and metal
                    pa = 0.25*g*((self.densitymetal*(hm+rt)**2 
                                  + self.densityslag*(rt*(4*hs-rt)
                                                      - hm*(hm+2*rt))) / rt)
                else:
                    # metal only
                    pa = g * (self.densityslag*(hs-hm) + self.densitymetal*hm)
            
            # calculate phase velocities
            converged = 1
            while converged > 1e-6:
                fdm = self.fdmodel(self.umetal, self.densitymetal, 
                                   self.viscositymetal,self.tapholediameter, 
                                   self.tapholeroughness)
                fds = self.fdmodel(self.uslag, self.densityslag, 
                                   self.viscosityslag, self.tapholediameter, 
                                   self.tapholeroughness)
                a_mc = a_m + (0.5 * self.densitymetal * fdm 
                              * self.tapholelength/self.tapholediameter)
                a_sc = a_s + (0.5 * self.densityslag * fds
                              * self.tapholelength/self.tapholediameter)
                nvm = (-b_m+numpy.sqrt(b_m**2 + 4*pa*a_mc)) / (2*a_mc)
                nvs = (-b_s+numpy.sqrt(b_s**2 + 4*pa*a_sc)) / (2*a_sc)
                converged = max(abs(nvm-self.umetal), abs(nvs-self.uslag))
                self.umetal, self.uslag = nvm, nvs
                
            # interface deformations
            h0_l = -rt-self.densityslag*self.uslag**2/(8*g*(self.densitymetal 
                                                            -self.densityslag))
            h0_h = rt+self.densitymetal*self.umetal**2/(8*g*(self.densitymetal 
                                                             -self.densityslag))
            
            if hs < -rt:
                # slag below taphole - no flow
                area_m = 0
                area_s = 0
            elif hs < rt:
                # taphole partially filled
                if hm < h0_l:
                    # slag only
                    theta = 2 * numpy.arccos(-hs/rt)
                    area_m = 0
                    area_s = 0.5 * rt**2 * (theta - numpy.sin(theta))
                else:
                    # slag and metal
                    hi = hs - (hs-hm) * (hs+rt) / (hs-h0_l)
                    theta = 2*numpy.arccos(-hi/rt)
                    area_m = 0.5 * rt**2 * (theta - numpy.sin(theta))
                    theta = 2*numpy.arccos(-hs/rt)
                    area_s = 0.5 * rt**2 * (theta - numpy.sin(theta)) - area_m
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
                    theta = 2 * numpy.arccos(-hi/rt)
                    area_m = 0.5 * rt**2 * (theta - numpy.sin(theta))
                    area_s = pi * rt**2 - area_m                
            self.vdotmetal_out = area_m * self.umetal
            self.vdotslag_out = area_s * self.uslag
            
            if not self.allownegativeheights_yn:
                if self.hslag <= 0:
                    self.uslag, self.vdotslag_out = 0, 0
                if self.hmetal <= 0:
                    self.umetal, self.vdotmetal_out = 0, 0
        else:
            self.uslag, self.vdotslag_out = 0, 0
            self.umetal, self.vdotmetal_out = 0, 0

    def calc_dt(self, dt):
        mdotmetal_in = (POWER_TIME_FACTOR * self.powerMVA * self.powerfactor 
                        / self.metalSER)
        vdotmetal_in = mdotmetal_in / self.densitymetal
        vdotslag_in = mdotmetal_in*self.slagmetalmassratio / self.densityslag
        dhmetal = dt * vdotmetal_in / (self.activearea * self.bedporosity)
        dhslag = dt * vdotslag_in / (self.activearea * self.bedporosity)
        self.hmetal += dhmetal
        self.hslag += dhmetal + dhslag
        self.calc_vdot_out()
        dhmetal = -dt*self.vdotmetal_out / (self.activearea*self.bedporosity)
        dhslag = -dt*self.vdotslag_out / (self.activearea*self.bedporosity)
        self.hmetal += dhmetal
        self.hslag += dhmetal + dhslag