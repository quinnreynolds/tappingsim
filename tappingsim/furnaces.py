import numpy
from scipy.constants import g, pi
POWER_TIME_FACTOR = 1000/3600


def frictionfactor_darcy(diameter, roughness, velocity, density, viscosity):
    d_over_e = diameter / roughness
    Re = diameter * velocity * density / viscosity
    a = 1 / (1 + (Re/2712) ** 8.4)
    b = 1 / (1 + (Re/(150*d_over_e)) ** 1.8)
    return ((64 / Re) ** a 
            * (0.75 * numpy.log(Re/5.37)) ** (2*(a-1)*b) 
            * (0.88 * numpy.log(3.41*d_over_e)) ** (2*(a-1)*(1-b)))


class SubmergedArcFurnace():
    def __init__(self, powerMW, metalSER, slagmetalmassratio, activearea, 
                 tapholediameter, tapholelength, tapholeheight, densitymetal, 
                 densityslag, viscositymetal, viscosityslag, particlediameter, 
                 particlesphericity, bedporosity, bedmaxradius, bedmodel, 
                 entrykl, channelfd, bedentryzone_yn, channellosses_yn, 
                 hmetal_init, hslag_init):
        self.activearea = activearea
        self.tapholediameter = tapholediameter
        self.tapholelength = tapholelength
        self.tapholeheight = tapholeheight
        self.densitymetal = densitymetal
        self.densityslag = densityslag
        self.viscositymetal = viscositymetal
        self.viscosityslag = viscosityslag
        self.particlediameter = particlediameter
        self.particlesphericity = particlesphericity
        self.bedporosity = bedporosity
        self.bedmaxradius = bedmaxradius
        self.bedmodel = bedmodel
        self.entrykl = entrykl
        self.channelfd = channelfd
        self.bedentryzone_yn = bedentryzone_yn
        self.channellosses_yn = channellosses_yn
        self.powerMW = powerMW
        self.metalSER = metalSER
        self.slagmetalmassratio = slagmetalmassratio
        self.hmetal = hmetal_init
        self.hslag = hslag_init
        self.tapholeopen_yn = False
        self.allownegativeheights_yn = False
        self.timetotaliser = 0
        self.powertotaliserkWh = 0
        self.tapmasstotaliser = 0
        self.tapvolumetotaliser = 0
    
    def toggle_tapholeopen(self):
        self.tapholeopen_yn = not self.tapholeopen_yn
        
    def toggle_allownegativeheights(self):
        self.allownegativeheights_yn = not self.allownegativeheights_yn

    def reset_time_totaliser(self):
        self.timetotaliser = 0
        
    def reset_power_totaliser(self):
        self.powertotaliserkWh = 0
        
    def reset_mass_totaliser(self):
        self.tapmasstotaliser = 0

    def reset_volume_totaliser(self):
        self.tapvolumetotaliser = 0

    def calc_vdot_out(self):
        """Return outlet flowrates as a function of the current furnace state. 
        Extended model with semi-empirical interface deformation near tap-hole 
        entry.
        """
        rt = 0.5 * self.tapholediameter
        relrad = rt / self.bedmaxradius
        eff_d = self.particlesphericity * self.particlediameter
        a_m = 0.5*(1+self.entrykl)*self.densitymetal
        a_s = 0.5*(1+self.entrykl)*self.densityslag        
        if self.channellosses_yn:
            a_m += (0.5*self.densitymetal*self.channelfd
                    * self.tapholelength/self.tapholediameter)
            a_s += (0.5*self.densityslag*self.channelfd
                    * self.tapholelength/self.tapholediameter)
            
        # Bernoulli coefficient calculation
        if self.bedmodel == 'ergun':
            if self.bedentryzone_yn:
                a_m += (1.75 * self.densitymetal * rt * (1-self.bedporosity) / 
                        (3*eff_d*self.bedporosity**3) * (1 - 0.25 * relrad**3))
                a_s += (1.75 * self.densityslag * rt * (1-self.bedporosity) / 
                        (3*eff_d*self.bedporosity**3) * (1 - 0.25 * relrad**3))
                b_m = (150*self.viscositymetal*rt*(1-self.bedporosity)**2 / 
                       (eff_d**2 * self.bedporosity**3) * (1 - 0.5 * relrad))
                b_s = (150*self.viscosityslag*rt*(1-self.bedporosity)**2 / 
                       (eff_d**2 * self.bedporosity**3) * (1 - 0.5 * relrad))
            else:
                a_m += (1.75 * self.densitymetal * rt * (1-self.bedporosity) / 
                        (12*eff_d*self.bedporosity**3) * (1 - relrad**3))
                a_s += (1.75 * self.densityslag * rt * (1-self.bedporosity) / 
                       (12*eff_d*self.bedporosity**3) * (1 - relrad**3))
                b_m = (150*self.viscositymetal*rt*(1-self.bedporosity)**2 / 
                       (2*eff_d**2 * self.bedporosity**3) * (1 - relrad))
                b_s = (150*self.viscosityslag*rt*(1-self.bedporosity)**2 / 
                       (2*eff_d**2 * self.bedporosity**3) * (1 - relrad))
        elif self.bedmodel == 'carmenkozeny':
            if self.bedentryzone_yn:
                b_m = (180*self.viscositymetal*rt*(1-self.bedporosity)**2 / 
                       (eff_d**2 * self.bedporosity**3) * (1 - 0.5 * relrad))
                b_s = (180*self.viscosityslag*rt*(1-self.bedporosity)**2 / 
                       (eff_d**2 * self.bedporosity**3) * (1 - 0.5 * relrad))
            else:
                b_m = (180*self.viscositymetal*rt*(1-self.bedporosity)**2 / 
                       (2*eff_d**2 * self.bedporosity**3) * (1 - relrad))
                b_s = (180*self.viscosityslag*rt*(1-self.bedporosity)**2 / 
                       (2*eff_d**2 * self.bedporosity**3) * (1 - relrad))
        else:
            raise NotImplementedError('Bed model ' + str(self.bedmodel) + 
                                      ' not implemented. Please use one of '
                                      '"ergun" or "carmenkozeny"')
            
        # Tapping pressure calculation (applied to both phases, based on 
        # undeformed interfaces). Integrated pressure field in z, averaged with 
        # total height in taphole.
        hs, hm = self.hslag-self.tapholeheight, self.hmetal-self.tapholeheight
        if hs < -rt:
            # liquid level below taphole - no flow
            pa = 0
        elif hs < rt:
            # liquid level between taphole limits
            if hm < -rt:
                # slag only
                pa = 0.5 * (hs+rt) * self.densityslag * g
            else:
                # slag and metal
                pa = 0.5*g*((self.densitymetal * (hm+rt)**2 + 
                             self.densityslag * (hs*(hs+2*rt) - hm*(hm+2*rt))) 
                             / (hs+rt))
        else:
            # liquid level above taphole
            if hm < -rt:
                # slag only
                pa = hs * self.densityslag * g
            elif hm < rt:
                # slag and metal
                pa = 0.25*g*((self.densitymetal * (hm+rt)**2 + 
                              self.densityslag * (rt*(4*hs-rt) - hm*(hm+2*rt))) 
                              / rt)
            else:
                # metal only
                pa = g * (self.densityslag * (hs-hm) + self.densitymetal * hm)
        
        # estimate single-phase-only velocities
        umetal = (-b_m+numpy.sqrt(b_m**2 + 4*pa*a_m)) / (2*a_m)
        uslag = (-b_s+numpy.sqrt(b_s**2 + 4*pa*a_s)) / (2*a_s)
        h0_l = -rt - self.densityslag*uslag**2 / (8*g*(self.densitymetal - 
                                                       self.densityslag))
        h0_h = rt + self.densitymetal*umetal**2 / (8*g*(self.densitymetal - 
                                                        self.densityslag))
        
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
        vdotmetal = area_m * umetal
        vdotslag = area_s * uslag
        
        if not self.allownegativeheights_yn:
            if self.hslag <= 0:
                vdotslag = 0
            if self.hmetal <= 0:
                vdotmetal = 0
        if not self.tapholeopen_yn:
            vdotmetal, vdotslag = 0, 0
        
        self.vdotmetal_out, self.vdotslag_out = vdotmetal, vdotslag
        self.umetal_out, self.uslag_out = umetal, uslag
    
    def calc_dt(self, dt):
        mdotmetal_in = POWER_TIME_FACTOR * self.powerMW / self.metalSER
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
        
        self.timetotaliser += dt
        self.powertotaliserkWh += POWER_TIME_FACTOR * dt * self.powerMW
        self.tapmasstotaliser += dt * (self.vdotmetal_out * self.densitymetal + 
                                       self.vdotslag_out * self.densityslag)
        self.tapvolumetotaliser += dt * (self.vdotmetal_out + self.vdotslag_out)