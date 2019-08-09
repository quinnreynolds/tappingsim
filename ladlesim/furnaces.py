# Submodule for simulation of FeMn ladle filling.
#
# Q Reynolds 2019

import numpy as np
from scipy.constants import g

class FeMnFurnace():
    def __init__(self, activearea, tapholediameter, kl, densitymetal, 
                 densityslag, viscositymetal, viscosityslag, particlediameter, 
                 particlesphericity, bedporosity, hmetal_init, hslag_init):
        self.activearea = activearea
        self.tapholediameter = tapholediameter
        self.kl = kl
        self.hmetal = hmetal_init
        self.hslag = hslag_init
        self.densitymetal = densitymetal
        self.densityslag = densityslag
        self.viscositymetal = viscositymetal
        self.viscosityslag = viscosityslag
        self.particlediameter = particlediameter
        self.particlesphericity = particlesphericity
        self.bedporosity = bedporosity
        self.tapholeopen_yn = False
    
    def calc_vdot_out(self):
        """Return outlet flowrates as a function of the current furnace state. 
        Extended model with semi-empirical interface deformation near tap-hole 
        entry.
        """
        rt = 0.5 * self.tapholediameter
        eff_d = self.particlesphericity * self.particlediameter
        a_m = (150 * self.viscositymetal * rt * (1-self.bedporosity)**2 / 
               (eff_d**2 * self.bedporosity**3) +
               0.5*(1+self.kl)*self.densitymetal)
        a_s = (150 * self.viscosityslag * rt * (1-self.bedporosity)**2 / 
               (eff_d**2 * self.bedporosity**3) +
               0.5*(1+self.kl)*self.densityslag)
        b_m = (1.75 * self.densitymetal * rt * (1-self.bedporosity) / 
               (3 * eff_d * self.bedporosity**3))
        b_s = (1.75 * self.densityslag * rt * (1-self.bedporosity) / 
               (3 * eff_d * self.bedporosity**3))
        
        # tapping pressure calculation (applied to both phases, based on 
        # undeformed interfaces). Integrated pressure field in z, averaged with 
        # total height in taphole.
        if self.hslag < -rt:
            # liquid level below taphole - no flow
            pa = 0
        elif self.hslag < rt:
            # liquid level between taphole limits
            if self.hmetal < -rt:
                # slag only
                pa = 0.5 * (self.hslag + rt) * self.densityslag * g
            else:
                # slag and metal
                pa = 0.5*g*(self.densitymetal*(self.hmetal+rt)**2 + 
                            self.densityslag*(self.hslag*(self.hslag+2*rt) - 
                            self.hmetal*(self.hmetal+2*rt))) / (self.hslag+rt)
        else:
            # liquid level above taphole
            if self.hmetal < -rt:
                # slag only
                pa = self.hslag * self.densityslag * g
            elif self.hmetal < rt:
                # slag and metal
                pa = 0.25*g*(self.densitymetal*(self.hmetal+rt)**2 + 
                             self.densityslag*(rt*(4*self.hslag-rt) - 
                                        self.hmetal*(self.hmetal+2*rt))) / rt
            else:
                # metal only
                pa = g * ((self.hslag-self.hmetal) * self.densityslag + 
                            self.hmetal * self.densitymetal)
        
        # estimate single-phase-only velocities
        umetal = (-b_m+np.sqrt(b_m**2 + 4*pa*a_m)) / (2*a_m)
        uslag = (-b_s+np.sqrt(b_s**2 + 4*pa*a_s)) / (2*a_s)
        h0_l = -rt - self.densityslag*uslag**2 / (8*g*(self.densitymetal - 
                                                       self.densityslag))
        h0_h = rt + self.densitymetal*umetal**2 / (8*g*(self.densitymetal - 
                                                        self.densityslag))
        
        if self.hslag < -rt:
            # slag below taphole - no flow
            area_m = 0
            area_s = 0
        elif self.hslag < rt:
            # taphole partially filled
            if self.hmetal < h0_l:
                # slag only
                theta = 2*np.arccos(-self.hslag/rt)
                area_m = 0
                area_s = 0.5*rt**2*(theta - np.sin(theta))
            else:
                # slag and metal
                hi = self.hslag - ((self.hslag-self.hmetal) * (self.hslag+rt) / 
                                   (self.hslag-h0_l))
                theta = 2*np.arccos(-hi/rt)
                area_m = 0.5*rt**2*(theta - np.sin(theta))
                theta = 2*np.arccos(-self.hslag/rt)
                area_s = 0.5*rt**2*(theta - np.sin(theta)) - area_m
        else:
            # taphole completely filled
            if self.hmetal < h0_l:
                # slag only
                area_m = 0
                area_s = np.pi*rt**2
            elif self.hmetal > h0_h:
                # metal only
                area_m = np.pi*rt**2
                area_s = 0
            else:
                # slag and metal
                hi = rt*(1 - 2*(h0_h-self.hmetal)/(h0_h-h0_l))
                theta = 2*np.arccos(-hi/rt)
                area_m = 0.5*rt**2*(theta - np.sin(theta))
                area_s = np.pi*rt**2 - area_m
                
        if not self.tapholeopen_yn:
            area_m = 0
            area_s = 0

        vdot_metal = area_m * umetal
        vdot_slag = area_s * uslag
        self.vdotmetal, self.vdotslag = vdot_metal, vdot_slag
        return vdot_metal, vdot_slag
    
    def calc_dt(self, dt, vdot_metal_in, vdot_slag_in):
        dhmetal = dt * vdot_metal_in / (self.activearea * self.bedporosity)
        dhslag = dt * vdot_slag_in / (self.activearea * self.bedporosity)
        self.hmetal += dhmetal
        self.hslag += dhmetal + dhslag
        vdm, vds = self.calc_vdot_out()
        dhmetal = -dt * (vdm / (self.activearea * self.bedporosity))
        dhslag = -dt * (vds / (self.activearea * self.bedporosity))
        self.hmetal += dhmetal
        self.hslag += dhmetal + dhslag
        return vdm, vds

class FeMnLadle():
    def __init__(self, diameter, depth, hmetal_init, hslag_init):
        self.diameter = diameter
        self.depth = depth
        self.xarea = 0.25 * np.pi * self.diameter**2
        self.volmetal = hmetal_init * self.xarea
        self.volslag = (hslag_init - hmetal_init) * self.xarea
        
    def calc_interfaces(self):
        """Return interface positions (heights) as a function of the current
        ladle state.
        """
        hmi = self.volmetal / self.xarea
        hsi = hmi + self.volslag / self.xarea
        self.hmetal, self.hslag = hmi, hsi
        return hmi, hsi
    
    def calc_vdot_out(self, dt):
        """Return outlet flowrates as a function of the current ladle state.
        Simplified model with no carry-over of metal until the interface reaches
        the ladle outlet.
        """
        hmi, hsi = self.calc_interfaces()
        if hsi < self.depth:
            vdot_slag = 0
            vdot_metal = 0
        else:
            if hmi < self.depth:
                vdot_slag = (hsi - self.depth) * self.xarea / dt
                vdot_metal = 0
            else:
                vdot_slag = (hsi - hmi) * self.xarea / dt
                vdot_metal = (hmi - self.depth) * self.xarea / dt
        self.vdotmetal, self.vdotslag = vdot_metal, vdot_slag
        return vdot_metal, vdot_slag
    
    def calc_dt(self, dt, vdot_metal_in, vdot_slag_in):
        self.volmetal += vdot_metal_in * dt
        self.volslag += vdot_slag_in * dt
        vdm, vds = self.calc_vdot_out(dt)
        self.volmetal -= vdm * dt
        self.volslag -= vds * dt
        self.calc_interfaces()
        return vdm, vds
