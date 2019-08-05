# Submodule for simulation of FeMn ladle filling.
#
# Q Reynolds 2019

import numpy as np

g = 9.81

class FeMnFurnace():
    def __init__(self, activearea, tapholediameter, kl, densitymetal, 
                 densityslag, viscositymetal, viscosityslag, particlediameter, 
                 bedporosity, hmetal_init, hslag_init):
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
        self.bedporosity = bedporosity
        self.tapholeopen_yn = False
    
    def calc_vdot_out(self):
        """Return outlet flowrates as a function of the current furnace state. 
        Simplified model without interface deformation near tap-hole entry.
        """
        rt = 0.5 * self.tapholediameter
        a_m = (150 * self.viscositymetal * rt * (1-self.bedporosity)**2 / 
               (self.particlediameter**2 * self.bedporosity**3) +
               0.5*(1+self.kl)*self.densitymetal)
        a_s = (150 * self.viscosityslag * rt * (1-self.bedporosity)**2 / 
               (self.particlediameter**2 * self.bedporosity**3) +
               0.5*(1+self.kl)*self.densityslag)
        b = (1.75 * rt * (1-self.bedporosity) / 
             (3 * self.particlediameter * self.bedporosity**3))
        
        if self.hslag < -rt:
            # slag below taphole - no flow
            pa_m = 0
            pa_s = 0
            area_m = 0
            area_s = 0
        elif self.hslag < rt:
            if self.hmetal < -rt:
                # only slag, partially filled taphole
                pa_m = 0
                pa_s = 0.5 * (self.hslag + rt) * self.densityslag * g
                theta = 2*np.arccos(1 - (self.hslag + rt)/rt)
                area_m = 0
                area_s = 0.5*rt**2 * (theta - np.sin(theta)) 
            else:
                # slag and metal, partially filled taphpole
                pa_m = 0.5 * g * (2*(self.hslag-self.hmetal)*self.densityslag +
                                  (self.hmetal+rt)*self.densitymetal)
                pa_s = 0.5 * (self.hslag-self.hmetal) * self.densityslag * g
                theta = 2*np.arccos(1 - (self.hslag + rt)/rt)
                area_tot = 0.5*rt**2 * (theta - np.sin(theta)) 
                theta = 2*np.arccos(1 - (self.hmetal + rt)/rt)                
                area_m = 0.5*rt**2 * (theta - np.sin(theta)) 
                area_s = area_tot - area_m
        else:
            if self.hmetal < -rt:
                # only slag, filled taphole
                pa_m = 0
                pa_s = self.hslag * self.densityslag * g
                area_m = 0
                area_s = 2*np.pi*rt**2
            elif self.hmetal < rt:
                # slag and metal, filled taphole
                pa_m = 0.5 * g * (2*(self.hslag-self.hmetal)*self.densityslag +
                                  (self.hmetal+rt)*self.densitymetal)
                pa_s = 0.5 * g * ((self.hslag-self.hmetal) * self.densityslag + 
                                  (self.hslag-rt) * self.densityslag)
                theta = 2*np.arccos(1 - (self.hmetal + rt)/rt)
                area_m = 0.5*rt**2 * (theta - np.sin(theta)) 
                area_s = np.pi*rt**2 - area_m
            else:
                # metal only, filled taphole
                pa_m = g * ((self.hslag-self.hmetal) * self.densityslag + 
                            self.hmetal * self.densitymetal)
                pa_s = 0
                area_m = np.pi*rt**2
                area_s = 0
        
        if not self.tapholeopen_yn:
            pa_m = 0
            pa_s = 0
            area_m = 0
            area_s = 0
            
        umetal = (-b+np.sqrt(b**2 + 4*pa_m*a_m)) / (2*a_m)
        uslag = (-b+np.sqrt(b**2 + 4*pa_s*a_s)) / (2*a_s)
        vdot_metal = area_m * umetal
        vdot_slag = area_s * uslag
        
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
    def __init__(self, diameter, height, hmetal_init, hslag_init):
        self.diameter = diameter
        self.height = height
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
        if hsi < self.height:
            vdot_slag = 0
            vdot_metal = 0
        else:
            if hmi < self.height:
                vdot_slag = (hsi - self.height) * self.xarea / dt
                vdot_metal = 0
            else:
                vdot_slag = (hsi - hmi) * self.xarea / dt
                vdot_metal = (hmi - self.height) * self.xarea / dt
        self.vdotmetal, self.vdotslag = vdot_metal, vdot_slag
        return vdot_metal, vdot_slag
    
    def calc_dt(self, dt, vdot_metal_in, vdot_slag_in):
        self.volmetal += vdot_metal_in * dt
        self.volslag += vdot_slag_in * dt
        vdm, vds = self.calc_vdot_out(dt)
        self.volmetal -= vdm * dt
        self.volslag -= vds * dt
        self.calc_interfaces()
    
