# Submodule for simulation of FeMn ladle filling.
#
# Q Reynolds 2019

import numpy as np

g = 9.81

class FeMnFurnace():
    def __init__(self, activearea, tapholediameter, kl, hmetal_init, hslag_init,
                 densitymetal, densityslag):
        self.activearea = activearea
        self.tapholediameter = tapholediameter
        self.kl = kl
        self.hmetal = hmetal_init
        self.hslag = hslag_init
        self.densitymetal = densitymetal
        self.densityslag = densityslag
    
    def calc_vdot_out(self):
        """Return outlet flowrates as a function of the current furnace state.
        """
        if self.hmetal > 0:
            delta_p = ((self.hslag-self.hmetal) * self.densityslag + 
                       self.hmetal * self.densitymetal) * g
        else:
            delta_p = self.hslag * self.densityslag * g
            
        if self.hmetal < -0.5*self.tapholediameter:
            #interface completely below taphole
            uslag = np.sqrt(2 * delta_p / (self.densityslag * (1 + self.kl)))
            vdot_metal = 0
            vdot_slag = 0.25 * np.pi * self.tapholediameter * uslag
        elif self.hmetal > 0.5 * self.tapholediameter:
            #interface completely above taphole
            umetal = np.sqrt(2 * delta_p / (self.densitymetal * (1 + self.kl)))
            vdot_metal = 0.25 * np.pi * self.tapholediameter * umetal
            vdot_slag = 0
        else:
            #interface in taphole
            umetal = np.sqrt(2 * delta_p / (self.densitymetal * (1 + self.kl)))
            uslag = np.sqrt(2 * delta_p / (self.densityslag * (1 + self.kl)))
            segheight = self.hmetal + 0.5 * self.tapholediameter
            theta = 2*np.arccos(1 - segheight/(0.5*self.tapholediameter))
            metalarea = 0.25*self.tapholediameter**2 * (theta - np.sin(theta)) 
            slagarea = 0.25*np.pi*self.tapholediameter**2 - metalarea
            vdot_metal = metalarea * umetal
            vdot_slag = slagarea * uslag
        return vdot_metal, vdot_slag
    
    def calc_dt(self, dt, vdot_metal_in, vdot_slag_in):
        self.volmetal += vdot_metal_in * dt
        self.volslag += vdot_slag_in * dt
        vdm, vds = self.calc_vdot_out(dt)
        self.volmetal -= vdm * dt
        self.volslag -= vds * dt
        self.calc_interfaces()

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
        """
        hmi, hsi = self.calc_interfaces()
        if hsi < self.height:
            vdot_slag = 0
            vdot_metal = 0
        else:
            # the following section can be modified to reflect more realisitic
            # carry-over effects.
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
    
