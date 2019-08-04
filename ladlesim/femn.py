# Submodule for simulation of FeMn ladle filling.
#
# Q Reynolds 2019

import numpy as np

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
    
