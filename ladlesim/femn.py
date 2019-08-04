# Submodule for simulation of FeMn ladle filling.
#
# Q Reynolds 2019

import numpy as np
from . import BaseModel

class FeMnLadle(BaseModel):
    def __init__(self, diameter, height, startvol_slag, startvol_metal):
        BaseModel.__init__(self, diameter, height)
        self.time = 0
        self.phasenames = ['metal', 'slag']
        self.volume = np.array([startvol_metal, startvol_slag])
        self.xarea = 0.25 * np.pi * self.diameter
        
    def calc_interfaces(self):
        """Return interface positions (heights) as a function of the current
        ladle state.
        """
        hi = self.volume / self.xarea
        hi[1] += hi[0]
        return hi
    
    def calc_vdot_out(self, dt):
        """Return outlet flowrates as a function of the current ladle state.
        """
        hi = self.calc_interfaces()
        if hi[1] < self.height:
            vdot_slag = 0
            vdot_metal = 0
        else:
            # the following section can be modified to reflect more realisitic
            # carry-over effects.
            if hi[0] < self.height:
                vdot_slag = (hi[1] - self.height) / dt
                vdot_metal = 0
            else:
                vdot_slag = (hi[1] - hi[0]) / dt
                vdot_metal = (hi[0] - self.height) / dt
        return np.array([vdot_metal, vdot_slag])
    
    def calc_dt(self, dt, vdot_in):
        self.volume += vdot_in * dt
        self.volume -= self.calc_vdot_out(dt) * dt
        self.time += dt
    
