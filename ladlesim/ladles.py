# Submodule for simulation of FeMn ladle filling.
#
# Q Reynolds 2019

from scipy.constants import pi

class FeMnLadle():
    def __init__(self, diameter, depth, hmetal_init, hslag_init):
        self.diameter = diameter
        self.depth = depth
        self.xarea = 0.25 * pi * self.diameter**2
        self.hmetal = hmetal_init
        self.hslag = hslag_init
        
    def empty_ladle(self):
        self.hmetal, self.hslag = 0, 0
        
    def calc_volumes(self):
        """Calculate slag and metal volume in ladle as a function of the current
        ladle state.
        """
        vm = self.hmetal * self.xarea
        vs = (self.hslag - self.hmetal) * self.xarea
        return vm, vs
    
    def calc_vdot_out(self, dt):
        """Calculate outlet flowrates as a function of the current ladle state.
        Simplified model with no carry-over of metal until the interface reaches
        the ladle outlet.
        """
        if self.hslag < self.depth:
            self.vdotmetal_out = 0
            self.vdotslag_out = 0
        else:
            if self.hmetal < self.depth:
                self.vdotmetal_out = 0
                self.vdotslag_out = (self.hslag-self.depth) * self.xarea / dt
            else:
                self.vdotmetal_out = (self.hmetal-self.depth) * self.xarea / dt
                self.vdotslag_out = (self.hslag-self.hmetal) * self.xarea / dt
    
    def calc_dt(self, dt, vdotmetal_in, vdotslag_in):
        dhmetal = dt*vdotmetal_in / self.xarea
        dhslag = dt*vdotslag_in / self.xarea
        self.hmetal += dhmetal
        self.hslag += dhmetal + dhslag
        self.calc_vdot_out(dt)
        dhmetal = -dt*self.vdotmetal_out / self.xarea
        dhslag = -dt*self.vdotslag_out / self.xarea
        self.hmetal += dhmetal
        self.hslag += dhmetal + dhslag
