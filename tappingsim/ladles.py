import numpy
from scipy.constants import pi

def overflowmodel_step(interfacedeltah, *consts):
    """Calculate metal volume fraction entrained in slag outflow stream as a 
    function of the location of the slag-metal interface (deltah = ladle 
    depth - interface position).
    
    fraction = 0
    """
    return 0

def overflowmodel_exp(interfacedeltah, *consts):
    """Calculate metal volume fraction entrained in slag outflow stream as a 
    function of the location of the slag-metal interface (deltah = ladle 
    depth - interface position).
    
    fraction = consts[0]*exp(-consts[1]*deltah)   if deltah >= 0
             = consts[0]                          if deltah < 0
    """
    if interfacedeltah > 0:
        return consts[0]*numpy.exp(-consts[1]*interfacedeltah)
    else:
        return consts[0]



class CylindricalLadle():
    def __init__(self, diameter, depth, hmetal_init, hslag_init,
                 overflowmodel, overflowconsts):
        self.diameter = diameter
        self.depth = depth
        self.hmetal = hmetal_init
        self.hslag = hslag_init
        self.overflowmodel = lambda dh: overflowmodel(dh, *overflowconsts)
        self.diameter = self.diameter
        
    def calc_vdot_out(self, dt):
        """Calculate outlet flowrates as a function of the current ladle state.
        """
        xarea = 0.25 * pi * self.diameter**2
        if self.hslag < self.depth:
            self.vdotmetal_out = 0
            self.vdotslag_out = 0
        else:
            vfrac = self.overflowmodel(self.depth - self.hmetal)
            if self.hmetal < self.depth:
                vslagout = (self.hslag-self.depth) * xarea
                self.vdotmetal_out = vslagout * vfrac / dt
                self.vdotslag_out = vslagout * (1-vfrac) / dt
            else:
                vslagout = (self.hslag-self.hmetal) * xarea
                self.vdotmetal_out = ((self.hmetal-self.depth) * xarea
                                      + vslagout * vfrac) / dt
                self.vdotslag_out = vslagout * (1-vfrac) / dt
    
    def calc_dt(self, dt, vdotmetal_in, vdotslag_in):
        xarea = 0.25 * pi * self.diameter**2
        dhmetal = dt*vdotmetal_in / xarea
        dhslag = dt*vdotslag_in / xarea
        self.hmetal += dhmetal
        self.hslag += dhmetal + dhslag
        self.calc_vdot_out(dt)
        dhmetal = -dt*self.vdotmetal_out / xarea
        dhslag = -dt*self.vdotslag_out / xarea
        self.hmetal += dhmetal
        self.hslag += dhmetal + dhslag
