import numpy
from scipy.constants import pi

def overflowmodel_step(interfacedeltah, *consts):
    """Calculate metal volume fraction at ladle outlet as a function of the 
    location of the slag-metal interface (deltah = ladle depth - interface 
    position).
    
    fraction = 0 if deltah < 0, else 1
    """
    if interfacedeltah > 0:
        return 0
    else:
        return 1

def overflowmodel_exp(interfacedeltah, *consts):
    """Calculate metal volume fraction at ladle outlet as a function of the 
    location of the slag-metal interface (deltah = ladle depth - interface 
    position).
    
    fraction = exp(-consts[0]*deltah)
    """
    if interfacedeltah > 0:
        return numpy.exp(-consts[0]*interfacedeltah)
    else:
        return 1

def overflowmodel_expoffset(interfacedeltah, *consts):
    """Calculate metal volume fraction at ladle outlet as a function of the 
    location of the slag-metal interface (deltah = ladle depth - interface 
    position).
    
    fraction = exp(-consts[0]*(deltah-consts[1])
    """
    if interfacedeltah > consts[1]:
        return numpy.exp(-consts[0]*(interfacedeltah-consts[1]))
    else:
        return 1


class CylindricalLadle():
    def __init__(self, diameter, depth, hmetal_init, hslag_init,
                 overflowmodel, overflowconsts):
        self.diameter = diameter
        self.depth = depth
        self.hmetal = hmetal_init
        self.hslag = hslag_init
        self.overflowmodel = lambda dh: overflowmodel(dh, *overflowconsts)
        self.diameter = self.diameter
        
    def empty_ladle(self):
        self.hmetal, self.hslag = 0, 0
        
    def calc_vdot_out(self, dt):
        """Calculate outlet flowrates as a function of the current ladle state.
        """
        xarea = 0.25 * pi * self.diameter**2
        if self.hslag < self.depth:
            self.vdotmetal_out = 0
            self.vdotslag_out = 0
        else:
            if self.hmetal < self.depth:
                vfrac = self.overflowmodel(self.depth - self.hmetal)
                vout = (self.hslag-self.depth) * xarea
                self.vdotmetal_out = vout * vfrac / dt
                self.vdotslag_out = vout * (1-vfrac) / dt
            else:
                self.vdotmetal_out = (self.hmetal-self.depth) * xarea / dt
                self.vdotslag_out = (self.hslag-self.hmetal) * xarea / dt
    
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
