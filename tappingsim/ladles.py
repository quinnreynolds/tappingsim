import numpy
from scipy.constants import pi

def overflowmodel_step(interfacedeltah, *consts):
    """Function to calculate metal volume fraction entrained in slag outflow 
    stream as a function of the location of the slag-metal interface.
    
    Parameters
    ----------
    interfacedeltah : float
        Not used in this overflow model.
    consts : list of float
        Not used in this overflow model.
        
    Returns
    -------
    float
        Metal volume fraction in slag outlet stream.
        
    Notes
    -----
    This model describes perfectly separated phase model with no metal 
    entrainment into the slag phase.

    Fraction = 0
    """
    return 0

def overflowmodel_exp(interfacedeltah, *consts):
    """Function to calculate metal volume fraction entrained in slag outflow 
    stream as a function of the location of the slag-metal interface.
    
    Parameters
    ----------
    interfacedeltah : float
        The distance of the slag-metal interface from the ladle outlet, ladle 
        depth minus interface position, in m.
    consts : list of float
        The exponential pre-multiplier and exponential factor, dimensionless 
        and 1/m respectively.

    Returns
    -------
    float
        Metal volume fraction in slag outlet stream.

    Notes
    -----
    This model describes an exponential decay of metal entrainment as a 
    function of the position of the interface relative to the outlet. [1]_
    
    Fraction = consts[0]*exp(-consts[1]*interfacedeltah) if interfacedeltah > 0
             = consts[0]                                 if interfacedeltah < 0
    
    References
    ----------
    .. [1] Q.G. Reynolds, J.E. Olsen, J.D. Steenkamp, Variability in Ferroalloy 
    Furnace Tapping - Insights from Modelling. Proceedings of the 16th 
    International Ferro-Alloys Congress (INFACON XVI) 2021, Available at 
    SSRN: https://ssrn.com/abstract=3926222 or 
    http://dx.doi.org/10.2139/ssrn.3926222.
    """
    if interfacedeltah > 0:
        return consts[0]*numpy.exp(-consts[1]*interfacedeltah)
    else:
        return consts[0]


class CylindricalLadle():
    def __init__(self, diameter, depth, hmetal_init, hslag_init,
                 overflowmodel, overflowconsts):
        """Tapping ladle class. This version is a cylindrical ladle with an 
        choice of empirical overflow models.
        
        Parameters
        ----------
        diameter : float
            The internal diameter of the ladle, m.
        depth : float
            The internal depth of the ladle, m.
        hmetal_init : float
            The initial level of metal in the ladle, m.
        hslag_init : float
            The initial level of slag in the ladle, m.
        overflowmodel : function
            Function to be used to model the phase overflow behaviour.
        overflowconsts : list of float
            Empirical parameters for the chosen overflow model.
                
        Attributes
        ----------
        diameter : float
            The internal diameter of the ladle, m.
        depth : float
            The internal depth of the ladle, m.
        hmetal : float
            The current level of metal in the ladle, m.
        hslag : float
            The current level of slag in the ladle, m.
        overflowmodel : function
            Internal lambda function created for overflow modelling.
        vdotmetal_out : float
            The current outlet volume flowrate of metal from the unit, m3/s.
        vdotslag_out : float
            The current outlet volume flowrate of slag from the unit, m3/s.
        """
        self.diameter = diameter
        self.depth = depth
        self.hmetal = hmetal_init
        self.hslag = hslag_init
        self.overflowmodel = lambda dh: overflowmodel(dh, *overflowconsts)
        self.diameter = self.diameter
        
    def _calc_vdot_out(self, dt):
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
        """
        Integrate object state forward over a single time step, and update 
        output flowrates.

        Parameters
        ----------
        dt : float
            Length of time step, in s.
        vdotmetal_in : float
            Flowrate of metal into the ladle, kg/s.
        vdotslag_in : float
            Flowrate of slag into the ladle, kg/s.

        Returns
        -------
        None.
        """
        xarea = 0.25 * pi * self.diameter**2
        dhmetal = dt*vdotmetal_in / xarea
        dhslag = dt*vdotslag_in / xarea
        self.hmetal += dhmetal
        self.hslag += dhmetal + dhslag
        self._calc_vdot_out(dt)
        dhmetal = -dt*self.vdotmetal_out / xarea
        dhslag = -dt*self.vdotslag_out / xarea
        self.hmetal += dhmetal
        self.hslag += dhmetal + dhslag
