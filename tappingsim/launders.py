"""Classes to create transfer launder objects. 

Launders are typically open-channel gravity flow conduits which redirect the 
flow of tapped material from the furnace into storage ladles or other 
post-treatment units.

"""

class SimpleSisoLaunder():
    """Tapping launder class. This version is a simple single-in-single-out
    launder channel with no accumulation or spillage effects.
    
    Parameters
    ----------
    None.
                
    Attributes
    ----------
    vdotmetal_out : float
        The output metal flowrate from the launder in the current state, 
        kg/s
    vdotslag_out : float
        The output slag flowrate from the launder in the current state, 
        kg/s
        
    """        
    def __init__(self):
        self.vdotmetal_out = 0
        self.vdotslag_out = 0
            
    def calc_dt(self, dt, vdotmetal_in, vdotslag_in):
        """
        Integrate object state forward over a single time step, and update 
        output flowrates.

        Parameters
        ----------
        dt : float
            Length of time step, in s.
        vdotmetal_in : float
            Flowrate of metal into the launder, kg/s.
        vdotslag_in : float
            Flowrate of slag into the launder, kg/s.

        Returns
        -------
        None.
        
        """
        self.vdotmetal_out = vdotmetal_in
        self.vdotslag_out = vdotslag_in
