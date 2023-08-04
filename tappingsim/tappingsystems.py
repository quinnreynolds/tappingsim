import numpy

POWER_TIME_FACTOR = 1000/3600


class SAF():
    def __init__(self, furnace):
        self.furnace = furnace
        self.timetotaliser = 0
        self.powertotaliserkWh = 0
        self.metalmasstotaliser = 0
        self.slagmasstotaliser = 0
    
    def open_taphole(self):
        self.furnace.tapholeopen_yn = True
        
    def close_taphole(self):
        self.furnace.tapholeopen_yn = False
    
    def reset_time_totaliser(self):
        self.timetotaliser = 0
        
    def reset_power_totaliser(self):
        self.powertotaliserkWh = 0
        
    def reset_mass_totaliser(self):
        self.metalmasstotaliser = 0
        self.slagmasstotaliser = 0
    
    def reset_all_totalisers(self):
        self.reset_time_totaliser()
        self.reset_power_totaliser()
        self.reset_mass_totaliser()
        
    def calc_dt(self, dt):
        """
        Integrate model over a single time step.

        Parameters
        ----------
        dt : float
            Length of time step, in s.

        Returns
        -------
        None.

        """
        self.furnace.calc_dt(dt)    
        self.timetotaliser += dt
        self.powertotaliserkWh += dt * (POWER_TIME_FACTOR 
                                        * self.furnace.powerMVA 
                                        * self.furnace.powerfactor)
        self.metalmasstotaliser += dt * (self.furnace.vdotmetal_out 
                                          * self.furnace.densitymetal)
        self.slagmasstotaliser += dt * (self.furnace.vdotslag_out 
                                          * self.furnace.densityslag)

    def calc_time_period(self, times):
        """
        Integrate model over a series of time steps. The model parameters are 
        assumed to remain fixed during this period.

        Parameters
        ----------
        times : 1D numpy array
            Time steps at which to perform calculation, in s.

        Returns
        -------
        None.

        """
        # TODO also initialise return arrays containing state vars during run
        mm, sm = self.furnace.calc_time_period(times)      
        deltat = times[-1] - times[0]
        self.timetotaliser += deltat
        self.powertotaliserkWh += deltat * (POWER_TIME_FACTOR 
                                            * self.furnace.powerMVA 
                                            * self.furnace.powerfactor)
        self.metalmasstotaliser += mm
        self.slagmasstotaliser += sm


class SAFWithLadles():
    def __init__(self, furnace, launder, ladles, taptotaptime, tappingtime):
        self.furnace = furnace
        self.launder = launder
        self.ladles = ladles
        self.taptotaptime = taptotaptime
        self.tappingtime = tappingtime
        self.timetotaliser = 0
        self.powertotaliserkWh = 0
        self.metalmasstotaliser = 0
        self.slagmasstotaliser = 0

    def open_taphole(self):
        self.furnace.tapholeopen_yn = True
        
    def close_taphole(self):
        self.furnace.tapholeopen_yn = False

    def empty_ladles(self):
        for ldl in self.ladles:
            ldl.hmetal, ldl.hslag = 0, 0
            
    def reset_time_totaliser(self):
        self.timetotaliser = 0
        
    def reset_power_totaliser(self):
        self.powertotaliserkWh = 0
        
    def reset_mass_totaliser(self):
        self.metalmasstotaliser = 0
        self.slagmasstotaliser = 0
        
    def ladle_masses(self):
        metalmasses, slagmasses = [], []
        for ldl in self.ladles:
            area = numpy.pi * 0.25 * ldl.diameter**2
            metalmasses.append(area * ldl.hmetal * self.furnace.densitymetal)
            slagmasses.append(area * (ldl.hslag-ldl.hmetal) 
                              * self.furnace.densityslag)
        return metalmasses, slagmasses
    
    def calc_dt(self, dt):
        self.furnace.calc_dt(dt)
        if self.furnace.tapholeopen_yn:
            self.launder.calc_dt(dt, self.furnace.vdotmetal_out, 
                                 self.furnace.vdotslag_out)
            self.ladles[0].calc_dt(dt, self.launder.vdotmetal_out, 
                                   self.launder.vdotslag_out)
            for lprev, lnext in zip(self.ladles[:-1], self.ladles[1:]):
                lnext.calc_dt(dt, lprev.vdotmetal_out, lprev.vdotslag_out)

        self.timetotaliser += dt
        self.powertotaliserkWh += dt * (POWER_TIME_FACTOR 
                                        * self.furnace.powerMVA 
                                        * self.furnace.powerfactor)
        self.metalmasstotaliser += dt * (self.furnace.vdotmetal_out 
                                         * self.furnace.densitymetal)
        self.slagmasstotaliser += dt * (self.furnace.vdotslag_out 
                                         * self.furnace.densityslag)
