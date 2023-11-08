import numpy


class SAF():
    def __init__(self, furnace):
        self.furnace = furnace
        self.timetotaliser = 0
        self.powertotaliserkWh = 0
        self.metalmasstotaliser = 0
        self.slagmasstotaliser = 0
    
    def open_taphole(self):
        """
        Set furnace tap-hole to open
        """
        self.furnace.tapholeopen_yn = True
        
    def close_taphole(self):
        """
        Set furnace tap-hole to closed.
        """
        self.furnace.tapholeopen_yn = False
    
    def reset_time_totaliser(self):
        """
        Reset time totaliser counter to zero.
        """
        self.timetotaliser = 0
        
    def reset_power_totaliser(self):
        """
        Reset power totaliser counter to zero.
        """
        self.powertotaliserkWh = 0
        
    def reset_mass_totaliser(self):
        """
        Reset mass totaliser counters to zero.
        """
        self.metalmasstotaliser = 0
        self.slagmasstotaliser = 0
    
    def reset_all_totalisers(self):
        """
        Reset all totaliser counters to zero.
        """
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
        self.powertotaliserkWh += dt * (self.furnace.powerMVA 
                                        * self.furnace.powerfactor / 3.6)
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
        mm, sm = self.furnace.calc_time_period(times)      
        deltat = times[-1] - times[0]
        self.timetotaliser += deltat
        self.powertotaliserkWh += deltat * (self.furnace.powerMVA 
                                            * self.furnace.powerfactor / 3.6)
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
        """
        Set furnace tap-hole to open.
        """
        self.furnace.tapholeopen_yn = True
        
    def close_taphole(self):
        """
        Set furnace tap-hole to closed.
        """
        self.furnace.tapholeopen_yn = False

    def empty_ladles(self):
        """
        Set slag and metal levels in all ladles to zero.
        """
        for ldl in self.ladles:
            ldl.hmetal, ldl.hslag = 0, 0
            
    def reset_time_totaliser(self):
        """
        Reset time totaliser counter to zero.
        """
        self.timetotaliser = 0
        
    def reset_power_totaliser(self):
        """
        Reset power totaliser counter to zero.
        """
        self.powertotaliserkWh = 0
        
    def reset_mass_totaliser(self):
        """
        Reset mass totaliser counters to zero.
        """
        self.metalmasstotaliser = 0
        self.slagmasstotaliser = 0

    def reset_all_totalisers(self):
        """
        Reset all totaliser counters to zero.
        """
        self.reset_time_totaliser()
        self.reset_power_totaliser()
        self.reset_mass_totaliser()
        
    def ladle_masses(self):
        """
        Calculate masses of metal and slag in each ladle.

        Parameters
        ----------
        None.
        
        Returns
        -------
        metalmasses : list of float
            Mass of metal in each ladle, in kg.
        slagmasses : list of float
            Mass of slag in each ladle, in kg.
        """
        metalmasses, slagmasses = [], []
        for ldl in self.ladles:
            area = numpy.pi * 0.25 * ldl.diameter**2
            metalmasses.append(area * ldl.hmetal * self.furnace.densitymetal)
            slagmasses.append(area * (ldl.hslag-ldl.hmetal) 
                              * self.furnace.densityslag)
        return metalmasses, slagmasses
    
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
        if self.furnace.tapholeopen_yn:
            self.launder.calc_dt(dt, self.furnace.vdotmetal_out, 
                                 self.furnace.vdotslag_out)
            self.ladles[0].calc_dt(dt, self.launder.vdotmetal_out, 
                                   self.launder.vdotslag_out)
            for lprev, lnext in zip(self.ladles[:-1], self.ladles[1:]):
                lnext.calc_dt(dt, lprev.vdotmetal_out, lprev.vdotslag_out)

        self.timetotaliser += dt
        self.powertotaliserkWh += dt * (self.furnace.powerMVA 
                                        * self.furnace.powerfactor) / 3.6
        self.metalmasstotaliser += dt * (self.furnace.vdotmetal_out 
                                         * self.furnace.densitymetal)
        self.slagmasstotaliser += dt * (self.furnace.vdotslag_out 
                                         * self.furnace.densityslag)
