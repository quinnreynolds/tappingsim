POWER_TIME_FACTOR = 1/3600

class FeMnSAF():
    def __init__(self, furnace, launder, ladle1, ladle2):
        self.furnace = furnace
        self.launder = launder
        self.ladles = [ladle1, ladle2]
        self.timetotaliser = 0
        self.powertotaliserMWh = 0
        self.metalmasstotaliser = 0
        self.slagmasstotaliser = 0
        self.metalmassfirstladle = 0
        self.slagmassfirstladle = 0

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
        self.powertotaliserMWh = 0
        
    def reset_mass_totaliser(self):
        self.metalmasstotaliser = 0
        self.slagmasstotaliser = 0
        self.metalmassfirstladle = 0
        self.slagmassfirstladle = 0
            
    def calc_dt(self, dt):
        self.furnace.calc_dt(dt)
        self.launder.calc_dt(dt, self.furnace.vdotmetal_out, 
                             self.furnace.vdotslag_out)
        self.ladles[0].calc_dt(dt, self.launder.vdotmetal_out, 
                               self.launder.vdotslag_out)
        for lprev, lnext in zip(self.ladles[:-1], self.ladles[1:]):
            lnext.calc_dt(dt, lprev.vdotmetal_out, lprev.vdotslag_out)

        self.timetotaliser += dt
        self.powertotaliserMWh += dt * (POWER_TIME_FACTOR 
                                        * self.furnace.powerMVA 
                                        * self.furnace.powerfactor)
        self.metalmasstotaliser += dt * (self.furnace.vdotmetal_out 
                                         * self.furnace.densitymetal)
        self.slagmasstotaliser += dt * (self.furnace.vdotslag_out 
                                         * self.furnace.densityslag)
        self.metalmassfirstladle += dt * ((self.furnace.vdotmetal_out
                                           - self.ladles[0].vdotmetal_out)
                                         * self.furnace.densitymetal)
        self.slagmassfirstladle += dt * ((self.furnace.vdotslag_out
                                           - self.ladles[0].vdotslag_out)
                                         * self.furnace.densityslag)
        
                