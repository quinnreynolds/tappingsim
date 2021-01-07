class SimpleSisoLaunder():
    def __init__(self):
        self.vdotmetal_out = 0
        self.vdotslag_out = 0
            
    def calc_dt(self, dt, vdotmetal_in, vdotslag_in):
        self.vdotmetal_out = vdotmetal_in
        self.vdotslag_out = vdotslag_in
