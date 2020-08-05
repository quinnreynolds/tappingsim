import numpy as np
from scipy.constants import g, pi

class TankWithPorousBed():
    def __init__(self, tankarea, tapholediameter, kl, densityfluid, 
                 viscosityfluid, particlediameter, particlesphericity, 
                 bedporosity, bedheight, hfluid_init):
        self.tankarea = tankarea
        self.tapholediameter = tapholediameter
        self.kl = kl
        self.densityfluid = densityfluid
        self.viscosityfluid = viscosityfluid
        self.particlediameter = particlediameter
        self.particlesphericity = particlesphericity
        self.bedporosity = bedporosity
        self.bedheight = bedheight
        self.hfluid = hfluid_init

    def calc_vdot_out(self):
        """Return outlet flowrates as a function of the current tank state. 
        Simplified model without interface deformation near tap-hole entry.
        """
        rt = 0.5 * self.tapholediameter
        eff_d = self.particlesphericity * self.particlediameter
        a = (1.75 * self.densityfluid * rt * (1-self.bedporosity) / 
             (3 * eff_d * self.bedporosity**3) +
             0.5*(1+self.kl)*self.densityfluid)
        b = (150 * self.viscosityfluid * rt * (1-self.bedporosity)**2 / 
             (eff_d**2 * self.bedporosity**3))
        if self.hfluid < self.tapholediameter:
            # partially filled taphole
            pa = 0.5 * self.hfluid * self.densityfluid * g
            theta = 2*np.arccos(1-self.hfluid/rt)
            xarea = 0.5 * rt**2 * (theta - np.sin(theta))
        else:
            pa = (self.hfluid - rt) * self.densityfluid * g
            xarea = pi * rt**2
        ufluid = (-b + np.sqrt(b**2 + 4*pa*a)) / (2*a)
        vdot = xarea * ufluid
        self.vdot = vdot
        return vdot
    
    def calc_dt(self, dt):
        vdf = self.calc_vdot_out()
        if self.hfluid > self.bedheight:
            dhfluid = -dt * (vdf / self.tankarea)
        else:
            dhfluid = -dt * (vdf / (self.tankarea * self.bedporosity))
        self.hfluid += dhfluid
        return vdf

def vango_parameters():
    return {'tankarea': 0.33*0.15, 'tapholediameter': 0.031, 'kl': 0.5, 
            'densityfluid': 1000, 'viscosityfluid': 0.001, 
            'particlediameter': 0.0065, 'particlesphericity': 0.75, 
            'bedporosity': 0.52, 'bedheight': 0.25, 'hfluid_init': 0.3}
    
def vango_experiment_data():
    return {'time': [7.87, 16.2, 27.8, 41.8, 63.8],
            'tapped_kg': [2.39, 4.7, 7.06, 8.46, 8.84]}
