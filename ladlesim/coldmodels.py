# Submodule for simulation of tank drainage in various cold-model simulation
# experiments.
#
# Q Reynolds 2019

import numpy as np

g = 9.81

class VangoExperiment():
    def __init__(self, tankarea, tapholearea, kl, densityfluid, viscosityfluid, 
                 particlediameter, bedporosity, bedheight, hfluid_init):
        self.tankarea = tankarea
        self.tapholearea = tapholearea
        self.kl = kl
        self.densityfluid = densityfluid
        self.viscosityfluid = viscosityfluid
        self.particlediameter = particlediameter
        self.bedporosity = bedporosity
        self.bedheight = bedheight
        self.hfluid = hfluid_init
    
    def calc_vdot_out(self):
        """Return outlet flowrates as a function of the current tank state. 
        Simplified model without interface deformation near tap-hole entry.
        """
        rt = np.sqrt(self.tapholearea / np.pi)
        a = (150 * self.viscosityfluid * rt * (1-self.bedporosity)**2 / 
               (self.particlediameter**2 * self.bedporosity**3) +
               0.5*(1+self.kl)*self.densitymetal)
        b = (1.75 * self.densityfluid * rt * (1-self.bedporosity) / 
               (3 * self.particlediameter * self.bedporosity**3))
        if self.hfluid < 2*rt:
            # partially filled taphole
            pa = 0.5 * (self.hfluid + 2*rt) * self.densityfluid * g
            theta = 2*np.arccos(1 - (self.hslag + rt)/rt)
            xarea = 0.5*rt**2 * (theta - np.sin(theta)) 
        else:
            pa = (self.hfluid - rt) * self.densityfluid * g
            xarea = 2*np.pi*rt**2
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
