#-----------------------------------------------------------------------
# Update soil temperatures
#-----------------------------------------------------------------------

import numpy as np
from pyFSM2_MODULES import Constants, Layers, Parameters, SoilProps
from utils import tridiag

import matplotlib.pyplot as plt

class SoilModel:
    def __init__(self):
    
        constants = Constants()
        layers = Layers()
        params = Parameters(SETPAR=2, DENSITY=0)
        soilprops = SoilProps()

        self.Dzsoil = layers.Dzsoil # Soil layer thicknesses (m)
        self.Nsoil = layers.Nsoil # Number of soil layers
        
        # in
        # dt,                # Timestep (s)
        # Gsoil,             # Heat flux into soil (W/m^2)
        # csoil(Nsoil),      # Areal heat capacity of soil layers (J/K/m^2)
        # ksoil(Nsoil) # Thermal conductivity of soil layers (W/m/K)

        # in/out
        # k # Soil layer counter

        # no in nor out
        self.a =  np.zeros(self.Nsoil) # Below-diagonal matrix elements
        self.b = np.zeros(self.Nsoil) # Diagonal matrix elements
        self.c = np.zeros(self.Nsoil) # Above-diagonal matrix elements
        self.dTs = np.zeros(self.Nsoil) # Temperature increments (k)
        self.gs = np.zeros(self.Nsoil) # Thermal conductivity between layers (W/m^2/k)
        self.rhs = np.zeros(self.Nsoil) # Matrix equation rhs

        self.eps = np.finfo(float).eps

    def run_timestep(self, dt, Gsoil, csoil, ksoil, Tsoil):

        for k in range(self.Nsoil-1):
            self.gs[k] = 2 / self.Dzsoil[k]/ksoil[k] + self.Dzsoil[k+1]/ksoil[k+1]

        self.a[0] = 0
        self.b[0] = csoil[0] + self.gs[0]*dt
        self.c[0] = - self.gs[0]*dt
        self.rhs[0] = (Gsoil - self.gs[0]*(Tsoil[0] - Tsoil[1]))*dt

        for k in range(1, self.Nsoil-1):
            self.a[k] = self.c[k-1]
            self.b[k] = csoil[k] + (self.gs[k-1] + self.gs[k])*dt
            self.c[k] = - self.gs[k]*dt
            self.rhs[k] = self.gs[k-1]*(Tsoil[k-1] - Tsoil[k])*dt + self.gs[k]*(Tsoil[k-1] - Tsoil[k])*dt 

        k = self.Nsoil-1
        self.gs[k] = ksoil[k]/self.Dzsoil[k]
        self.a[k] = self.c[k-1]
        self.b[k] = csoil[k] + (self.gs[k-1] + self.gs[k])*dt
        self.c[k] = 0
        self.rhs[k] = self.gs[k-1]*(Tsoil[k-1] - Tsoil[k])*dt

        self.dTs = tridiag(Nvec=self.Nsoil, Nmax=self.Nsoil, a=self.a, b=self.b, c=self.c, r=self.rhs)

        for k in range(self.Nsoil):
            Tsoil[k] = Tsoil[k] + self.dTs[k]

        return Tsoil





