#-----------------------------------------------------------------------
# Update soil temperatures
#-----------------------------------------------------------------------

import numpy as np
from pyFSM2_MODULES import Constants, Layers, Parameters, SoilProps
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
        self.Tsoil = np.zeros(self.Nsoil)+273.15 # Soil layer temperatures (K)
        # k # Soil layer counter

        # no in nor out
        self.a =  np.zeros(self.Nsoil) # Below-diagonal matrix elements
        self.b = np.zeros(self.Nsoil) # Diagonal matrix elements
        self.c = np.zeros(self.Nsoil) # Above-diagonal matrix elements
        self.dTs = np.zeros(self.Nsoil) # Temperature increments (k)
        self.gs = np.zeros(self.Nsoil) # Thermal conductivity between layers (W/m^2/k)
        self.rhs = np.zeros(self.Nsoil) # Matrix equation rhs

        self.eps = np.finfo(float).eps

        def run_timestep(dt, Gsoil, csoil, ksoil):
            ksoil = self.soil_thermal()
            for k in Nsoil:
                self.gs[k] = 2 / self.Dzsoil[k]/ksoil[k] + self.Dzsoil[k+1]/ksoil[k+1]

            self.a[0] = 0
            self.b[0] = csoil[0] + self.gs[0]*dt
            self.c[0] = - self.gs[0]*dt
            self.rhs[0] = (Gsoil - self.gs[0]*(Tsoil[0] - Tsoil[1]))*dt

            k = self.Nsoil
            self.gs[k] = ksoil[k]/self.Dzsoil[k]
            self.a[k] = self.c[k-1]
            self.b[k] = csoil[k] + (self.gs[k-1] + self.gs[k])*dt
            self.c[k] = 0
            self.rhs[k] = self.gs[k-1]*(Tsoil[k-1] - Tsoil[k])*dt

            self.dTs = self.tridiag(Nvec=self.Nsoil, Nmax=self.Nsoil)

            for k in Nsoil:
                self.Tsoil[k] = self.Tsoil[k] + self.dTs[k]

            Tsoil = self.Tsoil

            return Tsoil
    

    def tridiag(self, Nvec, Nmax):
        '''
        Input
        Nvec: Vector length
        Nmax: Maximum vector length
        a: Below-diagonal matrix elements
        b: Diagonal matrix elements
        c: Above-diagonal matrix elements
        r: Matrix equation rhs
        
        Output
        x: Solution vector
        '''

        x = np.zeros(Nmax)
        g = np.zeros(Nmax)
        r = self.rhs
        
        beta = self.b[0]
        x[0] = r[0] / beta

        for n in range(1, Nvec):
            g[n] = self.c[n-1] / beta
            beta = self.b[n] - self.a[n] * g[n]
            x[n] = (r[n] - self.a[n] * x[n-1]) / beta

        for n in range(Nvec - 2, 0, -1):
            x[n] = x[n] - g[n + 1] * x[n + 1]

        return x

    def soil_thermal(self):
        
        '''
        Heat capacity and thermal conductivity of soil
        '''
    
        dPsidT = -self.rho_ice*self.Lf/(self.rho_wat*self.g*self.Tm)
        for k in self.Nsoil:
            self.csoil[k] = self.hcap_soil*self.Dzsoil[k]
            self.ksoil[k] = self.hcon_soil
            if (Vsmc[k] > self.eps[k]):
                dthudT = 0
                sthu = Vsmc[k]
                sthf = 0
                Tc = self.Tsoil[k] - self.Tm
                Tmax = self.Tm + (self.sathh/dPsidT)*(self.Vsat/Vsmc[k])**self.bch
                if (self.Tsoil[k] < Tmax):
                    dthudT = (-dPsidT*self.Vsat/(self.bch*self.sathh)) * (dPsidT*Tc/self.sathh)**(-1/self.bch - 1)
                    sthu = self.Vsat*(dPsidT*Tc/self.sathh)**(-1/self.bch)
                    sthu = min(sthu, Vsmc[k])
                    sthf = (Vsmc[k] - sthu)*self.rho_wat/self.rho_ice
                Mf = self.rho_ice*self.Dzsoil[k]*sthf
                Mu = self.rho_wat*self.Dzsoil[k]*sthu
                self.csoil[k] = self.hcap_soil*self.Dzsoil[k] + self.hcap_ice*Mf + hcap_wat*Mu + rho_wat*Dzsoil[k]*((hcap_wat - hcap_ice)*Tc + Lf)*dthudT
                Smf = rho_ice*sthf/(rho_wat*self.Vsat)
                Smu = sthu/self.Vsat
                thice = 0
                if (Smf > 0):
                    thice = self.Vsat*Smf/(Smu + Smf)
                thwat = 0
                if (Smu > 0):
                    thwat = self.Vsat*Smu/(Smu + Smf)
                hcon_sat = self.hcon_soil*(self.hcon_wat**thwat)*(self.hcon_ice**thice) / (self.hcon_air**self.Vsat)
                ksoil[k] = (hcon_sat - self.hcon_soil)*(Smf + Smu) + self.hcon_soil
                if (k == 1):
                    gs1 = self.gsat*max((Smu*self.Vsat/self.Vcrit)**2, 1.)




