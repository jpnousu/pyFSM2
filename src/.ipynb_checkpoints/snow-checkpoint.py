import numpy as np
from pyFSM2_MODULES import Constants, Layers, Parameters
import matplotlib.pyplot as plt

class SnowModel:
    def __init__(self):

        constants = Constants()
        layers = Layers()
        params = Parameters(SETPAR=2, DENSITY=0)
        
        # from Constants
        self.g = constants.g # Acceleration due to gravity (m/s^2)
        self.hcap_ice = constants.hcap_ice # Specific heat capacity of ice (J/K/kg)
        self.hcap_wat = constants.hcap_wat # Specific heat capacity of water (J/K/kg)
        self.Lf = constants.Lf # Latent heat of fusion (J/kg)
        self.Ls = constants.Ls # Latent heat of sublimation (J/kg)
        self.mu_wat = constants.mu_wat # Dynamic viscosity of water (kg/m/s)
        self.rho_ice = constants.rho_ice # Density of ice (kg/m^3)
        self.rho_wat = constants.rho_wat # Density of water (kg/m^3)
        self.Tm = constants.Tm # Melting point (K)

        # from Layers
        self.Dzsnow = layers.Dzsnow # Minimum snow layer thicknesses (m)
        self.Dzsoil = layers.Dzsoil # Soil layer thicknesses (m)
        self.Nsmax = layers.Nsmax # Maximum number of snow layers
        self.Nsoil = layers.Nsoil # Number of soil layers

        # from Parameters
        self.eta0 = params.eta0 # Reference snow viscosity (Pa s)
        self.rcld = params.rcld # Maximum density for cold snow (kg/m^3)
        self.rfix = params.rfix # Fixed snow density (kg/m^3)
        self.rgr0 = params.rgr0 # Fresh snow grain radius (m)
        self.rhof = params.rhof # Fresh snow density (kg/m^3)
        self.rhow = params.rhow # Wind-packed snow density (kg/m^3)
        self.rmlt = params.rmlt # Maximum density for melting snow (kg/m^3)
        self.snda = params.snda # Thermal metamorphism parameter (1/s)
        self.trho = params.trho # Snow compaction timescale (s)
        self.Wirr = params.Wirr # Irreducible liquid water content of snow
        self.kfix = params.kfix # 

        self.HYDROL = 1 # NOTE THIS NEEDS TO COME FROM THE NAMELIST!
        self.CONDUCT = 1
        self.DENSITY = 1

        # dt should come in
        # dt # timestep (s)
        
        self.eps = np.finfo(float).eps
        
        # Model state variables (in/out)
        self.Nsnow = np.zeros(1) # Number of snow layers
        self.Dsnw = np.zeros(self.Nsmax) # Snow layer thicknesses (m)
        self.Rgrn = np.zeros(self.Nsmax) # Snow layer grain radius (m)
        self.Sice = np.zeros(self.Nsmax) # Ice content of snow layers (kg/m^2)
        self.Sliq = np.zeros(self.Nsmax) # Liquid content of snow layers (kg/m^2)
        self.Tsnow = np.ones(self.Nsmax)*273.15 # Snow layer temperatures (K)
        
        # Soil temp should come from elsewhere!
        # Tsoil = np.zeros(self.Nsmax) # Soil layer temperatures (K)

        # Only out!
        # Gsoil = np.zeros(1) # Heat flux into soil (W/m^2)
        # Roff = np.zeros(1) # Runoff from snow (kg/m^2/s)
        # hs = np.zeros(1) # Snow depth (m)
        # swe = np.zeros(1) # Total snow mass on ground (kg/m^2)
        self.Wflx = np.zeros(self.Nsmax) # Water flux into snow layer (kg/m^2/s)

        # no in nor out
        # i,j = np.zeros(1),np.zeros(1) # Hydrology iteration counters
        # k = np.zeros(1) # Snow layer counter
        # knew =  np.zeros(1) # New snow layer pointer
        # kold = np.zeros(1) # Old snow layer pointer
        # Nold = np.zeros(1) # Previous number of snow layers
        # coldcont = np.zeros(1) # Layer cold content (J/m^2)
        # dnew = np.zeros(1) # New snow layer thickness (m)
        # dSice = np.zeros(1) # Change in layer ice content (kg/m^2)
        # Esnow = np.zeros(1) # Snow sublimation rate (kg/m^2/s)
        # ggr = np.zeros(1) # Grain area growth rate (m^2/s)
        # mass = np.zeros(1) # Mass of overlying snow (kg/m^2)
        # rhos = np.zeros(1) # Density of snow layer (kg/m^3)
        # SliqMax = np.zeros(1) # Maximum liquid content for layer (kg/m^2)
        # wt = np.zeros(1) # Layer weighting
        self.a =  np.zeros(self.Nsmax) # Below-diagonal matrix elements
        self.b = np.zeros(self.Nsmax) # Diagonal matrix elements
        self.c = np.zeros(self.Nsmax) # Above-diagonal matrix elements

        # Below are variables that need to have shape rightaway!
        self.csnow = np.zeros(self.Nsmax) # Areal heat capacity of snow (J/K/m^2)
        self.dTs = np.zeros(self.Nsmax) # Temperature increments (k)
        self.D = np.zeros(self.Nsmax) # Layer thickness before adjustment (m)
        self.E = np.zeros(self.Nsmax) # Energy contents before adjustment (J/m^2)
        self.Gs = np.zeros(self.Nsmax) # Thermal conductivity between layers (W/m^2/k)
        self.phi = np.zeros(self.Nsmax) # Porosity of snow layers
        self.rhs = np.zeros(self.Nsmax) # Matrix equation rhs
        #self.R = np.zeros(self.Nsmax) # Snow grain radii before adjustment (kg/m^2)
        #self.S = np.zeros(self.Nsmax) # Ice contents before adjustment (kg/m^2)
        self.U = np.zeros(self.Nsmax) # Layer internal energy contents (J/m^2)
        #self.W = np.zeros(self.Nsmax) # Liquid contents before adjustment (kg/m^2)   
        # dth = np.zeros(1) # Hydrology timestep (s)

        # Below are variables that need to have shape rightaway!
        self.dtheta = np.zeros(self.Nsmax) # Change in liquid water content
        self.ksat = np.zeros(self.Nsmax) # Saturated hydraulic conductivity (m/s)
        self.thetar = np.zeros(self.Nsmax) # Irreducible water content
        self.thetaw = np.zeros(self.Nsmax) # Volumetric liquid water content
        self.theta0 = np.zeros(self.Nsmax) # Liquid water content at start of timestep
        self.Qw = np.zeros(self.Nsmax+1) # Water flux at snow layer boundaruess (m/s)

        # temp
        self.hs_list = []
    
    def run_timestep(self, dt, drip, Esrf, Gsrf, ksoil, Melt, Rf, Sf, Ta, trans, Tsrf, unload, Tsoil):
        '''
        '''
        print('initial self.Nsnow', self.Nsnow)
        print('initial self.Dsnw', self.Dsnw)
        print('initial self.Sice', self.Sice)
        print('initial self.Sliq', self.Sliq)


        # No snow
        Gsoil = Gsrf.copy()
        Roff = Rf + drip / dt
        self.Wflx[:] = 0

        # Existing snowpack
        if (self.Nsnow > 0):
            ksnow = self.snow_thermal()
            # Heat conduction
            for k in range(self.Nsnow):
                # Areal heat capacity
                self.csnow[k] = self.Sice[k]*self.hcap_ice + self.Sliq[k]*self.hcap_wat
            if (self.Nsnow == 1):
                self.Gs[0] = 2 / (self.Dsnw[0]/ksnow[0] + self.Dzsoil[0]/ksoil[0])
                self.dTs[0] = (Gsrf + self.Gs[0]*(Tsoil[0] - self.Tsnow[0]))* dt / (self.csnow[0] + self.Gs[0] * dt)
                
            else:
                for k in range(self.Nsnow-2):
                      self.Gs[k] = 2 / (self.Dsnw[k]/ksnow[k] + self.Dsnw[k+1]/ksnow[k+1])
                self.a[0] = 0
                self.b[0] = self.csnow[0] + self.Gs[0] * dt
                self.c[0] = - self.Gs[0] * dt

                self.rhs[0] = (Gsrf - self.Gs[0]*(self.Tsnow[0] - self.Tsnow[1]))*dt
                for k in range(1,self.Nsnow-1):
                    self.a[k] = self.c[k-1]
                    self.b[k] = self.csnow[k] + (self.Gs[k-1] + self.Gs[k])*dt
                    self.c[k] = - self.Gs[k]*dt
                    self.rhs[k] = self.Gs[k-1]*(self.Tsnow[k-1] - self.Tsnow[k])*dt + self.Gs[k]*(self.Tsnow[k+1] - self.Tsnow[k])*dt 
            k = self.Nsnow - 1
            print('initial ksnow', ksnow)
            print('initial ksoil', ksoil)
            self.Gs[k] = 2 / (self.Dsnw[k]/ksnow[k] + self.Dzsoil[0]/ksoil[0])
            self.a[k] = self.c[k-1]
            self.b[k] = self.csnow[k] + (self.Gs[k-1] + self.Gs[k])*dt
            self.c[k] = 0
            self.rhs[k] = self.Gs[k-1] * (self.Tsnow[k-1] - self.Tsnow[k])*dt + self.Gs[k]*(Tsoil[0] - self.Tsnow[k]) * dt
            self.dTs = self.tridiag(Nvec=self.Nsnow, Nmax=self.Nsmax)
            
            for k in range(self.Nsnow):
                self.Tsnow[k] = self.Tsnow[k] + self.dTs[k]
            k = self.Nsnow - 1
            Gsoil = self.Gs[k] * (self.Tsnow[k] - Tsoil[0])
            
            # Convert melting ice to liquid water
            dSice = Melt * dt
            for k in range(self.Nsnow):
                coldcont = self.csnow[k]*(self.Tm - self.Tsnow[k])
                if (coldcont < 0):
                    dSice = dSice - coldcont / self.Lf
                    self.Tsnow[k] = self.Tm
                if (dSice > 0):
                    if (dSice > self.Sice[k]):  # Layer melts completely
                        dSice = dSice - self.Sice[k]
                        self.Dsnw[k] = 0
                        self.Sliq[k] = self.Sliq[k] + self.Sice[k]
                        self.Sice[k] = 0
                    else:                       # Layer melts partially
                        self.Dsnw[k] = (1 - dSice/self.Sice[k])*self.Dsnw[k]
                        self.Sice[k] = self.Sice[k] - dSice
                        self.Sliq[k] = self.Sliq[k] + dSice
                        dSice = 0
                    
            # Remove snow by sublimation 
            dSice = Esrf * dt
            if (dSice > 0):
                for k in range(self.Nsnow):
                    if (dSice > self.Sice[k]):  # Layer sublimates completely
                        dSice = dSice - self.Sice[k]
                        self.Dsnw[k] = 0
                        self.Sice[k] = 0
                    else:                       # Layer sublimates partially
                        self.Dsnw[k] = (1 - dSice/self.Sice[k])*self.Dsnw[k]
                        self.Sice[k] = self.Sice[k] - dSice
                        dSice = 0

            # Remove wind-trasported snow 
            dSice = trans * dt
            if (dSice > 0):
                for k in range(self.Nsnow):
                    if (dSice > self.Sice[k]):  # Layer completely removed
                        dSice = dSice - self.Sice[k]
                        self.Dsnw[k] = 0
                        self.Sice[k] = 0
                    else:                       # Layer partially removed
                        self.Dsnw[k] = (1 - dSice/self.Sice[k])*self.Dsnw[k]
                        self.Sice[k] = self.Sice[k] - dSice
                        dSice = 0

            if self.DENSITY == 0:
                # Fixed snow density
                for k in range(self.Nsnow):
                    if (self.Dsnw[k] > self.eps):
                        self.Dsnw[k] = (self.Sice[k] + self.Sliq[k]) / rfix
            if self.DENSITY == 1:
                # Snow compaction with age
                for k in range(self.Nsnow):
                    if self.Dsnw[k] > self.eps: # epsillon different in FSM
                        self.rhos = (self.Sice[k] + self.Sliq[k]) / self.Dsnw[k]
                        if self.Tsnow[k] >= self.Tm:
                            if self.rhos < self.rmlt:
                                self.rhos = self.rmlt + (self.rhos - self.rmlt) * np.exp(-dt / self.trho)
                        else:
                            if self.rhos < self.rcld:
                                self.rhos = self.rcld + (self.rhos - self.rcld) * np.exp(-dt / self.trho)
                        self.Dsnw[k] = (self.Sice[k] + self.Sliq[k]) / self.rhos
                        
            if self.DENSITY == 2:
                # Snow compaction by overburden
                mass = 0
                for k in range(self.Nsnow):
                    mass = mass + 0.5*(self.Sice[k] + self.Sliq[k]) 
                    if (self.Dsnw[k] > np.finfo(float).eps):
                        self.rhos = (self.Sice[k] + self.Sliq[k]) / self.Dsnw[k]
                        self.rhos = self.rhos + (self.rhos*g*mass*dt/(eta0*np.exp(-(self.Tsnow[k] - self.Tm)/12.4 + self.rhos/55.6)) + dt * self.rhos*snda*np.exp((self.Tsnow[k] - self.Tm)/23.8 - max(self.rhos - 150, 0.)/21.7))
                        self.Dsnw[k] = (self.Sice[k] + self.Sliq[k]) / self.rhos
                    mass = mass + 0.5*(self.Sice[k] + self.Sliq[k])

            # Snow grain growth
            for k in range(self.Nsnow):
                ggr = 2e-13
                if (self.Tsnow[k] < self.Tm):
                    if (self.Rgrn[k] < 1.50e-4):
                        ggr = 2e-14
                    else:
                        ggr = 7.3e-8*np.exp(-4600/self.Tsnow[k])
                self.Rgrn[k] = self.Rgrn[k] + dt * ggr / self.Rgrn[k]

        # End if for existing snowpack

        # Add snowfall and frost to layer 1 with fresh snow density and grain size
        Esnow = 0
        if (Esrf < 0) & (Tsrf < self.Tm):
            Esnow = Esrf
        dSice = (Sf - Esnow)*dt
        self.Dsnw[0] = self.Dsnw[0] + dSice / self.rhof
        if (self.Sice[0] + dSice > self.eps):
            self.Rgrn[0] = (self.Sice[0]*self.Rgrn[0] + dSice*self.rgr0) / (self.Sice[0] + dSice)
        self.Sice[0] = self.Sice[0] + dSice
    
        # Add canopy unloading to layer 1 with bulk snow density and grain size
        self.rhos = self.rhof
        swe = sum(self.Sice[:]) + sum(self.Sliq[:])
        hs = sum(self.Dsnw[:])
        if (hs > self.eps):
            self.rhos = swe / hs
        self.Dsnw[0] = self.Dsnw[0] + unload / self.rhos
        if (self.Sice[0] + unload > self.eps):
            self.Rgrn[0] = (self.Sice[0]*self.Rgrn[0] + unload*self.rgr0) / (self.Sice[0] + unload)
        self.Sice[0] = self.Sice[0] + unload

        # Add wind-blown snow to layer 1 with wind-packed density and fresh grain size
        dSice = - trans*dt
        if (dSice > 0):
            self.Dsnw[0] = self.Dsnw[0] + dSice / rhow
            self.Rgrn[0] = (self.Sice[0]*self.Rgrn[0] + dSice*self.rgr0) / (self.Sice[0] + dSice)
            self.Sice[0] = self.Sice[0] + dSice

        # New snowpack
        if (self.Nsnow == 0) & (self.Sice[0] > 0):
            self.Nsnow = 1
            self.Rgrn[0] = self.rgr0
            self.Tsnow[0] = min(Ta, self.Tm)

        # Store state of old layers
        D = self.Dsnw[:].copy()
        R = self.Rgrn[:].copy()
        S = self.Sice[:].copy()
        W = self.Sliq[:].copy()
        if self.Nsnow > 0:
            for k in range(self.Nsnow):
                self.csnow[k] = self.Sice[k]*self.hcap_ice + self.Sliq[k]*self.hcap_wat
                self.E[k] = self.csnow[k]*(self.Tsnow[k] - self.Tm)
        Nold = self.Nsnow
        hs = sum(self.Dsnw[:])

        # Initialise new layers
        self.Dsnw[:] = 0
        self.Rgrn[:] = 0
        self.Sice[:] = 0
        self.Sliq[:] = 0
        self.Tsnow[:] = self.Tm
        self.U[:] = 0
        self.Nsnow = 0

        if (hs > 0):  # Existing or new snowpack
            # Re-assign and count snow layers
            dnew = hs
            self.Dsnw[0] = dnew
            if (self.Dsnw[0] > self.Dzsnow[0]):
                for k in range(self.Nsmax):
                    self.Dsnw[k] = self.Dzsnow[k]
                    dnew = dnew - self.Dzsnow[k]
                    if (dnew <= self.Dzsnow[k]) | (k == self.Nsmax):
                        self.Dsnw[k] = self.Dsnw[k] + dnew
                        break
            self.Nsnow = k + 1

            # Fill new layers from the top downwards
            knew = 0
            dnew = self.Dsnw[0]
            for kold in range(Nold):
                while True:
                    if (D[kold] < dnew):
                        print('kold', kold)
                        print('knew', knew)
                        # All snow from old layer partially fills new layer
                        self.Rgrn[knew] = self.Rgrn[knew] + S[kold] * R[kold]
                        self.Sice[knew] = self.Sice[knew] + S[kold]
                        self.Sliq[knew] = self.Sliq[knew] + W[kold]
                        self.U[knew] = self.U[knew] + self.E[kold]
                        dnew = dnew - D[kold]
                        break
                    else:
                        # Some snow from old layer fills new layer
                        wt = dnew / D[kold]
                        self.Rgrn[knew] = self.Rgrn[knew] + wt * S[kold] * R[kold]
                        self.Sice[knew] = self.Sice[knew] + wt * S[kold]
                        self.Sliq[knew] = self.Sliq[knew] + wt * W[kold]
                        self.U[knew] = self.U[knew] + wt * self.E[kold]
                        D[kold] = (1 - wt) * D[kold]
                        self.E[kold] = (1 - wt) * self.E[kold]
                        S[kold] = (1 - wt) * S[kold]
                        W[kold] = (1 - wt) * W[kold]
                        knew = knew + 1
                        if (knew > self.Nsnow-1):
                            break
                        dnew = self.Dsnw[knew]

        # Diagnose snow layer temperatures            
        for k in range(self.Nsnow):
            self.csnow[k] = self.Sice[k]*self.hcap_ice + self.Sliq[k]*self.hcap_wat
            self.Tsnow[k] = self.Tm + self.U[k] / self.csnow[k]
            self.Rgrn[k] = self.Rgrn[k] / self.Sice[k]

        # Drain, retain or freeze snow in layers
        if self.HYDROL == 0:
            # Free-draining snow, no retention or freezing 
            self.Wflx[0] = Roff
            for k in range(self.Nsnow):
                Roff = Roff + self.Sliq[k] / dt
                self.Sliq[k] = 0
                if (k < self.Nsnow):
                    self.Wflx[k+1] = Roff

        if self.HYDROL == 1:
            # Bucket storage 
            if (np.max(self.Sliq)) > 0 | (Rf > 0):
                for k in range(self.Nsnow):
                    self.phi[k] = 1 - self.Sice[k]/(self.rho_ice * self.Dsnw[k])
                    SliqMax = self.rho_wat * self.Dsnw[k] * self.phi[k] * self.Wirr
                    self.Sliq[k] = self.Sliq[k] + Roff * dt
                    self.Wflx[k] = Roff
                    Roff = 0
                if (self.Sliq[k] > SliqMax):       # Liquid capacity exceeded
                    Roff = (self.Sliq[k] - SliqMax)/dt   # so drainage to next layer
                    self.Sliq[k] = SliqMax
                self.csnow[k] = self.Sice[k]*self.hcap_ice + self.Sliq[k]*self.hcap_wat
                coldcont = self.csnow[k]*(self.Tm - self.Tsnow[k])
                if (coldcont > 0):            # Liquid can freeze
                    dSice = min(self.Sliq[k], coldcont / self.Lf)
                    self.Sliq[k] = self.Sliq[k] - dSice
                    self.Sice[k] = self.Sice[k] + dSice
                    self.Tsnow[k] = self.Tsnow[k] + self.Lf*dSice/self.csnow[k]

        '''
        if self.HYDROL == 2: # NOTE THIS NEEDS TESTING!
            # Gravitational drainage 
            if (np.max(self.Sliq) > 0 | Rf > 0):
                self.Qw[:] = 0
                self.Qw[0] = Rf/rho_wat
                Roff = 0
                for k in range(self.Nsnow):
                    self.ksat[k] = 0.31*(self.rho_wat*self.g/self.mu_wat) * self.Rgrn[k]**2 * np.exp(-7.8 * self.Sice[k]/(rho_wat * self.Dsnw[k]))
                    self.phi[k] = 1 - self.Sice[k]/(self.rho_ice*self.Dsnw[k])
                    self.thetar[k] = self.Wirr*self.phi[k]
                    self.thetaw[k] = self.Sliq[k]/(rho_wat*self.Dsnw[k])
                    if (self.thetaw[k]>self.phi[k]):
                        Roff = Roff + rho_wat * self.Dsnw[k]*(self.thetaw[k] - self.phi[k])/dt
                        self.thetaw[k] = self.phi[k]
                dth = 0.1*dt
                for i in range(10): # subdivide timestep NOTE CHECK THIS LATER!
                    self.theta0[:] = self.thetaw[:]
                    for j in range(10): # Newton-Raphson iteration
                        a[:] = 0
                        b[:] = 1/dth 
                        if (self.thetaw[0] > self.thetar[0]):
                            b[0] = 1/dth + 3*self.ksat[0]*(self.thetaw[0] - self.thetar[0])**2/(self.phi[0] - self.thetar[0])**3 / self.Dsnw[0]
                            self.Qw[1] = self.ksat[0]*((self.thetaw[0] - self.thetar[0])/(self.phi[0] - self.thetar[0]))**3
                        self.rhs[0] = (self.thetaw[0] - self.theta0[0])/dth + (self.Qw[1] - self.Qw[0])/self.Dsnw[0]
                        for k in range(1, self.Nsnow):
                            if (self.thetaw[k-1] > self.thetar[k-1]):
                                a[k] = - 3*self.ksat[k-1]*(self.thetaw[k-1] - self.thetar[k-1])**2/(self.phi[k-1] - self.thetar[k-1])**3 / self.Dsnw[k-1]
                            if (self.thetaw[k] > self.thetar[k]):
                                b[k] = 1/dth + 3*self.ksat[k]*(self.thetaw[k] - self.thetar[k])**2/(self.phi[k] - self.thetar[k])**3 / self.Dsnw[k]
                                self.Qw[k+1] = self.ksat[k]*((self.thetaw[k] - self.thetar[k])/(self.phi[k] - self.thetar[k]))**3
                            self.rhs[k] = (self.thetaw[k] - self.theta0[k])/dth + (self.Qw[k+1] - self.Qw[k]) / self.Dsnw[k]
                        dtheta[0] = - self.rhs[0]/b[0]
                        for k in range(1, self.Nsnow):
                            dtheta[k] = - (a[k]*dtheta[k-1] + self.rhs[k])/b[k]
                        for k in range(self.Nsnow):
                            self.thetaw[k] = self.thetaw[k] + dtheta[k]
                            if (self.thetaw[k] > self.phi[k]):
                                self.Qw[k+1] = self.Qw[k+1] + (self.thetaw[k] - self.phi[k]) * self.Dsnw[k]/dth
                                self.thetaw[k] = self.phi[k]
                    self.Wflx[:] = self.Wflx[:] + rho_wat*self.Qw[0:self.Nsmax]/10
                    Roff = Roff + rho_wat*self.Qw[self.Nsnow+1]/10
                self.Sliq[:] = rho_wat * self.Dsnw[:]*self.thetaw[:]
                for k in range(self.Nsnow):
                    self.csnow[k] = self.Sice[k]*self.hcap_ice + self.Sliq[k]*self.hcap_wat
                    coldcont = self.csnow[k]*(Tm - Tsnow[k])
                    if (coldcont > 0): # Liquid can freeze
                        dSice = min(self.Sliq[k], coldcont/Lf)
                        self.Sliq[k] = self.Sliq[k] - dSice
                        self.Sice[k] = self.Sice[k] + dSice
                        Tsnow[k] = Tsnow[k] + self.Lf*dSice/self.csnow[k]
        '''

        swe = sum(self.Sice[:]) + sum(self.Sliq[:])
        Wflx = self.Wflx # NOTE IS THIS REASONABLE?
        Sice = self.Sice
        # End if existing or new snowpack

        # NOTE SHOULD SAVE RESULTS BETTER!
        
        self.hs_list.append(hs)
        #x_temp = np.arange(0, len(self.hs_list), 1)
        #plt.plot(x_temp, self.hs_list)
        #plt.show(x_temp, self.hs_list)
        
        return Gsoil, Roff, hs, swe, Wflx, Sice, Sliq, Dsnw, Rgrn, Tsnow, Tsoil, Nsnow
    
    
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

    def snow_thermal(self):
        '''
        Thermal conductivity of snow
        '''
        # Here could add routine to create fixed ksnow
        ksnow = np.zeros(self.Nsnow)
        ksnow[:] = self.kfix
        if self.CONDUCT == 1:
            for k in range(self.Nsnow):
                self.rhos = self.rhof
            if self.DENSITY == 1:
                if (self.Dsnw[k] > self.eps):
                    self.rhos = (self.Sice[k] + self.Sliq[k]) / self.Dsnw[k]
                ksnow[k] = 2.224 * (self.rhos / self.rho_wat)**1.885
                
        return ksnow

        