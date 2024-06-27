#-----------------------------------------------------------------------
# Thermal properties of snow and soil
#-----------------------------------------------------------------------

import numpy as np
from pyFSM2_MODULES import Constants, Layers, Parameters, SoilProps
import matplotlib.pyplot as plt

class Thermal():
    def __init__(self):
        
        constants = Constants()
        layers = Layers()
        params = Parameters(SETPAR=2, DENSITY=0)
        soilprops = SoilProps()
        
        self.kfix = params.kfix
        self.fcly = params.fcly
        self.fsnd = params.fsnd
        self.rhof = params.rhof # Fresh snow density (kg/m^3)

        self.Nsoil = layers.Nsoil
        self.Nsmax = layers.Nsmax
        self.Dzsoil = layers.Dzsoil

        self.rho_ice = constants.rho_ice
        self.rho_wat = constants.rho_wat
        self.Lf = constants.Lf
        self.g = constants.g
        self.Tm = constants.Tm
        self.hcap_ice = constants.hcap_ice
        self.hcap_wat = constants.hcap_wat
        self.hcon_air = constants.hcon_air
        self.hcon_ice = constants.hcon_ice
        self.hcon_wat = constants.hcon_wat
        self.hcon_sand = constants.hcon_sand
        self.hcon_clay = constants.hcon_clay

        self.gsat = params.gsat

        self.bch = soilprops.bch
        self.hcap_soil = soilprops.hcap_soil
        self.hcon_soil = soilprops.hcon_soil
        self.sathh = soilprops.sathh
        self.Vcrit = soilprops.Vcrit
        self.Vsat = soilprops.Vsat

        # Soil properties
        self.bch = 3.1 + 15.7*self.fcly - 0.3*self.fsnd
        self.hcap_soil = (2.128*self.fcly + 2.385*self.fsnd)*1e6 / (self.fcly + self.fsnd)
        self.sathh = 10**(0.17 - 0.63*self.fcly - 1.58*self.fsnd)
        self.Vsat = 0.505 - 0.037*self.fcly - 0.142*self.fsnd
        self.Vcrit = self.Vsat*(self.sathh/3.364)**(1/self.bch)
        self.hcon_soil = (self.hcon_air**self.Vsat) * ((self.hcon_clay**self.fcly)*(self.hcon_sand**(1 - self.fcly))**(1 - self.Vsat))

        self.HYDROL = 1 # NOTE THIS NEEDS TO COME FROM THE NAMELIST!
        self.CONDUCT = 1
        self.DENSITY = 1

        self.eps = np.finfo(float).eps


    def run_timestep(self, Nsnow, Dsnw, Sice, Sliq, Tsnow, Tsoil, Vsmc):
        '''
        Thermal conductivity of snow
        '''
        # Here could add routine to create fixed ksnow
        ksoil = np.zeros(int(self.Nsoil))
        ksnow = np.zeros(int(self.Nsmax))
        csoil = np.zeros(int(self.Nsoil))

        ksnow[:] = self.kfix
        if self.CONDUCT == 1:
            for k in range(int(Nsnow)):
                self.rhos = self.rhof
                if self.DENSITY == 1:
                    if (Dsnw[k] > self.eps):
                        self.rhos = (Sice[k] + Sliq[k]) / Dsnw[k]
                ksnow[k] = 2.224 * (self.rhos / self.rho_wat)**1.885

        '''
        Heat capacity and thermal conductivity of soil
        '''
    
        dPsidT = -self.rho_ice*self.Lf/(self.rho_wat*self.g*self.Tm)
        for k in range(self.Nsoil):
            csoil[k] = self.hcap_soil*self.Dzsoil[k]
            ksoil[k] = self.hcon_soil
            if (Vsmc[k] > self.eps):
                dthudT = 0
                sthu = Vsmc[k]
                sthf = 0
                Tc = Tsoil[k] - self.Tm
                Tmax = self.Tm + (self.sathh/dPsidT)*(self.Vsat/Vsmc[k])**self.bch
                if (Tsoil[k] < Tmax):
                    dthudT = (-dPsidT*self.Vsat/(self.bch*self.sathh)) * (dPsidT*Tc/self.sathh)**(-1/self.bch - 1)
                    sthu = self.Vsat*(dPsidT*Tc/self.sathh)**(-1/self.bch)
                    sthu = min(sthu, Vsmc[k])
                    sthf = (Vsmc[k] - sthu)*self.rho_wat/self.rho_ice
                Mf = self.rho_ice*self.Dzsoil[k]*sthf
                Mu = self.rho_wat*self.Dzsoil[k]*sthu
                csoil[k] = self.hcap_soil*self.Dzsoil[k] + self.hcap_ice*Mf + self.hcap_wat*Mu + self.rho_wat*self.Dzsoil[k]*((self.hcap_wat - self.hcap_ice)*Tc + self.Lf)*dthudT
                Smf = self.rho_ice*sthf/(self.rho_wat*self.Vsat)
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

        '''
        Surface layer
        '''
        Ds1 = max(self.Dzsoil[0], Dsnw[0])
        Ts1 = Tsoil[0] + (Tsnow[0] - Tsoil[0])*Dsnw[0]/self.Dzsoil[0]
        print('ksnow', ksnow)
        ks1 = self.Dzsoil[0]/(2*Dsnw[0]/ksnow[0] + (self.Dzsoil[0] - 2*Dsnw[0])/ksoil[0])
        snd = sum(Dsnw)
        if (snd > 0.5*self.Dzsoil[0]):
            ks1 = ksnow[0]
        if (snd > self.Dzsoil[0]):
            Ts1 = Tsnow[0]
        
        return Ds1, gs1, ks1, Ts1, csoil, ksnow, ksoil

