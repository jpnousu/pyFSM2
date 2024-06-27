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

        self.gsat = params.gsat

        self.bch = soilprops.bch
        self.hcap_soil = soilprops.hcap_soil
        self.hcon_soil = soilprops.hcon_soil
        self.sathh = soilprops.sathh
        self.Vcrit = soilprops.Vcrit
        self.Vsat = soilprops.Vsat

    def run_timestep(Nsnow, Dsnw, Sice, Sliq, Tsnow, Tsoil, Vsmc):
        '''
        Thermal conductivity of snow
        '''
        # Here could add routine to create fixed ksnow
        ksnow = np.zeros(Nsnow)
        ksnow[:] = self.kfix
        if self.CONDUCT == 1:
            for k in range(Nsnow):
                self.rhos = self.rhof
            if self.DENSITY == 1:
                if (Dsnw[k] > self.eps):
                    self.rhos = (Sice[k] + Sliq[k]) / Dsnw[k]
                ksnow[k] = 2.224 * (self.rhos / self.rho_wat)**1.885

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

        '''
        Surface layer
        '''
        Ds1 = max(self.Dzsoil[0], self.Dsnw[0])
        Ts1 = Tsoil[0] + (Tsnow[0] - Tsoil[0])*Dsnw[0]/Dzsoil[0]
        ks1 = self.Dzsoil[0]/(2*Dsnw[0]/ksnow[0] + (Dzsoil[0] - 2*Dsnw[0])/ksoil[0])
        snd = sum(Dsnw)
        if (snd > 0.5*Dzsoil[0]):
            ks1 = ksnow[0]
        if (snd > Dzsoil[0]):
            Ts1 = Tsnow[0]
        
        return Ds1, gs1, ks1, Ts1, csoil, ksnow, ksoil

