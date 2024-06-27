import numpy as np
from pyFSM2_MODULES import Constants, Layers, Parameters

class SebModel:
    def __init__(self):

        constants = Constants()
        layers = Layers()
        params = Parameters(SETPAR=2, DENSITY=0)

        # From Constants
        self.cp = constants.cp # Specific heat capacity of air (J/K/kg)
        self.g = constants.g # Acceleration due to gravity (m/s^2)
        self.Lf = constants.Lf # Latent heat of fusion (J/kg)
        self.Ls = constants.Ls # Latent heat of sublimation (J/kg)
        self.Lv = constants.Lv # Latent heat of vapourisation (J/kg)
        self.Rair = constants.Rair # Gas constant for air (J/K/kg)
        self.Rwat = constants.Rwat # Gas constant for water vapour (J/K/kg)
        self.sb = constants.sb # Stefan-Boltzmann constant (W/m^2/K^4)
        self.Tm = constants.Tm # Melting point (K)
        self.vkman = constants.vkman # Von Karman constant
        self.e0 = constants.e0 # # Saturation vapour pressure at Tm (Pa)
                
        # From Layers
        self.Ncnpy = layers.Ncnpy # Number of canopy layers
        self.Nsmax = layers.Nsmax # Maximum number of snow layers
        self.fvg1 = layers.fvg1 # Fraction of vegetation in upper canopy layer
        self.zsub = layers.zsub # Subcanopy wind speed diagnostic height (m)

        # From Parameters
        self.gsnf = params.gsnf # Snow-free vegetation moisture conductance (m/s)
        self.hbas = params.hbas # Canopy base height (m)
        self.kext = params.kext # Vegetation light extinction coefficient
        self.leaf = params.leaf # Leaf boundary resistance (s/m)^(1/2)
        self.wcan = params.wcan # Canopy wind decay coefficient
        self.z0sf = params.z0sf # Snow-free surface roughness length (m)
        self.z0sn = params.z0sn # Snow roughness length (m)

        # Coming IN!
        #self.Ds1 = np.zeros(1) # Surface layer thickness (m)
        #self.dt = np.zeros(1) # Timestep (s)
        self.fsnow = np.zeros(1) # Ground snowcover fraction
        #self.gs1 = np.zeros(1) # Surface moisture conductance (m/s)
        #self.ks1 = np.zeros(1) # Surface layer thermal conductivity (W/m/K)
        #self.LW = np.zeros(1) # Incoming longwave radiation (W/m2)
        #self.Ps = np.zeros(1) # Surface pressure (Pa)
        #self.Qa = np.zeros(1) # Specific humidity (kg/kg)
        #self.SWsrf = np.zeros(1) # SW absorbed by snow/ground surface (W/m^2)
        #self.Ta = np.zeros(1) # Air temperature (K)
        #self.Ts1 = np.zeros(1) # Surface layer temperature (K)
        #self.Ua = np.zeros(1) # Wind speed (m/s)

        # These also coming IN but probably should be taken from parameter file etc. and read in initizaliation!
        self.VAI = np.zeros(1) # Vegetation area index
        self.vegh = np.zeros(1) # Canopy height (m)
        self.zT = np.zeros(1)+2 # Temperature and humidity measurement height (m)
        self.zU = np.zeros(1)+5 # Wind speed measurement height (m)

        # Canopy related (some coming only IN)
        #self.cveg = np.zeros(self.Ncnpy) # Vegetation heat capacities (J/K/m^2)
        #self.fcans = np.zeros(self.Ncnpy) # Canopy layer snowcover fractions
        #self.lveg =  = np.zeros(self.Ncnpy) # Canopy layer vegetation area indices
        #self.Sveg =  = np.zeros(self.Ncnpy) # Snow mass on vegetation layers (kg/m^2)
        #self.SWveg =  = np.zeros(self.Ncnpy) # SW absorbed by vegetation layers (W/m^2)
        #self.tdif =  = np.zeros(self.Ncnpy) # Canopy layer diffuse transmittances
        #self.Tveg0 =  = np.zeros(self.Ncnpy) # Vegetation temperatures at start of timestep (K)

        # next goes in and out
        #self.Tsrf = np.zeros(1) # Snow/ground surface temperature (K)
        #self.Qcan = np.zeros(Ncnpy) # Canopy air space humidities
        #self.Sice = np.zeros(Nsmax) # Ice content of snow layers (kg/m^2)
        #self.Tcan = np.zeros(Ncnpy) # Canopy air space temperatures (K)
        #self.Tveg = np.zeros(Ncnpy) # Vegetation layer temperatures (K)

        # These should also be returning from the run_timestep!
        #self.Esrf = np.zeros(1) # Moisture flux from the surface (kg/m^2/s)
        #self.Gsrf = np.zeros(1) # Heat flux into snow/ground surface (W/m^2)
        #self.H = np.zeros(1) # Sensible heat flux to the atmosphere (W/m^2)
        #self.LE = np.zeros(1) # Latent heat flux to the atmosphere (W/m^2)
        #self.LWout = np.zeros(1) # Outgoing LW radiation (W/m^2)
        #self.LWsub = np.zeros(1) # Subcanopy downward LW radiation (W/m^2)
        #self.Melt = np.zeros(1) # Surface melt rate (kg/m^2/s)
        #self.subl = np.zeros(1) # Sublimation rate (kg/m^2/s)
        #self.Usub = np.zeros(1) # Subcanopy wind speed (m/s)
        self.Eveg = np.zeros(self.Ncnpy) # Moisture flux from vegetation layers (kg/m^2/s)
        
        # Counters
        #self.k = np.zeros(1) # Canopy layer counter
        #self.ne = np.zeros(1) # Energy balance iteration counter
        self.eps = np.finfo(float).eps

        # State fluxes
        #self.B = np.zeros(1) # Kinematic bouyancy flux (Km/s)
        self.d = np.zeros(1) #  Displacement height (m)
        #self.Dsrf = np.zeros(1) # dQsat/dT at ground surface temperature (1/K)
        #self.dEs = np.zeros(1) # Change in surface moisture flux (kg/m^2/s)
        #self.dGs = np.zeros(1) # Change in surface heat flux (kg/m^2/s)
        #self.dHs = np.zeros(1) # Change in surface sensible heat flux (kg/m^2/s)
        #self.dTs = np.zeros(1) # Change in surface temperature (K)
        #self.E = np.zeros(1) # Moisture flux to the atmosphere (kg/m^2/s)
        #self.ebal = np.zeros(1) # Surface energy balance closure (W/m^2)
        #self.Ecan = np.zeros(1) # Within-canopy moisture flux (kg/m^2/s)
        self.fveg = np.zeros(1) # Vegetation weighting
        #self.ga = np.zeros(1) # Aerodynamic conductance to the atmosphere (m/s)
        #self.gc = np.zeros(1) # Conductance within canopy air space (m/s)
        #self.gs = np.zeros(1) # Surface to canopy air space conductance (m/s)
        #self.Hcan = np.zeros(1) # Within-canopy sensible heat flux (W/m^2)
        #self.Hsrf = np.zeros(1) # Sensible heat flux from the surface (W/m^2)
        #self.Kh = np.zeros(1) # Eddy diffusivity at canopy top (m^2/s)
        #self.Lsrf = np.zeros(1) # Latent heat for phase change on ground (J/kg)
        #self.psih = np.zeros(1) # Stability function for heat
        #self.psim = np.zeros(1) # Stability function for momentum
        #self.Qsrf = np.zeros(1) # Saturation humidity at surface temperature
        #self.rd = np.zeros(1) # Dense vegetation aerodynamic resistance (s/m)
        #self.rho = np.zeros(1) # Air density (kg/m^3)
        #self.rL = np.zeros(1) # Reciprocal of Obukhov length (1/m)
        #self.ro = np.zeros(1) # Open aerodynamic resistance (s/m)
        #self.Rsrf = np.zeros(1) # Net radiation absorbed by the surface (W/m^2)
        #self.Ssub = np.zeros(1) # Mass of snow available for sublimation (kg/m^2)
        #self.Uc = np.zeros(1) # Within-canopy wind speed (m/s)
        #self.Uh = np.zeros(1) # Wind speed at canopy top (m/s)
        #self.usd = np.zeros(1) # Dense canopy friction velocity (m/s)
        #self.uso = np.zeros(1) # Friction velocity (m/s)
        #self.ustar = np.zeros(1) # Open friction velocity (m/s)
        #self.wsrf = np.zeros(1) # Surface water availability factor
        
        # Below are defined before run timestep!
        #self.zT1 = np.zeros(1) # Temperature measurement height with offset (m)
        #self.zU1 = np.zeros(1) # Wind measurement height with offset (m)
        #self.z0g = np.zeros(1) # Snow/ground surface roughness length (m)
        #self.z0h = np.zeros(1) # Roughness length for heat (m)
        #self.z0v = np.zeros(1) # Vegetation roughness length (m)

        #self.dEv = np.array(Ncnpy) # Change in vegetation moisture flux (kg/m^2/s)
        #self.dHv = np.array(Ncnpy) # Change in veg sensible heat flux (kg/m^2/s)
        #self.dQc = np.array(Ncnpy) # Change in canopy air humidity (kg/kg)
        #self.dTv = np.array(Ncnpy) # Change in vegetation temperature (K)
        #self.dTc = np.array(Ncnpy) # Change in canopy air temperature (K)
        #self.Dveg = np.array(Ncnpy) # dQsat/dT at vegetation layer temperature (1/K)
        #self.gv = np.array(Ncnpy) # Vegetation to canopy air space conductance (m/s)
        self.Hveg = np.zeros(self.Ncnpy) # Sensible heat flux from vegetation (W/m^2)
        #self.Lcan = np.array(Ncnpy) # Latent heat for canopy water phase change (J/kg)
        #self.Qveg = np.array(Ncnpy) # Saturation humidity at vegetation temperature
        #self.Rveg = np.array(Ncnpy) # Net radiation absorbed by vegetation (W/m^2)
        #self.wveg = np.array(Ncnpy) # Vegetation water availability factor
        self.zh = np.array(self.Ncnpy) # Vegetation layer heights (m)

        #self.J = np.array(3*Ncnpy+1,3*Ncnpy+1) # Jacobian of energy and mass balance equations
        #self.f = np.array(3*Ncnpy+1) # Residuals of energy and mass balance equations
        #self.x = np.array(3*Ncnpy+1) # Temperature and humidity increments

        self.Tsrf = 273.15
        self.ZOFFST = 0
        self.CANMOD = 0
        self.EXCHNG = 1

        if self.ZOFFST == 0:
            # Heights specified above ground
            self.zU1 = self.zU
            self.zT1 = self.zT

        if self.ZOFFST == 1:
            # Heights specified above canopy top
            self.zU1 = self.zU + self.vegh
            self.zT1 = self.zT + self.vegh

        if self.CANMOD == 1:
            self.zh[0] = self.hbas + 0.5 * (self.vegh - self.hbas)

        if self.CANMOD == 2:
            self.zh[0] = (1 - 0.5 * self.fvg1) * self.vegh
            self.zh[1] = 0.5 * (1 - self.fvg1) * self.vegh    

    def run_timestep(self, cveg, Ds1, dt, fcans, fsnow, gs1, ks1, lveg, LW, Ps, Qa, SWsrf, 
                     Sveg, SWveg, Ta, tdif, Ts1, Tveg0, Ua, VAI, vegh, zT, zU, Sice):
        '''
        '''

        # Roughness lengths
        self.fveg = 1 - np.exp(-self.kext * self.VAI)
        self.d = 0.67 * self.fveg * self.vegh
        self.z0g = (self.z0sn**fsnow) * (self.z0sf**(1 - self.fsnow))
        self.z0h = 0.1 * self.z0g
        self.z0v = ((0.05 * self.vegh)**self.fveg) * (self.z0g**(1 - self.fveg))

        self.d = 0.67 * self.vegh
        self.z0v = 0.1 * self.vegh
        
        # Saturation humidity and air density
        Qsrf = self.qsat(Ps=Ps, T=self.Tsrf)
        Lsrf = self.Ls
        if (self.Tsrf > self.Tm):
            Lsrf = self.Lv
        Dsrf = Lsrf * Qsrf / (self.Rwat * self.Tsrf**2)
        rho = Ps / (self.Rair * Ta)

        if (VAI == 0):  # open
            self.Eveg[:] = 0
            self.Hveg[:] = 0
            ustar = self.vkman * Ua / np.log(self.zU1/self.z0g)
            print('ustar before', ustar)
            ga = self.vkman * ustar / np.log(self.zT1/self.z0h)

            for ne in range(20):
                if self.EXCHNG == 1:
                    if (ne < 10):
                        B = ga * (self.Tsrf - Ta)
                    if (Ta * ustar**3) != 0:
                        rL = -self.vkman * B / (Ta * ustar**3)
                        rL = max(min(rL,2.),-2.)
                    else:
                        rL = 0
                    ustar = self.vkman * Ua / (np.log(self.zU1/self.z0g) - self.psim(self.zU1, rL) + self.psim(self.z0g, rL))
                    ga = self.vkman * ustar / (np.log(self.zT1 / self.z0h) - self.psih(self.zT1, rL) + self.psih(self.z0h, rL)) # !+ 2/(rho * self.cp)
                    if ~np.isfinite(ga):
                        print('-self.vman', -self.vkman)
                        print('ustar', ustar)
                        print('rL', rL)
                        print('B', B)
                        print('Ua', Ua)
                        print('Ta', Ta)
                        print('psih zT1', self.psih(self.zT1, rL))
                        print('psih z0h', self.psih(self.z0h, rL))
                        break
                
                # Surface water availability
                if (Qa > Qsrf):
                    wsrf = 1
                else:
                    wsrf = fsnow + (1 - fsnow) * gs1 / (gs1 + ga)

                # Explicit fluxes
                Esrf = rho * wsrf * ga * (Qsrf - Qa)
                self.Eveg[:] = 0
                Gsrf = 2 * ks1 * (self.Tsrf - Ts1) / Ds1
                Hsrf = self.cp * rho * ga * (self.Tsrf - Ta)
                self.Hveg[:] = 0
                Melt = 0
                Rsrf = SWsrf + LW - self.sb * self.Tsrf**4

                # Surface energy balance increments without melt
                dTs = (Rsrf - Gsrf - Hsrf - Lsrf * Esrf) / \
                        (4 * self.sb * self.Tsrf**3 + 2 * ks1 / Ds1 + rho * (self.cp + Lsrf * Dsrf * wsrf) * ga) # NOTE CHECK THIS / &
                dEs = rho * wsrf * ga * Dsrf * dTs
                dGs = 2 * ks1 * dTs / Ds1 
                dHs = self.cp * rho * ga * dTs

                # Surface melting
                if (self.Tsrf + dTs > self.Tm) & (Sice[0] > 0):
                    Melt = sum(Sice) / dt
                    dTs = (Rsrf - Gsrf - Hsrf - Lsrf * Esrf - self.Lf * Melt) \
                            / (4 * self.sb * self.Tsrf**3 + 2 * ks1/Ds1 + rho*(self.cp + self.Ls * Dsrf * wsrf) * ga) # NOTE line change
                    dEs = rho * wsrf * ga * Dsrf * dTs
                    dGs = 2 * ks1 * dTs/Ds1
                    dHs = self.cp * rho * ga * dTs
                    if (self.Tsrf + dTs < self.Tm):
                        Qsrf = self.qsat(Ps=Ps,T=self.Tm)
                        Esrf = rho * wsrf * ga * (Qsrf - Qa)  
                        self.Gsrf = 2 * ks1 * (self.Tm - Ts1)/Ds1
                        Hsrf = self.cp * rho * ga * (self.Tm - Ta)
                        Rsrf = SWsrf + LW - self.sb * self.Tm**4 
                        Melt = (Rsrf - Gsrf - Hsrf - Lsrf * Esrf) / self.Lf
                        Melt = np.maximum(Melt, 0)
                        dEs = 0
                        dGs = 0
                        dHs = 0
                        dTs = self.Tm - self.Tsrf
                
                # Update surface temperature and fluxes
                Esrf = Esrf + dEs
                Gsrf = Gsrf + dGs
                Hsrf = Hsrf + dHs
                self.Tsrf = self.Tsrf + dTs
                # Diagnostics
                ebal = SWsrf + LW - self.sb * self.Tsrf**4 - Gsrf - Hsrf - Lsrf * Esrf - self.Lf * Melt
                LWout = self.sb * self.Tsrf**4
                LWsub = LW
                #print('self.zsub', self.zsub)
                #print('rL', rL)
                #print('self.z0g', self.z0g)
                Usub = (ustar / self.vkman) * (np.log(self.zsub/self.z0g) - self.psim(self.zsub,rL) + self.psim(self.z0g,rL))

                if (ne > 4) & (abs(ebal) < 0.01):
                    break
            print('rho', rho)
            print('wsrf', wsrf)
            print('ga', ga)
            print('Qsrf', Qsrf)
            print('Qa', Qa)                
            print('Explicit fluxes, Esrf:', Esrf)
            print('ne', ne)
            print('abs(ebal)', abs(ebal))
        '''
        else: # forest # NOTE BELOW THIS VARIABLES NEED TO BE CHECKED
            rL = 0
            usd = vkman * Ua / np.log((self.zU1-d)/self.z0v)
            Kh = vkman * usd * (self.vegh - d)
            rd = np.log((zT1-d)/(vegh-d)) / (vkman * usd) + vegh * (np.exp(wcan*(1 - zh[0])/self.vegh) - 1)/(wcan*Kh)
            uso = vkman * Ua / np.log(self.zU1/self.z0g)
            ro = np.log(self.zT1/self.zh[0]) / (vkman * uso)
            ga = fveg/rd + (1 - fveg)/ro
            for ne in range(20):
                # Aerodynamic resistance
                if EXCHNG == 1:
                    ustar = fveg * usd + (1 - fveg) * uso
                    if (ne<10): 
                        B = ga * (Tcan[0] - Ta)
                        rL = -vkman * B / (Ta * ustar**3)
                        rL = max(min(rL,2.),-2.)
                        usd = vkman * Ua / (np.log((self.zu1-d)/self.z0v) - self.psim(self.zU1-d,rL) + self.psim(self.z0v,rL))
                        if (rL > 0):
                            Kh = vkman * usd * (vegh - d) / (1 + 5 * (vegh - d) * rL)
                        else:
                            Kh = vkman * usd * (vegh - d) * sqrt(1 - 16 * (vegh - d) * rL)
                    rd = (np.log((self.zT1-d)/(self.vegh-d)) - self.psih(self.zT1-d,rL) + self.psih(self.vegh-d,rL))/(vkman*usd) \
                            + self.vegh * (np.exp(wcan*(1 - self.zh[0])/self.vegh) - 1)/(wcan*Kh)
                    uso = vkman * Ua / (np.log(self.zU1/self.z0g) - self.psim(self.zU1,rL) + self.psim(self.z0g,rL))
                    ro = (np.log(self.zT1/self.zh[0]) - self.psih(self.zT1,rL) + self.psih(self.zh[0]),rL)/(vkman*uso)
                    ga = self.fveg / rd + (1 - fveg)/ro # + 2/(rho*cp)
                Uh = (usd / vkman) * (np.log((vegh-d)/z0v) - self.psim(self.vegh-d,rl) + self.psim(self.z0v,rl))
                for k in range(self.Ncnpy):
                    Uc = fveg * np.exp(wcan*(zh[k]/vegh - 1))* Uh  +  \
                    (1 - fveg) * (uso / vkman) * (np.log(zh[k]/z0g) - self.psim(self.zh[k],rL) + self.psim(self.z0g,rL))
                    gv[k] = sqrt(Uc)*lveg[k]/leaf
                #if CANMOD == 2:
                #    rd = vegh * np.exp(wcan) * (np.exp(-wcan*zh[1])/vegh) - np.exp(-wcan*zh[0]/vegh)/(wcan*Kh)
                #    ro = (np.log(zh[0])/zh[1]) - self.psih(self.zh[1],rL) + self.psih(self.zh[2],rL))/(vkman*uso)
                #    gc = fveg/rd + (1 - fveg)/ro
                k = Ncnpy
                Uc = np.exp(wcan*(hbas/vegh - 1)) * Uh
                rd = np.log(hbas/z0g) * np.log(hbas/z0h) / (vkman**2 * Uc) \
                        + vegh * np.exp(wcan) * (np.exp(-wcan*hbas/vegh) - np.exp(-wcan * zh[k])) / (wcan * Kh)
                ro = (np.log(zh[k] / z0h) - self.psih(self.zh[k],rL) + self.psih(self.z0h,rL)) / (vkman * uso)
                gs = fveg / rd + (1 - fveg) / ro

                # rd = log((zT1-d)/z0v)/(vkman*usd)  !!!!!!!!!!!!!
                # rd = (log((zT1-d)/(vegh-d)) - psih(zT1-d,rL) + psih(vegh-d,rL))/(vkman*usd) + log(vegh/zh(1))/(vkman*usd)  !!!!!!!!!
                # Uc = fveg*(usd/vkman)*log(zh(k)/z0g) + (1 - fveg)*(uso/vkman)*log(zh(k)/z0g)  !!!!!!!!
                # rd = log(zh(1)/zh(2))/(vkman*usd)  !!!!!!!!!!!!!!!!!!
                # rd = 4*log(zh(k)/z0h)/(vkman*usd)   !!!!!!!!!!!!!!!!!!! 

                # ga = (vkman*usd)/(log((zT1-d)/z0v) - psih(zT1-d,rL) + psih(z0v,rL))
                # gv(:) = (vkman*usd)/log(1/0.999)
                # gc = 1
                # gs = (vkman*usd)/log(z0v/z0h)/4

                # Saturation humidity
                for k in range(self.Ncnpy):
                    self.qsat(Ps,Tveg[k],Qveg[k])
                    Lcan[k] = self.Ls
                    if (Tveg[k] > self.Tm):
                        Lcan[k] = Lv
                Dveg[:] = Lcan[:] * Qveg[:] / (Rwat * Tveg[:]**2)

                # Water availability
                if (Qcan[Ncnpy] > Qsrf):
                    wsrf = 1
                else:
                    wsrf = fsnow + (1 - fsnow) * gs1 / (gs1 + gs)

                for k in range(Ncnpy):
                    if (Qcan[k] > Qveg[k]):
                        wveg[k] = 1
                    else:
                        wveg[k] = fcans[k] + (1 - fcans[k]) * gsnf / (gsnf + gv[k])

                if CANMOD == 1:
                    # 1-layer canopy model
                    # Explicit fluxes
                    E = rho * ga * (Qcan[1] - Qa)
                    Esrf = rho * wsrf * gs * (Qsrf - Qcan[1])
                    Eveg[0] = rho * wveg[0] * gv[0] * (Qveg[0] - Qcan[0])
                    Gsrf = 2 * ks1 * (Tsrf - Ts1) / Ds1
                    H = rho * cp * ga * (Tcan[0] - Ta)
                    Hsrf = rho * cp * gs * (Tsrf - Tcan[0])
                    Hveg[0] = rho * cp * gv[0] * (Tveg[0] - Tcan[0])
                    Melt = 0
                    Rsrf = SWsrf + tdif[0] * LW - sb * Tsrf**4 + (1 - tdif[0]) * sb * Tveg[0]**4
                    Rveg[0] = SWveg[0] + (1 - tdif[0]) * (LW + sb * Tsrf**4 - 2 * sb * Tveg[0]**4)

                    # Surface energy balance increments without melt
                    J[0,0] = -rho * gs * (cp + Lsrf * Dsrf * wsrf) - 4 * sb * Tsrf**3 - 2 * ks1 / Ds1
                    J[0,1] = Lsrf * rho * wsrf * gs
                    J[0,2] = rho * cp * gs
                    J[0,3] = 4 * (1 - tdif[0]) * sb * Tveg[0]**3
                    J[1,0] = 4 * (1 - tdif[0]) * sb * Tsrf**3
                    J[1,1] = Lcan[0] * rho * wveg[0] * gv[0]
                    J[1,2] = rho * cp * gv[0]
                    J[1,3] = -rho * gv[0] * (cp + Lcan[0] * Dveg[0] * wveg[0]) \
                            -8 * (1 - tdif[0]) * sb * Tveg[0]**3 - cveg[0] / dt
                    J[2,0] = -gs
                    J[2,1] = 0
                    J[2,2] = ga + gs + gv[0]
                    J[2,3] = -gv[0]
                    J[3,0] = -Dsrf * wsrf * gs
                    J[3,1] = ga + wsrf * gs + wveg[0] * gv[0]
                    J[3,2] = 0
                    J[3,3] = -Dveg[0] * wveg[0] * gv[0]
                    f[1]   = -(Rsrf - Gsrf - Hsrf - Lsrf*Esrf)
                    f[2]   = -(Rveg[0] - Hveg[0] - Lcan[0] * Eveg[0] - \
                                cveg[0] * (Tveg[0] - Tveg0[0]) / dt)
                    f[3]   = -(H - Hveg[0] - Hsrf) / (rho * cp)
                    f[4]   = -(E - Eveg[0] - Esrf) / rho
                    x = self.ludcmp(4, J, f)
                    dTs = x[0]
                    dQc[0] = x[1]
                    dTc[0] = x[2]
                    dTv[0] = x[3]
                    dEs = rho * wsrf * gs * (Dsrf * dTs - dQc[0])
                    dEv[0] = rho * wveg[0] * gv[0] * (Dveg[0] * dTv[0] - dQc[0])
                    dGs = 2 * ks1 * dTs / Ds1
                    dHs = rho * cp * gs * (dTs - dTc[0])
                    dHv[0] = rho * cp * gv[0] * (dTv[0] - dTc[0])

                    # Surface melting
                    if (Tsrf + dTs > self.Tm) & (Sice[0] > 0):
                        Melt = sum(Sice) / dt
                        f[0] += Lf*Melt
                        x = self.ludcmp(4, J, f)
                        dTs = x[0]
                        dQc[0] = x[1]
                        dTc[0] = x[2]
                        dTv[0] = x[3]
                        dEs = rho * wsrf * gs * (Dsrf * dTs - dQc[0])
                        dEv[0] = rho * wveg[0] * gv[0] * (Dveg[0] * dTv[0] - dQc[0])
                        dGs = 2 * ks1 * dTs / Ds1
                        dHs = rho * cp * gs * (dTs - dTc[0])
                        dHv[0] = rho * cp * gv[0] * (dTv[0] - dTc[0])
                        if (Tsrf + dTs < self.Tm):
                            Qsrf = self.qsat(Ps=Ps, T=self.Tm)
                            Esrf = rho * wsrf * gs * (Qsrf - Qcan[0])
                            Gsrf = 2 * ks1 * (Tm - Ts1) / Ds1
                            Hsrf = rho * cp * gs * (Tm - Tcan[0])
                            Rsrf = SWsrf + tdif[0] * LW - sb * self.Tm**4 + (1 - tdif[0]) * sb * Tveg[0]**4
                            Rveg[0] = SWveg[0] + (1 - tdif[0]) * (LW + sb * self.Tm**4 - 2 * sb * Tveg[0]**4) 
                            J[0,0] = -1
                            J[1,0] = 0
                            J[2,0] = 0
                            J[3,0] = 0
                            f[0]   = -(Rsrf - Gsrf - Hsrf - Lsrf * Esrf)
                            f[1]   = -(Rveg[0] - Hveg[0] - Lcan[0] * Eveg[0] - cveg[0] * (Tveg[0] - Tveg0[0]) / dt)
                            f[2]   = -(H - Hveg[0] - Hsrf) / (rho * cp)
                            f[3]   = -(E - Eveg[0] - Esrf) / rho
                            x = self.ludcmp(4, J, f)
                            Melt = x[0] / Lf
                            dQc[0] = x[1]
                            dTc[0] = x[2]
                            dTv[0] = x[3]
                            dTs = self.Tm - Tsrf
                            dEs = 0
                            dEv[0] = rho * wveg[0] * gv[0] * (Dveg[0] * dTv[0] - dQc[0])
                            dGs = 0
                            dHs = 0
                            dHv[0] = rho * cp * gv[0] * (dTv[0] - dTc[0])
                            
                    LWout = (1 - tdif[0]) * sb * Tveg[0]**4 + tdif[0] * sb * Tsrf**4
                    LWsub = tdif[0] * LW + (1 - tdif[0]) * sb * Tveg[0]**4 

                # Update vegetation temperatures and fluxes
                Eveg[:] = Eveg[:] + dEv[:]
                Hveg[:] = Hveg[:] + dHv[:]
                Qcan[:] = Qcan[:] + dQc[:]
                Tcan[:] = Tcan[:] + dTc[:]
                Tveg[:] = Tveg[:] + dTv[:]

                # Update surface temperature and fluxes
                Esrf = Esrf + dEs
                Gsrf = Gsrf + dGs
                Hsrf = Hsrf + dHs
                Tsrf = Tsrf + dTs
                # Diagnostics
                ebal = SWsrf + LWsub - sb*Tsrf**4 - Gsrf - Hsrf - Lsrf * Esrf - Lf * Melt
                Uc = np.exp(wcan*(hbas/vegh - 1)) * Uh
                Usub = fveg * Uc * np.log(zsub/z0g) / np.log(hbas/z0g) \
                        + (1 - fveg) * Ua * (np.log(zsub/z0g) - self.psim(self.zsub,rL) + self.psim(self.z0g,rL)) \
                            / (np.log(self.zU/self.z0g) - self.psim(self.zU,rL) + self.psim(self.z0g,rL))

                if (ne>4) & (abs(ebal)<0.01):
                    break
            '''
        # Sublimation limited by available snow
        subl = 0
        Ssub = sum(Sice[:]) - Melt * dt
        if (Ssub > 0) | (self.Tsrf < self.Tm):
            Esrf = min(Esrf, Ssub / dt)
            subl = Esrf

        if (VAI > 0):
            for k in range(Ncnpy):
                if (Sveg[k] > 0) | (Tveg[k] < self.Tm):
                    Eveg[k] = min(Eveg[k], Sveg[k] / dt)
                    subl = subl + Eveg[k]
                
        # Fluxes to the atmosphere
        E = Esrf + np.sum(self.Eveg[:])
        H = Hsrf + np.sum(self.Hveg[:])
        LE = Lsrf * Esrf #+ sum(Lcan[:]*self.Eveg[:])
                
        Tsrf = self.Tsrf
        Eveg = self.Eveg
        Hveg = self.Hveg
                    
        return Esrf, Gsrf, H, LE, LWout, LWsub, Melt, subl, Usub, Eveg, Tsrf



    def qsat(self, Ps, T):
        '''
        '''
        Tc = T - self.Tm
        if (Tc > 0):
            es = self.e0 * np.exp(17.5043 * Tc / (241.3 + Tc))
        else:
            es = self.e0 * np.exp(22.4422 * Tc / (272.186 + Tc))

        Qsrf = self.eps * es / Ps

        return Qsrf


    def psim(self, z, rL):
        '''
        '''
        zeta = z * rL
        psim = np.where(zeta > 0, 
                        -5 * zeta, 
                        2 * np.log((1 + (1 - 16 * zeta) ** 0.25) / 2) + 
                        np.log((1 + (1 - 16 * zeta) ** 0.5) / 2) - 
                        2 * np.arctan((1 - 16 * zeta) ** 0.25) + 
                        np.pi / 2)
        
        return psim


    def psih(self, z, rL):
        '''
        '''
        zeta = z * rL
        psih = np.where(zeta > 0, 
                        -5 * zeta, 
                        2 * np.log((1 + (1 - 16 * zeta) ** 0.5) / 2))
        
        return psih


    def swrad():
        
        # acn0,              &! Snow-free dense canopy albedo
        # acns,              &! Snow-covered dense canopy albedo
        # asmx,              &! Maximum albedo for fresh snow
        # asmn,              &! Minimum albedo for melting snow
        # hfsn,              &! Snowcover fraction depth scale (m)
        # kext,              &! Vegetation light extinction coefficient
        # Salb,              &! Snowfall to refresh albedo (kg/m^2)
        # Talb,              &! Snow albedo decay temperature threshold (C)
        # tcld,              &! Cold snow albedo decay time scale (s)
        # tmlt                ! Melting snow albedo decay time scale (s)
        # alb0,              &! Snow-free ground albedo
        # dt,                &! Timestep (s)
        # elev,              &! Solar elevation (radians)
        # Sdif,              &! Diffuse shortwave radiation (W/m^2)
        # Sdir,              &! Direct-beam shortwave radiation (W/m^2)
        # Sf,                &! Snowfall rate (kg/m2/s)
        # Tsrf,              &! Snow/ground surface temperature (K)
        # Dsnw(Nsmax),       &! Snow layer thicknesses (m)
        # fcans(Ncnpy),      &! Canopy layer snowcover fractions
        # lveg(Ncnpy)         ! Canopy layer vegetation area indices
        # albs                ! Snow albedo
        # fsnow,             &! Ground snowcover fraction
        # SWout,             &! Outgoing SW radiation (W/m^2)
        # SWsrf,             &! SW absorbed by snow/ground surface (W/m^2)
        # SWsub,             &! Subcanopy downward SW radiation (W/m^2)
        # SWveg(Ncnpy),      &! SW absorbed by vegetation layers (W/m^2)
        # tdif(Ncnpy)         ! Canopy layer diffuse transmittances
        # alim,              &! Limiting snow albedo
        # asrf,              &! Snow/ground surface albedo
        # snd,               &! Snow depth (m)
        # tdec                ! Snow albedo decay time scale (s)

        # A(2*Ncnpy+1,2*Ncnpy+1),   &! Canopy radiative transfer matrix
        # b(2*Ncnpy+1),      &! Canopy layer boundary SW fluxes (W/m^2)
        # x(2*Ncnpy+1),      &! Canopy SW sources (W/m^2)
        # acan(Ncnpy),       &! Dense canopy albedo
        # rdif(Ncnpy),       &! Canopy layer diffuse reflectance
        # rdir(Ncnpy),       &! Canopy layer direct-beam reflectance
        # tdir(Ncnpy)         ! Canopy layer direct-beam transmittance

        if ALBEDO == 1:
            # Diagnostic snow albedo
            albs = asmn + (asmx - asmn)*(Tsrf - Tm) / Talb
        
        if ALBEDO == 2:
            # Prognostic snow albedo
            tdec = tcld
            if (Tsrf >= Tm):
                tdec = tmlt
            alim = (asmn/tdec + asmx*Sf/Salb)/(1/tdec + Sf/Salb)
            albs = alim + (albs - alim)*exp(-(1/tdec + Sf/Salb)*dt)
        albs = max(min(albs,asmx),asmn)

        # Partial snowcover on ground
        snd = sum(Dsnw[:])
        if SNFRAC == 1:
            fsnow = min(snd/hfsn, 1.)
        if SNFRAC == 2:
            fsnow = snd / (snd + hfsn)

        # Surface and vegetation net shortwave radiation
        asrf = (1 - fsnow)*alb0 + fsnow*albs
        SWsrf = (1 - asrf)*(Sdif + Sdir)
        SWveg[:] = 0
        SWout = asrf*(Sdif + Sdir)
        SWsub = Sdif + Sdir
        tdif[:] = 0
        tdir[:] = 0
        '''
        if (lveg(1) > 0):
            if CANRAD == 1:
                acan[:] = (1 - fcans[:])*acn0 + fcans[:]*acns
                tdif[:] = exp(-1.6*kext*lveg[:])
                tdir[:] = tdif[:]
                if (elev > 0) tdir[:] = exp(-kext*lveg[:]/sin(elev))
                    rdif[:] = (1 - tdif[:])*acan[:]
                    rdir[:] = (1 - tdir[:])*acan[:]
            if CANRAD == 2:
                for k in range(Ncnpy):
                    self.twostream(elev,fcans(k),lveg(k),rdif(k),rdir(k),tdif(k),tdir(k))
            A[:,:] = 0
            for k in range(2*Ncnpy + 1):
                A[k,k] = 1
            if CANMOD == 1:
                A(1,2) = -rdif(1)
                A(2,1) = -asrf
                A(3,2) = -tdif(1)
                b(1) = tdif(1)*Sdif
                b(2) = asrf*tdir(1)*Sdir
                b(3) = rdif(1)*Sdif + rdir(1)*Sdir
                LUDCMP(3,A,b,x)
                SWout = x(3)
                SWveg(1) = Sdif - x(1) + x(2) - x(3) + (1 - tdir(1))*Sdir
                SWsub = x(1) + tdir(1)*Sdir
                SWsrf = (1 - asrf)*SWsub
            if CANMOD == 2:
                A(1,4) = -rdif(1)
                A(2,1) = -tdif(2)
                A(2,3) = -rdif(2)
                A(3,2) = -asrf
                A(4,1) = -rdif(2)
                A(4,3) = -tdif(2)
                A(5,4) = -tdif(1)
                b(1) = tdif(1)*Sdif
                b(2) = 0
                b(3) = asrf*tdir(1)*tdir(2)*Sdir
                b(4) = rdir(2)*tdir(1)*Sdir
                b(5) = rdif(1)*Sdif + rdir(1)*Sdir
                LUDCMP(5,A,b,x)
                SWout = x(5)
                SWveg(1) = Sdif - x(1) + x(4) - x(5) + (1 - tdir(1))*Sdir
                SWveg(2) = x(1) - x(2) + x(3) - x(4) + tdir(1)*(1 - tdir(2))*Sdir
                SWsub = x(2) + tdir(1)*tdir(2)*Sdir
                SWsrf = (1 - asrf)*SWsub
        '''
    
    
    
    def ludcmp(N, A, b):
        '''
        #
        Solve matrix equation Ax = b for x by LU decomposition
        #
        
        Args:
        N # Number of equations to solve
        A(N,N) # Matrix
        b(N) # RHS of matrix equation
        Out:
        x(N) # Solution of matrix equation

        integer :: i,ii,imax,j,k,ll,indx(N)

        real :: Acp(N,N),aamax,dum,sum,vv(N)
        '''

        Acp = A[:,:]
        x = b[:]

        vv = np.zeros(N)
        indx = np.zeros(N, dtype=int)

        for i in range(N):
            aamax = 0
            for j in range(N):
                if (abs(Acp[i,j]) > aamax):
                    aamax = abs(Acp[i,j])
            vv[i] = 1/aamax

        for j in range(N):
            for i in range(j):
                sum = Acp[i,j]
                if (i > 1):
                    for k in range(i):
                        sum -= Acp[i,k] * Acp[k,j]
                    Acp[i,j] = sum
                    
            aamax = 0
            for i in range(j, N):
                sum = Acp[i,j]
                for k in range(j):
                    sum -= Acp[i,k] * Acp[k,j]
                Acp[i,j] = sum

                dum = vv[i] * abs(sum)
                if dum >= aamax:
                    imax = i
                    aamax = dum
            if j != imax:
                for k in range(N):
                    dum = Acp[imax, k]
                    Acp[[imax, j], :] = Acp[[j, imax], :]
                vv[imax] = vv[j]

            indx[j] = imax
            if (Acp[j,j] == 0):
                Acp[j,j] = 1e-20
            if j != N-1:
                dum = 1 / Acp[j,j]
                for i in range(j+1, N):
                    Acp[i,j] *= dum

        ii = 0
        for i in range(N):
            ll = indx[i]
            sum = x[ll]
            x[ll] = x[i]

            if ii != 0:
                for j in range(ii, i):
                    sum -= Acp[i,j] * x[j]
            elif sum != 0:
                ii = i

            x[i] = sum

        for i in range(N-1, 0, -1):
            sum = x[i]
            for j in range(i+1, N):
                x[i] = sum / Acp[i,i]

        return x


