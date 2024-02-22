import numpy as np
from pyFSM2_MODULES import Constants, Layers, Parameters

class SrfEbal:
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
            ga = self.vkman * ustar / np.log(self.zT1/self.z0h)

            for ne in range(20):
                if self.EXCHNG == 1:
                    if (ne < 10):
                        B = ga * (self.Tsrf - Ta)
                    rL = -self.vkman*B/(Ta*ustar**3)
                    rL = max(min(rL,2.),-2.)
                    ustar = self.vkman * Ua / (np.log(self.zU1/self.z0g) - self.psim(self.zU1, rL) + self.psim(self.z0g, rL))
                    ga = self.vkman * ustar / (np.log(self.zT1 / self.z0h) - self.psih(self.zT1, rL) + self.psih(self.z0h, rL)) # !+ 2/(rho * self.cp) # NOTE WHAT IS !+

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
                    dTs = (self.Rsrf - Gsrf - Hsrf - Lsrf * Esrf - self.Lf * Melt) \
                            / (4 * self.sb * self.Tsrf**3 + 2 * ks1/Ds1 + rho*(self.cp + self.Ls * Dsrf * wsrf) * ga) # NOTE line change
                    dEs = rho * wsrf * ga * Dsrf * dTs
                    dGs = 2 * ks1 * dTs/Ds1
                    dHs = self.cp * rho * ga * dTs
                    if (self.Tsrf + dTs < Tm):
                        Qsrf = qsat(Ps=Ps,T=Tm)
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
                Usub = (ustar / self.vkman) * (np.log(self.zsub/self.z0g) - self.psim(self.zsub,rL) + self.psim(self.z0g,rL))

                if (ne > 4) & (abs(ebal) < 0.01):
                    break

        else: # forest # NOTE BELOW THIS VARIABLES NEED TO BE CHECKED
            rL = 0
            usd = vkman * Ua / np.log((self.zU1-d)/self.z0v)
            Kh = vkman * usd * (self.vegh - d)
            rd = np.log((zT1-d)/(vegh-d)) / (vkman * usd) + vegh * (np.exp(wcan*(1 - zh[0])/self.vegh) - 1))/(wcan*Kh)
            uso = vkman * Ua / np.log(self.zU1/self.z0g)
            ro = np.log(self.zT1/self.zh[0]) / (vkman * uso)
            ga = fveg/rd + (1 - fveg)/ro
            for ne in range(20):
                # Aerodynamic resistance
                if EXCHNG == 1:
                    ustar = fveg * usd + (1 - fveg) * uso
                    if (ne<10): 
                        B = ga * (Tcan[0]) - Ta)
                        rL = -vkman * B / (Ta * ustar**3)
                        rL = max(min(rL,2.),-2.)
                        usd = vkman * Ua / (np.log((self.zu1-d)/self.z0v) - self.psim(self.zU1-d,rL) + self.psim(self.z0v,rL))
                        if (rL > 0):
                            Kh = vkman * usd * (vegh - d) / (1 + 5 * (vegh - d) * rL)
                        else:
                            Kh = vkman * usd * (vegh - d) * sqrt(1 - 16 * (vegh - d) * rL)
                    rd = (np.log((self.zT1-d)/(self.vegh-d)) - self.psih(self.zT1-d,rL) + self.psih(self.vegh-d,rL))/(vkman*usd) + \ 
                    self.vegh * (np.exp(wcan*(1 - self.zh[0])/self.vegh) - 1))/(wcan*Kh)
                    uso = vkman * Ua / (np.log(self.zU1/self.z0g) - self.psim(self.zU1,rL) + self.psim(self.z0g,rL))
                    ro = (np.log(self.zT1/self.zh[0]) - self.psih(self.zT1,rL) + self.psih(self.zh[0]),rL))/(vkman*uso)
                    ga = self.fveg / rd + (1 - fveg)/ro # + 2/(rho*cp)
                Uh = (usd / vkman) * (np.log((vegh-d)/z0v) - self.psim(self.vegh-d,rl) + self.psim(self.z0v,rl))
                for k in range(self.Ncnpy):
                    Uc = fveg * np.exp(wcan*(zh[k]/vegh - 1))* Uh  +  \
                    (1 - fveg) * (uso / vkman) * (np.log(zh[k]/z0g) - self.psim(self.zh[k],rL) + self.psim(self.z0g,rL))
                    gv[k] = sqrt(Uc)*lveg[k]/leaf
                if CANMOD == 2:
                    rd = vegh * np.exp(wcan) * (np.exp(-wcan*zh[1])/vegh) - np.exp(-wcan*zh[0])/vegh))/(wcan*Kh)
                    ro = (np.log(zh[0])/zh[1]) - self.psih(self.zh[1],rL) + self.psih(self.zh[2],rL))/(vkman*uso)
                    gc = fveg/rd + (1 - fveg)/ro
                k = Ncnpy
                Uc = np.exp(wcan*(hbas/vegh - 1)) * Uh
                rd = np.log(hbas/z0g) * np.log(hbas/z0h) / (vkman**2 * Uc) + \ 
                vegh * np.exp(wcan) * (np.exp(-wcan*hbas/vegh) - np.exp(-wcan * zh[k])) / (wcan * Kh)
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
                    call LUDCMP(4,J,f,x)
                    dTs = x(1)
                    dQc(1) = x(2)
                    dTc(1) = x(3)
                    dTv(1) = x(4)
                    dEs = rho*wsrf*gs*(Dsrf*dTs - dQc(1))
                    dEv(1) = rho*wveg(1)*gv(1)*(Dveg(1)*dTv(1) - dQc(1))
                    dGs = 2*ks1*dTs/Ds1
                    dHs = rho*cp*gs*(dTs - dTc(1))
                    dHv(1) = rho*cp*gv(1)*(dTv(1) - dTc(1))

                    # Surface melting
                    if (Tsrf + dTs > Tm .and. Sice(1) > 0) then
    Melt = sum(Sice) / dt
    f(1) = f(1) + Lf*Melt
    call LUDCMP(4,J,f,x)
    dTs = x(1)
    dQc(1) = x(2)
    dTc(1) = x(3)
    dTv(1) = x(4)
    dEs = rho*wsrf*gs*(Dsrf*dTs - dQc(1))
    dEv(1) = rho*wveg(1)*gv(1)*(Dveg(1)*dTv(1) - dQc(1))
    dGs = 2*ks1*dTs/Ds1
    dHs = rho*cp*gs*(dTs - dTc(1))
    dHv(1) = rho*cp*gv(1)*(dTv(1) - dTc(1))
    if (Tsrf + dTs < Tm) then
      call QSAT(Ps,Tm,Qsrf)
      Esrf = rho*wsrf*gs*(Qsrf - Qcan(1))
      Gsrf = 2*ks1*(Tm - Ts1)/Ds1
      Hsrf = rho*cp*gs*(Tm - Tcan(1))
      Rsrf = SWsrf + tdif(1)*LW - sb*Tm**4 + (1 - tdif(1))*sb*Tveg(1)**4
      Rveg(1) = SWveg(1) + (1 - tdif(1))*(LW + sb*Tm**4 - 2*sb*Tveg(1)**4) 
      J(1,1) = -1
      J(2,1) = 0
      J(3,1) = 0
      J(4,1) = 0
      f(1)   = -(Rsrf - Gsrf - Hsrf - Lsrf*Esrf)
      f(2)   = -(Rveg(1) - Hveg(1) - Lcan(1)*Eveg(1) -  &
                 cveg(1)*(Tveg(1) - Tveg0(1))/dt)
      f(3)   = -(H - Hveg(1) - Hsrf)/(rho*cp)
      f(4)   = -(E - Eveg(1) - Esrf)/rho
      call LUDCMP(4,J,f,x)
      Melt = x(1)/Lf
      dQc(1) = x(2)
      dTc(1) = x(3)
      dTv(1) = x(4)
      dTs = Tm - Tsrf
      dEs = 0
      dEv(1) = rho*wveg(1)*gv(1)*(Dveg(1)*dTv(1) - dQc(1))
      dGs = 0
      dHs = 0
      dHv(1) = rho*cp*gv(1)*(dTv(1) - dTc(1))
    end if
  end if
  LWout = (1 - tdif(1))*sb*Tveg(1)**4 + tdif(1)*sb*Tsrf**4
  LWsub = tdif(1)*LW + (1 - tdif(1))*sb*Tveg(1)**4 
#endif

'''
#if CANMOD == 2
! 2-layer canopy model

  ! Explicit fluxes
  E = rho*ga*(Qcan(1) - Qa)
  Ecan = rho*gc*(Qcan(2) - Qcan(1))
  Esrf = rho*wsrf*gs*(Qsrf - Qcan(2))
  Eveg(1) = rho*wveg(1)*gv(1)*(Qveg(1) - Qcan(1))
  Eveg(2) = rho*wveg(2)*gv(2)*(Qveg(2) - Qcan(2))
  Gsrf = 2*ks1*(Tsrf - Ts1)/Ds1
  H = rho*cp*ga*(Tcan(1) - Ta)
  Hcan = rho*cp*gc*(Tcan(2) - Tcan(1))
  Hsrf = rho*cp*gs*(Tsrf - Tcan(2))
  Hveg(1) = rho*cp*gv(1)*(Tveg(1) - Tcan(1))
  Hveg(2) = rho*cp*gv(2)*(Tveg(2) - Tcan(2))
  Melt = 0
  Rsrf = SWsrf + tdif(1)*tdif(2)*LW +              &
         (1 - tdif(1))*tdif(2)*sb*Tveg(1)**4 +     &
         (1 - tdif(2))*sb*Tveg(2)**4 - sb*Tsrf**4 
  Rveg(1) = SWveg(1) + (1 - tdif(1))*(LW - 2*sb*Tveg(1)**4    & 
            + (1 - tdif(2))*sb*Tveg(2)**4 + tdif(2)*sb*Tsrf**4) 
  Rveg(2) = SWveg(2) +                                               & 
            (1 - tdif(2))*(tdif(1)*LW + (1 - tdif(1))*sb*Tveg(1)**4  & 
               - 2*sb*Tveg(2)**4 + sb*Tsrf**4) 

! Surface energy balance increments without melt
  J(1,1) = - rho*gs*(cp + Lsrf*Dsrf*wsrf) - 4*sb*Tsrf**3 - 2*ks1/Ds1
  J(1,2) = 0
  J(1,3) = 0
  J(1,4) = 4*(1 - tdif(1))*tdif(2)*sb*Tveg(1)**3
  J(1,5) = Lsrf*rho*wsrf*gs
  J(1,6) = rho*cp*gs
  J(1,7) = 4*(1 - tdif(2))*sb*Tveg(2)**3
  J(2,1) = 4*(1 - tdif(1))*tdif(2)*sb*Tveg(2)**3 
  J(2,2) = Lcan(1)*rho*wveg(1)*gv(1)
  J(2,3) = rho*cp*gv(1)
  J(2,4) = - rho*gv(1)*(cp + Lcan(1)*Dveg(1)*wveg(1))    & 
           - 8*(1 - tdif(1))*sb*Tveg(1)**3 - cveg(1)/dt
  J(2,5) = 0
  J(2,6) = 0
  J(2,7) = 4*(1 - tdif(1))*(1 - tdif(2))*sb*Tveg(2)**3
  J(3,1) = 4*(1 - tdif(2))*sb*Tsrf**3 
  J(3,2) = 0
  J(3,3) = 0
  J(3,4) = 4*(1 - tdif(1))*(1 - tdif(2))*sb*Tveg(1)**3
  J(3,5) = Lcan(2)*rho*wveg(2)*gv(2)
  J(3,6) = rho*cp*gv(2)
  J(3,7) = - rho*gv(2)*(cp + Lcan(2)*Dveg(2)*wveg(2))    &
           - 8*(1 - tdif(2))*sb*Tveg(2)**3 - cveg(2)/dt
  J(4,1) = 0
  J(4,2) = 0
  J(4,3) = ga + gc + gv(1)
  J(4,4) = -gv(1)
  J(4,5) = 0
  J(4,6) = -gc
  J(4,7) = 0
  J(5,1) = -gs
  J(5,2) = 0
  J(5,3) = -gc
  J(5,4) = 0
  J(5,5) = 0
  J(5,6) = gc + gs + gv(2)
  J(5,7) = -gv(2)
  J(6,1) = 0
  J(6,2) = ga + gc + wveg(1)*gv(1)
  J(6,3) = 0
  J(6,4) = -Dveg(1)*wveg(1)*gv(1)
  J(6,5) = -gc
  J(6,6) = 0
  J(6,7) = 0
  J(7,1) = -Dsrf*wsrf*gs
  J(7,2) = -gc
  J(7,3) = 0
  J(7,4) = 0
  J(7,5) = gc + wsrf*gs + wveg(2)*gv(2)
  J(7,6) = 0
  J(7,7) = -Dveg(2)*wveg(2)*gv(2)
  f(1)   = -(Rsrf - Gsrf - Hsrf - Lsrf*Esrf)
  f(2)   = -(Rveg(1) - Hveg(1) - Lcan(1)*Eveg(1) -  & 
             cveg(1)*(Tveg(1) - Tveg0(1))/dt)
  f(3)   = -(Rveg(2) - Hveg(2) - Lcan(2)*Eveg(2) -  &
             cveg(2)*(Tveg(2) - Tveg0(2))/dt)
  f(4)   = -(H - Hcan - Hveg(1))/(rho*cp)
  f(5)   = -(Hcan - Hsrf - Hveg(2))/(rho*cp)
  f(6)   = -(E - Ecan - Eveg(1))/rho
  f(7)   = -(Ecan - Esrf - Eveg(2))/rho
  call LUDCMP(7,J,f,x)
  dTs    = x(1)
  dQc(1) = x(2)
  dTc(1) = x(3)
  dTv(1) = x(4)
  dQc(2) = x(5)
  dTc(2) = x(6)
  dTv(2) = x(7)
  dEs = rho*wsrf*gs*(Dsrf*dTs - dQc(2))
  dEv(1) = rho*wveg(1)*gv(1)*(Dveg(1)*dTv(1) - dQc(1))
  dEv(2) = rho*wveg(2)*gv(2)*(Dveg(2)*dTv(2) - dQc(2))
  dGs = 2*ks1*dTs/Ds1
  dHs = rho*cp*gs*(dTs - dTc(2))
  dHv(1) = rho*cp*gv(1)*(dTv(1) - dTc(1))
  dHv(2) = rho*cp*gv(2)*(dTv(2) - dTc(2))

  ! Surface melting
  if (Tsrf + dTs > Tm .and. Sice(1) > 0) then
    Melt = sum(Sice)/dt
    f(1) = f(1) + Lf*Melt
    call LUDCMP(7,J,f,x)
    dTs    = x(1)
    dQc(1) = x(2)
    dTc(1) = x(3)
    dTv(1) = x(4)
    dQc(2) = x(5)
    dTc(2) = x(6)
    dTv(2) = x(7)
    dEs = rho*wsrf*gs*(Dsrf*dTs - dQc(2))
    dEv(1) = rho*wveg(1)*gv(1)*(Dveg(1)*dTv(1) - dQc(1))
    dEv(2) = rho*wveg(2)*gv(2)*(Dveg(2)*dTv(2) - dQc(2))
    dGs = 2*ks1*dTs/Ds1
    dHs = rho*cp*gs*(dTs - dTc(2))
    dHv(1) = rho*cp*gv(1)*(dTv(1) - dTc(1))
    dHv(2) = rho*cp*gv(2)*(dTv(2) - dTc(2))
    if (Tsrf + dTs < Tm) then
      call QSAT(Ps,Tm,Qsrf)
      Esrf = rho*wsrf*gs*(Qsrf - Qcan(2))
      Hsrf = rho*cp*gs*(Tm - Tcan(2))
      Rsrf = Rsrf + sb*Tsrf**4 - sb*Tm**4
      Rveg(1) = Rveg(1) + (1 - tdif(1))*tdif(2)*sb*(Tm**4 - Tsrf**4) 
      Rveg(2) = Rveg(2) + (1 - tdif(2))*sb*(Tm**4 - Tsrf**4)
      J(1,1) = -1
      J(2,1) = 0
      J(3,1) = 0
      J(4,1) = 0
      J(5,1) = 0
      J(6,1) = 0
      J(7,1) = 0
      f(1) = -(Rsrf - Gsrf - Hsrf - Lsrf*Esrf)
      f(2) = -(Rveg(1) - Hveg(1) - Lcan(1)*Eveg(1) -  & 
               cveg(1)*(Tveg(1) - Tveg0(1))/dt)
      f(3) = -(Rveg(2) - Hveg(2) - Lcan(2)*Eveg(2) -  &
               cveg(2)*(Tveg(2) - Tveg0(2))/dt)
      f(4) = -(H - Hcan - Hveg(1))/(rho*cp)
      f(5) = -(Hcan - Hsrf - Hveg(2))/(rho*cp)
      f(6) = -(E - Ecan - Eveg(1))/rho
      f(7) = -(Ecan - Esrf - Eveg(2))/rho
      call LUDCMP(7,J,f,x)
      Melt = x(1)/Lf
      dQc(1) = x(2)
      dTc(1) = x(3)
      dTv(1) = x(4)
      dQc(2) = x(5)
      dTc(2) = x(6)
      dTv(2) = x(7)
      dTs = Tm - Tsrf
      dEs = 0
      dEv(1) = rho*wveg(1)*gv(1)*(Dveg(1)*dTv(1) - dQc(1))
      dEv(2) = rho*wveg(2)*gv(2)*(Dveg(2)*dTv(2) - dQc(2))
      dGs = 0
      dHs = 0
      dHv(1) = rho*cp*gv(1)*(dTv(1) - dTc(1))
      dHv(2) = rho*cp*gv(2)*(dTv(2) - dTc(2))
    end if
  end if
  LWout = (1 - tdif(1))*sb*Tveg(1)**4 +          &
          (1 - tdif(2))*tdif(1)*sb*Tveg(1)**4 +  &
          tdif(1)*tdif(2)*sb*Tsrf**4
  LWsub = tdif(1)*tdif(2)*LW +                   &
          (1 - tdif(1))*tdif(2)*sb*Tveg(1)**4 +  &
          (1 - tdif(2))*sb*Tveg(2)**4
#endif
'''
'''
  ! Update vegetation temperatures and fluxes
  Eveg(:) = Eveg(:) + dEv(:)
  Hveg(:) = Hveg(:) + dHv(:)
  Qcan(:) = Qcan(:) + dQc(:)
  Tcan(:) = Tcan(:) + dTc(:)
  Tveg(:) = Tveg(:) + dTv(:)

  ! Update surface temperature and fluxes
  Esrf = Esrf + dEs
  Gsrf = Gsrf + dGs
  Hsrf = Hsrf + dHs
  Tsrf = Tsrf + dTs
  ! Diagnostics
  ebal = SWsrf + LWsub - sb*Tsrf**4 - Gsrf - Hsrf - Lsrf*Esrf - Lf*Melt
  Uc = exp(wcan*(hbas/vegh - 1))*Uh
  Usub = fveg*Uc*log(zsub/z0g)/log(hbas/z0g) +  &
         (1 - fveg)*Ua*(log(zsub/z0g) - psim(zsub,rL) + psim(z0g,rL)) / &
                       (log(zU/z0g) - psim(zU,rL) + psim(z0g,rL))

if (ne>4 .and. abs(ebal)<0.01) exit
end do
end if  ! forest
!print*,ne,ebal
!write(31,*) SWsrf,LWsub - sb*Tsrf**4,Gsrf,Hsrf,Lsrf*Esrf,Lf*Melt

! Sublimation limited by available snow
subl = 0
Ssub = sum(Sice(:)) - Melt*dt
if (Ssub > 0 .or. Tsrf<Tm) then
  Esrf = min(Esrf, Ssub/dt)
  subl = Esrf
end if
if (VAI>0) then
  do k = 1, Ncnpy
    if (Sveg(k)>0 .or. Tveg(k)<Tm) then
      Eveg(k) = min(Eveg(k), Sveg(k)/dt)
      subl = subl + Eveg(k)
    end if
  end do
end if                
'''
                # Fluxes to the atmosphere
                E = Esrf + np.sum(self.Eveg[:])
                H = Hsrf + np.sum(self.Hveg[:])
                LE = Lsrf * Esrf #+ sum(Lcan[:]*self.Eveg[:])
                
                Tsrf = self.Tsrf.copy()
                Eveg = self.Eveg.copy()
                Hveg = self.Hveg.copy()
                    
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
        zeta = z * rL
        psim = np.where(zeta > 0, -5 * zeta, 2 * np.log((1 + (1 - 16 * zeta) ** 0.25) / 2) + np.log((1 + (1 - 16 * zeta) ** 0.5) / 2) - 2 * np.arctan((1 - 16 * zeta) ** 0.25) + np.pi / 2)
        
        return psim

    def psih(self, z, rL):
        zeta = z * rL
        psih = np.where(zeta > 0, -5 * zeta, 2 * np.log((1 + (1 - 16 * zeta) ** 0.5) / 2))
        
        return psih
