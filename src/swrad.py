#-----------------------------------------------------------------------
# Surface and vegetation net shortwave radiation
#-----------------------------------------------------------------------

from pyFSM2_MODULES import Constants, Layers, Parameters
import numpy as np

class SWrad:
    def __init__(self):

        constants = Constants()
        layers = Layers()
        params = Parameters(SETPAR=2, DENSITY=0)

        # Constants
        self.Tm = constants.Tm          # Melting point (K)

        # Layers
        self.Ncnpy = layers.Ncnpy       # Number of canopy layers
        self.Nsmax = layers.Nsmax       # Maximum number of snow layers

        # Parameters
        self.acn0 = params.acn0         # Snow-free dense canopy albedo
        self.acns = params.acns         # Snow-covered dense canopy albedo
        self.asmx = params.asmx         # Maximum albedo for fresh snow
        self.asmn = params.asmn         # Minimum albedo for melting snow
        self.hfsn = params.hfsn         # Snowcover fraction depth scale (m)
        self.kext = params.kext         # Vegetation light extinction coefficient
        self.Salb = params.Salb         # Snowfall to refresh albedo (kg/m^2)
        self.Talb = params.Talb         # Snow albedo decay temperature threshold (C)
        self.tcld = params.tcld         # Cold snow albedo decay time scale (s)
        self.tmlt = params.tmlt         # Melting snow albedo decay time scale (s)
 

        #alb0,              &! Snow-free ground albedo
        #dt,                &! Timestep (s)
        #elev,              &! Solar elevation (radians)
        #Sdif,              &! Diffuse shortwave radiation (W/m^2)
        #Sdir,              &! Direct-beam shortwave radiation (W/m^2)
        #Sf,                &! Snowfall rate (kg/m2/s)
        #Tsrf,              &! Snow/ground surface temperature (K)
        #Dsnw(Nsmax),       &! Snow layer thicknesses (m)
        #fcans(Ncnpy),      &! Canopy layer snowcover fractions
        #lveg(Ncnpy)         ! Canopy layer vegetation area indices

        #real, intent(inout) :: &
        #self.albs = np.zeros(1)                # Snow albedo

        #real, intent(out) :: &
        #fsnow,             &! Ground snowcover fraction
        #SWout,             &! Outgoing SW radiation (W/m^2)
        #SWsrf,             &! SW absorbed by snow/ground surface (W/m^2)
        #SWsub,             &! Subcanopy downward SW radiation (W/m^2)
        self.SWveg = np.zeros(self.Ncnpy)     # SW absorbed by vegetation layers (W/m^2)
        self.tdif = np.zeros(self.Ncnpy)        # Canopy layer diffuse transmittances

        #k          ! Canopy layer counter

        #alim,              &! Limiting snow albedo
        #asrf,              &! Snow/ground surface albedo
        #snd,               &! Snow depth (m)
        #tdec                ! Snow albedo decay time scale (s)

        self.A  = np.zeros([2*self.Ncnpy+1,2*self.Ncnpy+1])   # Canopy radiative transfer matrix
        self.b = np.zeros(2*self.Ncnpy+1)                   # Canopy layer boundary SW fluxes (W/m^2)
        self.x = np.zeros(2*self.Ncnpy+1)                        # Canopy SW sources (W/m^2)
        self.acan = np.zeros(self.Ncnpy)                         # Dense canopy albedo
        self.rdif = np.zeros(self.Ncnpy)                         # Canopy layer diffuse reflectance
        self.rdir = np.zeros(self.Ncnpy)                         # Canopy layer direct-beam reflectance
        self.tdir = np.zeros(self.Ncnpy)                         # Canopy layer direct-beam transmittance

        self.ALBEDO = 2
        self.SNFRAC = 1

    def run_timestep(self, albs, alb0, dt, elev, Sdif, Sdir, Sf, Tsrf, Dsnw, fcans, lveg):

        if self.ALBEDO == 1:
            # Diagnostic snow albedo
            albs = self.asmn + (self.asmx - self.asmn)*(Tsrf - self.Tm) / self.Talb
        
        if self.ALBEDO == 2:
            # Prognostic snow albedo
            tdec = self.tcld
            if (Tsrf >= self.Tm):
                tdec = self.tmlt
            alim = (self.asmn/tdec + self.asmx*Sf/self.Salb)/(1/tdec + Sf/self.Salb)
            albs = alim + (albs - alim)*np.exp(-(1/tdec + Sf/self.Salb)*dt)
        albs = max(min(albs,self.asmx),self.asmn)

        # Partial snowcover on ground
        hs = sum(Dsnw[:])
        if self.SNFRAC == 1:
            fsnow = min(hs/self.hfsn, 1.)
        if self.SNFRAC == 2:
            fsnow = hs / (hs + self.hfsn)

        # Surface and vegetation net shortwave radiation
        asrf = (1 - fsnow)*alb0 + fsnow * albs
        SWsrf = (1 - asrf)*(Sdif + Sdir)
        self.SWveg[:] = 0
        SWout = asrf*(Sdif + Sdir)
        SWsub = Sdif + Sdir
        self.tdif[:] = 0
        self.tdir[:] = 0
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
    
        return albs, fsnow, SWout, SWsrf, SWsub, self.SWveg, self.tdif