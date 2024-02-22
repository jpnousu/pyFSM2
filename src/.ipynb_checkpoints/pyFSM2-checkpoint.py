#----------------------------------------------------------------------!
# Flexible Snow Model (FSM version 2.1.0)                              !
#                                                                      !
# Richard Essery                                                       !
# School of GeoSciences                                                !
# University of Edinburgh                                              !
#----------------------------------------------------------------------!

from pyFSM2_MODULES import Constants, IOUnits, Layers, Parameters, Soilprops

import numpy as np

class pyFSM2():
    def __init__(self, Constants, IOUnits, Layers, Parameters, Soilprops):
        
        # Grid dimensions
        Ncols = []  # Number of columns in grid
        Nrows = []  # Number of rows in grid

        # Site characteristics
        lat = []    # Latitude (radians)
        noon = []   # Time of solar noon (hour)

        # Meteorological driving data
        met_file = ""   # Meteorological driving file name
        year = []       # Year
        month = []      # Month of year
        day = []        # Day of month
        EoF = []        # End-of-file flag
        dt = []         # Timestep (s)
        elev = []       # Solar elevation (radians)
        hour = []       # Hour of day
        zT = []         # Temperature and humidity measurement height (m)
        zU = []         # Wind speed measurement height (m)

        # Meteorological driving data arrays
        LW = []      # Incoming longwave radiation (W/m^2)
        Ps = []      # Surface pressure (Pa)
        Qa = []      # Specific humidity (kg/kg)
        Rf = []      # Rainfall rate (kg/m^2/s)
        Sdif = []    # Diffuse shortwave radiation (W/m^2)
        Sdir = []    # Direct-beam shortwave radiation (W/m^2)
        Sf = []      # Snowfall rate (kg/m^2/s)
        Ta = []      # Air temperature (K)
        trans = []   # Wind-blown snow transport rate (kg/m^2/s)
        Ua = []      # Wind speed (m/s)

        # Model state variables
        Nsnow = []   # Number of snow layers
        albs = []    # Snow albedo
        Tsrf = []    # Snow/ground surface temperature (K)
        Dsnw = []    # Snow layer thicknesses (m)
        Qcan = []    # Canopy air space humidities
        Rgrn = []    # Snow layer grain radii (m)
        Sice = []    # Ice content of snow layers (kg/m^2)
        Sliq = []    # Liquid content of snow layers (kg/m^2)
        Sveg = []    # Snow mass on vegetation layers (kg/m^2)
        Tcan = []    # Canopy air space temperatures (K)
        Tsnow = []   # Snow layer temperatures (K)
        Tsoil = []   # Soil layer temperatures (K)
        Tveg = []    # Vegetation layer temperatures (K)
        Vsmc = []    # Volumetric moisture content of soil layers
        fsat = []    # Initial soil layer moisture/saturation
        Tprf = []    # Initial soil layer temperatures (K)

        # Diagnostics
        H = []       # Sensible heat flux to the atmosphere (W/m^2)
        LE = []      # Latent heat flux to the atmosphere (W/m^2)
        LWout = []   # Outgoing LW radiation (W/m^2)
        LWsub = []   # Subcanopy downward LW radiation (W/m^2)
        Melt = []    # Surface melt rate (kg/m^2/s)
        Roff = []    # Runoff from snow (kg/m^2/s)
        snd = []     # Snow depth (m)
        snw = []     # Total snow mass on ground (kg/m^2)
        subl = []    # Sublimation rate (kg/m^2/s)
        svg = []     # Total snow mass on vegetation (kg/m^2)
        SWout = []   # Outgoing SW radiation (W/m^2)
        SWsub = []   # Subcanopy downward SW radiation (W/m^2)
        Usub = []    # Subcanopy wind speed (m/s)
        Wflx = []    # Water flux into snow layer (kg/m^2/s)

        # Vegetation characteristics
        alb0_file = ""  # Snow-free ground albedo map file name
        fsky_file = ""  # Skyview fraction map file name
        vegh_file = ""  # Canopy height map file name
        VAI_file = ""   # Vegetation area index map file name
        alb0 = []       # Snow-free ground albedo
        fsky = []       # Skyview not obstructed by remote vegetation
        vegh = []       # Canopy height (m)
        VAI = []        # Vegetation area index

        # Start and dump file names
        dump_file = ""  # End dump file name
        runid = ""      # Run identifier
        start_file = ""  # Start file name

        # NetCDF variables
        ncid = []          # Dataset ID
        rec = []           # Record number
        status = []        # Error status
        varid = [0] * 17  # Variable IDs

        #namelist    /drive/ met_file,dt,lat,noon,zT,zU
        #namelist /gridpnts/ Ncols,Nrows,Nsmax,Nsoil
        #namelist /gridlevs/ Dzsnow,Dzsoil,fvg1,zsub
        #namelist  /initial/ fsat,Tprf,start_file
        #namelist  /outputs/ dump_file,runid
        #namelist      /veg/ alb0,fsky,vegh,VAI,  &
        #                    alb0_file,fsky_file,vegh_file,VAI_file

        if SETPAR == 1:
            import FSM2_PARAMS

        # Grid dimensions
        Ncols = 1
        Nrows = 1
        Nsmax = 3
        Nsoil = 4
        #read(5,gridpnts) # I think this is supposed to overwrite above

        # Canopy, snow and soil layers
        if CANMOD == 1:
            Ncnpy = 1
        if CANMOD == 2:
            Ncnpy = 2
        fvg1 = 0.5
        zsub = 1.5
        allocate(Dzsnow(Nsmax))
        allocate(Dzsoil(Nsoil))
        if (Nsmax == 3) Dzsnow = (/0.1, 0.2, 0.4/)
        if (Nsoil == 4) Dzsoil = (/0.1, 0.2, 0.4, 0.8/)
        read(5, gridlevs)

        # Site and driving data characteristics
        met_file = 'met'
        dt = 3600
        lat = 0
        noon = 12
        zT = 2
        zU = 10
        #read(5,drive) # Overwrites above
        open(umet, file = met_file)
        read(umet,*) year,month,day,hour
        rewind(umet)
        lat = (3.14159/180)*lat

        # Allocate driving data arrays
        LW = np.empty((Ncols, Nrows))
        Ps = np.empty((Ncols, Nrows))
        Qa = np.empty((Ncols, Nrows))
        Rf = np.empty((Ncols, Nrows))
        Sdif = np.empty((Ncols, Nrows))
        Sdir = np.empty((Ncols, Nrows))
        Sf = np.empty((Ncols, Nrows))
        Ta = np.empty((Ncols, Nrows))
        trans = np.empty((Ncols, Nrows))
        Ua = np.empty((Ncols, Nrows))
        trans[:,:] = 0

        # Vegetation characteristics from defaults, namelist or named map files
        alb0 = np.empty((Ncols, Nrows))
        fsky = np.empty((Ncols, Nrows))
        vegh = np.empty((Ncols, Nrows))
        VAI = np.empty((Ncols, Nrows))

        alb0_file = 'none'
        fsky_file = 'none'
        vegh_file = 'none'
        VAI_file  = 'none'
        alb0[:,:] = 0.2
        fsky[:,:] = 1
        vegh[:,:] = 0
        VAI[:,:]  = 0

        #read(5,veg) 
        if (alb0_file /= 'none') call FSM2_MAP(alb0_file,Ncols,Nrows,alb0)
        if (fsky_file /= 'none') call FSM2_MAP(fsky_file,Ncols,Nrows,fsky)
        if (vegh_file /= 'none') call FSM2_MAP(vegh_file,Ncols,Nrows,vegh)
        if (VAI_file  /= 'none') call FSM2_MAP(VAI_file,Ncols,Nrows,VAI)

        # Soil properties
        b = 3.1 + 15.7*fcly - 0.3*fsnd
        hcap_soil = (2.128*fcly + 2.385*fsnd)*1e6 / (fcly + fsnd)
        sathh = 10**(0.17 - 0.63*fcly - 1.58*fsnd)
        Vsat = 0.505 - 0.037*fcly - 0.142*fsnd
        Vcrit = Vsat*(sathh/3.364)**(1/b)
        hcon_soil = (hcon_air**Vsat) * ((hcon_clay**fcly)*(hcon_sand**(1 - fcly))**(1 - Vsat))

        # Allocate state variable arrays
        albs = np.empty((Ncols, Nrows))
        Nsnow = np.empty((Ncols, Nrows))
        Tsrf = np.empty((Ncols, Nrows))
        Dsnw = np.empty((Nsmax, Ncols, Nrows))
        Qcan = np.empty((Nsnpy, Ncols, Nrows))
        Rgrn = np.empty((Nsmax, Ncols, Nrows))
        Sice = np.empty((Nsmax, Ncols, Nrows))
        Sliq = np.empty((Nsmax, Ncols, Nrows))
        Sveg = np.empty((Ncnpy, Ncols, Nrows))
        Tcan = np.empty((Ncnpy, Ncols, Nrows))
        Tsnow = np.empty((Nsmax, Ncols, Nrows))
        Tsoil = np.empty((Nsoil, Ncols, Nrows))
        Tveg = np.empty((Ncnpy, Ncols, Nrows))
        Vsmc = np.empty((Nsoil, Ncols, Nrows))

        # Default initialization of state variables
        albs[:,:]  = 0.8
        Dsnw[:,:]  = 0
        Nsnow[:,:] = 0
        Qcan[:,:]  = 0
        Rgrn[:,:]  = rgr0
        Sice[:,:]  = 0
        Sliq[:,:]  = 0
        Sveg[:,:]  = 0
        Tcan[:,:]  = 285
        Tsnow[:,:] = 273
        Tveg[:,:]  = 285


        # Missing values for vegetation at non-forest points
        for k in range(Ncnpy):
            mask = VAI == 0
            Sveg[k, mask, :] = -999. / Ncnpy
            Tveg[k, mask, :] = -999

        # Initial soil profiles from namelist
        fsat = np.empty(Nsoil)
        Tprf = np.empty(Nsoil)

        fsat[:] = 0.5
        Tprf[:] = 285
        start_file = 'none'
        read(5,initial)

        for k in range(Nsoil):
            Tsoil[k, :, :] = Tprf[k]
            Vsmc[k, :, :] = fsat[k] * Vsat
    
        Tsrf = Tsoil[1,:,:]

        # Initialize state variables from a named start file
        if (start_file != 'none'):
          open(udmp,file = start_file)
          read(udmp,*) albs
          read(udmp,*) Dsnw
          read(udmp,*) Nsnow
          read(udmp,*) Qcan # still fortran
          read(udmp,*) Rgrn
          read(udmp,*) Sice
          read(udmp,*) Sliq
          read(udmp,*) Sveg
          read(udmp,*) Tcan
          read(udmp,*) Tsnow
          read(udmp,*) Tsoil
          read(udmp,*) Tsrf
          read(udmp,*) Tveg
          read(udmp,*) Vsmc
          close(udmp)

        # Allocate diagnostic output arrays
        H = np.empty((Ncols, Nrows))
        LE = np.empty((Ncols, Nrows))
        LWout = np.empty((Ncols, Nrows))
        LWsub = np.empty((Ncols, Nrows))
        Melt = np.empty((Ncols, Nrows))
        Roff = np.empty((Ncols, Nrows))
        snd = np.empty((Ncols, Nrows))
        snw = np.empty((Ncols, Nrows))
        subl = np.empty((Ncols, Nrows))
        svg = np.empty((Ncols, Nrows))
        SWout = np.empty((Ncols, Nrows))
        SWsub = np.empty((Ncols, Nrows))
        Usub = np.empty((Ncols, Nrows))
        Wflx = np.empty((Nsmax, Ncols, Nrows))

        # Output files
        dump_file = 'dump'
        runid = 'none'
        read(5,outputs)
        if (runid == 'none'):
            runid = ''
        # Check whether NetCDF output is enabled
        if PROFNC == 1:
            # Verify that the grid has only one column and one row; otherwise, raise an error
            if (Ncols*Nrows) > 1:
                raise ValueError('NetCDF output only available for Nrows = Ncols = 1')

            # Call a function to prepare NetCDF output
            #call_FSM2_PREPNC(runid, year, month, day, hour, ncid, rec, varid)
        else:
            # If NetCDF output is not enabled, check if the maximum value of VAI is greater than zero
            if np.max(VAI) > 0:
                # Open output files with filenames based on the runid string
                ucan = open(runid + 'subc.txt', 'w')
                uflx = open(runid + 'flux.txt', 'w')
                usta = open(runid + 'stat.txt', 'w')

# here 16.1. 
# Run the model
EoF = .false.
do
  call FSM2_DRIVE(Ncols,Nrows,fsky,lat,noon,                           &
                  year,month,day,hour,elev,EoF,                        &
                  LW,Ps,Qa,Rf,Sdif,Sdir,Sf,Ta,Ua)  
  if (EoF) goto 1
  do i = 1, Nrows
  do j = 1, Ncols
    call FSM2_TIMESTEP(                                                &
                       ! Driving variables                             &
                       dt,elev,zT,zU,LW(j,i),Ps(j,i),Qa(j,i),          &
                       Rf(j,i),Sdif(j,i),Sdir(j,i),Sf(j,i),            &
                       Ta(j,i),trans(j,i),Ua(j,i),                     &
                       ! Vegetation characteristics                    &
                       alb0(j,i),vegh(j,i),VAI(j,i),                   &
                       ! State variables                               &
                       albs(j,i),Tsrf(j,i),Dsnw(:,j,i),Nsnow(j,i),     &
                       Qcan(:,j,i),Rgrn(:,j,i),Sice(:,j,i),            &
                       Sliq(:,j,i),Sveg(:,j,i),Tcan(:,j,i),            &
                       Tsnow(:,j,i),Tsoil(:,j,i),Tveg(:,j,i),          &
                       Vsmc(:,j,i),                                    &
                       ! Diagnostics                                   &
                       H(j,i),LE(j,i),LWout(j,i),LWsub(j,i),           &
                       Melt(j,i),Roff(j,i),snd(j,i),snw(j,i),          &
                       subl(j,i),svg(j,i),SWout(j,i),SWsub(j,i),       &
                       Usub(j,i),Wflx(:,j,i)                           )
  end do
  end do
#if PROFNC == 1
  call FSM2_WRITENC(Dsnw(:,1,1),dt,H(1,1),LE(1,1),LWout(1,1),Melt(1,1),&
                    ncid,Nsnow(1,1),Rgrn(:,1,1),Roff(1,1),Sice(:,1,1), &
                    Sliq(:,1,1),snd(1,1),snw(1,1),SWout(1,1),          &
                    Tsnow(:,1,1),Tsoil(:,1,1),Tsrf(1,1),varid,         &
                    Wflx(:,1,1),rec) 
#else
  call FSM2_OUTPUT(Ncols,Nrows,year,month,day,hour,                    &
                   H,LE,LWout,LWsub,Melt,Roff,snd,snw,subl,svg,SWout,  &
                   SWsub,Tsoil,Tsrf,Tveg,Usub,VAI)
#endif
end do
1 continue

# Write out state variables at end of run
open(udmp,file = trim(runid) // trim(dump_file))
write(udmp,*) albs
write(udmp,*) Dsnw
write(udmp,*) Nsnow
write(udmp,*) Qcan
write(udmp,*) Rgrn
write(udmp,*) Sice
write(udmp,*) Sliq
write(udmp,*) Sveg
write(udmp,*) Tcan
write(udmp,*) Tsnow
write(udmp,*) Tsoil
write(udmp,*) Tsrf
write(udmp,*) Tveg
write(udmp,*) Vsmc
close(udmp)

#if PROFNC == 1
status = nf90_close(ncid) 
#endif
















