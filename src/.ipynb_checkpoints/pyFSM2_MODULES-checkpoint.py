import numpy as np
#SETPAR = 2
#DENSTY = 0

# Physical constants
class Constants:
    def __init__(self):
        self.cp = 1005  # Specific heat capacity of air (J/K/kg)
        self.eps = 0.622  # Ratio of molecular weights of water and dry air
        self.e0 = 611.213  # Saturation vapour pressure at Tm (Pa)
        self.g = 9.81  # Acceleration due to gravity (m/s^2)
        self.hcap_ice = 2100  # Specific heat capacity of ice (J/K/kg)
        self.hcap_wat = 4180  # Specific heat capacity of water (J/K/kg)
        self.hcon_air = 0.025  # Thermal conductivity of air (W/m/K)
        self.hcon_clay = 1.16  # Thermal conductivity of clay (W/m/K)
        self.hcon_ice = 2.24  # Thermal conductivity of ice (W/m/K)
        self.hcon_sand = 1.57  # Thermal conductivity of sand (W/m/K)
        self.hcon_wat = 0.56  # Thermal conductivity of water (W/m/K)
        self.I0 = 1367  # Solar constant (W/m^2)
        self.Lf = 0.334e6  # Latent heat of fusion (J/kg)
        self.Lv = 2.501e6  # Latent heat of vapourisation (J/kg)
        self.Ls = self.Lf + self.Lv  # Latent heat of sublimation (J/kg)
        self.mu_wat = 1.78e-3  # Dynamic viscosity of water (kg/m/s)
        self.pi = 3.14159  # pi
        self.Rair = 287  # Gas constant for air (J/K/kg)
        self.Rwat = 462  # Gas constant for water vapour (J/K/kg)
        self.rho_ice = 917  # Density of ice (kg/m^3)
        self.rho_wat = 1000  # Density of water (kg/m^3)
        self.sb = 5.67e-8  # Stefan-Boltzmann constant (W/m^2/K^4)
        self.Tm = 273.15  # Melting point (K)
        self.vkman = 0.4  # von Karman constant

# Input/Output file unit numbers
class IOUnits:
    def __init__(self):    
        self.ucan = 11  # Subcanopy diagnostics file unit number
        self.udmp = 12  # Start / dump file unit number
        self.uflx = 13  # Flux output file unit number
        self.umap = 14  # Map input file unit number
        self.umet = 15  # Meteorological driving file unit number
        self.usta = 16  # State output file unit number

# Canopy, snow, and soil layers
class Layers:
    def __init__(self):    
        self.Ncnpy = 1  # Number of canopy layers
        self.Nsmax = 3  # Maximum number of snow layers
        self.Nsoil = 3  # Number of soil layers
        self.Dzsnow = np.array([0.1, 0.2, 0.4]) # Minimum snow layer thicknesses (m)
        self.Dzsoil = np.array([0.1, 0.2, 0.4])  # Soil layer thicknesses (m)
        self.fvg1 = []  # Fraction of vegetation in the upper canopy layer
        self.zsub = 2.0  # Subcanopy wind speed diagnostic height (m)

# Parameters
class Parameters:
    def __init__(self, SETPAR, DENSITY):    
        if SETPAR == 1:
            # Vegetation parameters
            self.acn0 = []  # Snow-free dense canopy albedo
            self.acns = []  # Snow-covered dense canopy albedo
            self.avg0 = []  # Canopy element reflectivity
            self.avgs = []  # Canopy snow reflectivity
            self.cvai = []  # Vegetation heat capacity per unit VAI (J/K/m^2)
            self.gsnf = []  # Snow-free vegetation moisture conductance (m/s)
            self.hbas = []  # Canopy base height (m)
            self.kext = []  # Vegetation light extinction coefficient
            self.leaf = []  # Leaf boundary resistance (s/m)^(1/2)
            self.svai = []  # Intercepted snow capacity per unit VAI (kg/m^2)
            self.tunl = []  # Canopy snow unloading time scale (s)
            self.wcan = []  # Canopy wind decay coefficient
        
            # Snow parameters
            self.asmn = []  # Minimum albedo for melting snow
            self.asmx = []  # Maximum albedo for fresh snow
            self.eta0 = []  # Reference snow viscosity (Pa s)
            self.hfsn = []  # Snowcover fraction depth scale (m)
            self.kfix = []  # Fixed thermal conductivity of snow (W/m/K)
            self.rcld = []  # Maximum density for cold snow (kg/m^3)
            self.rfix = []  # Fixed snow density (kg/m^3)
            self.rgr0 = []  # Fresh snow grain radius (m)
            self.rhof = []  # Fresh snow density (kg/m^3)
            self.rhow = []  # Wind-packed snow density (kg/m^3)
            self.rmlt = []  # Maximum density for melting snow (kg/m^3)
            self.Salb = []  # Snowfall to refresh albedo (kg/m^2)
            self.snda = []  # Thermal metamorphism parameter (1/s)
            self.Talb = []  # Snow albedo decay temperature threshold (C)
            self.tcld = []  # Cold snow albedo decay time scale (s)
            self.tmlt = []  # Melting snow albedo decay time scale (s)
            self.trho = []  # Snow compaction timescale (s)
            self.Wirr = []  # Irreducible liquid water content of snow
            self.z0sn = []  # Snow roughness length (m)
        
            # Ground surface and soil parameters
            self.fcly = []  # Soil clay fraction
            self.fsnd = []  # Soil sand fraction
            self.gsat = []  # Surface conductance for saturated soil (m/s)
            self.z0sf = []  # Snow-free surface roughness length (m)
        
        else:
            # Vegetation parameters
            self.acn0 = 0.1        # Snow-free dense canopy albedo
            self.acns = 0.4        # Snow-covered dense canopy albedo
            self.avg0 = 0.21       # Canopy element reflectivity
            self.avgs = 0.6        # Canopy snow reflectivity
            self.cvai = 3.6e4      # Vegetation heat capacity per unit VAI (J/K/m^2)
            self.gsnf = 0.01       # Snow-free vegetation moisture conductance (m/s)
            self.hbas = 2          # Canopy base height (m)
            self.kext = 0.5        # Vegetation light extinction coefficient
            self.leaf = 20         # Leaf boundary resistance (s/m)^(1/2)
            self.svai = 4.4        # Intercepted snow capacity per unit VAI (kg/m^2)
            self.tunl = 240*3600   # Canopy snow unloading time scale (s)
            self.wcan = 2.5        # Canopy wind decay coefficient
        
            # Snow parameters
            self.asmn = 0.5        # Minimum albedo for melting snow
            self.asmx = 0.85       # Maximum albedo for fresh snow
            self.eta0 = 3.7e7      # Reference snow viscosity (Pa s)
            self.hfsn = 0.1        # Snowcover fraction depth scale (m)
            self.kfix = 0.24       # Fixed thermal conductivity of snow (W/m/K)
            self.rcld = 300        # Maximum density for cold snow (kg/m^3)
            self.rfix = 300        # Fixed snow density (kg/m^3)
            self.rgr0 = 5e-5       # Fresh snow grain radius (m)
            if DENSITY == 0:
                self.rhof = self.rfix   # Fresh snow density (kg/m^3)
            else:
                self.rhof = 100    # Fresh snow density (kg/m^3)
            self.rhow = 300        # Wind-packed snow density (kg/m^3)
            self.rmlt = 500        # Maximum density for melting snow (kg/m^3)
            self.Salb = 10         # Snowfall to refresh albedo (kg/m^2)
            self.snda = 2.8e-6     # Thermal metamorphism parameter (1/s)
            self.Talb = -2         # Snow albedo decay temperature threshold (C)
            self.tcld = 3.6e6      # Cold snow albedo decay time scale (s)
            self.tmlt = 3.6e5      # Melting snow albedo decay time scale (s)
            self.trho = 200*3600   # Snow compaction timescale (s)
            self.Wirr = 0.03       # Irreducible liquid water content of snow
            self.z0sn = 0.001      # Snow roughness length (m)
        
            # Ground surface and soil parameters
            self.fcly = 0.3        # Soil clay fraction
            self.fsnd = 0.6        # Soil sand fraction
            self.gsat = 0.01       # Surface conductance for saturated soil (m/s)
            self.z0sf = 0.1        # Snow-free surface roughness length (m)

# Soil properties
class SoilProps:
    def __init__(self):    
        self.b = []  # Clapp-Hornberger exponent
        self.hcap_soil = []  # Volumetric heat capacity of dry soil (J/K/m^3)
        self.hcon_soil = []  # Thermal conductivity of dry soil (W/m/K)
        self.sathh = []  # Saturated soil water pressure (m)
        self.Vcrit = []  # Volumetric soil moisture concentration at critical point
        self.Vsat = []  # Volumetric soil moisture concentration at saturation


