#-----------------------------------------------------------------------
# Set default parameter values or read from namelist
#-----------------------------------------------------------------------

class Parameters:
    def __init__(self):
        # Vegetation parameters
        self.acn0 = 0.1            # Snow-free dense canopy albedo
        self.acns = 0.4            # Snow-covered dense canopy albedo
        self.avg0 = 0.21           # Canopy element reflectivity
        self.avgs = 0.6            # Canopy snow reflectivity
        self.cvai = 3.6e4          # Vegetation heat capacity per unit VAI (J/K/m^2)
        self.gsnf = 0.01           # Snow-free vegetation moisture conductance (m/s)
        self.hbas = 2              # Canopy base height (m)
        self.kext = 0.5            # Vegetation light extinction coefficient
        self.leaf = 20             # Leaf boundary resistance (s/m)^(1/2)
        self.svai = 4.4            # Intercepted snow capacity per unit VAI (kg/m^2)
        self.tunl = 240 * 3600     # Canopy snow unloading time scale (s)
        self.wcan = 2.5            # Canopy wind decay coefficient

        # Snow parameters
        self.asmn = 0.5            # Minimum albedo for melting snow
        self.asmx = 0.85           # Maximum albedo for fresh snow
        self.eta0 = 3.7e7          # Reference snow viscosity (Pa s)
        self.hfsn = 0.1            # Snowcover fraction depth scale (m)
        self.kfix = 0.24           # Fixed thermal conductivity of snow (W/m/K)
        self.rcld = 300            # Maximum density for cold snow (kg/m^3)
        self.rfix = 300            # Fixed snow density (kg/m^3)
        self.rgr0 = 5e-5           # Fresh snow grain radius (m)
        self.rhof = 100            # Fresh snow density (kg/m^3)
        self.rhow = 300            # Wind-packed snow density (kg/m^3)
        self.rmlt = 500            # Maximum density for melting snow (kg/m^3)
        self.Salb = 10             # Snowfall to refresh albedo (kg/m^2)
        self.snda = 2.8e-6         # Thermal metamorphism parameter (1/s)
        self.Talb = -2             # Snow albedo decay temperature threshold (C)
        self.tcld = 3.6e6          # Cold snow albedo decay time scale (s)
        self.tmlt = 3.6e5          # Melting snow albedo decay time scale (s)
        self.trho = 200 * 3600     # Snow compaction timescale (s)
        self.Wirr = 0.03           # Irreducible liquid water content of snow
        self.z0sn = 0.001          # Snow roughness length (m)

        # Ground surface and soil parameters
        self.fcly = 0.3            # Soil clay fraction
        self.fsnd = 0.6            # Soil sand fraction
        self.gsat = 0.01           # Surface conductance for saturated soil (m/s)
        self.z0sf = 0.1            # Snow-free surface roughness length (m)

