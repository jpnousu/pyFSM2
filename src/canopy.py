#-----------------------------------------------------------------------
# Properties of vegetation canopy layers
#-----------------------------------------------------------------------

import numpy as np
from pyFSM2_MODULES import Constants, Layers, Parameters

class Canopy:
    def __init__(self):

        constants = Constants()
        layers = Layers()
        params = Parameters(SETPAR=2, DENSITY=1)

        # from Constants
        self.hcap_ice = constants.hcap_ice # Specific heat capacity of ice (J/K/kg)

        # from Layers
        self.Ncnpy = layers.Ncnpy
        self.fvg1 = layers.fvg1

        # from Parameters
        self.cvai = params.cvai # Vegetation heat capacity per unit VAI (J/K/m^2)
        self.svai = params.svai # Intercepted snow capacity per unit VAI (kg/m^2)

        # Model state variables
        self.cveg = np.zeros(self.Ncnpy) # Vegetation heat capacities (J/K/m^2)
        self.fcans = np.zeros(self.Ncnpy) # Canopy layer snowcover fractions
        self.lveg = np.zeros(self.Ncnpy) # Canopy layer vegetation area indices
        self.Scap = np.zeros(self.Ncnpy) # Canopy layer snow capacities (kg/m^2)
        self.Tveg0 = np.zeros(self.Ncnpy) # Vegetation temperatures at start of timestep (K)

        self.eps = np.finfo(float).eps

  def run_timestep(Sveg, Tveg, VAI):
    '''
    '''

    if CANMOD == 1:
      self.lveg[0] = VAI

    if CANMOD == 2:
       self.lveg[0] = self.fvg1 * VAI
       self.lveg[1] = (1 - self.fvg1) * VAI

    self.cveg[:] = self.cvai * self.lveg[:] + self.hcap_ice * Sveg[:]
    self.fcans[:] = 0
    self.Scap[:] = 0
    if (VAI > 0):
       self.Scap[:] = self.svai * self.lveg[:]
       if (self.svai > 0):
          self.fcans[:] = (Sveg[:]/Scap[:])**0.67

    self.fcans[self.fcans > 1] = 1

    return cveg, fcans, lveg, Scap, Tveg0