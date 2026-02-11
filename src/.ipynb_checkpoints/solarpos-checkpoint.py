#-----------------------------------------------------------------------
# Azimuth and elevation angles of the sun
#-----------------------------------------------------------------------

from pyFSM2_MODULES import Constants
import numpy as np
from math import sin, asin, cos, acos

class SolarPos():
    def __init__(self):
         
        constants = Constants()

        self.pi = constants.pi  # pi

    def run_timestep(self, year, month, day, hour, lat, noon):

        DoY = (7*year)/4 - 7*(year+(month+9)/12)/4 + (275*month)/9 + day - 30

        dangle = 2*self.pi*(DoY - 1)/365

        declin = (0.006918 - 0.399912*cos(dangle)   + 0.070257*sin(dangle)
                        - 0.006758*cos(2*dangle) + 0.000907*sin(2*dangle)
                        - 0.002697*cos(3*dangle) + 0.001480*sin(3*dangle))
        eqtime = ((0.000075 + 0.001868*cos(dangle)   - 0.032077*sin(dangle)
                        - 0.014615*cos(2*dangle) - 0.04089*sin(2*dangle))
                *(12/self.pi))
        
        hangle = (self.pi/12)*(noon - hour - eqtime)

        elev = asin(sin(declin)*sin(lat) + cos(declin)*cos(lat)*cos(hangle))

        azim = acos((sin(elev)*sin(lat) - sin(declin))/(cos(elev)*cos(lat)))

        if (hangle < 0):
               azim = - azim

        return azim, elev