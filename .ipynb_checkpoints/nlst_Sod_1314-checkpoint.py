"""
PARAMETERS
@author: jpnousu
"""

def params(folder=''):
    gridpnts = {'Ncols': 2}
    gridlevs = {}
    drive = {'met_file': 'met_Sod_1314.txt',
             'zT': 18,
             'zU': 18,
             'lat': 67.37,
             'noon': 10
            }
    veg = {'vegh': [0, 15],
           'VAI': [0, 2],
          }
    initial = {}
    outputs = {'runid': 'Sod_1314_'}