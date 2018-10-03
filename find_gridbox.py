''' Function to find the gridbox closest to the inputted coordinates.

======================================================================================================================

    Inputs:

    - 1-D arrays of 'real' (i.e. already rotated if the data is modelled on a rotated pole grid) latitudes and longitudes

    - coordinates of the lats and lons you want to find

    Outputs:

    - indices for the inputted coordinates in your model domain.

'''

import numpy as np

def find_gridbox(x, y, real_lat, real_lon):
    global lon_index, lat_index
    lat_index = np.argmin((real_lat - x) ** 2)  # take whole array and subtract lat you want from
    lon_index = np.argmin((real_lon - y) ** 2)  # each point, then find the smallest difference
    return lon_index, lat_index