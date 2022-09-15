#!/usr/bin/env python
from __future__ import print_function

# Copyright 2022 Juliane Mai - juliane.mai(at)uwaterloo.ca
#
# License
# This file is part of the HFE code library for NRCan for "Improving the
# characterization of the flood-producing precipitation events in the
# Historic Flood Event (HFE) database by using the CaPA reanalysis/analysis".
#
# The NRCan-HFE code library is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 2.1 of the License, or
# (at your option) any later version.
#
# The NRCan-HFE code library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with The NRCan-HFE code library.
# If not, see <https://github.com/julemai/nrcan-hfe/blob/main/LICENSE>.

import warnings
import datetime as datetime
import pygrib as pg
from pathlib import Path
import numpy as np

__all__ = ['read_geomet_grib2']

def read_geomet_grib2(filenames=None,lintransform={'a':1.0,'b':0.0},silent=True):
    """
        Read GRIB2 retrieved from Geomet.

        Definition
        ----------
        def read_geomet_grib2(filename=None,lintransform={'a':1.0,'b':0.0},silent=True)


        Input           Format          Description
        -----           -----           -----------
        filenames       dict            Dictionary that provides the name of the file containing
                                        each of the requested time steps:
                                        filenames = { date_1: [ filename_1, filename_2 ],
                                                      date_2: [ filename_2 ],
                                                      date_3: [ ],
                                                      ... }
                                        --> date_1 is available in two files.
                                        --> date_2 is available in one file.
                                        --> date_3 is available in no file.
                                        This variable is the output of "request_geomet_grib2()".
                                        Default: None

        lintransform    dict            Dictionary containing slope (a) and intercept (b) of
                                        linear transform applied to data read for variable:
                                        var_new = var_old * a + b
                                        Can be used to convert units, e.g.
                                        mm/h --> mm/d  :: {'a'=24.,'b'=0.0}
                                        m/d  --> mm/h  :: {'a'=1000./24.,'b'=0.0}
                                        ft/d --> mm/d  :: {'a'=304.8, 'b'=0.0}
                                        Default: {'a'=1.0, 'b'=0.0}

        silent          Boolean         If set to True, nothing will be printed to terminal.
                                        Default: True

        Output          Format          Description
        -----           -----           -----------
        {var,time,      dict            Dictionary {"var":var,"time":time,"lat":lat,"lon":lon}
        lat,lon}                        containing first variable of file and latitudes
                                        and longitudes of each grid cell. If several files
                                        are read, "var" is returned as list of 2-dimensional
                                        fields while "time", "lat", and "lon" are expected to be
                                        consistent across the files and will be only returned
                                        as 2-dimensional fields.


        Description
        -----------
        Reads a GRIB2 file(s) and returns dictionary containing gridded variable ("var")
        and latitudes ("lat") and longitudes ("lon") of each grid cell as well
        as the time steps ("time").


        Restrictions
        ------------
        - Assumes that each file contains one variable and one time step.
        - There is no checking which variable is read and returned.
        - If several files are read, it is checked that all have consistent
          lat/lon fields. If this is not given, an error is returned.

        Examples
        --------

        >>> # Read data (1 file)

        >>> filenames = {datetime.datetime(2022, 8, 24, 12, 0): ['test-data/rdpa-6h_2022082412.grib2']}
        >>> data_geomet = read_geomet_grib2(filenames)
        >>> print('var[0,0,0:4] = '+str(data_geomet["var"][0,0,0:4]))
        var[0,0,0:4] = [0.09765625 0.0625     0.09375    0.1015625 ]
        >>> print('lat[0,0:4] = '+str(data_geomet["lat"][0,0:4]))
        lat[0,0:4] = [45.359908   45.31104666 45.2620931  45.21304786]
        >>> print('lon[0,0:4] = '+str(data_geomet["lon"][0,0:4]))
        lon[0,0:4] = [-74.714051   -74.61954144 -74.52526103 -74.43120944]

        >>> # Read data (4 files)

        >>> filenames = {datetime.datetime(2022, 8, 24, 0, 0): ['test-data/rdpa-6h_2022082400.grib2'], datetime.datetime(2022, 8, 24, 6, 0): ['test-data/rdpa-6h_2022082406.grib2'], datetime.datetime(2022, 8, 24, 12, 0): ['test-data/rdpa-6h_2022082412.grib2'], datetime.datetime(2022, 8, 24, 18, 0): ['test-data/rdpa-6h_2022082418.grib2']}
        >>> data_geomet = read_geomet_grib2(filenames)
        >>> print('var[0,0,0:4] = '+str(data_geomet["var"][0,0,0:4]))
        var[0,0,0:4] = [0.1875  0.125   0.34375 0.21875]
        >>> print('lat[0,0:4] = '+str(data_geomet["lat"][0,0:4]))
        lat[0,0:4] = [45.359908   45.31104666 45.2620931  45.21304786]
        >>> print('lon[0,0:4] = '+str(data_geomet["lon"][0,0:4]))
        lon[0,0:4] = [-74.714051   -74.61954144 -74.52526103 -74.43120944]


        License
        -------
        This file is part of the HFE code library for NRCan for "Improving the
        characterization of the flood-producing precipitation events in the
        Historic Flood Event (HFE) database by using the CaPA reanalysis/analysis".

        The NRCan-HFE code library is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 2.1 of the License, or
        (at your option) any later version.

        The NRCan-HFE code library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with The NRCan-HFE code library.
        If not, see <https://github.com/julemai/nrcan-hfe/blob/main/LICENSE>.

        Copyright 2022 Juliane Mai - juliane.mai@uwaterloo.ca


        History
        -------
        Written,  Juliane Mai, August 2022
    """

    # checking inputs
    if filenames is None:
        raise ValueError("read_geomet_grib2: filename needs to be specified")

    result = {}
    for idate in filenames:

        if len(filenames[idate]) > 1:
            warnings.warn("read_geomet_grib2: Date '{}' is available in several files. Only first file will be used.".format(idate))
            filename = filenames[idate][0]
        elif len(filenames[idate]) == 0:
            warnings.warn("read_geomet_grib2: Date '{}' is not available in any file. Will be skipped.".format(idate))
            continue
        else:
            filename = filenames[idate][0]

        # reading first variable of file
        if Path(filename).is_file():
            ff = pg.open(filename)
            vv = ff[1]  # numbering starts with 1 in GRIB2
            (var,lat,lon) = vv.data()
            ff.close()
        else:
            raise ValueError("read_geomet_grib2: File '{}' not found.".format(filename))

        # unit conversion
        var = var * lintransform['a'] + lintransform['b']

        # save time
        if "time" in result.keys():
            # append to existing time
            result["time"] = np.append(result["time"], idate )
        else:
            # initialize time
            result["time"] = np.array([ idate ])

        # save latitude (or check that it is consistent)
        if not( "lat" in result.keys() ):
            result["lat"] = lat
        else:
            if np.any(result["lat"] != lat):
                raise ValueError("read_geomet_grib2: latitude field in file '{}' is not consistent with the previous files.".format(filename))

        # save longitude (or check that it is consistent)
        if not( "lon" in result.keys() ):
            result["lon"] = lon
        else:
            if np.any(result["lon"] != lon):
                raise ValueError("read_geomet_grib2: longitude field in file '{}' is not consistent with the previous files.".format(filename))

        # save variable
        if "var" in result.keys():
            # append to existing var
            result["var"] = np.append(result["var"], [var],axis=0)
        else:
            # initialize var
            result["var"] = np.array([ var ])

    # return dictionary
    return result



if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # filenames = {datetime.datetime(2022, 8, 24, 0, 0): ['test-data/rdpa-6h_2022082400.grib2'], datetime.datetime(2022, 8, 24, 6, 0): ['test-data/rdpa-6h_2022082406.grib2'], datetime.datetime(2022, 8, 24, 12, 0): ['test-data/rdpa-6h_2022082412.grib2'], datetime.datetime(2022, 8, 24, 18, 0): ['test-data/rdpa-6h_2022082418.grib2']}
