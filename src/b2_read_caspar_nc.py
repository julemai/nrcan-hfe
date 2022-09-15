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
import numpy as np
import glob as glob
import datetime as datetime
from pathlib import Path
import xarray as xr

__all__ = ['read_caspar_nc']

def read_caspar_nc(variable=None,filenames=None,bbox=None,lintransform={'a':1.0,'b':0.0},silent=True):
    """
        Read NetCDF files retrieved from CaSPAr.

        Definition
        ----------
        def read_caspar_nc(product=None,variable=None,date=None,bbox=None,foldername='/tmp/test/',lintransform={'a':1.0,'b':0.0},silent=True)


        Input           Format          Description
        -----           -----           -----------
        variable        string          Name of variable in CaSPAr and in files, for example:
                                            "RDRS_v2.1_A_PR0_SFC" ... precipitation analysis
                                                                      in RDRS v2.1
                                        Will be used to check that files found in folder
                                        (foldername) are for this product.
                                        Default: None

        filenames       dict            Dictionary that provides the name of the file containing
                                        each of the requested time steps:
                                        filenames = { date_1: [ filename_1, filename_2 ],
                                                      date_2: [ filename_2 ],
                                                      date_3: [ ],
                                                      ... }
                                        --> date_1 is available in two files.
                                        --> date_2 is available in one file.
                                        --> date_3 is available in no file.
                                        This variable is the output of "request_caspar_nc()".
                                        Default: None

        bbox            dict            Dictionary specifying the bounding box, i.e.,
                                          {"lat":{"min":min_lat,"max":max_lat},
                                           "lon":{"min":min_lon,"max":max_lon}}
                                        If not specified, all data are returned.
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
        Reads NetCDF file(s) retrieved from CaSPAr and returns dictionary containing gridded
        variable ("var") and latitudes ("lat") and longitudes ("lon") of each grid cell as well
        as the time steps ("time").


        Restrictions
        ------------
        - Assumes that files are retrieved from CaSPAr ensuring a standard file format.
        - No handling of ensemble products (yet).
        - If several files are read, it is checked that all have consistent
          lat/lon fields. If this is not given, an error is returned.

        Examples
        --------

        >>> variable='RDRS_v2.1_A_PR0_SFC'
        >>> bbox={"lat":{"min":45.0,"max":46.0},"lon":{"min":-74.0,"max":-73.0}}
        >>> lintransform={'a':1000.0,'b':0.0} # convert from m/h to mm/h

        >>> # Read data (1 file)

        >>> filenames={datetime.datetime(2018, 8, 9, 7, 0): ['test-data/2018080812.nc']}

        >>> data_caspar = read_caspar_nc(variable=variable,filenames=filenames,bbox=bbox,lintransform=lintransform,silent=True)
        >>> print('shape(var) = '+str(np.shape(data_caspar["var"])))
        shape(var) = (1, 14, 12)
        >>> print('var[0,8,0:4] = '+str(data_caspar["var"][0,8,0:4]))
        var[0,8,0:4] = [0.12983558 0.0376055  0.         0.        ]
        >>> print('lat[8,0:4] = '+str(data_caspar["lat"][8,0:4]))
        lat[8,0:4] = [45.78684 45.76281 45.73859 45.71424]
        >>> print('lon[8,0:4] = '+str(data_caspar["lon"][8,0:4]))
        lon[8,0:4] = [-74.150696 -74.02869  -73.9068   -73.785   ]

        >>> # Read data (2 files)

        >>> filenames={datetime.datetime(2018, 8, 9, 7, 0): ['test-data/2018080812.nc'], datetime.datetime(2018, 8, 9, 14, 0): ['test-data/2018080912.nc']}

        >>> data_caspar = read_caspar_nc(variable=variable,filenames=filenames,bbox=bbox,lintransform=lintransform,silent=True)
        >>> print('shape(var) = '+str(np.shape(data_caspar["var"])))
        shape(var) = (2, 14, 12)
        >>> print('var[0,8,0:4] = '+str(data_caspar["var"][0,8,0:4]))
        var[0,8,0:4] = [0.12983558 0.0376055  0.         0.        ]
        >>> print('lat[8,0:4] = '+str(data_caspar["lat"][8,0:4]))
        lat[8,0:4] = [45.78684 45.76281 45.73859 45.71424]
        >>> print('lon[8,0:4] = '+str(data_caspar["lon"][8,0:4]))
        lon[8,0:4] = [-74.150696 -74.02869  -73.9068   -73.785   ]

        >>> # Read data (2 files; no file for one date available)

        >>> filenames={datetime.datetime(2018, 8, 9, 7, 0): ['test-data/2018080812.nc'], datetime.datetime(2018, 8, 9, 14, 0): ['test-data/2018080912.nc'], datetime.datetime(2018, 8, 10, 13, 0): []}

        >>> data_caspar = read_caspar_nc(variable=variable,filenames=filenames,bbox=bbox,lintransform=lintransform,silent=True)
        >>> print('shape(var) = '+str(np.shape(data_caspar["var"])))
        shape(var) = (2, 14, 12)
        >>> print('var[0,8,0:4] = '+str(data_caspar["var"][0,8,0:4]))
        var[0,8,0:4] = [0.12983558 0.0376055  0.         0.        ]
        >>> print('lat[8,0:4] = '+str(data_caspar["lat"][8,0:4]))
        lat[8,0:4] = [45.78684 45.76281 45.73859 45.71424]
        >>> print('lon[8,0:4] = '+str(data_caspar["lon"][8,0:4]))
        lon[8,0:4] = [-74.150696 -74.02869  -73.9068   -73.785   ]


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
    if variable is None:
        raise ValueError("read_caspar_nc: variable needs to be specified")
    if filenames is None:
        raise ValueError("read_caspar_nc: dictiornary of filenames per date (filenames) to be specified")

    # go through all dates and extract data
    result = {}

    # find files that contain this date
    for idate in filenames:

        if len(filenames[idate]) > 1:
            warnings.warn("read_caspar_nc: Date '{}' is available in several files. Only first file will be used.".format(idate))
            filename = filenames[idate][0]
        elif len(filenames[idate]) == 0:
            warnings.warn("read_caspar_nc: Date '{}' is not available in any file. Will be skipped.".format(idate))
            continue
        else:
            filename = filenames[idate][0]

        # check that required variables exist
        ds = xr.open_dataset(filename)
        if not('time' in list(ds.variables)):
            raise ValueError("read_caspar_nc: file {} does not contain variable 'time'.".format(filename))
        if not('lat' in list(ds.variables)):
            raise ValueError("read_caspar_nc: file {} does not contain variable 'lat'.".format(filename))
        if not('lon' in list(ds.variables)):
            raise ValueError("read_caspar_nc: file {} does not contain variable 'lon'.".format(filename))
        if not(variable in list(ds.variables)):
            raise ValueError("read_caspar_nc: file {} does not contain variable '{}'.".format(filename,variable))

        # get all available time steps of current file
        time_datetime64 = ds['time'].data   # type is "datetime64[ns]"
        time_datetime   = [ datetime.datetime.strptime(str(dd).split('.')[0],'%Y-%m-%dT%H:%M:%S') for dd in time_datetime64 ]   # type is "datetime.datetime"
        avail_dates = time_datetime

        # time index to read
        tidx = avail_dates.index(idate)

        # find grid cells within bounding box and read data
        lat = ds['lat']
        lon = ds['lon']
        nlat = np.shape(lat)[0]
        nlon = np.shape(lon)[1]
        if bbox is None:
            latidx = slice(0,nlat)
            lonidx = slice(0,nlon)
        else:
            # this is where the actual magic happens
            idx = np.where((lat >= bbox["lat"]["min"]) & (lat <= bbox["lat"]["max"]) & (lon >= bbox["lon"]["min"]) & (lon <= bbox["lon"]["max"]))
            # this is just returning a rectangle around the identified bounding box (same as Geomet)
            # -1 and +1 is to have one grid-cell buffer
            latidx = slice(np.min(idx[0])-1,np.max(idx[0])+2)
            lonidx = slice(np.min(idx[1])-1,np.max(idx[1])+2)
        slices = (tidx, latidx, lonidx)
        var = ds[variable][slices[0]][slices[1],slices[2]].data * lintransform['a'] + lintransform['b']
        lat = lat[slices[1],slices[2]].data
        lon = lon[slices[1],slices[2]].data
        ds.close()

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
                raise ValueError("read_caspar_nc: latitude field in file '{}' is not consistent with the previous files.".format(filename))

        # save longitude (or check that it is consistent)
        if not( "lon" in result.keys() ):
            result["lon"] = lon
        else:
            if np.any(result["lon"] != lon):
                raise ValueError("read_caspar_nc: longitude field in file '{}' is not consistent with the previous files.".format(filename))

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
