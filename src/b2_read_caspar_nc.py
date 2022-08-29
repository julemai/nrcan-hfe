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

import numpy as np
import glob as glob
import datetime as datetime
from pathlib import Path
import xarray as xr

__all__ = ['read_caspar_nc']

def read_caspar_nc(product=None,variable=None,date=None,bbox=None,foldername='/tmp/test/',lintransform={'a':1.0,'b':0.0},silent=True):
    """
        Read NetCDF files retrieved from CaSPAr.

        Definition
        ----------
        def read_caspar_nc(product=None,variable=None,date=None,bbox=None,foldername='/tmp/test/',lintransform={'a':1.0,'b':0.0},silent=True)


        Input           Format          Description
        -----           -----           -----------
        product         string          Name of product in CaSPAr, for example:
                                            "RDRS_v2.1" ... RDRS v2.1
                                        Will be used to check that files found in folder
                                        (foldername) are for this product.
                                        Default: None

        variable        string          Name of variable in CaSPAr and in files, for example:
                                            "RDRS_v2.1_A_PR0_SFC" ... precipitation analysis
                                                                      in RDRS v2.1
                                        Will be used to check that files found in folder
                                        (foldername) are for this product.
                                        Default: None

        date            datetime or     Datetime object specifying time to request data for.
                        list(datetime)  The CaSPAr files are expected to contain part of the
                                        date (YYYYMMDD). So, do not rename them when you received
                                        them from CaSPAr. A list of dates can be read at the
                                        same time. The returned variable (var) will then be 3D
                                        with the first dimension being the time dimension.
                                        Default: None

        bbox            dict            Dictionary specifying the bounding box, i.e.,
                                          {"lat":{"min":min_lat,"max":max_lat},
                                           "lon":{"min":min_lon,"max":max_lon}}
                                        If not specified all data are returned.
                                        Default: None

        foldername      string          Name of the folder containing the files requested from
                                        CaSPAr. Filenames are expected to have the following
                                        pattern YYYYMMDDHH.nc which is the filename pattern of
                                        CaSPAr files. So, do not rename them. At the moment,
                                        it is not possible to read data from ensemble products
                                        with file pattern YYYYMMDDHH_EEE.nc.
                                        Default: '/tmp/test/'

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
        {var,lat,lon}   dict            Dictionary {"var":var,"lat":lat,"lon":lon}
                                        containing first variable of file and latitudes
                                        and longitudes of each grid cell. If several files
                                        are read, "var" is returned as list of 2-dimensional
                                        fields while "lat" and "lon" are expected to be
                                        consistent across the files and will be only returned
                                        as 2-dimensional fields.


        Description
        -----------
        Reads NetCDF file(s) retrieved from CaSPAr and returns dictionary containing gridded
        variable ("var") and latitudes ("lat") and longitudes ("lon") of each grid cell.


        Restrictions
        ------------
        - Assumes that files are retrieved from CaSPAr ensuring a standard file format.
        - No handling of ensemble products (yet).
        - If several files are read, it is checked that all have consistent
          lat/lon fields. If this is not given, an error is returned.

        Examples
        --------

        >>> product='RDRS_v2.1'
        >>> variable='RDRS_v2.1_A_PR0_SFC'
        >>> bbox={"lat":{"min":45.0,"max":46.0},"lon":{"min":-74.0,"max":-73.0}}
        >>> foldername='test-data/'
        >>> lintransform={'a':1000.0,'b':0.0} # convert from m/h to mm/h

        Read data (1 file provided as string --> returns var as 2D array)

        >>> date=datetime.datetime(1980,1,14,18,0)

        >>> data_caspar = read_caspar_nc(product=product,variable=variable,date=date,bbox=bbox,foldername=foldername,lintransform=lintransform,silent=True)
        >>> print('var[8,0:4] = '+str(data_caspar["var"][8,0:4]))
        var[8,0:4] = [0.53641945 0.46970782 0.48578462 0.50169617]
        >>> print('lat[8,0:4] = '+str(data_caspar["lat"][8,0:4]))
        lat[8,0:4] = [45.84939  45.825138 45.800735 45.776184]
        >>> print('lon[8,0:4] = '+str(data_caspar["lon"][8,0:4]))
        lon[8,0:4] = [-73.99344 -73.87134 -73.74933 -73.62747]

        Read data (1 file provided as list --> returns var as 3D array)

        >>> date=[datetime.datetime(1980,1,14,18,0)]

        >>> data_caspar = read_caspar_nc(product=product,variable=variable,date=date,bbox=bbox,foldername=foldername,lintransform=lintransform,silent=True)
        >>> print('var[0,8,0:4] = '+str(data_caspar["var"][0,8,0:4]))
        var[0,8,0:4] = [0.53641945 0.46970782 0.48578462 0.50169617]
        >>> print('lat[8,0:4] = '+str(data_caspar["lat"][8,0:4]))
        lat[8,0:4] = [45.84939  45.825138 45.800735 45.776184]
        >>> print('lon[8,0:4] = '+str(data_caspar["lon"][8,0:4]))
        lon[8,0:4] = [-73.99344 -73.87134 -73.74933 -73.62747]

        Read data (2 files --> returns var as 3D array)

        >>> date=[datetime.datetime(1980,1,14,18,0),datetime.datetime(1980,1,14,20,0)]

        >>> data_caspar = read_caspar_nc(product=product,variable=variable,date=date,bbox=bbox,foldername=foldername,lintransform=lintransform,silent=True)
        >>> print('var[0,8,0:4] = '+str(data_caspar["var"][0,8,0:4]))
        var[0,8,0:4] = [0.53641945 0.46970782 0.48578462 0.50169617]
        >>> print('lat[8,0:4] = '+str(data_caspar["lat"][8,0:4]))
        lat[8,0:4] = [45.84939  45.825138 45.800735 45.776184]
        >>> print('lon[8,0:4] = '+str(data_caspar["lon"][8,0:4]))
        lon[8,0:4] = [-73.99344 -73.87134 -73.74933 -73.62747]


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
    if product is None:
        raise ValueError("read_caspar_nc: product need to be specified")
    if variable is None:
        raise ValueError("read_caspar_nc: variable need to be specified")
    if date is None:
        raise ValueError("read_caspar_nc: date(s) need(s) to be specified")

    # make sure date input is always list
    date_was_string = False
    if type(date) == datetime.datetime:
        date = [ date ]
        date_was_string = True
    date = np.sort(date)

    # gather all time steps available in all files
    filenames = np.sort( glob.glob(foldername+'*.nc') )
    avail_dates = {}
    for filename in filenames:

        ds = xr.open_dataset(filename)
        time_datetime64 = ds['time'].data   # type is "datetime64[ns]"
        time_datetime   = [ datetime.datetime.strptime(str(dd).split('.')[0],'%Y-%m-%dT%H:%M:%S') for dd in time_datetime64 ]   # type is "datetime.datetime"

        avail_dates[filename] = time_datetime
        ds.close()

    # go through all dates and extract data
    result = {}
    for idate in date:

        found_date = False
        if not(silent): print("Search for date '{}' ...".format(idate))

        # find files that contain this date
        for filename in filenames:

            # if date exists in multiple files, it will be overwritten until last time found
            # does not matter for RDRS because each time step exists only once
            if idate in avail_dates[filename]:

                found_date = True
                if not(silent): print(">>> Found in file '{}' ...".format(filename))

                # time index to read
                tidx = avail_dates[filename].index(idate)

                # find grid cells within bounding box
                ds = xr.open_dataset(filename)
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

        if not(found_date):
            raise ValueError("Date '{}' does not seem to exist in any NetCDF file in folder '{}'!".format(idate,foldername))

    # make sure a 2D array is returned if only one filename (not in list) was requested
    if date_was_string:
        result["var"] = np.array(result["var"][0])

    # return dictionary
    return result


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # product='RDRS_v2.1'
    # variable='RDRS_v2.1_A_PR0_SFC'
    # date=datetime.datetime(1980,1,13,16,0)
    # bbox=None
    # bbox={"lat":{"min":45.0,"max":46.0},"lon":{"min":-74.0,"max":-73.0}}
    # foldername='test-data/'
    # lintransform={'a':1.0,'b':0.0}
    # silent=False

    # # checking inputs
    # if product is None:
    #     raise ValueError("read_caspar_nc: product need to be specified")
    # if variable is None:
    #     raise ValueError("read_caspar_nc: variable need to be specified")
    # if date is None:
    #     raise ValueError("read_caspar_nc: date(s) need(s) to be specified")

    # # make sure date input is always list
    # date_was_string = False
    # if type(date) == datetime.datetime:
    #     date = [ date ]
    #     date_was_string = True
    # date = np.sort(date)

    # # gather all time steps available in all files
    # filenames = np.sort( glob.glob(foldername+'*.nc') )
    # avail_dates = {}
    # for filename in filenames:

    #     ds = xr.open_dataset(filename)
    #     time_datetime64 = ds['time'].data   # type is "datetime64[ns]"
    #     time_datetime   = [ datetime.datetime.strptime(str(dd).split('.')[0],'%Y-%m-%dT%H:%M:%S') for dd in time_datetime64 ]   # type is "datetime.datetime"

    #     avail_dates[filename] = time_datetime
    #     ds.close()

    # # go through all dates and extract data
    # result = {}
    # for idate in date:

    #     # find files that contain this date
    #     for filename in filenames:

    #         # if date exists in multiple files, it will be overwritten until last time found
    #         # does not matter for RDRS because each time step exists only once
    #         if idate in avail_dates[filename]:

    #             # time index to read
    #             tidx = avail_dates[filename].index(idate)

    #             # find grid cells within bounding box
    #             ds = xr.open_dataset(filename)
    #             lat = ds['lat']
    #             lon = ds['lon']
    #             nlat = np.shape(lat)[0]
    #             nlon = np.shape(lon)[1]
    #             if bbox is None:
    #                 latidx = slice(0,nlat)
    #                 lonidx = slice(0,nlon)
    #             else:
    #                 # this is where the actual magic happens
    #                 idx = np.where((lat >= bbox["lat"]["min"]) & (lat <= bbox["lat"]["max"]) & (lon >= bbox["lon"]["min"]) & (lon <= bbox["lon"]["max"]))
    #                 latidx = slice(np.min(idx[0]),np.max(idx[0])+1)  # this is just returning a rectangle around the identified bounding box (same as Geomet)
    #                 lonidx = slice(np.min(idx[1]),np.max(idx[1])+1)  # this is just returning a rectangle around the identified bounding box (same as Geomet)
    #             slices = (tidx, latidx, lonidx)
    #             var = ds[variable][slices[0]][slices[1],slices[2]].data
    #             lat = lat[slices[1],slices[2]].data
    #             lon = lon[slices[1],slices[2]].data
    #             ds.close()

    #             # save latitude (or check that it is consistent)
    #             if not( "lat" in result.keys() ):
    #                 result["lat"] = lat
    #             else:
    #                 if np.any(result["lat"] != lat):
    #                     raise ValueError("read_caspar_nc: latitude field in file '{}' is not consistent with the previous files.".format(filename))

    #             # save longitude (or check that it is consistent)
    #             if not( "lon" in result.keys() ):
    #                 result["lon"] = lon
    #             else:
    #                 if np.any(result["lon"] != lon):
    #                     raise ValueError("read_caspar_nc: longitude field in file '{}' is not consistent with the previous files.".format(filename))

    #             # save variable
    #             if "var" in result.keys():
    #                 # append to existing var
    #                 result["var"] = np.append(result["var"], [var],axis=0)
    #             else:
    #                 # initialize var
    #                 result["var"] = np.array([ var ])

    # # make sure a 2D array is returned if only one filename (not in list) was requested
    # if date_was_string:
    #     result["var"] = np.array(result["var"][0])
