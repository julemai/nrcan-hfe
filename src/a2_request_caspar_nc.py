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
import json as json

__all__ = ['request_caspar_nc']

def request_caspar_nc(product=None,variable=None,date=None,foldername='/tmp/test/',silent=True):
    """
        Request NetCDF from CaSPAr.

        This is a dummy function. Requests of data from CaSPAr have to be made manually.
        A description is provided here:
        https://github.com/julemai/nrcan-hfe/blob/main/data/caspar/README.md

        Or directly under CaSPAr:
        https://github.com/julemai/CaSPAr/wiki/How-to-get-started-and-download-your-first-data

        This function just checks if the files for the product specified and the (list of)
        date(s) provided are available under the specified foldername.


        Definition
        ----------
        def request_caspar_nc(product=None,variable=None,date=None,foldername='/tmp/test/',silent=True)


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

        foldername      string          Name of the folder containing the files requested from
                                        CaSPAr. Filenames are expected to have the following
                                        pattern YYYYMMDDHH.nc which is the filename pattern of
                                        CaSPAr files. So, do not rename them. At the moment,
                                        it is not possible to read data from ensemble products
                                        with file pattern YYYYMMDDHH_EEE.nc.
                                        Default: '/tmp/test/'

        silent          Boolean         If set to True, nothing will be printed to terminal.
                                        Default: True

        Output          Format          Description
        -----           -----           -----------
        filenames       dict            Dictionary that provides the name of the file containing
                                        each of the requested time setps:
                                        filenames = { date_1: [ filename_1, filename_2 ],
                                                      date_2: [ filename_2 ],
                                                      date_3: [ ],
                                                      ... }
                                        date_1 is available in two files.
                                        date_2 is available in one file.
                                        date_3 is available in no file.


        Description
        -----------
        A dummy function to provide some information on how to request data from CaSPAr. The function
        is actually only checking whether data are already available.


        Restrictions
        ------------
        Requests need to be made manually.


        Examples
        --------

        Request data

        >>> product='RDRS_v2.1'
        >>> variable='RDRS_v2.1_A_PR0_SFC'
        >>> date=[ datetime.datetime(2018,8,9,7,0), datetime.datetime(2018,8,9,8,0), datetime.datetime(2018,8,10,7,0), datetime.datetime(2018,8,10,13,0) ]
        >>> foldername='test-data/'

        >>> file_caspar = request_caspar_nc(product=product,variable=variable,date=date,foldername=foldername,silent=False)
        >>> print('file_caspar = '+str(file_caspar))
        file_caspar = {datetime.datetime(2018, 8, 9, 7, 0): ['test-data/2018080812.nc'], datetime.datetime(2018, 8, 9, 8, 0): ['test-data/2018080812.nc'], datetime.datetime(2018, 8, 10, 7, 0): ['test-data/2018080912.nc'], datetime.datetime(2018, 8, 10, 13, 0): []}


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

    warnings.warn("request_caspar_nc: CaSPAr data need to be requested manually. Find information under 'https://github.com/julemai/nrcan-hfe/blob/main/data/caspar/README.md'.")


    # checking inputs
    if product is None:
        raise ValueError("request_caspar_nc: product needs to be specified")
    if variable is None:
        raise ValueError("request_caspar_nc: variable needs to be specified")
    if date is None:
        raise ValueError("request_caspar_nc: date(s) need(s) to be specified")

    # initialize return
    filenames = {}

    # make sure date input is always list
    date_was_string = False
    if type(date) == datetime.datetime:
        date = [ date ]
        date_was_string = True
    date = np.sort(date)

    # check content of NetCDFs in specified folder
    # --> since this can take ages when there is lots of files in this folder, this is done once and results saved in a JSON
    # --> if this JSON is found only read results

    # all NetCDF files in specified folder
    ncfiles = np.sort( glob.glob(foldername+'/*.nc') )

    # check-file with filenames and available time steps in each file (only created when not existing --> saving lots of time)
    jsonfile = foldername+'/check_summary.json'

    if not( Path(jsonfile).is_file() ):

        # make sure all files contain specified product
        for ncfile in ncfiles:

            ds = xr.open_dataset(ncfile)
            iproduct = ds.attrs['product']
            if iproduct != product:
                raise ValueError("request_caspar_nc: file {} does not contain specified product {}.".format(ncfile,product))

            ds.close()

        # gather all time steps available in all files
        avail_dates = {}
        for ncfile in ncfiles:

            ds = xr.open_dataset(ncfile)
            time_datetime64 = ds['time'].data   # type is "datetime64[ns]"
            time_datetime   = [ datetime.datetime.strptime(str(dd).split('.')[0],'%Y-%m-%dT%H:%M:%S') for dd in time_datetime64 ]   # type is "datetime.datetime"

            avail_dates[ncfile] = time_datetime
            ds.close()

        # save in JSON (no datetime objects)
        avail_dates_JSON = { ii:[str(iii) for iii in avail_dates[ii]] for ii in avail_dates }
        json_dump = json.dumps(avail_dates_JSON)
        ff = open(jsonfile, "w")
        ff.write(json_dump)
        ff.close()

    with open(jsonfile, 'r') as ff:
        avail_dates_JSON = json.load(ff)

    # make str into datetime again
    avail_dates = { ii:[datetime.datetime.strptime(iii,'%Y-%m-%d %H:%M:%S') for iii in avail_dates_JSON[ii]] for ii in avail_dates_JSON }

    # checking if information for all files in folder exists
    # if not, the check-file might be outdated
    all_files_in_checkfile = list(avail_dates.keys())
    if not( np.all([ ifile in ncfiles for ifile in all_files_in_checkfile ]) ):
        raise ValueError("request_caspar_nc: files in check-file '{}' is not consistent with files found in folder '{}'.\nYou should update the check-file by deleting '{}' and re-running this script. This way the file will be recreated. Be aware that this might take a few minutes if there are lots of NetCDF files located in folder '{}'.".format(jsonfile,foldername,jsonfile,foldername))


    # return only filenames for requested dates
    filenames = { idate : [ avail_date for avail_date in avail_dates if idate in avail_dates[avail_date] ] for idate in date }

    # check if all lists of filenames are empty
    if np.sum([ len(filenames[idate]) for idate in filenames ]) == 0:
        raise ValueError("request_caspar_nc: No file was found. The foldername '{}' is probably wrong or data have not been manually requested from CaSPAr yet.".format(foldername))

    return filenames



if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # product='RDRS_v2.1'
    # variable='RDRS_v2.1_A_PR0_SFC'
    # date=[ datetime.datetime(2018,8,9,7,0), datetime.datetime(2018,8,9,8,0), datetime.datetime(2018,8,10,7,0), datetime.datetime(2018,8,10,13,0) ]
    # foldername='test-data/'

    # file_caspar = request_caspar_nc(product=product,variable=variable,date=date,foldername=foldername,silent=False)



    # # checking inputs
    # if product is None:
    #     raise ValueError("request_caspar_nc: product needs to be specified")
    # if variable is None:
    #     raise ValueError("request_caspar_nc: variable needs to be specified")
    # if date is None:
    #     raise ValueError("request_caspar_nc: date(s) need(s) to be specified")


    # # make sure date input is always list
    # date_was_string = False
    # if type(date) == datetime.datetime:
    #     date = [ date ]
    #     date_was_string = True
    # date = np.sort(date)

    # # make sure all files contain specified product
    # filenames = np.sort( glob.glob(foldername+'*.nc') )
    # for filename in filenames:

    #     ds = xr.open_dataset(filename)
    #     iproduct = ds.attrs['product']
    #     if iproduct != product:
    #         raise ValueError("request_caspar_nc: file {} does not contain specified product {}.".format(filename,product))

    #     ds.close()

    # # gather all time steps available in all files
    # filenames = np.sort( glob.glob(foldername+'*.nc') )
    # avail_dates = {}
    # for filename in filenames:

    #     ds = xr.open_dataset(filename)
    #     time_datetime64 = ds['time'].data   # type is "datetime64[ns]"
    #     time_datetime   = [ datetime.datetime.strptime(str(dd).split('.')[0],'%Y-%m-%dT%H:%M:%S') for dd in time_datetime64 ]   # type is "datetime.datetime"

    #     avail_dates[filename] = time_datetime
    #     ds.close()

    # # go through all dates and collect file containing required time step and variable
    # filenames = { idate : [ avail_date for avail_date in avail_dates if idate in avail_dates[avail_date] ] for idate in date }
