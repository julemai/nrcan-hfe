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

import pygrib as pg
import datetime as datetime
import urllib.request
from pathlib import Path
import warnings

__all__ = ['request_geomet_grib2']

def request_geomet_grib2(product=None,date=None,bbox=None,crs='EPSG:4326',filename='/tmp/test',overwrite=False,silent=True):
    """
        Request GRIB2 from Geomet.

        The URL to request RDPA 24-h accumulations for noon Aug 24, 2022 of a
        bounding box around Montreal (lat=45.508888, lon=-73.561668) from
        Geometdata is the following:
                https://api.weather.gc.ca/collections/weather:rdpa:10km:24f/
                coverage?f=GRIB&datetime=2022-08-24T12Z&CRS=EPSG:4326&bbox=-74,45,-73,46
        See documentation for GeoMet-OGC-API here:
                https://api.weather.gc.ca/openapi?f=html#/
        Specifically the documentation on how to retrieve RDPA 6h accumulations:
                https://api.weather.gc.ca/openapi?f=html#/weather%3Ardpa%3A10km%3A6f

        This code requests this file and stores it under a specified path and name.


        Definition
        ----------
        def request_geomet_grib2(product=None,date=None,bbox=None,crs='EPSG:4326',
        filename='/tmp/test',overwrite=False)


        Input           Format          Description
        -----           -----           -----------
        product         string          Name of product to retrieve currently available in
                                        GeoMet-OGC-API:
                                            "rdpa:10km:24f" ... RDPA 10km 24h accumulation
                                            "rdpa:10km:6f"  ... RDPA 10km  6h accumulation
                                        Default: None

        date            datetime or     Datetime object specifying time to request data for
                        list(datetime)  (i.e., end of accumulation period). The date will
                                        be used in the filename. Hence only one filename
                                        (i.e., basename) needs to be provided.
                                        Default: None

        bbox            dict            Dictionary specifying the bounding box, i.e.,
                                          {"lat":{"min":min_lat,"max":max_lat},
                                           "lon":{"min":min_lon,"max":max_lon}}
                                        Default: None

        crs             string          String specifying coordinate reference system (CRS).
                                        Default: "EPSG:4326"

        filename        string or       Filename of where to store retrieved file. The filename
                                        can be a path but is expected to not contain the file ending.
                                        The filename provided will be appended by the date in format
                                        YYYYMMDDHH.
                                        Default: "/tmp/test"

        overwrite       Boolean         If file already exists at specified location (filename) it
                                        will not be overwritten; unless overwrite is set to True.
                                        Default: False

        silent          Boolean         If set to True, nothing will be printed to terminal.
                                        Default: True

        Output          Format          Description
        -----           -----           -----------
        {var,lat,lon}   dict            Dictionary {"var":var,"lat":lat,"lon":lon}
                                        containing first variable of file and latitudes
                                        and longitudes of each grid cell


        Description
        -----------
        Requests GRIB2 file from Geomet for specified product ("product"), region ("bbox"), and time ("date").
        It returns a string specifying where the retrieved file has been stored.


        Restrictions
        ------------
        Assumes folder for output file exists. Only works for products that
        are supported by MSC GeoMet - GeoMet-OGC-API.

        Examples
        --------

        Request data (1 file)

        >>> date = datetime.datetime(2022,8,24,12,0)
        >>> bbox = {"lat":{"min":45.0,"max":46.0},"lon":{"min":-74.0,"max":-73.0}}
        >>> filename = 'test-data/test'
        >>> file_geomet = request_geomet_grib2(product='rdpa:10km:6f',date=date,bbox=bbox,crs='EPSG:4326',filename=filename,overwrite=True)
        >>> print('file_geomet = '+str(file_geomet))
        file_geomet = ['test-data/test_2022082412.grib2']

        Request data (4 files)

        >>> dates = [ datetime.datetime(2022,8,24,0,0) + datetime.timedelta(hours=6*ii) for ii in range(4) ]
        >>> bbox = {"lat":{"min":45.0,"max":46.0},"lon":{"min":-74.0,"max":-73.0}}
        >>> filename = 'test-data/test'
        >>> file_geomet = request_geomet_grib2(product='rdpa:10km:6f',date=dates,bbox=bbox,crs='EPSG:4326',filename=filename,overwrite=False)
        >>> print('file_geomet = ',file_geomet)
        file_geomet = ['test-data/test_2022082400.grib2', 'test-data/test_2022082406.grib2', 'test-data/test_2022082412.grib2', 'test-data/test_2022082418.grib2']


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
        raise ValueError("request_geomet_grib2: product needs to be specified")
    if date is None:
        raise ValueError("request_geomet_grib2: date needs to be specified")
    if bbox is None:
        raise ValueError("request_geomet_grib2: bounding box needs to be specified")

    # make sure date input is always list
    if type(date) == datetime.datetime:
        date = [ date ]

    # initialize return
    filenames = []

    for idate in date:

        # filename including file extension and date
        ifilename = filename+"_"+idate.strftime("%Y%m%d%H")+".grib2"

        # check if product is available at requested HH
        hour = int(idate.strftime("%H"))
        if product == "rdpa:10km:24f":
            if not( hour in [12]):
                raise ValueError("request_geomet_grib2: Analysis for product '{}' not available at hour '{}'.".format(product, hour))
        elif  product == "rdpa:10km:6f":
            if not( hour in [0,6,12,18]):
                raise ValueError("request_geomet_grib2: Analysis for product '{}' not available at hour '{}'.".format(product, hour))
        else:
            raise ValueError("request_geomet_grib2: Product '{}' not known. It is either not available in Geomet or you need to add it here with the available issue hours.".format(product))

        # make sure folder to store file exists; otherwise create
        Path(ifilename).parent.mkdir(parents=True, exist_ok=True)

        # create request string
        if overwrite or not( Path(ifilename).is_file() ):
            api_str = "https://api.weather.gc.ca/collections/weather:{}/coverage?f=GRIB&datetime={}&CRS={}&bbox={},{},{},{}".format(
                product,
                idate.strftime("%Y-%m-%dT%HZ"),   # 2022-08-24T12Z
                crs,
                bbox["lon"]["min"],
                bbox["lat"]["min"],
                bbox["lon"]["max"],
                bbox["lat"]["max"])
            if not(silent): print("Downloading: {}".format(api_str))
            urllib.request.urlretrieve(api_str, ifilename)
        else:
            warnings.warn("File '{}' already exists. Will not be downloaded again.".format(ifilename))

        # append file to return
        filenames.append(ifilename)

    return filenames



if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)


    # date = datetime.datetime(2022,8,24,12,0)
    # bbox = {"lat":{"min":45.0,"max":46.0},"lon":{"min":-74.0,"max":-73.0}}
    # filename = '/tmp/tmf/test'
    # file_geomet = request_geomet_grib2(product='rdpa:10km:6f',date=date,bbox=bbox,crs='EPSG:4326',filename=filename)
