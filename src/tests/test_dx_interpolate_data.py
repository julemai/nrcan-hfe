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


# -----------------------
# add subolder scripts/lib to search path
# -----------------------
import sys
import os
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path+'/../../src')


from dx_interpolate_data import interpolate_data
import numpy as np
import datetime as datetime
import pytest


@pytest.mark.filterwarnings("ignore:request_geomet_grib2")
def test_interpolate_geomet():

    # --------------------------------------
    # Request and read data from Geomet
    # --------------------------------------

    from a1_request_geomet_grib2 import request_geomet_grib2
    from b1_read_geomet_grib2 import read_geomet_grib2

    product = 'rdpa:10km:24f'
    crs = 'EPSG:4326'
    bbox = {"lat":{"min":45.0,"max":46.0},"lon":{"min":-74.0,"max":-73.0}}

    # four dates 6h apart from each other
    dates = [ datetime.datetime(2022,8,24,12,0) + datetime.timedelta(hours=24*ii) for ii in range(2) ]

    # request data
    filename = '/tmp/test'
    files_geomet = request_geomet_grib2(product=product,date=dates,bbox=bbox,crs=crs,filename=filename,silent=True)

    # read data
    data_geomet = read_geomet_grib2(files_geomet,silent=True)

    # --------------------------------------
    # Interpolate data from Geomet
    # --------------------------------------

    # 4 locations ; w/ post-process (negative values set to zero)
    locations = {"lat":[45.25,45.5,45.75,45.75], "lon":[-73.5,-73.5,-73.5,-73.75]}
    interpolate_geomet = interpolate_data(var=data_geomet["var"],lat=data_geomet["lat"],lon=data_geomet["lon"],locations=locations,
                                              bbox=bbox,post_process=True,silent=True)
    np.testing.assert_almost_equal( interpolate_geomet['var'][:,0], [0.06768329, 2.35446854] )
    np.testing.assert_almost_equal( interpolate_geomet['var'][:,1], [0.03124098, 0.65609319] )
    np.testing.assert_almost_equal( interpolate_geomet['var'][:,2], [0.00268678, 0.04813891] )
    np.testing.assert_almost_equal( interpolate_geomet['var'][:,3], [0.00427919, 0.02936855] )

    # 1 location (as list) ; w/ post-process (negative values set to zero)
    locations = {"lat":[45.25], "lon":[-73.5]}
    interpolate_geomet = interpolate_data(var=data_geomet["var"],lat=data_geomet["lat"],lon=data_geomet["lon"],locations=locations,
                                              bbox=bbox,post_process=True,silent=True)
    np.testing.assert_almost_equal( interpolate_geomet['var'][:,0], [0.06768329, 2.35446854] )

    # 1 location (as scalar) ; w/ post-process (negative values set to zero)
    locations = {"lat":45.25, "lon":-73.5}
    interpolate_geomet = interpolate_data(var=data_geomet["var"],lat=data_geomet["lat"],lon=data_geomet["lon"],locations=locations,
                                              bbox=bbox,post_process=True,silent=True)
    np.testing.assert_almost_equal( interpolate_geomet['var'][:,0], [0.06768329, 2.35446854] )


@pytest.mark.filterwarnings("ignore:request_caspar_nc")
def test_interpolate_caspar():

    # --------------------------------------
    # Request and read data from CaSPAr
    # --------------------------------------

    from b2_read_caspar_nc import read_caspar_nc
    product='RDRS_v2.1'
    variable='RDRS_v2.1_A_PR0_SFC'
    dates=[ datetime.datetime(2018,8,8,15,0) + datetime.timedelta(hours=ii) for ii in range(10) ]
    bbox={"lat":{"min":45.0,"max":46.0},"lon":{"min":-74.0,"max":-73.0}}
    foldername='test-data/'
    lintransform={'a':1000.0,'b':0.0}  # convert from m/h to mm/h
    silent=False

    data_caspar = read_caspar_nc( product=product, variable=variable, date=dates, bbox=bbox, foldername=foldername, lintransform=lintransform, silent=True)

    # --------------------------------------
    # Interpolate data from CaSPAr
    # --------------------------------------

    # 4 locations; w/o post-process (leads to negative values)
    locations = {"lat":[45.25,45.5,45.75,45.75], "lon":[-73.5,-73.5,-73.5,-73.75]}
    interpolate_caspar = interpolate_data(var=data_caspar["var"],lat=data_caspar["lat"],lon=data_caspar["lon"],locations=locations,
                                              bbox=bbox,post_process=False,silent=True)

    np.testing.assert_almost_equal( interpolate_caspar['var'][:4,0], [-0.00077367, -0.00034989, -0.00029589, -0.00056612] )
    np.testing.assert_almost_equal( interpolate_caspar['var'][:4,1], [0.04931655, 0.32463379, 0.36022639, 0.0334888 ] )
    np.testing.assert_almost_equal( interpolate_caspar['var'][:4,2], [ 0.00068256,  0.03492004,  0.20845831, -0.0208557 ] )
    np.testing.assert_almost_equal( interpolate_caspar['var'][:4,3], [7.03383407e-05, 2.11351033e-02, 1.49056495e-01, 8.51443785e-02] )

    # 1 location (as list); w/o post-process (leads to negative values)
    locations = {"lat":[45.25], "lon":[-73.5]}
    interpolate_caspar = interpolate_data(var=data_caspar["var"],lat=data_caspar["lat"],lon=data_caspar["lon"],locations=locations,
                                              bbox=bbox,post_process=False,silent=True)
    np.testing.assert_almost_equal( interpolate_caspar['var'][:4,0], [-0.00077367, -0.00034989, -0.00029589, -0.00056612] )

    # 1 location (as scalar); w/o post-process (leads to negative values)
    locations = {"lat":45.25, "lon":-73.5}
    interpolate_caspar = interpolate_data(var=data_caspar["var"],lat=data_caspar["lat"],lon=data_caspar["lon"],locations=locations,
                                              bbox=bbox,return_tmp=False,post_process=False,silent=True)
    np.testing.assert_almost_equal( interpolate_caspar['var'][:4,0], [-0.00077367, -0.00034989, -0.00029589, -0.00056612] )

    # 1 location (as scalar) ; w/ post-process (negative values set to zero)
    locations = {"lat":45.25, "lon":-73.5}
    interpolate_caspar = interpolate_data(var=data_caspar["var"],lat=data_caspar["lat"],lon=data_caspar["lon"],locations=locations,
                                               bbox=bbox,return_tmp=False,post_process=True,silent=True)
    np.testing.assert_almost_equal( interpolate_caspar['var'][:4,0], [0., 0., 0., 0.] )
