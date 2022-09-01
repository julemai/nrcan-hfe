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


from cx_plot_data import plot_data
import numpy as np
import datetime as datetime
import pytest


@pytest.mark.filterwarnings("ignore:request_geomet_grib2")
def test_plot_geomet_2ts():

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
    # Plot data from Geomet
    # --------------------------------------

    # plot (2 time steps)
    plot_geomet = plot_data(var=data_geomet["var"],lat=data_geomet["lat"],lon=data_geomet["lon"],date=dates,png=True,gif=True,legend=True,cities=True,bbox=bbox,basefilename='/tmp/test_4dates',silent=True )
    assert np.all( plot_geomet['png'] == ['/tmp/test_4dates_2022082412.png', '/tmp/test_4dates_2022082512.png'] )

@pytest.mark.filterwarnings("ignore:request_geomet_grib2")
def test_plot_geomet_1ts():

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
    # Plot data from Geomet
    # --------------------------------------

    # plot (1 time step)
    plot_geomet = plot_data(var=data_geomet["var"][0],lat=data_geomet["lat"],lon=data_geomet["lon"],date=dates[0],png=True,gif=False,legend=False,cities=True,bbox=bbox,basefilename='/tmp/test_1date',silent=True )
    assert np.all( plot_geomet['png'] == ['/tmp/test_1date_2022082412.png'] )

def test_plot_caspar_2ts():

    # --------------------------------------
    # Request and read data from CaSPAr
    # --------------------------------------

    from b2_read_caspar_nc import read_caspar_nc
    product='RDRS_v2.1'
    variable='RDRS_v2.1_A_PR0_SFC'
    dates=[ datetime.datetime(2018,8,8,15,0) + datetime.timedelta(days=ii) for ii in range(2) ]
    bbox={"lat":{"min":45.0,"max":46.0},"lon":{"min":-74.0,"max":-73.0}}
    foldername='test-data/'
    lintransform={'a':1000.0,'b':0.0}  # convert from m/h to mm/h
    silent=False

    data_caspar = read_caspar_nc( product=product, variable=variable, date=dates, bbox=bbox, foldername=foldername, lintransform=lintransform, silent=True)

    # --------------------------------------
    # Plot data from CaSPAr
    # --------------------------------------

    plot_caspar = plot_data( var=data_caspar["var"], lat=data_caspar["lat"], lon=data_caspar["lon"], date=dates, png=True, gif=True, legend=True, cities=True,bbox=bbox, basefilename='/tmp/test_2dates_caspar', silent=True )
    assert np.all( plot_caspar['png']    == ['/tmp/test_2dates_caspar_2018080815.png', '/tmp/test_2dates_caspar_2018080915.png'] )
    assert np.all( plot_caspar['gif']    == ['/tmp/test_2dates_caspar.gif'] )
    assert np.all( plot_caspar['legend'] == ['/tmp/test_2dates_caspar_legend.png'] )
