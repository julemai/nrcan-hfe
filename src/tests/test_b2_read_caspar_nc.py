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


from b2_read_caspar_nc import read_caspar_nc
import numpy as np
import datetime as datetime


product='RDRS_v2.1'
variable='RDRS_v2.1_A_PR0_SFC'
bbox={"lat":{"min":45.0,"max":46.0},"lon":{"min":-74.0,"max":-73.0}}
foldername='test-data/'
lintransform={'a':1000.0,'b':0.0} # convert from m/h to mm/h


def test_read_caspar_1file_string():

    # Read data (1 file provided as string --> returns var as 2D array)

    date=datetime.datetime(2018,8,9,7,0)

    data_caspar = read_caspar_nc(product=product,variable=variable,date=date,bbox=bbox,foldername=foldername,lintransform=lintransform,silent=True)

    np.testing.assert_almost_equal( data_caspar["var"][8,0:4], [0.12983558, 0.0376055, 0., 0. ] )
    np.testing.assert_almost_equal( data_caspar["lat"][8,0:4], [45.78684, 45.76281, 45.73859, 45.71424], 5 )
    np.testing.assert_almost_equal( data_caspar["lon"][8,0:4], [-74.150696, -74.02869, -73.9068, -73.785 ], 5 )

def test_read_caspar_1file_list():

    # Read data (1 file provided as list --> returns var as 3D array)

    date=[datetime.datetime(2018,8,9,7,0)]

    data_caspar = read_caspar_nc(product=product,variable=variable,date=date,bbox=bbox,foldername=foldername,lintransform=lintransform,silent=True)

    np.testing.assert_almost_equal( data_caspar["var"][0,8,0:4], [0.12983558, 0.0376055, 0., 0. ] )
    np.testing.assert_almost_equal( data_caspar["lat"][8,0:4],   [45.78684, 45.76281, 45.73859, 45.71424], 5 )
    np.testing.assert_almost_equal( data_caspar["lon"][8,0:4],   [-74.150696, -74.02869, -73.9068, -73.785 ], 5 )

def test_read_caspar_2files_list():

    # Read data (2 files --> returns var as 3D array)

    date=[datetime.datetime(2018,8,9,7,0),datetime.datetime(2018,8,9,14,0)]

    data_caspar = read_caspar_nc(product=product,variable=variable,date=date,bbox=bbox,foldername=foldername,lintransform=lintransform,silent=True)

    np.testing.assert_almost_equal( data_caspar["var"][0,8,0:4], [0.12983558, 0.0376055, 0., 0. ] )
    np.testing.assert_almost_equal( data_caspar["lat"][8,0:4],   [45.78684, 45.76281, 45.73859, 45.71424], 5 )
    np.testing.assert_almost_equal( data_caspar["lon"][8,0:4],   [-74.150696, -74.02869, -73.9068, -73.785 ], 5 )
