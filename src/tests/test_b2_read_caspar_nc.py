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
import pytest


variable='RDRS_v2.1_A_PR0_SFC'
bbox={"lat":{"min":45.0,"max":46.0},"lon":{"min":-74.0,"max":-73.0}}
lintransform={'a':1000.0,'b':0.0} # convert from m/h to mm/h


def test_read_caspar_1file():

    # Read data (1 file)

    filenames={datetime.datetime(2018, 8, 9, 7, 0): ['test-data/2018080812.nc']}

    data_caspar = read_caspar_nc(variable=variable,filenames=filenames,bbox=bbox,lintransform=lintransform,silent=True)

    np.testing.assert_almost_equal( np.shape(data_caspar["var"]), (1, 14, 12) )
    np.testing.assert_almost_equal( data_caspar["var"][0,8,0:4], [0.12983558, 0.0376055, 0., 0. ] )
    np.testing.assert_almost_equal( data_caspar["lat"][8,0:4],   [45.78684, 45.76281, 45.73859, 45.71424], 5 )
    np.testing.assert_almost_equal( data_caspar["lon"][8,0:4],   [-74.150696, -74.02869, -73.9068, -73.785 ], 5 )
    assert( np.all( data_caspar["time"] == [datetime.datetime(2018, 8, 9, 7, 0)]) )

def test_read_caspar_2files():

    # Read data (2 files)

    filenames={datetime.datetime(2018, 8, 9, 7, 0): ['test-data/2018080812.nc'], datetime.datetime(2018, 8, 9, 14, 0): ['test-data/2018080912.nc']}

    data_caspar = read_caspar_nc(variable=variable,filenames=filenames,bbox=bbox,lintransform=lintransform,silent=True)

    np.testing.assert_almost_equal( np.shape(data_caspar["var"]), (2, 14, 12) )
    np.testing.assert_almost_equal( data_caspar["var"][0,8,0:4], [0.12983558, 0.0376055, 0., 0. ] )
    np.testing.assert_almost_equal( data_caspar["lat"][8,0:4],   [45.78684, 45.76281, 45.73859, 45.71424], 5 )
    np.testing.assert_almost_equal( data_caspar["lon"][8,0:4],   [-74.150696, -74.02869, -73.9068, -73.785 ], 5 )
    assert( np.all( data_caspar["time"] == [datetime.datetime(2018, 8, 9, 7, 0),datetime.datetime(2018, 8, 9, 14, 0)]) )

@pytest.mark.filterwarnings("ignore:read_caspar_nc")
def test_read_caspar_3files():

    # Read data (2 files; no file for third date found)

    filenames={datetime.datetime(2018, 8, 9, 7, 0): ['test-data/2018080812.nc'], datetime.datetime(2018, 8, 9, 14, 0): ['test-data/2018080912.nc'], datetime.datetime(2018, 8, 10, 13, 0): []}

    data_caspar = read_caspar_nc(variable=variable,filenames=filenames,bbox=bbox,lintransform=lintransform,silent=True)

    np.testing.assert_almost_equal( np.shape(data_caspar["var"]), (2, 14, 12) )
    np.testing.assert_almost_equal( data_caspar["var"][0,8,0:4], [0.12983558, 0.0376055, 0., 0. ] )
    np.testing.assert_almost_equal( data_caspar["lat"][8,0:4],   [45.78684, 45.76281, 45.73859, 45.71424], 5 )
    np.testing.assert_almost_equal( data_caspar["lon"][8,0:4],   [-74.150696, -74.02869, -73.9068, -73.785 ], 5 )
    assert( np.all( data_caspar["time"] == [datetime.datetime(2018, 8, 9, 7, 0),datetime.datetime(2018, 8, 9, 14, 0)]) )
