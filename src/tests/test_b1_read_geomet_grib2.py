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


from b1_read_geomet_grib2 import read_geomet_grib2
import datetime as datetime
import numpy as np


def test_read_geomet_1file():

    # Read data (1 file)

    filenames = {datetime.datetime(2022, 8, 24, 12, 0): ['test-data/rdpa-6h_2022082412.grib2']}
    data_geomet = read_geomet_grib2(filenames)

    np.testing.assert_almost_equal( data_geomet["var"][0,0,0:4], [0.09765625, 0.0625, 0.09375, 0.1015625 ] )
    np.testing.assert_almost_equal( data_geomet["lat"][0,0:4],   [45.359908, 45.31104666, 45.2620931, 45.21304786] )
    np.testing.assert_almost_equal( data_geomet["lon"][0,0:4],   [-74.714051, -74.61954144, -74.52526103, -74.43120944] )

def test_read_geomet_4files():

    # Read data (4 files)

    filenames = {datetime.datetime(2022, 8, 24, 0, 0): ['test-data/rdpa-6h_2022082400.grib2'], datetime.datetime(2022, 8, 24, 6, 0): ['test-data/rdpa-6h_2022082406.grib2'], datetime.datetime(2022, 8, 24, 12, 0): ['test-data/rdpa-6h_2022082412.grib2'], datetime.datetime(2022, 8, 24, 18, 0): ['test-data/rdpa-6h_2022082418.grib2']}
    data_geomet = read_geomet_grib2(filenames)

    np.testing.assert_almost_equal( data_geomet["var"][0,0,0:4], [0.1875, 0.125, 0.34375, 0.21875] )
    np.testing.assert_almost_equal( data_geomet["lat"][0,0:4],   [45.359908, 45.31104666, 45.2620931, 45.21304786] )
    np.testing.assert_almost_equal( data_geomet["lon"][0,0:4],   [-74.714051, -74.61954144, -74.52526103, -74.43120944] )
