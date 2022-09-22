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


from a1_request_geomet_grib2 import request_geomet_grib2
import datetime as datetime
import pytest



@pytest.mark.filterwarnings("ignore:request_geomet_grib2")
def test_request_geomet_1file():

    # Request data (1 file)

    date = datetime.datetime(2022,8,24,12,0)
    bbox = {"lat":{"min":45.0,"max":46.0},"lon":{"min":-74.0,"max":-73.0}}
    filename = dir_path+'/../test-data/rdpa-6h'
    file_geomet = request_geomet_grib2(product='rdpa:10km:6f',date=date,bbox=bbox,crs='EPSG:4326',filename=filename,overwrite=False)

    assert file_geomet == {datetime.datetime(2022, 8, 24, 12, 0): [dir_path+'/../test-data/rdpa-6h_2022082412.grib2']}

@pytest.mark.filterwarnings("ignore:request_geomet_grib2")
def test_request_geomet_4files():

    # Request data (4 files)

    dates = [ datetime.datetime(2022,8,24,0,0) + datetime.timedelta(hours=6*ii) for ii in range(4) ]
    bbox = {"lat":{"min":45.0,"max":46.0},"lon":{"min":-74.0,"max":-73.0}}
    filename = dir_path+'/../test-data/rdpa-6h'
    file_geomet = request_geomet_grib2(product='rdpa:10km:6f',date=dates,bbox=bbox,crs='EPSG:4326',filename=filename,overwrite=False)

    assert file_geomet == {datetime.datetime(2022, 8, 24, 0, 0): [dir_path+'/../test-data/rdpa-6h_2022082400.grib2'], datetime.datetime(2022, 8, 24, 6, 0): [dir_path+'/../test-data/rdpa-6h_2022082406.grib2'], datetime.datetime(2022, 8, 24, 12, 0): [dir_path+'/../test-data/rdpa-6h_2022082412.grib2'], datetime.datetime(2022, 8, 24, 18, 0): [dir_path+'/../test-data/rdpa-6h_2022082418.grib2']}
