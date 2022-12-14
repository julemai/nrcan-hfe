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


from b3_read_hfe_json import read_hfe_json
import numpy as np
import pytest



# polygon to check if points are within
polygon = [[-30,30], [-130,30], [-130,80], [-30,80]]


@pytest.mark.filterwarnings("ignore:read_hfe_json")
def test_read_hfe_events():

    # ------------------------
    # Read flood events (MULTIPOINTS)
    # ------------------------

    filename=dir_path+'/../../data/hfe/historical_flood_event.json'

    # no filtering; just reading
    data_hfe = read_hfe_json(filename=filename,filtering=False,polygon=None,silent=True)
    np.testing.assert_almost_equal( len(data_hfe['data']['features']) , 376 )

    # filtering w/o checking if points are in polygon
    data_hfe = read_hfe_json(filename=filename,filtering=True,polygon=None,silent=True)
    np.testing.assert_almost_equal( len(data_hfe['data']['features']) , 363 )

    # filtering w/  checking if points are in polygon
    data_hfe = read_hfe_json(filename=filename,filtering=True,polygon=polygon,silent=True)
    np.testing.assert_almost_equal( len(data_hfe['data']['features']) , 361 )

@pytest.mark.filterwarnings("ignore:read_hfe_json")
def test_read_hfe_occurrences():

    # ------------------------
    # Read flood occurences (POINTS)
    # ------------------------

    filename=dir_path+'/../../data/hfe/historical_flood.json'

    # no filtering; just reading
    data_hfe = read_hfe_json(filename=filename,filtering=False,polygon=None,silent=True)
    np.testing.assert_almost_equal( len(data_hfe['data']['features']) , 1999 )

    # filtering w/o checking if points are in polygon
    data_hfe = read_hfe_json(filename=filename,filtering=True,polygon=None,silent=True)
    np.testing.assert_almost_equal( len(data_hfe['data']['features']) , 1949 )

    # filtering w/  checking if points are in polygon
    data_hfe = read_hfe_json(filename=filename,filtering=True,polygon=polygon,silent=True)
    np.testing.assert_almost_equal( len(data_hfe['data']['features']) , 1947 )
