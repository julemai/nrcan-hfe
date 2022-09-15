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


from a2_request_caspar_nc import request_caspar_nc
import datetime as datetime
import pytest


@pytest.mark.filterwarnings("ignore:request_caspar_nc")
def test_request_caspar():

    product    = 'RDRS_v2.1'
    variable   = 'RDRS_v2.1_A_PR0_SFC'
    date       = [ datetime.datetime(2018,8,9,7,0), datetime.datetime(2018,8,9,8,0), datetime.datetime(2018,8,10,7,0), datetime.datetime(2018,8,10,13,0) ]
    foldername = 'test-data/'

    file_caspar = request_caspar_nc(product=product,variable=variable,date=date,foldername=foldername,silent=True)

    assert file_caspar[date[0]] == ['test-data/2018080812.nc']
    assert file_caspar[date[1]] == ['test-data/2018080812.nc']
    assert file_caspar[date[2]] == ['test-data/2018080912.nc']
    assert file_caspar[date[3]] == []
