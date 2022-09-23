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


from a3_request_hfe_json import request_hfe_json
import pytest
from pathlib import Path


@pytest.mark.filterwarnings("ignore:request_hfe_json")
def test_request_hfe():

    jsonfilebase = dir_path+'/../../data/hfe/historical_flood'
    file_hfe = request_hfe_json(filename=None, jsonfilebase=jsonfilebase,silent=False)

    assert str(Path(file_hfe['json'][0])) == str(Path(dir_path+'/../../data/hfe/historical_flood.json'))
    assert str(Path(file_hfe['json'][1])) == str(Path(dir_path+'/../../data/hfe/historical_flood_event.json'))
