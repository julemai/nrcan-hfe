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
import pytest


@pytest.mark.filterwarnings("ignore:request_caspar_nc")
def test_request_caspar():

    # just a warning
    file_caspar = request_caspar_nc()

    assert file_caspar is None
