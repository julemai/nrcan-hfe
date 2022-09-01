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


from fx_determine_dates import determine_dates
import numpy as np
import datetime as datetime



period = [datetime.datetime(2000,5,28,3,0),datetime.datetime(2000,5,29,20,59)]
feature_w_end  = {"type" : "Feature", "geometry" : {"type" : "Point","coordinates" : [ -73.5, 45.25 ]}, "properties" : { "start_date" : "2000-05-28", "end_date" : "2000-05-29" }}
feature_wo_end = {"type" : "Feature", "geometry" : {"type" : "Point","coordinates" : [ -73.5, 45.25 ]}, "properties" : { "start_date" : "2000-05-28", "end_date" : None }}


def test_derive_tsteps_rdrs_period_wobuffer():
    # --------------------------------------
    # Derive time steps for RDRS_v2.1 w/ start_date and end-date specified; w/o buffer
    # --------------------------------------
    tsteps = determine_dates(start_date=period[0],end_date=period[1],feature=None,product='RDRS_v2.1',dates_buffer=[0.0,0.0],silent=True)

    assert tsteps[0:1]   == [datetime.datetime(2000, 5, 28, 3, 0)]
    assert tsteps[10:11] == [datetime.datetime(2000, 5, 28, 13, 0)]
    assert tsteps[-1:]   == [datetime.datetime(2000, 5, 29, 20, 0)]
    assert len(tsteps)   == 42

def test_derive_tsteps_rdrs_period_wbuffer():
    # --------------------------------------
    # Derive time steps for RDRS_v2.1 w/ start_date and end-date specified; w/ buffer
    # --------------------------------------
    tsteps = determine_dates(start_date=period[0],end_date=period[1],feature=None,product='RDRS_v2.1',dates_buffer=[3.0,1.0],silent=True)

    assert tsteps[0:1]   == [datetime.datetime(2000, 5, 25, 3, 0)]
    assert tsteps[10:11] == [datetime.datetime(2000, 5, 25, 13, 0)]
    assert tsteps[-1:]   == [datetime.datetime(2000, 5, 30, 20, 0)]
    assert len(tsteps)   == 138

def test_derive_tsteps_rdrs_start_wbuffer():
    # --------------------------------------
    # Derive time steps for RDRS_v2.1 w/ start_date; w/ buffer
    # --------------------------------------
    tsteps = determine_dates(start_date=period[0],end_date=None,feature=None,product='RDRS_v2.1',dates_buffer=[3.0,1.0],silent=True)

    assert tsteps[0:1]   == [datetime.datetime(2000, 5, 25, 3, 0)]
    assert tsteps[10:11] == [datetime.datetime(2000, 5, 25, 13, 0)]
    assert tsteps[-1:]   == [datetime.datetime(2000, 5, 30, 3, 0)]
    assert len(tsteps)   == 121

def test_derive_tsteps_rdrs_feature_wend_wbuffer():
    # --------------------------------------
    # Derive time steps for RDRS_v2.1 w/ end-date specified in feature; w/ buffer
    # --------------------------------------
    tsteps = determine_dates(start_date=None,end_date=None,feature=feature_w_end,product='RDRS_v2.1',dates_buffer=[3.0,1.0],silent=True)

    assert tsteps[0:1]   == [datetime.datetime(2000, 5, 25, 0, 0)]
    assert tsteps[10:11] == [datetime.datetime(2000, 5, 25, 10, 0)]
    assert tsteps[-1:]   == [datetime.datetime(2000, 5, 30, 23, 0)]
    assert len(tsteps)   == 144

def test_derive_tsteps_rdrs_feature_woend_wbuffer():
    # --------------------------------------
    # Derive time steps for RDRS_v2.1 w/o end-date specified in feature; w/ buffer
    # --------------------------------------
    tsteps = determine_dates(start_date=None,end_date=None,feature=feature_wo_end,product='RDRS_v2.1',dates_buffer=[3.0,1.0],silent=True)

    assert tsteps[0:1]   == [datetime.datetime(2000, 5, 25, 0, 0)]
    assert tsteps[10:11] == [datetime.datetime(2000, 5, 25, 10, 0)]
    assert tsteps[-1:]   == [datetime.datetime(2000, 5, 30, 0, 0)]
    assert len(tsteps)   == 121

def test_derive_tsteps_rdpa24_period_wobuffer():
    # --------------------------------------
    # Derive time steps for RDPA-24h w/ start_date and end-date specified; w/o buffer
    # --------------------------------------
    tsteps = determine_dates(start_date=period[0],end_date=period[1],feature=None,product='rdpa:10km:24f',dates_buffer=[0.0,0.0],silent=True)

    assert np.all(tsteps == [datetime.datetime(2000, 5, 28, 12, 0), datetime.datetime(2000, 5, 29, 12, 0)])

def test_derive_tsteps_rdpa24_period_wbuffer():
    # --------------------------------------
    # Derive time steps for RDPA-24h w/ start_date and end-date specified; w/ buffer
    # --------------------------------------
    tsteps = determine_dates(start_date=period[0],end_date=period[1],feature=None,product='rdpa:10km:24f',dates_buffer=[3.0,1.0],silent=True)

    assert np.all(tsteps == [    datetime.datetime(2000, 5, 25, 12, 0), datetime.datetime(2000, 5, 26, 12, 0),
                                 datetime.datetime(2000, 5, 27, 12, 0), datetime.datetime(2000, 5, 28, 12, 0),
                                 datetime.datetime(2000, 5, 29, 12, 0), datetime.datetime(2000, 5, 30, 12, 0)])

def test_derive_tsteps_rdpa6_period_wobuffer():
    # --------------------------------------
    # Derive time steps for RDPA-6h w/ start_date and end-date specified; w/o buffer
    # --------------------------------------
    tsteps = determine_dates(start_date=period[0],end_date=period[1],feature=None,product='rdpa:10km:6f',dates_buffer=[0.0,0.0],silent=True)

    assert np.all(tsteps == [    datetime.datetime(2000, 5, 28, 6, 0), datetime.datetime(2000, 5, 28, 12, 0),
                                 datetime.datetime(2000, 5, 28, 18, 0), datetime.datetime(2000, 5, 29, 0, 0),
                                 datetime.datetime(2000, 5, 29, 6, 0), datetime.datetime(2000, 5, 29, 12, 0),
                                 datetime.datetime(2000, 5, 29, 18, 0)])

def test_derive_tsteps_rdpa6_period_wbuffer():
    # --------------------------------------
    # Derive time steps for RDPA-6h w/ start_date and end-date specified; w/ buffer
    # --------------------------------------
    tsteps = determine_dates(start_date=period[0],end_date=period[1],feature=None,product='rdpa:10km:6f',dates_buffer=[3.0,1.0],silent=True)

    assert np.all(tsteps == [    datetime.datetime(2000, 5, 25, 6, 0), datetime.datetime(2000, 5, 25, 12, 0),
                                 datetime.datetime(2000, 5, 25, 18, 0), datetime.datetime(2000, 5, 26, 0, 0),
                                 datetime.datetime(2000, 5, 26, 6, 0), datetime.datetime(2000, 5, 26, 12, 0),
                                 datetime.datetime(2000, 5, 26, 18, 0), datetime.datetime(2000, 5, 27, 0, 0),
                                 datetime.datetime(2000, 5, 27, 6, 0), datetime.datetime(2000, 5, 27, 12, 0),
                                 datetime.datetime(2000, 5, 27, 18, 0), datetime.datetime(2000, 5, 28, 0, 0),
                                 datetime.datetime(2000, 5, 28, 6, 0), datetime.datetime(2000, 5, 28, 12, 0),
                                 datetime.datetime(2000, 5, 28, 18, 0), datetime.datetime(2000, 5, 29, 0, 0),
                                 datetime.datetime(2000, 5, 29, 6, 0), datetime.datetime(2000, 5, 29, 12, 0),
                                 datetime.datetime(2000, 5, 29, 18, 0), datetime.datetime(2000, 5, 30, 0, 0),
                                 datetime.datetime(2000, 5, 30, 6, 0), datetime.datetime(2000, 5, 30, 12, 0),
                                 datetime.datetime(2000, 5, 30, 18, 0)])
