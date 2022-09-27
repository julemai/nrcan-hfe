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


from gx_identify_precipitation_event import identify_precipitation_event
import numpy as np
import datetime as datetime
import pytest





@pytest.mark.filterwarnings("ignore:request_geomet_grib2")
@pytest.mark.filterwarnings("ignore:read_hfe_json")
def test_identify_precipitation_event():

    # --------------------------------------
    # Read HFE database and select one feature
    # --------------------------------------

    from b3_read_hfe_json import read_hfe_json

    filename        = dir_path+'/../../data/hfe/historical_flood.json'
    filtering       = True
    polygon         = None
    return_filtered = False
    silent          = True

    data_hfe = read_hfe_json(filename=filename,filtering=filtering,polygon=polygon,return_filtered=return_filtered,silent=silent)

    ifeature = 1277
    feature = data_hfe['data']['features'][ifeature]

    # --------------------------------------
    # Determine dates and bbox to request
    # --------------------------------------

    from ex_determine_bbox import determine_bbox
    from fx_determine_dates import determine_dates

    product = 'rdpa:10km:6f'
    crs = 'EPSG:4326'
    bbox_buffer  = 0.5
    dates_buffer = [5.0,5.0]

    bbox = determine_bbox(feature=feature,bbox_buffer=bbox_buffer,silent=True)
    dates = determine_dates(feature=feature,product=product,dates_buffer=dates_buffer,silent=True)

    # --------------------------------------
    # Request and read data from Geomet
    # --------------------------------------

    from a1_request_geomet_grib2 import request_geomet_grib2
    from b1_read_geomet_grib2 import read_geomet_grib2

    # request data
    filename = '/tmp/pytest_rdpa_10km_6f_feature_'+str(ifeature)
    files_geomet = request_geomet_grib2(product=product,date=dates,bbox=bbox,crs=crs,filename=filename,silent=True)

    # read data
    data_geomet = read_geomet_grib2(filenames=files_geomet,silent=True)

    # --------------------------------------
    # Interpolate data from Geomet
    # --------------------------------------

    from dx_interpolate_data import interpolate_data, plot_interpolated

    var       = data_geomet["var"]
    lat       = data_geomet["lat"]
    lon       = data_geomet["lon"]
    dates     = data_geomet["time"]   # does not contain time steps with missing data anymore

    locations = {'lon':np.array([feature['geometry']['coordinates'][0]]),'lat':np.array([feature['geometry']['coordinates'][1]])}
    interpolate_geomet = interpolate_data(var=var,lat=lat,lon=lon,locations=locations,bbox=bbox,post_process=True,silent=True)

    # --------------------------------------
    # Derive time steps for RDRS_v2.1 w/ start_date and end-date specified; w/o buffer
    # --------------------------------------
    idx_event = identify_precipitation_event(feature=feature,product=product,dates=dates,data=interpolate_geomet,length_window_d=2,min_prec_window=3.0,min_prec=0.001,silent=True)

    nlocations = 1
    sum_prec = [ np.sum(interpolate_geomet['var'][idx_event[ilocation],ilocation]) for ilocation in range(nlocations) ]

    assert idx_event         == [[15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28]]
    assert len(idx_event[0]) == 14
    np.testing.assert_almost_equal( sum_prec[0] , 121.89938715300597 )
