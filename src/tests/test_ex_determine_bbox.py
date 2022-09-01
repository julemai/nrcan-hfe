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


from ex_determine_bbox import determine_bbox



locations = {"lat":[45.25,45.5,45.75,45.75], "lon":[-73.5,-73.5,-73.5,-73.75]}
feature_point = {"type" : "Feature","geometry" : {"type" : "Point","coordinates" : [ -73.5, 45.25 ]}, "properties" : { }}
feature_multipoint = {"type" : "Feature","geometry" : {"type" : "MultiPoint","coordinates" : [[-73.5,45.25],[-73.5,45.5],[-73.5,45.75],[-73.75,45.75]]}, "properties" : { }}


def test_derive_bbox_scalar():

    # --------------------------------------
    # Derive bounding box lat/lon (scalar)
    # --------------------------------------
    bbox = determine_bbox(lat=locations["lat"][0],lon=locations["lon"][0],bbox_buffer=0.5,silent=True)
    assert bbox == {'lat': {'min': 44.75, 'max': 45.75}, 'lon': {'min': -74.0, 'max': -73.0}}

def test_derive_bbox_list():

    # --------------------------------------
    # Derive bounding box lat/lon (list)
    # --------------------------------------
    bbox = determine_bbox(lat=locations["lat"],lon=locations["lon"],bbox_buffer=0.5,silent=True)
    assert bbox == {'lat': {'min': 44.75, 'max': 46.25}, 'lon': {'min': -74.25, 'max': -73.0}}

def test_derive_bbox_point():

    # --------------------------------------
    # Derive bounding box feature (point)
    # --------------------------------------
    bbox = determine_bbox(feature=feature_point,bbox_buffer=0.5,silent=True)
    assert bbox == {'lat': {'min': 44.75, 'max': 45.75}, 'lon': {'min': -74.0, 'max': -73.0}}

def test_derive_bbox_multipoint():

    # --------------------------------------
    # Derive bounding box feature (multipoint)
    # --------------------------------------
    bbox = determine_bbox(feature=feature_multipoint,bbox_buffer=0.5,silent=True)
    assert bbox == {'lat': {'min': 44.75, 'max': 46.25}, 'lon': {'min': -74.25, 'max': -73.0}}
