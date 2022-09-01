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

import numpy as np
import datetime as datetime
import scipy as scipy


__all__ = ['determine_bbox']


def determine_bbox(lat=None,lon=None,feature=None,bbox_buffer=0.5,silent=True):
    """
        Derives the bounding box around a set of locations.


        Definition
        ----------
        def determine_bbox(lat=None,lon=None,bbox_buffer=0.5,silent=True)


        Input           Format         Description
        -----           -----          -----------
        lat             scalar or      Latitude(s) of location(s) the bounding box is requested for. Number of
                        1D array       latitudes must be the same as number of longitudes. If "lat" is not specified,
                                       "feature" needs to be provided.
                                       Default: None

        lon             scalar or      Longitudes(s) of location(s) the bounding box is requested for. Number of
                        1D array       longitudes must be the same as number of latitudes. If "lon" is not specified,
                                       "feature" needs to be provided.
                                       Default: None

        feature         dict           Dictionary containing a feature with either a point or multipoints,e.g.,
                                            {    "type" : "Feature",
                                                 "geometry" : {
                                                    "type" : "Point",
                                                    "coordinates" : [ -126.0, 52.0 ]
                                                 },
                                                 "properties" : { ... }
                                            }
                                       If "feature" is not provided, "lat" and "lon" needs to be specified.
                                       Default: None

        bbox_buffer     float          Buffer around points specified (in degree). If set to zero, points will be on
                                       edge/vertex of bounding box. If set to positive value, points will be truly
                                       located within the bounding box with at least "bbox_buffer" distance to an
                                       edge/vertex. Negative values for bbox_buffer are not allowed.
                                       Values are in [degree].
                                       Default: 0.5

        silent          Boolean        If set to True, nothing will be printed to terminal.
                                       Default: True

        Output          Format         Description
        -----           -----          -----------
        {lat,lon}       dict           Dictionary specifying the bounding box, i.e.,
                                         {"lat":{"min":min_lat,"max":max_lat},
                                          "lon":{"min":min_lon,"max":max_lon}}


        Description
        -----------
        Derives the bounding box around a set of locations. (lat,lon given as lists)
        or for a feature which can contain either a single point or multipoints. The
        buffer area around the points can be specified or the deafult value of 0.5
        [degrees] can be used. The function returns a dictionary conatining the information
        of the bounding box, i.e.,
             {"lat":{"min":min_lat,"max":max_lat},
              "lon":{"min":min_lon,"max":max_lon}}


        Restrictions
        ------------
        Only individual features can be handled; no list of features.


        Examples
        --------

        >>> locations = {"lat":[45.25,45.5,45.75,45.75], "lon":[-73.5,-73.5,-73.5,-73.75]}
        >>> feature_point = {"type" : "Feature","geometry" : {"type" : "Point","coordinates" : [ -73.5, 45.25 ]}, "properties" : { }}
        >>> feature_multipoint = {"type" : "Feature","geometry" : {"type" : "MultiPoint","coordinates" : [[-73.5,45.25],[-73.5,45.5],[-73.5,45.75],[-73.75,45.75]]}, "properties" : { }}

        >>> # --------------------------------------
        >>> # Derive bounding box lat/lon (scalar)
        >>> # --------------------------------------
        >>> bbox = determine_bbox(lat=locations["lat"][0],lon=locations["lon"][0],bbox_buffer=0.5,silent=True)
        >>> print("bbox = ",bbox)
        bbox =  {'lat': {'min': 44.75, 'max': 45.75}, 'lon': {'min': -74.0, 'max': -73.0}}

        >>> # --------------------------------------
        >>> # Derive bounding box lat/lon (list)
        >>> # --------------------------------------
        >>> bbox = determine_bbox(lat=locations["lat"],lon=locations["lon"],bbox_buffer=0.5,silent=True)
        >>> print("bbox = ",bbox)
        bbox =  {'lat': {'min': 44.75, 'max': 46.25}, 'lon': {'min': -74.25, 'max': -73.0}}

        >>> # --------------------------------------
        >>> # Derive bounding box feature (point)
        >>> # --------------------------------------
        >>> bbox = determine_bbox(feature=feature_point,bbox_buffer=0.5,silent=True)
        >>> print("bbox = ",bbox)
        bbox =  {'lat': {'min': 44.75, 'max': 45.75}, 'lon': {'min': -74.0, 'max': -73.0}}

        >>> # --------------------------------------
        >>> # Derive bounding box feature (multipoint)
        >>> # --------------------------------------
        >>> bbox = determine_bbox(feature=feature_multipoint,bbox_buffer=0.5,silent=True)
        >>> print("bbox = ",bbox)
        bbox =  {'lat': {'min': 44.75, 'max': 46.25}, 'lon': {'min': -74.25, 'max': -73.0}}



        License
        -------
        This file is part of the HFE code library for NRCan for "Improving the
        characterization of the flood-producing precipitation events in the
        Historic Flood Event (HFE) database by using the CaPA reanalysis/analysis".

        The NRCan-HFE code library is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 2.1 of the License, or
        (at your option) any later version.

        The NRCan-HFE code library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with The NRCan-HFE code library.
        If not, see <https://github.com/julemai/nrcan-hfe/blob/main/LICENSE>.

        Copyright 2022 Juliane Mai - juliane.mai@uwaterloo.ca


        History
        -------
        Written,  Juliane Mai, August 2022
    """

    # checking inputs
    if ((lat is None) and (lon is None) and (feature is None)):
        raise ValueError("determine_bbox: either lat/lon or feature needs to be specified")
    if (not(lat is None) and not(lon is None) and not(feature is None)):
        raise ValueError("determine_bbox: either lat/lon or feature needs to be specified")
    if ((lat is None) and not(lon is None)):
        raise ValueError("determine_bbox: if lon is specified, lat also needs to be specified")
    if ((lon is None) and not(lat is None)):
        raise ValueError("determine_bbox: if lat is specified, lon also needs to be specified")
    if (bbox_buffer < 0.0):
        raise ValueError("determine_bbox: value for bbox_buffer needs to be non-negative")
    if (not(lat is None) and not(lon is None)) and (type(lat) != float) and (type(lat) != int):
        if (len(lat) != len(lon)):
            raise ValueError("determine_bbox: length of lat and lon needs to be the same")


    # initialize return
    result = {}

    # get coordinates and make sure they are always a list of coordinates
    if (not(lat is None) and not(lon is None)):
        if type(lat) == float or type(lat) == int:
            coords = [ [lon,lat] ]
        else:
            coords = [ [lon[iii],lat[iii]] for iii,ii in enumerate(lat) ]
    if (not(feature is None)):
        coords = feature["geometry"]["coordinates"]
        if (feature["geometry"]["type"] == "Point"):
            coords = [ coords ]
    coords = np.array(coords)

    # find min and max in lat and lon
    min_lon = np.min( coords[:,0] )
    max_lon = np.max( coords[:,0] )
    min_lat = np.min( coords[:,1] )
    max_lat = np.max( coords[:,1] )

    # derive bounding box
    bbox = {"lat":{"min":min_lat-bbox_buffer,"max":max_lat+bbox_buffer},
            "lon":{"min":min_lon-bbox_buffer,"max":max_lon+bbox_buffer}}

    result = bbox

    # ----------------------------------------------------
    # Done.
    # ----------------------------------------------------
    return result



if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # locations = {"lat":[45.25,45.5,45.75,45.75], "lon":[-73.5,-73.5,-73.5,-73.75]}
    # feature_point = {"type" : "Feature","geometry" : {"type" : "Point","coordinates" : [ -73.5, 45.25 ]}, "properties" : { }}
    # feature_multipoint = {"type" : "Feature","geometry" : {"type" : "MultiPoint","coordinates" : [[-73.5,45.25],[-73.5,45.5],[-73.5,45.75],[-73.75,45.75]]}, "properties" : { }}

    # feature=None
    # lat=locations["lat"]
    # lon=locations["lon"]
    # bbox_buffer=0.5

    # # get coordinates and make sure they are always a list of coordinates
    # if (not(lat is None) and not(lon is None)):
    #     if type(lat) == float or type(lat) == int:
    #         coords = [ [lon,lat] ]
    #     else:
    #         coords = [ [lon[iii],lat[iii]] for iii,ii in enumerate(lat) ]
    # if (not(feature is None)):
    #     coords = feature["geometry"]["coordinates"]
    #     if (feature["geometry"]["type"] == "Point"):
    #         coords = [ coords ]
    # coords = np.array(coords)

    # # find min and max in lat and lon
    # min_lon = np.min( coords[:,0] )
    # max_lon = np.max( coords[:,0] )
    # min_lat = np.min( coords[:,1] )
    # max_lat = np.max( coords[:,1] )

    # # derive bounding box
    # bbox = {"lat":{"min":min_lat-bbox_buffer,"max":max_lat+bbox_buffer},
    #         "lon":{"min":min_lon-bbox_buffer,"max":max_lon+bbox_buffer}}

    # result = bbox



    # # --------------------------------------
    # # Derive bounding box lat/lon (scalar)
    # # --------------------------------------
    # bbox = determine_bbox(lat=locations["lat"][0],lon=locations["lon"][0],bbox_buffer=0.5,silent=True)
    # print("bbox = ",bbox)

    # # --------------------------------------
    # # Derive bounding box lat/lon (list)
    # # --------------------------------------
    # bbox = determine_bbox(lat=locations["lat"],lon=locations["lon"],bbox_buffer=0.5,silent=True)
    # print("bbox = ",bbox)

    # # --------------------------------------
    # # Derive bounding box feature (point)
    # # --------------------------------------
    # bbox = determine_bbox(feature=feature_point,bbox_buffer=0.5,silent=True)
    # print("bbox = ",bbox)

    # # --------------------------------------
    # # Derive bounding box feature (multipoint)
    # # --------------------------------------
    # bbox = determine_bbox(feature=feature_multipoint,bbox_buffer=0.5,silent=True)
    # print("bbox = ",bbox)
