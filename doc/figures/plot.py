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


# pyenv activate env-3.8.5-nrcan


import datetime as datetime
import urllib.request


# -----------------------
# add subolder scripts/lib to search path
# -----------------------
import sys
import os
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path+'/../../src')


# ---------------------------------
# create:
#   'test-map-geomet-geo-weather-rdpa24_2018080912.png'
#   'test-map-geomet-geo-weather-rdpa24_legend.png'
# ---------------------------------

plot_api = {'png': [], 'gif': [], 'legend': []}

api_str = "https://geo.weather.gc.ca/geomet?SERVICE=WMS&VERSION=1.3.0&REQUEST=GetMap&LAYERS=RDPA.24F_PR&STYLES=RDPA-WXO&CRS=EPSG:4326&BBOX=45,-74,46,-73&WIDTH=400&HEIGHT=400&FORMAT=image/png&TIME=2018-08-09T12:00:00Z&DIM_REFERENCE_TIME=2018-08-09T12:00:00Z"
ifilename = 'test-map-geomet-geo-weather-rdpa24_2018080912.png'
urllib.request.urlretrieve(api_str, ifilename)
plot_api['png'] = [ifilename]

api_str = "https://geo.weather.gc.ca/geomet?SERVICE=WMS&VERSION=1.3.0&REQUEST=GetLegendGraphic&LAYERS=RDPA.24F_PR&STYLES=RDPA-WXO&CRS=EPSG:4326&BBOX=45,-74,46,-73&WIDTH=400&HEIGHT=400&FORMAT=image/png&TIME=2018-08-09T12:00:00Z&DIM_REFERENCE_TIME=2018-08-09T12:00:00Z"
ifilename = 'test-map-geomet-geo-weather-rdpa24_legend.png'
urllib.request.urlretrieve(api_str, ifilename)
plot_api['legend'] = [ifilename]

print("All files created: ", plot_api)


# ---------------------------------
# create:
#   'test-map-geomet-geo-weather-rdpa6_2018080906.png'
#   'test-map-geomet-geo-weather-rdpa6_legend.png'
# ---------------------------------

plot_api = {'png': [], 'gif': [], 'legend': []}

api_str = "https://geo.weather.gc.ca/geomet?SERVICE=WMS&VERSION=1.3.0&REQUEST=GetMap&LAYERS=RDPA.6F_PR&STYLES=RDPA-WXO&CRS=EPSG:4326&BBOX=45,-74,46,-73&WIDTH=400&HEIGHT=400&FORMAT=image/png&TIME=2018-08-09T06:00:00Z&DIM_REFERENCE_TIME=2018-08-09T06:00:00Z"
ifilename = 'test-map-geomet-geo-weather-rdpa6_2018080906.png'
urllib.request.urlretrieve(api_str, ifilename)
plot_api['png'] = [ifilename]

api_str = "https://geo.weather.gc.ca/geomet?SERVICE=WMS&VERSION=1.3.0&REQUEST=GetLegendGraphic&LAYERS=RDPA.6F_PR&STYLES=RDPA-WXO&CRS=EPSG:4326&BBOX=45,-74,46,-73&WIDTH=400&HEIGHT=400&FORMAT=image/png&TIME=2018-08-09T06:00:00Z&DIM_REFERENCE_TIME=2018-08-09T06:00:00Z"
ifilename = 'test-map-geomet-geo-weather-rdpa6_legend.png'
urllib.request.urlretrieve(api_str, ifilename)
plot_api['legend'] = [ifilename]

print("All files created: ", plot_api)


# ---------------------------------
# create:
#   'test-map-geomet-nrcan-hfe-rdpa24_2018080912.png'
#   'test-map-geomet-nrcan-hfe-rdpa24_legend.png'
# ---------------------------------

from a1_request_geomet_grib2 import request_geomet_grib2
from b1_read_geomet_grib2 import read_geomet_grib2
from cx_plot_data import plot_data

product = 'rdpa:10km:24f'
crs = 'EPSG:4326'
bbox = {"lat":{"min":45.0,"max":46.0},"lon":{"min":-74.0,"max":-73.0}}
dates = [ datetime.datetime(2018,8,9,12,0) ]
lintransform={'a':1.0,'b':0.0}  # no convert of units

# request data
filename = '/tmp/test-map-geomet-nrcan-hfe-rdpa24'
files_geomet = request_geomet_grib2(product=product,date=dates,bbox=bbox,crs=crs,overwrite=True,filename=filename)

# read data
data_geomet = read_geomet_grib2(
    files_geomet,
    lintransform=lintransform)

# plot (1 time step)
plot_geomet = plot_data(
    var=data_geomet["var"],
    lat=data_geomet["lat"],
    lon=data_geomet["lon"],
    date=dates,
    png=True,
    gif=False,
    legend=True,
    cities=True,
    bbox=bbox,
    basefilename='test-map-geomet-nrcan-hfe-rdpa24',
    silent=True)
print("All files created: ", plot_geomet)


# ---------------------------------
# create:
#   'test-map-geomet-nrcan-hfe-rdpa6_2018080906.png'
#   'test-map-geomet-nrcan-hfe-rdpa6_legend.png'
# ---------------------------------

from a1_request_geomet_grib2 import request_geomet_grib2
from b1_read_geomet_grib2 import read_geomet_grib2
from cx_plot_data import plot_data

product = 'rdpa:10km:6f'
crs = 'EPSG:4326'
bbox = {"lat":{"min":45.0,"max":46.0},"lon":{"min":-74.0,"max":-73.0}}
dates = [ datetime.datetime(2018,8,9,6,0) ]
lintransform={'a':1.0,'b':0.0}  # no convert of units

# request data
filename = '/tmp/test-map-geomet-nrcan-hfe-rdpa6'
files_geomet = request_geomet_grib2(product=product,date=dates,bbox=bbox,crs=crs,overwrite=True,filename=filename)

# read data
data_geomet = read_geomet_grib2(
    files_geomet,
    lintransform=lintransform)

# plot (1 time step)
plot_geomet = plot_data(
    var=data_geomet["var"],
    lat=data_geomet["lat"],
    lon=data_geomet["lon"],
    date=dates,
    png=True,
    gif=False,
    legend=True,
    cities=True,
    bbox=bbox,
    basefilename='test-map-geomet-nrcan-hfe-rdpa6',
    silent=True)
print("All files created: ", plot_geomet)


# ---------------------------------
# create:
#     'test-map-caspar-nrcan-hfe_2018080907.png'
#     'test-map-caspar-nrcan-hfe_legend.png'
# ---------------------------------

from b2_read_caspar_nc import read_caspar_nc
from cx_plot_data import plot_data

product='RDRS_v2.1'
variable='RDRS_v2.1_A_PR0_SFC'
dates=[ datetime.datetime(2018,8,9,7,0) ]
bbox={"lat":{"min":45.0,"max":46.0},"lon":{"min":-74.0,"max":-73.0}}
foldername='../../src/test-data/'
lintransform={'a':1000.0,'b':0.0}  # convert from m/h to mm/h
silent=False

data_caspar = read_caspar_nc(
    product=product,
    variable=variable,
    date=dates,
    bbox=bbox,
    foldername=foldername,
    lintransform=lintransform,
    silent=True)

plot_caspar = plot_data(
    var=data_caspar["var"],
    lat=data_caspar["lat"],
    lon=data_caspar["lon"],
    date=dates,
    png=True,
    gif=False,
    legend=True,
    cities=True,
    bbox=bbox,
    basefilename='test-map-caspar-nrcan-hfe',
    silent=True)
print("All files created: ", plot_caspar)


# ---------------------------------
# create:
#     'test-bilinear-interpolation'
# ---------------------------------

os.system('python plot_bilinear_interpolation_results.py > /dev/null 2>&1')
print("All files created:  {'png': ['test-bilinear-interpolation.png']}")


# ---------------------------------
# create:
#     'test-all-coordinates-in-hfe'
# ---------------------------------

os.system('python plot_all-coordinates-in-hfe.py > /dev/null 2>&1')
print("All files created:  {'png': ['test-all-coordinates-in-hfe.png']}")
