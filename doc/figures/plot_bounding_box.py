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
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import json as json

# -----------------------
# add subolder scripts/lib to search path
# -----------------------
import sys
import os
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path+'/../../src')


from a1_request_geomet_grib2 import request_geomet_grib2
from b1_read_geomet_grib2 import read_geomet_grib2
from ex_determine_bbox import determine_bbox



# locations shown in plot
locations_singlepoint = {"lat":[45.25], "lon":[-73.5]}
locations_multipoint  = {"lat":[45.25,45.5,45.75,45.75], "lon":[-73.5,-73.5,-73.5,-73.75]}

# bounding box buffers
bbox_buffer = [0.1,0.5,1.0]

# determine bounding boxes
bbox = {}
data = {}

# single-point
bbox_tmp = {}
data_tmp = {}
for ibbox_buffer in bbox_buffer:

    # get bounding box
    ibbox = determine_bbox(lat=locations_singlepoint["lat"],lon=locations_singlepoint["lon"],bbox_buffer=ibbox_buffer,silent=True)
    bbox_tmp[ibbox_buffer] = ibbox

    # get data for bounding box
    product = 'rdpa:10km:6f'
    crs = 'EPSG:4326'
    dates = [ datetime.datetime(2018,8,9,6,0) ]
    lintransform={'a':1.0,'b':0.0}  # no convert of units

    # request data
    filename = '/tmp/test-map-geomet-nrcan-hfe-rdpa6-bbox-spoint-buf'+str(ibbox_buffer)
    files_geomet = request_geomet_grib2(product=product,date=dates,bbox=ibbox,crs=crs,overwrite=False,filename=filename)

    # read data
    data_geomet = read_geomet_grib2(
        files_geomet,
        lintransform=lintransform)
    data_tmp[ibbox_buffer] = data_geomet

bbox['single'] = bbox_tmp
data['single'] = data_tmp

# multi-point
bbox_tmp = {}
data_tmp = {}
for ibbox_buffer in bbox_buffer:

    # get bounding box
    ibbox = determine_bbox(lat=locations_multipoint["lat"],lon=locations_multipoint["lon"],bbox_buffer=ibbox_buffer,silent=True)
    bbox_tmp[ibbox_buffer] = ibbox

    # get data for bounding box
    product = 'rdpa:10km:6f'
    crs = 'EPSG:4326'
    dates = [ datetime.datetime(2018,8,9,6,0) ]
    lintransform={'a':1.0,'b':0.0}  # no convert of units

    # request data
    filename = '/tmp/test-map-geomet-nrcan-hfe-rdpa6-bbox-mpoint-buf'+str(ibbox_buffer)
    files_geomet = request_geomet_grib2(product=product,date=dates,bbox=ibbox,crs=crs,overwrite=False,filename=filename)

    # read data
    data_geomet = read_geomet_grib2(
        files_geomet,
        lintransform=lintransform)
    data_tmp[ibbox_buffer] = data_geomet

bbox['multi'] = bbox_tmp
data['multi'] = data_tmp








# -------------------------------------------------------------------------
# create plot
# -------------------------------------------------------------------------

# plot settings
textsize    = 14          # standard text size
lwidth      = 1.5         # linewidth
alwidth     = 1.0         # axis line width
dpi         = 100         # dpi=100 --> filesize 220kB, dpi=600 --> filesize 1.8MB
ifig        = 0           # initialize counter
transparent = False
bbox_inches = 'tight'
pad_inches  = 0.035
mpl.use('Agg') # set directly after import matplotlib
mpl.rc('figure', figsize=(8.27,11.69)) # a4 portrait
mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
mpl.rc('text.latex') #, unicode=True)
mpl.rc('savefig', dpi=dpi, format='png')
mpl.rc('font', size=textsize)
mpl.rc('lines', linewidth=lwidth, color='black')
mpl.rc('axes', linewidth=alwidth, labelcolor='black')
mpl.rc('path', simplify=False) # do not remove

llxbbox     = 0.0         # x-anchor legend bounding box
llybbox     = 1.0         # y-anchor legend bounding box
llrspace    = 0.25          # spacing between rows in legend
llcspace    = 1.0         # spacing between columns in legend
llhtextpad  = 0.8         # the pad between the legend handle and text
llhlength   = 1.5         # the length of the legend handles
frameon     = False       # if True, draw a frame around the legend. If None, use rc

# colors (gathered from reference legend created by Geomet)
# e.g., run the following request and extract colors:
# "https://geo.weather.gc.ca/geomet?SERVICE=WMS&VERSION=1.3.0&REQUEST=GetLegendGraphic&LAYERS=RDPA.6F_PR&STYLES=RDPA-WXO&CRS=EPSG:4326&BBOX=45,-74,46,-73&WIDTH=400&HEIGHT=400&FORMAT=image/png"
ocean_color = (151/256., 183/256., 224/256.)
colors = [
    '#ffffff', #   0.00           - white
    '#98cbfe', #   0.10 -    0.50 - light blue
    '#0098fe', #   0.50 -    1.00 - medium blue
    '#222cff', #   1.00 -    2.50 - dark blue
    '#00fe65', #   2.50 -    5.00 - light green
    '#00cb00', #   5.00 -    7.50 - medium green
    '#009800', #   7.50 -   10.00 - green
    '#006500', #  10.00 -   15.00 - dark green
    '#fefe32', #  15.00 -   20.00 - yellow
    '#fecb00', #  20.00 -   25.00 - light orange
    '#fe9800', #  25.00 -   30.00 - orange
    '#fe6500', #  30.00 -   40.00 - dark orange
    '#fe0000', #  40.00 -   50.00 - red
    '#fe0098', #  50.00 -   75.00 - pink
    '#9832cb', #  75.00 -  100.00 - light purple
    '#650098', # 100.00 -  150.00 - purple
    '#989898'  # 150.00 -  250.00 - gray
    ]
cmap = mpl.colors.ListedColormap(colors)
bounds = [ 0.00, 0.10, 0.50, 1.00, 2.50, 5.00, 7.50, 10.00, 15.00, 20.00, 25.00, 30.00, 40.00, 50.00, 75.00, 100.00, 150.00, 250.00]
norm = mpl.colors.BoundaryNorm(bounds, len(colors))



# -------------------------------------------------------------------------
# Create figure object
# -------------------------------------------------------------------------
pngfile = 'test-bounding-box.png'
ifig += 1
iplot = 0
print('     Plot - Fig ', ifig, ' ::  ',pngfile)
fig = plt.figure(ifig)


# -------------------------------------------------------------------------
# Create map on figure
# -------------------------------------------------------------------------

point_types = ['single','multi']
for ipoint_type,point_type in enumerate(point_types):
    for iibbox_buffer,ibbox_buffer in enumerate(bbox_buffer):

        var = data[point_type][ibbox_buffer]['var']
        lat = data[point_type][ibbox_buffer]['lat']
        lon = data[point_type][ibbox_buffer]['lon']

        ntime = np.shape(var)[0]
        nlat = np.shape(var)[1]
        nlon = np.shape(var)[2]

        # boundaries between lats and lons
        lonh = np.empty((nlat+1,nlon+1), dtype=float)
        lath = np.empty((nlat+1,nlon+1), dtype=float)

        tmp = [ [ (lat[ii+1,jj+1]-lat[ii,jj])/2 for jj in range(nlon-1) ] + [ (lat[ii+1,nlon-1]-lat[ii,nlon-2])/2 ] for ii in range(nlat-1) ]
        dlat = np.array(tmp + [ tmp[-1] ])

        tmp = [ [ (lon[ii+1,jj+1]-lon[ii,jj])/2 for jj in range(nlon-1) ] + [ (lon[ii+1,nlon-1]-lon[ii,nlon-2])/2 ] for ii in range(nlat-1) ]
        dlon = np.array(tmp + [ tmp[-1] ])

        lonh[0:nlat,0:nlon] = lon - dlon
        lath[0:nlat,0:nlon] = lat - dlat

        # make lat and lon one column and row wider such that all
        lonh[nlat,0:nlon] = lonh[nlat-1,0:nlon] + (lonh[nlat-1,0:nlon] - lonh[nlat-2,0:nlon])
        lath[nlat,0:nlon] = lath[nlat-1,0:nlon] + (lath[nlat-1,0:nlon] - lath[nlat-2,0:nlon])
        lonh[0:nlat,nlon] = lonh[0:nlat,nlon-1] + (lonh[0:nlat,nlon-1] - lonh[0:nlat,nlon-2])
        lath[0:nlat,nlon] = lath[0:nlat,nlon-1] + (lath[0:nlat,nlon-1] - lath[0:nlat,nlon-2])
        lonh[nlat,nlon]   = lonh[nlat-1,nlon-1] + (lonh[nlat-1,nlon-1] - lonh[nlat-2,nlon-2])
        lath[nlat,nlon]   = lath[nlat-1,nlon-1] + (lath[nlat-1,nlon-1] - lath[nlat-2,nlon-2])

        # bounding box
        llcrnrlon_raw =  np.min(lonh)
        urcrnrlon_raw =  np.max(lonh)
        llcrnrlat_raw =  np.min(lath)
        urcrnrlat_raw =  np.max(lath)
        llcrnrlon_buf =  np.min(lonh) - (np.max(lonh) - np.min(lonh))*.05
        urcrnrlon_buf =  np.max(lonh) + (np.max(lonh) - np.min(lonh))*.05
        llcrnrlat_buf =  np.min(lath) - (np.max(lath) - np.min(lath))*.05
        urcrnrlat_buf =  np.max(lath) + (np.max(lath) - np.min(lath))*.05

        llcrnrlon_buf = -76.78610120487753
        urcrnrlon_buf = -70.37920020067847
        llcrnrlat_buf =  43.243040983915165
        urcrnrlat_buf =  47.787848958385084

        idate = dates[0]

        timestep_title = idate.strftime('%d %h %Y %H:%M:%S')+' UTC'   # '02 Oct 2017 18:00:00 UTC'

        # -------------------------------------------------------------------------
        # Create map on figure
        # -------------------------------------------------------------------------
        iplot += 1

        pos = [0.92+0.35*(iibbox_buffer),1.0-(ipoint_type)*0.25-0.1,0.3,0.27]
        # print("pos = ",pos)
        sub    = fig.add_axes( pos )  # [left, bottom, width, height]

        # Basemap
        lat_1     =  (llcrnrlat_buf+urcrnrlat_buf)/2  # first  "equator"
        lat_2     =  (llcrnrlat_buf+urcrnrlat_buf)/2  # second "equator"
        lat_0     =  (llcrnrlat_buf+urcrnrlat_buf)/2  # center of the map
        lon_0     =  (llcrnrlon_buf+urcrnrlon_buf)/2  # center of the map
        bmap = Basemap(projection='merc',area_thresh=2000.,
                llcrnrlon=llcrnrlon_buf, urcrnrlon=urcrnrlon_buf, llcrnrlat=llcrnrlat_buf, urcrnrlat=urcrnrlat_buf,
                lat_1=lat_1, lat_2=lat_2, lat_0=lat_0, lon_0=lon_0,
                resolution='i') # Lambert conformal

        # plot coastlines
        bmap.drawcoastlines(linewidth=0.3)

        #bmap.drawmapboundary(color='black', fill_color=ocean_color, linewidth=0.3)
        #bmap.drawcountries(color='black', linewidth=0.3)
        bmap.fillcontinents(color='white', lake_color=ocean_color)

        # latitudes
        iparallel_minor = [ 43.0+ii*0.5 for ii in range(12) ]
        iparallel = [44.0,45.0,46.0,47,0]
        if iibbox_buffer == len(bbox_buffer)-1:
            bmap.drawparallels(iparallel_minor,
                               labels=[0,0,0,0], dashes=[1,1], linewidth=0.25, color='0.5')
            bmap.drawparallels(iparallel,
                               labels=[0,1,0,0], dashes=[1,1], linewidth=0.25, color='0.5')
        else:
            bmap.drawparallels(iparallel_minor,
                               labels=[0,0,0,0], dashes=[1,1], linewidth=0.25, color='0.5')
            bmap.drawparallels(iparallel,
                               labels=[0,0,0,0], dashes=[1,1], linewidth=0.25, color='0.5')


        # longitudes
        imeridian_minor = [ -77.0+ii*0.5 for ii in range(14) ]
        imeridian = [-75.5,-73.5,-71.5]
        if ipoint_type == len(point_types)-1:
            bmap.drawmeridians(imeridian_minor,
                               labels=[0,0,0,0], dashes=[1,1], linewidth=0.25, color='0.5')
            bmap.drawmeridians(imeridian,
                               labels=[0,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5')
        else:
            bmap.drawmeridians(imeridian_minor,
                               labels=[0,0,0,0], dashes=[1,1], linewidth=0.25, color='0.5')
            bmap.drawmeridians(imeridian,
                               labels=[0,0,0,0], dashes=[1,1], linewidth=0.25, color='0.5')

        # geo-referenced
        # xx, yy = bmap(lon,lat)
        xxh, yyh = bmap(lonh,lath)
        zz = var[0,:,:]
        variable_plot = bmap.pcolor(xxh, yyh, zz, cmap=cmap, norm=norm, zorder = 100, alpha=0.7)

        # plot locations
        labelshiftx = 0.0
        labelshifty = (urcrnrlat_raw - llcrnrlat_raw)*.005
        if (point_type == 'single'):
            locations = locations_singlepoint
        elif (point_type == 'multi'):
            locations = locations_multipoint
        else:
            raise ValueError('Dont know point type.')
        nlocation = len(locations["lon"])
        for ilocation in range(nlocation):
            locx, locy   = bmap(locations["lon"][ilocation],locations["lat"][ilocation])
            xshift, yshift = bmap(locations["lon"][ilocation]+labelshiftx,locations["lat"][ilocation]-labelshifty)
            sub.plot(locx, locy, 'ok',      markersize=5, zorder = 200)
            #sub.text(xshift, yshift, "Loc #"+str(ilocation+1) , fontsize=textsize-1, zorder = 200, horizontalalignment="center", verticalalignment="top")

        if not(bbox is None):

            xbbox = []
            ybbox = []
            xpt, ypt = bmap(bbox[point_type][ibbox_buffer]["lon"]["min"], bbox[point_type][ibbox_buffer]["lat"]["min"])
            xbbox.append(xpt)
            ybbox.append(ypt)
            xpt, ypt = bmap(bbox[point_type][ibbox_buffer]["lon"]["max"], bbox[point_type][ibbox_buffer]["lat"]["min"])
            xbbox.append(xpt)
            ybbox.append(ypt)
            xpt, ypt = bmap(bbox[point_type][ibbox_buffer]["lon"]["max"], bbox[point_type][ibbox_buffer]["lat"]["max"])
            xbbox.append(xpt)
            ybbox.append(ypt)
            xpt, ypt = bmap(bbox[point_type][ibbox_buffer]["lon"]["min"], bbox[point_type][ibbox_buffer]["lat"]["max"])
            xbbox.append(xpt)
            ybbox.append(ypt)
            xpt, ypt = bmap(bbox[point_type][ibbox_buffer]["lon"]["min"], bbox[point_type][ibbox_buffer]["lat"]["min"])
            xbbox.append(xpt)
            ybbox.append(ypt)
            bmap.plot(xbbox, ybbox, color='black', linewidth=1., zorder = 400)

        #
        if iibbox_buffer == 0:
            if point_type == 'single':
                lab = 'point (occurrence)'
            elif point_type == 'multi':
                lab = 'multipoint (event)'
            else:
                raise ValueError('Dont know point type.')

            label = sub.text(-0.02,0.5, lab,
                 transform=sub.transAxes,
                 rotation=90, fontsize=textsize-1,
                 horizontalalignment='right', verticalalignment='center')

        # set title as time step
        if ipoint_type == 0:
            sub.set_title("bbox_buffer = "+str(ibbox_buffer),fontsize=textsize-1)





# -------------------------------------------------------------------------
# Done.
# -------------------------------------------------------------------------

fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
plt.close(fig)
