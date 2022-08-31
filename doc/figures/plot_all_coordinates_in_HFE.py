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


import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
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


# range of plot
llcrnrlat =   35.
urcrnrlat =   65.
llcrnrlon = -151.
urcrnrlon =  -45.

# settings
cities = True


from b3_read_hfe_json import read_hfe_json
from in_poly import in_poly     # in lib/


coords = {}

# -----------------------
# read HFE events from "data/hfe/historical_flood_event.json" (filtered using checks 1-3) and extract all coordinates
# -----------------------
filename='../../data/hfe/historical_flood_event.json'

# filtering w/o checking if points are in polygon
data_hfe = read_hfe_json(filename=filename,filtering=True,polygon=None,silent=True)

coords_event = []
for ifeature in data_hfe['data']['features']:

    icoord = ifeature["geometry"]["coordinates"]
    coords_event += icoord

coords["event"] = coords_event


# -----------------------
# read HFE occurrences from "data/hfe/historical_flood.json" (filtered using checks 1-3) and extract all coordinates
# -----------------------
filename='../../data/hfe/historical_flood.json'

# filtering w/o checking if points are in polygon
data_hfe = read_hfe_json(filename=filename,filtering=True,polygon=None,silent=True)

coords_occurrence = []
for ifeature in data_hfe['data']['features']:

    icoord = ifeature["geometry"]["coordinates"]
    coords_occurrence += [ icoord ]

coords["occurrence"] = coords_occurrence


# -----------------------
# load cities if needed
# -----------------------
if cities:

    with open(dir_path+'/../../src/lib/canadacities.json', 'r') as ff:    # w/  population
        dict_cities = json.load(ff)

    # find all cities in bounding box
    city_inregion = {}
    for cc in dict_cities:
        if (    (dict_cities[cc]["lat"] > llcrnrlat) and
                (dict_cities[cc]["lat"] < urcrnrlat) and
                (dict_cities[cc]["lon"] > llcrnrlon) and
                (dict_cities[cc]["lon"] < urcrnrlon)):
            city_inregion[cc] = dict_cities[cc]
    nmaxcities = 10
    if len(city_inregion) > nmaxcities:
        icities = np.array([ cc for cc in city_inregion ])
        ipopulation = np.array([ dict_cities[cc]["population"] for cc in city_inregion ])
        idx_largest = np.argsort(ipopulation)[::-1][0:nmaxcities]

        # thinning out cities to plot
        city_inregion = { cc:city_inregion[cc] for cc in icities[idx_largest] }

    print("Found {} cities in region.".format(len(city_inregion)))



# -------------------------------------------------------------------------
# create plot
# -------------------------------------------------------------------------

# plot settings
textsize    = 16          # standard text size
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

msize       = 4.0         # marker size
mwidth      = 1.5         # marker edge width

llxbbox     = 0.0         # x-anchor legend bounding box
llybbox     = 0.0         # y-anchor legend bounding box
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
pngfile = 'test-all-coordinates-in-hfe.png'
ifig += 1
iplot = 0
print('     ')
print('     Plot - Fig ', ifig, ' ::  ',pngfile)
fig = plt.figure(ifig)

# -------------------------------------------------------------------------
# Create plot with all locations
# -------------------------------------------------------------------------
nlocation_event = len(coords["event"])
nlocation_occurrence = len(coords["occurrence"])

iplot += 1
sub    = fig.add_axes( [0.0,0.0,1.0,1.0] )  # [left, bottom, width, height]

# Basemap
lat_1     =  (llcrnrlat+urcrnrlat)/2  # first  "equator"
lat_2     =  (llcrnrlat+urcrnrlat)/2  # second "equator"
lat_0     =  (llcrnrlat+urcrnrlat)/2  # center of the map
lon_0     =  (llcrnrlon+urcrnrlon)/2  # center of the map
bmap = Basemap(projection='merc',area_thresh=5000.,
        llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon, llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat,
        lat_1=lat_1, lat_2=lat_2, lat_0=lat_0, lon_0=lon_0,
        resolution='i') # Lambert conformal

# plot coastlines
bmap.drawcoastlines(linewidth=0.3)

#bmap.drawmapboundary(color='black', fill_color=ocean_color, linewidth=0.3)
#bmap.drawcountries(color='black', linewidth=0.3)
bmap.fillcontinents(color='white', lake_color=ocean_color)

# latitudes - draw 3 parallels, labels = [left, right, top, bottom]
iparallel = np.arange(0,90,10)
bmap.drawparallels(iparallel,
                   labels=[0,1,0,0], dashes=[1,1], linewidth=0.25, color='0.5')


# longitudes - draw 3 meridians, labels = [left, right, top, bottom]
imeridian = np.arange(-180,180,20)
bmap.drawmeridians(imeridian,
                    labels=[0,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5')


# plot polygon used to request data from caspar
coords_caspar = np.array([[-61.43692000000001, 42.107902], [-69.679413, 43.733895], [-75.802917, 41.948767], [-81.803513, 40.336600000000004], [-85.34317, 40.905729], [-86.717148, 43.926583], [-90.00412, 46.520603], [-98.191681, 46.843286], [-109.592743, 46.635294], [-121.171646, 46.26819], [-127.60894799999998, 47.161241], [-133.259354, 50.607646], [-136.426163, 54.110943], [-135.518814, 56.883876], [-119.99999999999999, 59.99999999999999], [-110, 59.99999999999999], [-90, 52], [-70, 52], [-52.67807, 53.008173], [-50.820694, 50.59152], [-49.156265, 47.78086600000001], [-52.49267600000001, 44.657909], [-57.761306999999995, 45.09049], [-61.43692000000001, 42.107902]])

sub.add_patch(Polygon([ bmap(ii[0],ii[1]) for ii in coords_caspar ], facecolor=colors[2], edgecolor=colors[2], linewidth=0.5, alpha=0.2, zorder=400,label='Domain requested from CaSPAr'))


# plot cities
city_to_plot = ['Toronto', 'Montreal', 'Vancouver', 'Calgary', 'Edmonton', 'Winnipeg', 'Quebec City']
if cities:

    labelshiftx = 0.0
    labelshifty = (urcrnrlat - llcrnrlat)*.005
    for icity in city_inregion:
        if icity in city_to_plot:
            cityx, cityy           = bmap(city_inregion[icity]["lon"],city_inregion[icity]["lat"])
            xshift, yshift = bmap(city_inregion[icity]["lon"]+labelshiftx,city_inregion[icity]["lat"]-labelshifty)
            sub.plot(cityx, cityy, 'ok',      markersize=3, zorder = 200)
            sub.text(xshift, yshift, icity , fontsize=textsize-5, zorder = 200, horizontalalignment="center", verticalalignment="top")

# plot locations in "event"
for iilocation,ilocation in enumerate(coords["event"]):

    ilon = ilocation[0]
    ilat = ilocation[1]
    locx, locy   = bmap(ilon,ilat)
    if iilocation == 0:
        lab = "Location of flood event (multipoint)"
        sub.plot(locx, locy, linestyle='None',
                 marker='s', markeredgecolor=colors[2], markerfacecolor='None',
                 markersize=msize, markeredgewidth=mwidth/2, label=lab)
    else:
        sub.plot(locx, locy, linestyle='None',
                 marker='s', markeredgecolor=colors[2], markerfacecolor='None',
                 markersize=msize, markeredgewidth=mwidth/2)


# plot locations in "occurrence"
for iilocation,ilocation in enumerate(coords["occurrence"]):

    ilon = ilocation[0]
    ilat = ilocation[1]
    locx, locy   = bmap(ilon,ilat)
    if iilocation == 0:
        lab = "Location of flood occurrence (point)"
        sub.plot(locx, locy, linestyle='None',
                 marker='x', markeredgecolor=colors[2], markerfacecolor='None',
                 markersize=msize, markeredgewidth=mwidth/2, label=lab)
    else:
        sub.plot(locx, locy, linestyle='None',
                 marker='x', markeredgecolor=colors[2], markerfacecolor='None',
                 markersize=msize, markeredgewidth=mwidth/2)

# plot legend
ll = sub.legend(frameon=True, ncol=1,
                            labelspacing=llrspace, handletextpad=llhtextpad, handlelength=llhlength,
                            loc='lower left', bbox_to_anchor=(llxbbox,llybbox), scatterpoints=1, numpoints=1,
                            fontsize = 'small')


# check if all points are in polygon
detected_point_outside = False

for coord in coords["event"]:
    # in_poly returns ::  1 if point inside polygon
    #                    -1 if point outside polygon
    #                     0 if point on vertex/edge of polygon
    inpoly = in_poly( coord, coords_caspar[:,0], coords_caspar[:,1] )
    if inpoly < 1:
        print("coordinate not in polygon: {}".format(coord))
        detected_point_outside = True

for coord in coords["occurrence"]:
    # in_poly returns ::  1 if point inside polygon
    #                    -1 if point outside polygon
    #                     0 if point on vertex/edge of polygon
    inpoly = in_poly( coord, coords_caspar[:,0], coords_caspar[:,1] )
    if inpoly < 1:
        print("coordinate not in polygon: {}".format(coord))
        detected_point_outside = True

if not(detected_point_outside):
    print("     >>>> Checked all points and none is outside the polygon!! Yay!!")





# -------------------------------------------------------------------------
# Done.
# -------------------------------------------------------------------------

fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
plt.close(fig)
