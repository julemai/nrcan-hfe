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
from dx_interpolate_data import interpolate_data


# -----------------------
# data and locations shown in plot
# -----------------------

product = 'rdpa:10km:6f'
crs = 'EPSG:4326'
bbox = {"lat":{"min":45.0,"max":46.0},"lon":{"min":-74.0,"max":-73.0}}
dates = [ datetime.datetime(2018,8,7,12,0) + datetime.timedelta(hours=ii) for ii in range(0,49,6) ]
lintransform={'a':1.0,'b':0.0}  # no convert of units

# request data
filename = '/tmp/test-map-geomet-nrcan-hfe-rdpa6'
files_geomet = request_geomet_grib2(product=product,date=dates,bbox=bbox,crs=crs,overwrite=False,filename=filename)

# read data
data_geomet = read_geomet_grib2(
    files_geomet,
    lintransform=lintransform)

# interpolate at three locations
locations = {"lat":[46.1,45.4], "lon":[-73.5,-72.8]}
interpolate_geomet = interpolate_data(var=data_geomet["var"],lat=data_geomet["lat"],lon=data_geomet["lon"],locations=locations,return_tmp=True,silent=True) # bbox=bbox,








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
pngfile = 'test-bilinear-interpolation.png'
ifig += 1
iplot = 0
print('     Plot - Fig ', ifig, ' ::  ',pngfile)
fig = plt.figure(ifig)

# -------------------------------------------------------------------------
# Create line plots on figure
# -------------------------------------------------------------------------
nlocation = len(locations["lat"])
for ilocation in range(nlocation):

    iplot += 1
    sub    = fig.add_axes( [0.0,1.0-0.415*ilocation-0.44,0.85,0.31] )  # [left, bottom, width, height]

    # 6h accumlations --> for stepwise plot divide by this value
    timedelta = (dates[1]-dates[0]).seconds/60/60

    # plot time series
    ntime = len(dates)
    for ii in range(3):
        for jj in range(3):

            if ii == 1 and jj == 1:

                # plot time series of center cell
                ts = np.array([ interpolate_geomet['tmp']['var'][ilocation][itime][ii][jj] for itime in range(ntime) ])
                label = 'center grid cell (lat={},lon={})'.format(round(interpolate_geomet['tmp']['lat'][ilocation][ii][jj],2),round(interpolate_geomet['tmp']['lon'][ilocation][ii][jj],2))
                sub.step(dates,ts/timedelta,color='0.4',linewidth=1.5*lwidth,linestyle='--',label=label,zorder=100)

            else:

                # eight neighboring cells
                ts = np.array([ interpolate_geomet['tmp']['var'][ilocation][itime][ii][jj] for itime in range(ntime) ])
                if ii==0 and jj==0:
                    label = 'eight neighboring grid cells'
                else:
                    label = ''
                sub.step(dates,ts/timedelta,color='0.8',label=label,zorder=40)

    # plot time series of interpolated location
    ts = interpolate_geomet['var'][:,ilocation]
    label = 'interpolated at location #{}'.format(ilocation+1)
    sub.step(dates,ts/timedelta,color='0.1',linewidth=1.5*lwidth,label=label,zorder=120)

    # set title as time step
    sub.set_title('Location #{}: (lat={},lon={})'.format(ilocation+1,locations["lat"][ilocation],locations["lon"][ilocation]),fontsize=textsize+2)

    # format date axis
    import matplotlib.dates as mdates
    from matplotlib.ticker import FormatStrFormatter
    myFmt = mdates.DateFormatter("%d %h\n%H:%M")
    sub.xaxis.set_major_formatter(myFmt)
    sub.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    #sub.xaxis.set_major_locator(mdates.HourLocator(interval=6))

    # plt.setp(sub,  xlabel="Date in 2018 [MM-DD HH:MM]") # axis labels
    plt.setp(sub,  ylabel="Precipitation Rate [mm h$^{-1}$]")   # [mm/6h]

    # set range of plot
    xmin, xmax = sub.get_xlim()
    ymin, ymax = sub.get_ylim()
    if ilocation ==0:
        plt.setp(sub, xlim=[xmin,xmax], ylim=[ymin,2.0])
    if ilocation ==1:
        plt.setp(sub, xlim=[xmin,xmax], ylim=[ymin,5.2])

    # plot legend
    ll = sub.legend(frameon=frameon, ncol=1,
                            labelspacing=llrspace, handletextpad=llhtextpad, handlelength=llhlength,
                            loc='upper left', bbox_to_anchor=(llxbbox,llybbox), scatterpoints=1, numpoints=1,
                            fontsize = 'small')


# -------------------------------------------------------------------------
# Create map on figure
# -------------------------------------------------------------------------

cities = False

var = data_geomet['var']
lat = data_geomet['lat']
lon = data_geomet['lon']

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

# load cities if needed
if cities:

    with open(dir_path+'/../../src/lib/canadacities.json', 'r') as ff:    # w/  population
        dict_cities = json.load(ff)

    # find all cities in bounding box
    city_inregion = {}
    for cc in dict_cities:
        if (    (dict_cities[cc]["lat"] > llcrnrlat_raw) and
                (dict_cities[cc]["lat"] < urcrnrlat_raw) and
                (dict_cities[cc]["lon"] > llcrnrlon_raw) and
                (dict_cities[cc]["lon"] < urcrnrlon_raw)):
            city_inregion[cc] = dict_cities[cc]
    nmaxcities = 10
    if len(city_inregion) > nmaxcities:
        icities = np.array([ cc for cc in city_inregion ])
        ipopulation = np.array([ dict_cities[cc]["population"] for cc in city_inregion ])
        idx_largest = np.argsort(ipopulation)[::-1][0:nmaxcities]

        # thinning out cities to plot
        city_inregion = { cc:city_inregion[cc] for cc in icities[idx_largest] }

    print("Found {} cities in region.".format(len(city_inregion)))

date = dates[0:]
for iidate,idate in enumerate(date):

    timestep_title = idate.strftime('%d %h %Y %H:%M:%S')+' UTC'   # '02 Oct 2017 18:00:00 UTC'

    # -------------------------------------------------------------------------
    # Create map on figure
    # -------------------------------------------------------------------------
    iplot += 1

    pos = [0.92+0.35*(iidate//3),1.0-(iidate%3+1)*0.27-0.1,0.3,0.27]
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

    # latitudes - draw 3 parallels, labels = [left, right, top, bottom]
    nparallels = 3.
    if iidate//3+1 == 3:
        # iparallel = np.arange(np.round(np.min(lath),1),np.round(np.max(lath),1)+0.1,np.round((np.max(lath)+0.1-np.min(lath))/nparallels,1))
        iparallel = [44.7,45.5, 46.3]
        bmap.drawparallels(iparallel,
                           labels=[0,1,0,0], dashes=[1,1], linewidth=0.25, color='0.5')
    else:
        bmap.drawparallels(np.arange(np.round(np.min(lath),1),np.round(np.max(lath),1)+0.1,np.round((np.max(lath)+0.1-np.min(lath))/nparallels,1)),
                           labels=[0,0,0,0], dashes=[1,1], linewidth=0.25, color='0.5')


    # longitudes - draw 3 meridians, labels = [left, right, top, bottom]
    nmeridians = 2.
    if iidate%3+1 == 3:
        # imeridian = np.arange(np.round(np.min(lonh),1),np.round(np.max(lonh),1)+0.1,np.round((np.max(lonh)+0.1-np.min(lonh))/nmeridians,1))
        imeridian = [-74.5,-73.5,-72.5]
        bmap.drawmeridians(imeridian,
                           labels=[0,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5')
    else:
        bmap.drawmeridians(np.arange(np.round(np.min(lonh),1),np.round(np.max(lonh),1)+0.1,np.round((np.max(lonh)+0.1-np.min(lonh))/nmeridians,1)),
                           labels=[0,0,0,0], dashes=[1,1], linewidth=0.25, color='0.5')

    # geo-referenced
    # xx, yy = bmap(lon,lat)
    xxh, yyh = bmap(lonh,lath)
    zz = var[iidate,:,:]
    variable_plot = bmap.pcolor(xxh, yyh, zz, cmap=cmap, norm=norm, zorder = 100, alpha=0.7)

    # plot cities
    if cities:

        labelshiftx = 0.0 #(urcrnrlon_raw - llcrnrlon_raw)*.01
        labelshifty = (urcrnrlat_raw - llcrnrlat_raw)*.005
        for icity in city_inregion:
            cityx, cityy           = bmap(city_inregion[icity]["lon"],city_inregion[icity]["lat"])
            xshift, yshift = bmap(city_inregion[icity]["lon"]+labelshiftx,city_inregion[icity]["lat"]-labelshifty)
            sub.plot(cityx, cityy, 'ok',      markersize=3, zorder = 200)
            sub.text(xshift, yshift, icity , fontsize=textsize-1, zorder = 200, horizontalalignment="center", verticalalignment="top")

    # plot locations
    labelshiftx = 0.0 #(urcrnrlon_raw - llcrnrlon_raw)*.01
    labelshifty = (urcrnrlat_raw - llcrnrlat_raw)*.005
    for ilocation in range(nlocation):
        locx, locy   = bmap(locations["lon"][ilocation],locations["lat"][ilocation])
        xshift, yshift = bmap(locations["lon"][ilocation]+labelshiftx,locations["lat"][ilocation]-labelshifty)
        sub.plot(locx, locy, 'ok',      markersize=3, zorder = 200)
        sub.text(xshift, yshift, "Loc #"+str(ilocation+1) , fontsize=textsize-1, zorder = 200, horizontalalignment="center", verticalalignment="top")

    if not(bbox is None):

        xbbox = []
        ybbox = []
        xpt, ypt = bmap(bbox["lon"]["min"], bbox["lat"]["min"])
        xbbox.append(xpt)
        ybbox.append(ypt)
        xpt, ypt = bmap(bbox["lon"]["max"], bbox["lat"]["min"])
        xbbox.append(xpt)
        ybbox.append(ypt)
        xpt, ypt = bmap(bbox["lon"]["max"], bbox["lat"]["max"])
        xbbox.append(xpt)
        ybbox.append(ypt)
        xpt, ypt = bmap(bbox["lon"]["min"], bbox["lat"]["max"])
        xbbox.append(xpt)
        ybbox.append(ypt)
        xpt, ypt = bmap(bbox["lon"]["min"], bbox["lat"]["min"])
        xbbox.append(xpt)
        ybbox.append(ypt)
        bmap.plot(xbbox, ybbox, color='black', linewidth=1., zorder = 400)

    # set title as time step
    sub.set_title(timestep_title,fontsize=textsize-1)





# -------------------------------------------------------------------------
# Done.
# -------------------------------------------------------------------------

fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
plt.close(fig)
