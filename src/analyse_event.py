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
import warnings

from a1_request_geomet_grib2 import request_geomet_grib2
from a2_request_caspar_nc import request_caspar_nc
from b1_read_geomet_grib2 import read_geomet_grib2
from b3_read_hfe_json import read_hfe_json

from cx_plot_data import plot_data
from dx_interpolate_data import interpolate_data
from ex_determine_bbox import determine_bbox
from fx_determine_dates import determine_dates


def plot_interpolated(locations=None,dates=None,data=None,pngfile=None,start_date=None,end_date=None):

    # checking inputs
    if (locations is None):
        raise ValueError("plot_interpolated: location needs to be specified")
    if (dates is None):
        raise ValueError("plot_interpolated: dates needs to be specified")
    if (data is None):
        raise ValueError("plot_interpolated: data needs to be specified")
    if (pngfile is None):
        raise ValueError("plot_interpolated: pngfile needs to be specified")

    import matplotlib as mpl
    import matplotlib.pyplot as plt


    # -------------------------------------------------------------------------
    # create plot
    # -------------------------------------------------------------------------

    # plot settings
    textsize    = 12          # standard text size
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

    llxbbox     = 0.5         # x-anchor legend bounding box
    llybbox     = 1.0        # y-anchor legend bounding box
    llrspace    = 0.25        # spacing between rows in legend
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
    ifig += 1
    iplot = 0
    print('     Plot - Fig ', ifig, ' ::  ',pngfile)
    fig = plt.figure(ifig)

    # -------------------------------------------------------------------------
    # Create line plots on figure
    # -------------------------------------------------------------------------

    sub    = fig.add_axes( [0.0,0.0,1.0,0.25] )  # [left, bottom, width, height]

    # 6h accumlations --> for stepwise plot divide by this value
    timedelta = (dates[1]-dates[0]).seconds/60/60

    label = [ "Loc #"+str(iloc+1)+"  ("+str(locations["lat"][iloc])+","+str(locations["lon"][iloc])+")" for iloc in range(len(locations["lat"])) ]
    sub.step(dates,data['var']/timedelta,color='0.4',linewidth=1.0*lwidth,linestyle='-',zorder=100) #label=label,

    xmin, xmax = sub.get_xlim()
    ymin, ymax = sub.get_ylim()

    # -------------------------------------------------------------------------
    # Plot some vertical lines to visualize time period requested (incl. buffer) and start and end date of event
    # -------------------------------------------------------------------------
    sub.plot([dates[0],dates[0]],[ymin,ymax],label='Start date incl. buffer',color='0.8',linewidth=1.0*lwidth,linestyle='--',zorder=50)

    if not(start_date is None):
        sub.plot([start_date,start_date],[ymin,ymax],label='Event start date',color='0.8',linewidth=1.0*lwidth,linestyle='dotted',zorder=50)
    if not(end_date is None):
        sub.plot([end_date,end_date],[ymin,ymax],label='Event end date',color='0.8',linewidth=1.0*lwidth,linestyle='dotted',zorder=50)

    sub.plot([dates[-1],dates[-1]],[ymin,ymax],label='End date incl. buffer',color='0.8',linewidth=1.0*lwidth,linestyle='--',zorder=50)

    # -------------------------------------------------------------------------
    # Format date axis
    # -------------------------------------------------------------------------
    import matplotlib.dates as mdates
    from matplotlib.ticker import FormatStrFormatter
    myFmt = mdates.DateFormatter("%d %h\n%H:%M")
    sub.xaxis.set_major_formatter(myFmt)
    sub.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    # -------------------------------------------------------------------------
    # Label on axes
    # -------------------------------------------------------------------------
    plt.setp(sub,  ylabel="Precipitation Rate [mm h$^{-1}$]")   # [mm/6h]

    # -------------------------------------------------------------------------
    # plot legend
    # -------------------------------------------------------------------------
    ll = sub.legend(frameon=frameon, ncol=4,
                            labelspacing=llrspace, handletextpad=llhtextpad, handlelength=llhlength,
                            loc='lower center', bbox_to_anchor=(llxbbox,llybbox), scatterpoints=1, numpoints=1,
                            fontsize = 'small')

    # -------------------------------------------------------------------------
    # Done.
    # -------------------------------------------------------------------------

    fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
    plt.close(fig)


# --------------------
# SPECIFICS OF WHAT TO ANALYSE
# --------------------
ifeature     = 320 # starts 2020-11-30 (using Geomet data)
bbox_buffer  = 0.5
dates_buffer = [1.0,2.0]

# --------------------
# Load HFE database
# --------------------
print("\n\nReading HFE database")

# ignore warnings produced
warnings.filterwarnings("ignore", category=UserWarning, message='b3_read_hfe_json')

filename        = '../data/hfe/historical_flood_event.json'
filtering       = True
polygon         = None
return_filtered = False
silent          = True

data_hfe = read_hfe_json(filename=filename,filtering=filtering,polygon=polygon,return_filtered=return_filtered,silent=silent)
nfeatures = len(data_hfe['data']['features'])
print("   Number of flood occurrences/events found: {}".format(nfeatures))

# for iidate,idata in enumerate(data_hfe['data']['features']):
#     print(iidate, idata['properties']['start_date'])

# --------------------
# Pick one event/occurrence
# --------------------
if ifeature >= nfeatures or ifeature < 0:
    raise ValueError("Feature needs to be between {} and {}.".format(0,nfeatures))
print("\n\nAnalysing event #{}:".format(ifeature))

feature = data_hfe['data']['features'][ifeature]
start_date = datetime.datetime(int(feature['properties']['start_date'][0:4]),int(feature['properties']['start_date'][5:7]),int(feature['properties']['start_date'][8:10]),0,0)
if not(feature['properties']['end_date'] is None):
    end_date = datetime.datetime(int(feature['properties']['end_date'][0:4]),int(feature['properties']['end_date'][5:7]),int(feature['properties']['end_date'][8:10]),0,0)
else:
    end_date = None

print("   Type: {}\n   Number of locations: {}\n   OBJECTID: {}\n   Start: {}\n   End: {}\n   Cause: {}".format(
    feature['geometry']['type'],
    len(feature['geometry']['coordinates']),
    feature['properties']['OBJECTID'],
    feature['properties']['start_date'],
    feature['properties']['end_date'],
    feature['properties']['flood_cause'],
    ))


# --------------------
# Determine bounding box
# --------------------
print("\n\nDetermine parameters for data requests:")

silent      = True

bbox = determine_bbox(feature=feature,bbox_buffer=bbox_buffer,silent=silent)
print("   bbox : {}".format(bbox))


# --------------------------------------
# Determine time steps for rdpa:10km:6f
# --------------------------------------
if (start_date - datetime.datetime(1980,1,1)).days < 0:
    raise ValueError("analyse_event: start date is '{}' which is before 1980-01-01 which is when precip data will become available".format(
        feature['properties']['start_date']))
elif (start_date - datetime.datetime(2018,1,1)).days < 0:
    product  = 'RDRS_v2.1'
else:
    product  = 'rdpa:10km:6f'
silent       = True

date = determine_dates(feature=feature,product=product,dates_buffer=dates_buffer,silent=silent)
print("   date : [ {}, {}, ..., {}, {} ] (in total {} time steps)".format(date[0], date[1],date[-2],date[-1],len(date)))


# --------------------
# Request data
# --------------------
print("\n\nRequest data:")

if (product == 'rdpa:10km:6f') or (product == 'rdpa:10km:24f'):

    # ignore warnings produced
    warnings.filterwarnings("ignore", category=UserWarning, message='request_geomet_grib2')

    filename  = '/tmp/tmp/analyse_event_'+str(ifeature)+'/geomet'
    crs       = 'EPSG:4326'
    overwrite = False
    silent    = True
    file_geomet = request_geomet_grib2(product=product,date=date,bbox=bbox,crs=crs,filename=filename,overwrite=overwrite,silent=silent)

    nfiles = len(file_geomet)
    print("   Number of Geomet files downloaded: {}".format(nfiles))

elif (product == 'RDRS_v2.1'):

    file_caspar = request_caspar_nc( )

    nfiles = 0
    print("   Number of CaSPAr files downloaded: {}".format(nfiles))

else:

    raise ValueError("analyse_event: product not known")


# --------------------
# Read data
# --------------------
print("\n\nRead data:")

if (product == 'rdpa:10km:6f'):

    lintransform = {'a':1.0,'b':0.0}
    silent       = True
    data = read_geomet_grib2(filenames=file_geomet,lintransform=lintransform,silent=silent)

    ntime = np.shape(data['var'])[0]
    nlat = np.shape(data['lat'])[0]
    nlon = np.shape(data['lon'])[1]
    print("   Number of time steps read: {}".format(ntime))
    print("   Number of latitudes  read: {}".format(nlat))
    print("   Number of longitudes read: {}".format(nlon))

elif (product == 'RDRS_v2.1'):

    file_caspar = read_caspar_nc( )   # TODO

    nfiles = 0
    print("   Number of CaSPAr files downloaded: {}".format(nfiles))

else:

    raise ValueError("analyse_event: product not known")



# --------------------
# Interpolate data at all locations of event
# --------------------
print("\n\nInterpolate data:")


locations = {'lon':np.array(feature['geometry']['coordinates'])[:,0],'lat':np.array(feature['geometry']['coordinates'])[:,1]}
var       = data["var"]
lat       = data["lat"]
lon       = data["lon"]
interpolated_data = interpolate_data(var=var,lat=lat,lon=lon,locations=locations,bbox=bbox,post_process=True,silent=True)
print("   Sum of precipitation at all {} locations over the time period evaluated: {}".format(len(feature['geometry']['coordinates']),np.sum(interpolated_data['var'],axis=0)))

plot_interpolated(locations,date[0:-1],interpolated_data,pngfile='/tmp/tmp/analyse_event_'+str(ifeature)+'/interpolated_at_stations.png',start_date=start_date,end_date=end_date)




# --------------------
# Plot data
# --------------------
print("\n\nPlot data:")

var          = data["var"]
lat          = data["lat"]
lon          = data["lon"]
date         = data["time"]
png          = True
gif          = True
legend       = True
cities       = True
basefilename = '/tmp/tmp/analyse_event_'+str(ifeature)+'/plot_event_'+str(ifeature)
overwrite    = False
silent       = True

plots_data = plot_data(var=var,lat=lat,lon=lon,date=date,
                           png=png,
                           gif=gif,
                           legend=legend,
                           cities=cities,
                           bbox=bbox,
                           basefilename=basefilename,
                           overwrite=overwrite,
                           silent=silent)

print("   Number of PNGs    plotted: {}".format(len(plots_data['png'])))
print("   Number of GIFs    plotted: {}".format(len(plots_data['gif'])))
print("   Number of Legends plotted: {}".format(len(plots_data['legend'])))
