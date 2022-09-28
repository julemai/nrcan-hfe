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
import glob as glob
import shutil
from pathlib import Path  # Pathlib docu: https://docs.python.org/3/library/pathlib.html

# -----------------------
# add subolder scripts/lib to search path
# -----------------------
import sys
import os
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path+'/../../src')


zipfiles = np.sort(glob.glob(str(Path(dir_path+'/../../data/output/analyse_event_*.zip'))))

results = {}
nfiles = 0
for zipfile in zipfiles:

    unpacked = False

    extract_dir = Path(zipfile).parent                                          # dir_path+'/../../data/output
    unzippedfoldername = Path(zipfile).parent.joinpath(Path(zipfile).stem)      # dir_path+'/../../data/output/analyse_event_123

    # unzip file
    if not( unzippedfoldername.exists() ):
        shutil.unpack_archive(zipfile, extract_dir)
        unpacked = True

    # find JSON file with results (there should be only one)
    jsonfile = glob.glob(str(unzippedfoldername)+'/*.json')
    # print("jsonfile = ",jsonfile)

    if len(jsonfile) != 1:
        exceptionfile = glob.glob(str(unzippedfoldername)+'/exception.token')
        if len(jsonfile) == 0:
            print("Event {} probably still processing. Skip for now.".format(str(unzippedfoldername)))
            # remove unzipped files and folder if they were created here
            if unpacked:
                shutil.rmtree(unzippedfoldername)
            continue
        elif len(exceptionfile) == 1:
            print("Event {} not analysed because too long.".format(str(unzippedfoldername)))
            # remove unzipped files and folder if they were created here
            if unpacked:
                shutil.rmtree(unzippedfoldername)
            continue
        else:
            raise ValueError("Folder: {}\nNumber of JSON files found is {} which is not the expected number of 1.\n".format(str(unzippedfoldername),len(jsonfile)))
    else:
        # collect results
        ff = jsonfile[0]
        ifeature_event_id = ff.split('/')[-1].split('_')[1]     # 'e76ba3f9-9d8e-40f3-9199-6da1579e00d6'
        ifeature_idx  = ff.split('/')[-2].split('_')[2]         # '13'

        tmp = {}
        tmp['event_id'] = ifeature_event_id

        with open(ff, 'r') as ifile:
            tmp['results'] = json.load(ifile)

        results[int(ifeature_idx)] = tmp

        # count files evaluated
        nfiles += 1

    # remove unzipped files and folder if they were created here
    if unpacked:
        shutil.rmtree(unzippedfoldername)

print("Results for {} multi-point event zip-files found and processed.".format(nfiles))

features_idx = np.sort(list(results.keys()))





# # --------------
# # read all results files
# # --------------

# files = glob.glob(dir_path+'/../../data/output/analyse_event_*/*.json')
# nfiles = len(files)
# print("Results for {} multi-point features found.".format(nfiles))

# results = {}
# for ff in files:

#     ifeature_event_id = ff.split('/')[-1].split('_')[1]
#     ifeature_idx  = ff.split('/')[-2].split('_')[2]

#     tmp = {}
#     tmp['uuid'] = ifeature_event_id

#     with open(ff, 'r') as ifile:
#         tmp['results'] = json.load(ifile)

#     results[int(ifeature_idx)] = tmp

# features_idx = np.sort(list(results.keys()))






# --------------
# derive some stats
# --------------
missing_timesteps_n = []
missing_timesteps_perc = []
available_timesteps_n = []
available_timesteps_perc = []
accumulated_mm = []
length_event = []
for ifeature_idx in features_idx:

    nlocations = len(results[ifeature_idx]['results']['missing_timesteps_n'])

    missing_timesteps_n.append(     results[ifeature_idx]['results']['missing_timesteps_n'])
    missing_timesteps_perc.append(  results[ifeature_idx]['results']['missing_timesteps_%'])
    available_timesteps_n.append(   results[ifeature_idx]['results']['available_timesteps_n'])
    available_timesteps_perc.append(results[ifeature_idx]['results']['available_timesteps_%'])
    accumulated_mm.append(          results[ifeature_idx]['results']['accumulated_mm'])

    dates = [ np.sort(list(results[ifeature_idx]['results']['available_timesteps_precip_mm/dt'][iloc].keys())) for iloc in range(nlocations) ]
    length_event_iloc = []
    for iloc in range(nlocations):
        if len(dates[iloc]) > 0:
            start_event = datetime.datetime.strptime(dates[iloc][0],"%Y-%m-%d %H:%M:%S")
            end_event   = datetime.datetime.strptime(dates[iloc][-1],"%Y-%m-%d %H:%M:%S")
            length_event_iloc.append((end_event-start_event).days+(end_event-start_event).seconds/60./60./24.)
        else:
            length_event_iloc.append(0.0)
    length_event.append( length_event_iloc )

# missing_timesteps_n       = np.array(missing_timesteps_n)
# missing_timesteps_perc    = np.array(missing_timesteps_perc)
# available_timesteps_n     = np.array(available_timesteps_n)
# available_timesteps_perc  = np.array(available_timesteps_perc)
# accumulated_mm            = np.array(accumulated_mm)
# length_event              = np.array(length_event)

miss="??"

no_precip_event_found = [ np.where((np.array(available_timesteps_n[iloc]) == 0) & (np.array(missing_timesteps_n[iloc]) == 0))[0] for iloc in range(nlocations) ]
precip_small = [ np.where((np.array(accumulated_mm[iloc]) < 10.0))[0] for iloc in range(nlocations) ]
precip_large = [ np.where((np.array(accumulated_mm[iloc]) > 1000.0))[0] for iloc in range(nlocations) ]
nlocations_all = np.sum([ len(missing_timesteps_n[iloc]) for iloc in range(nlocations) ])

print("")
print("\\item number of available time steps: {} of {} ({:.4f} \%)".format(
    np.sum([ np.sum(available_timesteps_n[iloc]) for iloc in range(nlocations) ]),
    np.sum([ np.sum(available_timesteps_n[iloc]) for iloc in range(nlocations) ])+np.sum([ np.sum(missing_timesteps_n[iloc]) for iloc in range(nlocations) ]),
    np.sum([ np.sum(available_timesteps_n[iloc]) for iloc in range(nlocations) ])/(np.sum([ np.sum(available_timesteps_n[iloc]) for iloc in range(nlocations) ])+np.sum([ np.sum(missing_timesteps_n[iloc]) for iloc in range(nlocations) ]))*100.))
print("\\item number of missing time steps: {} of {} ({:.4f} \%)".format(
    np.sum([ np.sum(missing_timesteps_n[iloc]) for iloc in range(nlocations) ]),
    np.sum([ np.sum(available_timesteps_n[iloc]) for iloc in range(nlocations) ])+np.sum([ np.sum(missing_timesteps_n[iloc]) for iloc in range(nlocations) ]),
    np.sum([ np.sum(missing_timesteps_n[iloc]) for iloc in range(nlocations) ])/(np.sum([ np.sum(available_timesteps_n[iloc]) for iloc in range(nlocations) ])+np.sum([ np.sum(missing_timesteps_n[iloc]) for iloc in range(nlocations) ]))*100.))
print("\\item precipitation sum below 10~mm for {} of {} locations across {} features".format(np.sum([len(precip_small[iloc]) for iloc in range(nlocations)]),nlocations_all,len(features_idx)))
print("\\item precipitation sum above 1000~mm for {} of {} locations across {} features \\\\(all have multi-year period specified in HFE database)\\\\".format(np.sum([len(precip_large[iloc]) for iloc in range(nlocations)]),nlocations_all,len(features_idx)))
for iii,ii in enumerate(precip_large):
    if len(precip_large[iii]) > 1:
        print("   {{\\scriptsize Event ID: {} $\curvearrowright$ Start and end date = [{},{}]}}\\\\[-4pt]".format(
            #features_idx[ii],
            results[features_idx[iii]]['event_id'],
            results[features_idx[iii]]['results']['start_date_w_buffer'],
            results[features_idx[iii]]['results']['end_date_w_buffer']))
print("\\item no precipitation event found for {} of {} locations across {} features\\\\".format(np.sum([len(no_precip_event_found[iloc]) for iloc in range(nlocations)]),nlocations_all,len(features_idx)))
for iii,ii in enumerate(no_precip_event_found):
    if len(no_precip_event_found[iii]) > 1:
        print("   {{\\scriptsize Event ID: {} $\curvearrowright$ Start and end date = [{},{}]}}\\\\[-4pt]".format(
            #features_idx[ii],
            results[features_idx[iii]]['event_id'],
            results[features_idx[iii]]['results']['start_date_w_buffer'],
            results[features_idx[iii]]['results']['end_date_w_buffer']))
print("")
print("")









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
pngfile = 'test-analyse-event-all-accum-precip.png'
ifig += 1
iplot = 0
print('     Plot - Fig ', ifig, ' ::  ',pngfile)
fig = plt.figure(ifig)


# -------------------------------------------------------------------------
# Create map on figure
# -------------------------------------------------------------------------

iplot += 1

pos = [0.0,0.0,1.0,0.3]  # [left, bottom, width, height]
# print("pos = ",pos)
sub    = fig.add_axes( pos )  # [left, bottom, width, height]

#sub.plot(np.cum(),color='0.4',linewidth=1.5*lwidth,linestyle='--',label=label,zorder=100)

flatten_accumulated_mm = np.array([item for sublist in accumulated_mm for item in sublist])

min_val = 0.0 #np.percentile(flatten_accumulated_mm,1.0)
max_val = np.percentile(flatten_accumulated_mm,99.3)
filtered = flatten_accumulated_mm[ (flatten_accumulated_mm>min_val) & (flatten_accumulated_mm<max_val) ]
sub.hist(flatten_accumulated_mm, 50000, density=True,
             cumulative=True, label='CDF',  #normed=True,
             histtype='step',
             alpha=0.8, color=ocean_color,
             zorder = 200,
             linewidth=0.8*lwidth
             )

xmin, xmax = sub.get_xlim()
ymin, ymax = sub.get_ylim()
plt.setp(sub, xlim=[min_val,max_val*0.99], ylim=[ymin,ymax])

# axis lables
sub.set_xlabel("Accumulated precipitation of identified precipitation event [mm]",fontsize=textsize)
sub.set_ylabel("Cumulative density function [-]",fontsize=textsize)

# percentiles
percentiles = [5,50,95,99]
for ipercentile in percentiles:
    sub.text(np.percentile(flatten_accumulated_mm,ipercentile)+5., ipercentile/100., "$p_{{{:d}}} = {:.2f}$ mm".format(ipercentile,np.percentile(flatten_accumulated_mm,ipercentile)) ,
             #transform=sub.transAxes,
             fontsize=textsize-5, zorder = 50,
             horizontalalignment="left",
             verticalalignment="top")
    sub.plot([min_val,np.percentile(flatten_accumulated_mm,ipercentile),np.percentile(flatten_accumulated_mm,ipercentile)],[ipercentile/100.,ipercentile/100.,0.0],
                 color='0.8',
                 linestyle='--',
                 linewidth=0.6*lwidth)

# nfeatures
sub.text(0.98, 0.02, "Features displayed: {}\nFeatures analysed: {}".format(len(filtered),len(flatten_accumulated_mm)) ,
             transform=sub.transAxes,
             fontsize=textsize-1, zorder = 200,
             horizontalalignment="right",
             verticalalignment="bottom")

# set title as time step
sub.set_title("Results of analyzing {} locations across {} valid events of HFE database".format(nlocations_all,len(features_idx)),fontsize=textsize+2)





# -------------------------------------------------------------------------
# Done.
# -------------------------------------------------------------------------

fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
plt.close(fig)





# -------------------------------------------------------------------------
# Create figure object
# -------------------------------------------------------------------------
pngfile = 'test-analyse-event-all-len-precip-event.png'
ifig += 1
iplot = 0
print('     Plot - Fig ', ifig, ' ::  ',pngfile)
fig = plt.figure(ifig)


# -------------------------------------------------------------------------
# Create map on figure
# -------------------------------------------------------------------------

iplot += 1

pos = [0.0,0.0,1.0,0.3]  # [left, bottom, width, height]
# print("pos = ",pos)
sub    = fig.add_axes( pos )  # [left, bottom, width, height]

#sub.plot(np.cum(),color='0.4',linewidth=1.5*lwidth,linestyle='--',label=label,zorder=100)

flatten_length_event = np.array([item for sublist in length_event for item in sublist])

min_val = 0.0 #np.percentile(flatten_length_event,1.0)
max_val = np.percentile(flatten_length_event,99.3)
filtered = flatten_length_event[ (flatten_length_event>min_val) & (flatten_length_event<max_val) ]
sub.hist(flatten_length_event, 50000, density=True,
             cumulative=True, label='CDF',  #normed=True,
             histtype='step',
             alpha=0.8, color=ocean_color,
             zorder = 200,
             linewidth=0.8*lwidth
             )

xmin, xmax = sub.get_xlim()
ymin, ymax = sub.get_ylim()
plt.setp(sub, xlim=[min_val,max_val*0.99], ylim=[ymin,ymax])

# axis lables
sub.set_xlabel("Length identified precipitation event [days]",fontsize=textsize)
sub.set_ylabel("Cumulative density function [-]",fontsize=textsize)

# percentiles
percentiles = [5,50,95,99]
for ipercentile in percentiles:
    sub.text(np.percentile(flatten_length_event,ipercentile)+0.5, ipercentile/100., "$p_{{{:d}}} = {:.1f}$ days".format(ipercentile,np.percentile(flatten_length_event,ipercentile)) ,
             #transform=sub.transAxes,
             fontsize=textsize-5, zorder = 50,
             horizontalalignment="left",
             verticalalignment="top")
    sub.plot([min_val,np.percentile(flatten_length_event,ipercentile),np.percentile(flatten_length_event,ipercentile)],[ipercentile/100.,ipercentile/100.,0.0],
                 color='0.8',
                 linestyle='--',
                 linewidth=0.6*lwidth)

# nfeatures
sub.text(0.98, 0.02, "Features displayed: {}\nFeatures analysed: {}".format(len(filtered),len(flatten_length_event)) ,
             transform=sub.transAxes,
             fontsize=textsize-1, zorder = 200,
             horizontalalignment="right",
             verticalalignment="bottom")

# set title as time step
sub.set_title("Results of analyzing {} locations across {} valid events of HFE database".format(nlocations_all,len(features_idx)),fontsize=textsize+2)





# -------------------------------------------------------------------------
# Done.
# -------------------------------------------------------------------------

fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
plt.close(fig)
