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
from b2_read_caspar_nc import read_caspar_nc
from b3_read_hfe_json import read_hfe_json

from cx_plot_data import plot_data
from dx_interpolate_data import interpolate_data, plot_interpolated
from ex_determine_bbox import determine_bbox
from fx_determine_dates import determine_dates
from gx_identify_precipitation_event import identify_precipitation_event





# --------------------
# SPECIFICS OF WHAT TO ANALYSE
# --------------------
ifeatures    = [     873, 1092, 1098, 1139, 1140, 1174, 1200, 1211, 1247, 1248,
                    1249, 1270, 1277, 1278, 1293, 1298, 1434, 1446, 1457, 1465,
                    1512, 1558, 1564, 1646, 1647, 1648, 1651, 1658, 1759, 1794  ] # starting after 2018-01-01 (using Geomet data)
#ifeatures    = [    1794  ]
#ifeatures    = [    1140  ]
#ifeatures    = [     875  ]   # needs addition of time step one before
#ifeatures    = [    1200  ]   # has missing time steps (2019-07-03T06 and 2019-07-06T12) but not in highlighted domain
ifeatures    = [    1277  ]   # has missing time steps (2019-07-03T06 and 2019-07-06T12) but     in highlighted domain
#ifeatures    = [    1512  ]   # has missing time steps (2019-07-03T06 and 2019-07-06T12) but all at the end

# ifeatures    = [     794, 883, 932, 936, 963, 965, 1008, 1009, 1036, 1043, 1054,
#                     1057, 1058, 1059, 1157, 1303, 1362, 1464, 1478, 1534, 1686,
#                     1737, 1756 ]   # 2017 features

#ifeatures    = range(1854)
#ifeatures    = range(3)

bbox_buffer  = 0.5
dates_buffer = [5.0,5.0]

# --------------------
# Load HFE database
# --------------------
print("\n\nReading HFE database")

# ignore warnings produced
warnings.filterwarnings("ignore", category=UserWarning, message='read_hfe_json')

filename        = '../data/hfe/historical_flood.json'
filtering       = True
polygon         = None
return_filtered = False
silent          = True

data_hfe = read_hfe_json(filename=filename,filtering=filtering,polygon=polygon,return_filtered=return_filtered,silent=silent)
nfeatures = len(data_hfe['data']['features'])
print("   Number of flood occurrences/events found: {}".format(nfeatures))

all_files_nc = {}
for iifeature,ifeature in enumerate(ifeatures):

    print("\n\n")
    print("--------------------------------")
    print("Working on feature {} ({} of {})".format(ifeature,iifeature+1,len(ifeatures)))
    print("--------------------------------")
    print("\n\n")

    # --------------------
    # Pick one event/occurrence
    # --------------------
    if ifeature >= nfeatures or ifeature < 0:
        raise ValueError("Feature needs to be between {} and {}.".format(0,nfeatures))

    feature = data_hfe['data']['features'][ifeature]
    start_date = datetime.datetime(int(feature['properties']['start_date'][0:4]),int(feature['properties']['start_date'][5:7]),int(feature['properties']['start_date'][8:10]),0,0)
    if not(feature['properties']['end_date'] is None):
        end_date = datetime.datetime(int(feature['properties']['end_date'][0:4]),int(feature['properties']['end_date'][5:7]),int(feature['properties']['end_date'][8:10]),0,0)
    else:
        end_date = None

    if feature['geometry']['type'] == 'Point':
        nlocations = 1
    else:
        nlocations = len(feature['geometry']['coordinates'])

    print("   Type: {}\n   Number of locations: {}\n   OBJECTID: {}\n   Start: {}\n   End: {}\n   Cause: {}".format(
        feature['geometry']['type'],
        nlocations,
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
        # print(">>>> CASPAR: ",ifeature)
    else:
        product  = 'rdpa:10km:6f'
        # print(">>>> GEOMET: ",ifeature)

    silent       = True

    if feature['properties']['end_date'] is None:
        # if there is no end-date make end buffer a bit larger
        date = determine_dates(feature=feature,product=product,dates_buffer=list(np.array(dates_buffer)+[0.0,3.0]),silent=silent)
    else:
        date = determine_dates(feature=feature,product=product,dates_buffer=dates_buffer,silent=silent)
    print("   date : [ {}, {}, ..., {}, {} ] (in total {} time steps)".format(date[0], date[1],date[-2],date[-1],len(date)))


    # # --------------------------------------
    # # Collect files that will be needed from CaSPAr
    # # --------------------------------------
    # # which files do we need? not used but helpful to know which data to request/copy from remote
    # if product == 'RDRS_v2.1':

    #     # full list of dates and files
    #     tmp = [ [tt,datetime.datetime.strftime(tt-datetime.timedelta(hours=+12),'%Y%m%d12.nc')] for tt in np.array(date) ]

    #     # unique files
    #     ifiles = np.unique(np.array(tmp)[:,1])

    #     all_files_nc[ifeature] = ifiles

    # continue


    # --------------------
    # Request data
    # --------------------
    print("\n\nRequest data:")

    if (product == 'rdpa:10km:6f') or (product == 'rdpa:10km:24f'):

        # ignore warnings produced
        warnings.filterwarnings("ignore", category=UserWarning, message="request_geomet_grib2: File '.*' already exists. Will not be downloaded again.")

        filename  = '/tmp/tmp/analyse_occurrence_'+str(ifeature)+'/geomet'
        crs       = 'EPSG:4326'
        overwrite = False
        silent    = True
        file_geomet = request_geomet_grib2(product=product,date=date,bbox=bbox,crs=crs,filename=filename,overwrite=overwrite,silent=silent)

        nfiles = len(file_geomet)
        print("   Number of Geomet files downloaded: {}".format(nfiles))
        print("   Number of Geomet files missing:    {}".format(len(date) - nfiles))

    elif (product == 'RDRS_v2.1'):

        variable = 'RDRS_v2.1_A_PR0_SFC'
        file_caspar = request_caspar_nc(product=product,variable=variable,date=date,foldername='../data/caspar/rdrs_v2.1/',silent=False)

        #nfiles = len(file_caspar)
        nfiles = len(np.unique([item for sublist in [ list(file_caspar[ff]) for ff in file_caspar ] for item in sublist ]))
        print("   Number of CaSPAr files required: {}".format(nfiles))

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

        lintransform = {'a':1000.0,'b':0.0}
        silent       = True
        data = read_caspar_nc(variable=variable,filenames=file_caspar,bbox=bbox,lintransform=lintransform,silent=silent)

        ntime = np.shape(data['var'])[0]
        nlat = np.shape(data['lat'])[0]
        nlon = np.shape(data['lon'])[1]
        print("   Number of time steps read: {}".format(ntime))
        print("   Number of latitudes  read: {}".format(nlat))
        print("   Number of longitudes read: {}".format(nlon))

    else:

        raise ValueError("analyse_event: product not known")



    # --------------------
    # Interpolate data at all locations of event
    # --------------------
    print("\n\nInterpolate data:")

    if feature['geometry']['type'] == 'Point':
        locations = {'lon':np.array([feature['geometry']['coordinates'][0]]),'lat':np.array([feature['geometry']['coordinates'][1]])}
    else:
        locations = {'lon':np.array(feature['geometry']['coordinates'])[:,0],'lat':np.array(feature['geometry']['coordinates'])[:,1]}
    nlocations = len(locations['lon'])
    var       = data["var"]
    lat       = data["lat"]
    lon       = data["lon"]
    dates     = data["time"]
    interpolated_data = interpolate_data(var=var,lat=lat,lon=lon,locations=locations,bbox=bbox,post_process=True,silent=True)

    sum_prec = np.sum(interpolated_data['var'],axis=0)
    print("   Sum of precipitation [mm] at all {} locations over the time period evaluated: {}".format(nlocations,sum_prec))

    plot_interpolated(locations=locations,
                          dates=dates,
                          data=interpolated_data,
                          pngfile='/tmp/tmp/analyse_occurrence_'+str(ifeature)+'/interpolated_at_stations_occurrence_'+str(ifeature)+'_all-timesteps.png',
                          start_date=start_date,
                          end_date=end_date,
                          start_date_buffer=date[0],   # start date with buffer (no matter if avail or not)
                          end_date_buffer=date[-1],    # end   date with buffer (no matter if avail or not)
                          )


    # --------------------
    # Find cluster of large precipitation for event (within start to end incl. buffer)
    # --------------------
    print("\n\nIdentify large precipitation cluster:")

    highlight_dates_idx = identify_precipitation_event(feature=feature,product=product,dates=dates,data=interpolated_data,length_window_d=2,min_prec_window=3.0,min_prec=0.001,silent=True)

    sum_prec = [ np.sum(interpolated_data['var'][highlight_dates_idx[ilocation],ilocation]) for ilocation in range(nlocations) ]
    print("   Sum of precipitation [mm] at all {} locations over the time period identified: {}".format(nlocations,sum_prec))

    pngfile = '/tmp/tmp/analyse_occurrence_'+str(ifeature)+'/interpolated_at_stations_occurrence_'+str(ifeature)+'_identified-timesteps_'+product+'.png'
    file_interpolated = plot_interpolated(locations=locations,
                          dates=dates,
                          data=interpolated_data,
                          highlight_dates_idx=highlight_dates_idx,
                          pngfile=pngfile,
                          start_date=start_date,
                          end_date=end_date,
                          start_date_buffer=date[0],   # start date with buffer (no matter if avail or not)
                          end_date_buffer=date[-1],    # end   date with buffer (no matter if avail or not)
                          label="uuid = '{}'\nevent precip. considered = {:.2f} mm".format(feature['properties']['uuid'],sum_prec[0]),
                          )
    print("\n\nPlotted: \n  ",file_interpolated['png'])



    # # --------------------
    # # Plot data
    # # --------------------
    # print("\n\nPlot data:")

    # var          = data["var"]
    # lat          = data["lat"]
    # lon          = data["lon"]
    # date         = data["time"]
    # png          = True
    # gif          = True
    # legend       = True
    # cities       = True
    # basefilename = '/tmp/tmp/analyse_occurrence_'+str(ifeature)+'/plot_occurrence_'+str(ifeature)
    # overwrite    = False
    # silent       = True

    # plots_data = plot_data(var=var,lat=lat,lon=lon,date=date,
    #                            png=png,
    #                            gif=gif,
    #                            legend=legend,
    #                            cities=cities,
    #                            bbox=bbox,
    #                            basefilename=basefilename,
    #                            overwrite=overwrite,
    #                            silent=silent)

    # print("   Number of PNGs    plotted: {}".format(len(plots_data['png'])))
    # print("   Number of GIFs    plotted: {}".format(len(plots_data['gif'])))
    # print("   Number of Legends plotted: {}".format(len(plots_data['legend'])))



# all NC files from CaSPAr
#    total: 9507
files_needed_from_caspar = np.unique([item for sublist in [ list(all_files_nc[ff]) for ff in all_files_nc ] for item in sublist ])
nfiles_needed_from_caspar = len(files_needed_from_caspar)






# In total 30 out of 1854 features need Geomet

# Feature 875  starts 2018-08-08.    >>>> GEOMET:   873
# Feature 1096 starts 2018-07-24.    >>>> GEOMET:  1092
# Feature 1102 starts 2018-07-24.    >>>> GEOMET:  1098
# Feature 1143 starts 2019-10-01.    >>>> GEOMET:  1139
# Feature 1144 starts 2019-01-24.    >>>> GEOMET:  1140
# Feature 1178 starts 2018-07-24.    >>>> GEOMET:  1174
# Feature 1204 starts 2019-06-29.    >>>> GEOMET:  1200
# Feature 1215 starts 2018-07-27.    >>>> GEOMET:  1211
# Feature 1251 starts 2018-07-24.    >>>> GEOMET:  1247
# Feature 1252 starts 2018-07-24.    >>>> GEOMET:  1248
# Feature 1253 starts 2018-07-24.    >>>> GEOMET:  1249
# Feature 1274 starts 2018-08-04.    >>>> GEOMET:  1270
# Feature 1281 starts 2018-07-24.    >>>> GEOMET:  1277
# Feature 1282 starts 2018-07-24.    >>>> GEOMET:  1278
# Feature 1297 starts 2018-08-29.    >>>> GEOMET:  1293
# Feature 1302 starts 2018-08-29.    >>>> GEOMET:  1298
# Feature 1443 starts 2019-09-04.    >>>> GEOMET:  1434
# Feature 1455 starts 2019-09-04.    >>>> GEOMET:  1446
# Feature 1467 starts 2018-07-24.    >>>> GEOMET:  1457
# Feature 1476 starts 2018-06-30.    >>>> GEOMET:  1465
# Feature 1524 starts 2018-08-29.    >>>> GEOMET:  1512
# Feature 1573 starts 2018-07-27.    >>>> GEOMET:  1558
# Feature 1582 starts 2018-07-27.    >>>> GEOMET:  1564
# Feature 1667 starts 2019-08-17.    >>>> GEOMET:  1646
# Feature 1668 starts 2019-08-17.    >>>> GEOMET:  1647
# Feature 1669 starts 2018-08-29.    >>>> GEOMET:  1648
# Feature 1673 starts 2018-08-29.    >>>> GEOMET:  1651
# Feature 1681 starts 2018-08-29.    >>>> GEOMET:  1658
# Feature 1785 starts 2018-07-24.    >>>> GEOMET:  1759
# Feature 1822 starts 2018-08-29.    >>>> GEOMET:  1794


# In total 23 out of 1854 are 2017 occurrences (compare RDRS and GEOMET)

# >>>> 2017 feature:  794
# >>>> 2017 feature:  883
# >>>> 2017 feature:  932
# >>>> 2017 feature:  936
# >>>> 2017 feature:  963
# >>>> 2017 feature:  965
# >>>> 2017 feature:  1008
# >>>> 2017 feature:  1009
# >>>> 2017 feature:  1036
# >>>> 2017 feature:  1043
# >>>> 2017 feature:  1054
# >>>> 2017 feature:  1057
# >>>> 2017 feature:  1058
# >>>> 2017 feature:  1059
# >>>> 2017 feature:  1157
# >>>> 2017 feature:  1303
# >>>> 2017 feature:  1362
# >>>> 2017 feature:  1464
# >>>> 2017 feature:  1478
# >>>> 2017 feature:  1534
# >>>> 2017 feature:  1686
# >>>> 2017 feature:  1737
# >>>> 2017 feature:  1756