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



import numpy as np
import datetime as datetime
import warnings
import json as json
import argparse
from pathlib import Path

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




__all__ = ['analyse_occurrence']

def analyse_occurrence(ifeatures=None,tmpdir='/tmp/',bbox_buffer=0.5,dates_buffer=[5.0,0.0],silent=True):
    """
        Analyses single-point feature (occurrence) of HFE database.


        Definition
        ----------
        analyse_occurrence(ifeatures=None,tmpdir='/tmp/',bbox_buffer=0.5,dates_buffer=[5.0,0.0])


        Input           Format          Description
        -----           -----           -----------
        ifeatures       list(int)       List of indexes of feature(s) to analyse.
                                        Numbering starts with 0.
                                        Default: None

        tmpdir          string          Name of folder where outputs will be stored.
                                        Default: /tmp/

        bbox_buffer     float           Buffer around points specified (in degree). If set to zero, points will be on
                                        edge/vertex of bounding box. If set to positive value, points will be truly
                                        located within the bounding box with at least "bbox_buffer" distance to an
                                        edge/vertex. Negative values for bbox_buffer are not allowed.
                                        Values are in [degree].
                                        Default: 0.5

        dates_buffer   [float,float]    Buffer around period [start_date,end_date] specified (in days). The earliest timesteps
                                        that will be returned will be "start_date-dates_buffer[0]" while the latest
                                        timestep returned will be "end_date+dates_buffer[1]".
                                        Values are in [days].
                                        Negative values of dates_buffer are not allowed.
                                        Default: [5.0,0.0]

        silent          Boolean         If set to True, nothing will be printed to terminal.
                                        Default: True


        Output          Format          Description
        -----           -----           -----------
        {json}          dict            Dictionary of files that are produced. Has keys 'png' and 'json'.


        Description
        -----------
        The function summarizes the workflow to analyse single-point events (occurrences) for the HFE database, i.e.,
        "historical_flood.json". The user needs to specify a feature index (or several). The other arguments are optional and
        set to default values.


        Restrictions
        ------------
        None.


        Examples
        --------

        >>> # Analyse events

        >>> files_occur = analyse_occurrence(ifeatures=[873,1092],tmpdir='/tmp/',bbox_buffer=0.5,dates_buffer=[5.0,0.0],silent=True)
        >>> print("files_occur['png'] = "+str(files_occur['png'][0]))
        files_occur['png'] = ['/tmp/analyse_occurrence_873/occurrence_99843381-1cc2-483b-a5f7-baf4a1a0102c_rdpa-10km-6f.png']
        >>> print("files_occur['png'] = "+str(files_occur['png'][1]))
        files_occur['png'] = ['/tmp/analyse_occurrence_1092/occurrence_73c4cd94-22c6-4859-bb74-2221e0c32ca1_rdpa-10km-6f.png']
        >>> print("files_occur['json'] = "+str(files_occur['json'][0]))
        files_occur['json'] = ['/tmp/analyse_occurrence_873/occurrence_99843381-1cc2-483b-a5f7-baf4a1a0102c_rdpa-10km-6f.json']
        >>> print("files_occur['json'] = "+str(files_occur['json'][1]))
        files_occur['json'] = ['/tmp/analyse_occurrence_1092/occurrence_73c4cd94-22c6-4859-bb74-2221e0c32ca1_rdpa-10km-6f.json']



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
        Written,  Juliane Mai, September 2022
    """

    if ifeatures is None:
        raise ValueError("analyse_occurrence: List of feature index(es), i.e., ifeatures, must be specified.")

    # initialize return
    result = {}
    result['png-ts'] = []
    result['png'] = []
    result['legend'] = []
    result['json'] = []
    result['gif'] = []

    # # --------------------
    # # SPECIFICS OF WHAT TO ANALYSE
    # # --------------------
    # ifeatures    = [     873, 1092, 1098, 1139, 1140, 1174, 1200, 1211, 1247, 1248,
    #                     1249, 1270, 1277, 1278, 1293, 1298, 1434, 1446, 1457, 1465,
    #                     1512, 1558, 1564, 1646, 1647, 1648, 1651, 1658, 1759, 1794  ] # starting after 2018-01-01 (using Geomet data)
    # #ifeatures    = [    1794  ]
    # #ifeatures    = [    1140  ]
    # #ifeatures    = [     875  ]   # needs addition of time step one before
    # #ifeatures    = [    1200  ]   # has missing time steps (2019-07-03T06 and 2019-07-06T12) but not in highlighted domain
    # #ifeatures    = [    1277  ]   # has missing time steps (2019-07-03T06 and 2019-07-06T12) but     in highlighted domain
    # #ifeatures    = [    1512  ]   # has missing time steps (2019-07-03T06 and 2019-07-06T12) but all at the end
    # #ifeatures    = [    1759  ]   # has missing time steps in highlighted period; has start and end date set in HFE

    # # ifeatures    = [     794, 883, 932, 936, 963, 965, 1008, 1009, 1036, 1043, 1054,
    # #                     1057, 1058, 1059, 1157, 1303, 1362, 1464, 1478, 1534, 1686,
    # #                     1737, 1756 ]   # 2017 features

    # #ifeatures    = range(1854)
    # ifeatures    = range(3)

    # bbox_buffer  = 0.5
    # dates_buffer = [5.0,0.0]

    # --------------------
    # Load HFE database (occurrences)
    # --------------------
    if not(silent): print("\n\nReading HFE database")

    # ignore warnings produced
    warnings.filterwarnings("ignore", category=UserWarning, message='read_hfe_json')

    filename        = str(Path(dir_path).parent)+'/data/hfe/historical_flood.json'
    filtering       = True
    polygon         = None
    return_filtered = False

    data_hfe = read_hfe_json(filename=filename,filtering=filtering,polygon=polygon,return_filtered=return_filtered,silent=True)
    nfeatures = len(data_hfe['data']['features'])
    if not(silent): print("   Number of flood occurrences/events found: {}".format(nfeatures))

    # --------------------
    # Process each occurrence
    # --------------------

    all_files_nc = {}
    for iifeature,ifeature in enumerate(ifeatures):

        if not(silent): print("\n\n")
        if not(silent): print("--------------------------------")
        if not(silent): print("Working on feature {} ({} of {})".format(ifeature,iifeature+1,len(ifeatures)))
        if not(silent): print("--------------------------------")
        if not(silent): print("\n\n")

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

        if not(silent):
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
        if not(silent): print("\n\nDetermine parameters for data requests:")

        bbox = determine_bbox(feature=feature,bbox_buffer=bbox_buffer,silent=True)
        if not(silent): print("   bbox : {}".format(bbox))


        # --------------------------------------
        # Determine time steps for rdpa:10km:6f
        # --------------------------------------
        if (start_date - datetime.datetime(1980,1,1)).days < 0:
            raise ValueError("analyse_occurrence: start date is '{}' which is before 1980-01-01 which is when precip data will become available".format(
                feature['properties']['start_date']))
        elif (start_date - datetime.datetime(2018,1,1)).days < 0:
            product  = 'RDRS_v2.1'
            # if not(silent): print(">>>> CASPAR: ",ifeature)
        else:
            product  = 'rdpa:10km:6f'
            # if not(silent): print(">>>> GEOMET: ",ifeature)

        if feature['properties']['end_date'] is None:
            # if there is no end-date make end buffer a bit larger
            date = determine_dates(feature=feature,product=product,dates_buffer=list(np.array(dates_buffer)+[0.0,8.0]),silent=True)
        else:
            date = determine_dates(feature=feature,product=product,dates_buffer=dates_buffer,silent=True)
        if not(silent): print("   date : [ {}, {}, ..., {}, {} ] (in total {} time steps)".format(date[0], date[1],date[-2],date[-1],len(date)))


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
        if not(silent): print("\n\nRequest data:")

        if (product == 'rdpa:10km:6f') or (product == 'rdpa:10km:24f'):

            # ignore warnings produced
            warnings.filterwarnings("ignore", category=UserWarning, message="request_geomet_grib2: File '.*' already exists. Will not be downloaded again.")

            filename  = str(Path(tmpdir+'/analyse_occurrence_'+str(ifeature)+'/geomet'))
            crs       = 'EPSG:4326'
            overwrite = False
            file_geomet = request_geomet_grib2(product=product,date=date,bbox=bbox,crs=crs,filename=filename,overwrite=overwrite,silent=True)

            nfiles = len(file_geomet)
            missing_dates_perc   = (len(date) - nfiles) *100. / len(date)
            missing_dates_n      = (len(date) - nfiles)

            if not(silent): print("   Number of Geomet files downloaded: {}".format(nfiles))
            if not(silent): print("   Number of Geomet files missing:    {} ({} %)".format(missing_dates_n,missing_dates_perc))

        elif (product == 'RDRS_v2.1'):

            variable = 'RDRS_v2.1_A_PR0_SFC'
            foldername = str(Path(dir_path).parent)+'/data/caspar/rdrs_v2.1/'
            file_caspar = request_caspar_nc(product=product,variable=variable,date=date,foldername=foldername,silent=False)

            #nfiles = len(file_caspar)
            nfiles = len(np.unique([item for sublist in [ list(file_caspar[ff]) for ff in file_caspar ] for item in sublist ]))
            missing_dates_perc = (len(date) - nfiles) *100. / len(date)
            missing_dates_n    = (len(date) - nfiles)
            if not(silent): print("   Number of CaSPAr files required: {}".format(nfiles))

        else:

            raise ValueError("analyse_occurrence: product not known")

        # --------------------
        # Read data
        # --------------------
        if not(silent): print("\n\nRead data:")

        if (product == 'rdpa:10km:6f'):

            lintransform = {'a':1.0,'b':0.0}
            data = read_geomet_grib2(filenames=file_geomet,lintransform=lintransform,silent=True)

            ntime = np.shape(data['var'])[0]
            nlat = np.shape(data['lat'])[0]
            nlon = np.shape(data['lon'])[1]
            if not(silent): print("   Number of time steps read: {}".format(ntime))
            if not(silent): print("   Number of latitudes  read: {}".format(nlat))
            if not(silent): print("   Number of longitudes read: {}".format(nlon))

        elif (product == 'RDRS_v2.1'):

            lintransform = {'a':1000.0,'b':0.0}
            data = read_caspar_nc(variable=variable,filenames=file_caspar,bbox=bbox,lintransform=lintransform,silent=True)

            ntime = np.shape(data['var'])[0]
            nlat = np.shape(data['lat'])[0]
            nlon = np.shape(data['lon'])[1]
            if not(silent): print("   Number of time steps read: {}".format(ntime))
            if not(silent): print("   Number of latitudes  read: {}".format(nlat))
            if not(silent): print("   Number of longitudes read: {}".format(nlon))

        else:

            raise ValueError("analyse_occurrence: product not known")



        # --------------------
        # Interpolate data at all locations of event
        # --------------------
        if not(silent): print("\n\nInterpolate data:")

        if feature['geometry']['type'] == 'Point':
            locations = {
                'lon':np.array([feature['geometry']['coordinates'][0]]),
                'lat':np.array([feature['geometry']['coordinates'][1]]),
                'name':np.array([feature['properties']['locality']])}
        else:
            locations = {'lon':np.array(feature['geometry']['coordinates'])[:,0],'lat':np.array(feature['geometry']['coordinates'])[:,1]}
        nlocations = len(locations['lon'])
        var       = data["var"]
        lat       = data["lat"]
        lon       = data["lon"]
        dates     = data["time"]
        interpolated_data = interpolate_data(var=var,lat=lat,lon=lon,locations=locations,bbox=bbox,post_process=True,silent=True)

        sum_prec = np.sum(interpolated_data['var'],axis=0)
        if not(silent): print("   Sum of precipitation [mm] at all {} locations over the time period evaluated: {}".format(nlocations,np.round(sum_prec,2)))

        # # plot interpolated data w/o highlighted time steps
        # pngfile = str(Path(tmpdir+'/analyse_occurrence_'+str(ifeature)+'/occurrence_'+str(ifeature)+'_'+product.replace(":","-")+'_no-highlight.png'))
        # pngfile = str(Path(tmpdir+'/analyse_occurrence_'+str(ifeature)+'/occurrence_'+feature['properties']['uuid']+'_'+product.replace(":","-")+'_no-highlight.png'))
        # file_interpolated = plot_interpolated(locations=locations,
        #                       dates=dates,
        #                       data=interpolated_data,
        #                       pngfile=pngfile,
        #                       start_date=start_date,
        #                       end_date=end_date,
        #                       start_date_buffer=date[0],   # start date with buffer (no matter if avail or not)
        #                       end_date_buffer=date[-1],    # end   date with buffer (no matter if avail or not)
        #                       )
        # print("\n\nPlotted: \n  ",file_interpolated['png'])
        # result['png-ts'].append(file_interpolated['png'])


        # --------------------
        # Find cluster of large precipitation for event (within start to end incl. buffer)
        # --------------------
        if not(silent): print("\n\nIdentify large precipitation cluster:")

        highlight_dates_idx = identify_precipitation_event(feature=feature,product=product,dates=dates,data=interpolated_data,length_window_d=2,min_prec_window=3.0,min_prec=0.001,silent=True)

        sum_prec = [ np.sum(interpolated_data['var'][highlight_dates_idx[ilocation],ilocation]) for ilocation in range(nlocations) ]
        if not(silent): print("   Sum of precipitation [mm] at all {} locations over the time period identified: {}".format(nlocations,np.round(sum_prec,2)))

        # plot interpolated data w/ highlighted time steps
        languages = ['en_CA','fr_CA']
        for ilanguage in languages:

            if (ilanguage.startswith('en_')):
                ilabel = "uuid = '{}'\nevent precip. considered = {:.2f} mm".format(feature['properties']['uuid'],sum_prec[0])
            elif (ilanguage.startswith('fr_')):
                ilabel = "uuid = '{}'\nprécip. de l’événement considérée = {:.2f} mm".format(feature['properties']['uuid'],sum_prec[0])
            else:
                raise ValueError("analyse_occurrence: language not known")


            pngfile = str(Path(tmpdir+'/analyse_occurrence_'+str(ifeature)+'/occurrence_'+ilanguage+'_'+str(ifeature)+'_'+product.replace(":","-")+'.png'))
            pngfile = str(Path(tmpdir+'/analyse_occurrence_'+str(ifeature)+'/occurrence_'+ilanguage+'_'+feature['properties']['uuid']+'_'+product.replace(":","-")+'.png'))
            file_interpolated = plot_interpolated(locations=locations,
                                  dates=dates,
                                  data=interpolated_data,
                                  highlight_dates_idx=highlight_dates_idx,
                                  pngfile=pngfile,
                                  start_date=start_date,
                                  end_date=end_date,
                                  start_date_buffer=date[0],   # start date with buffer (no matter if avail or not)
                                  end_date_buffer=date[-1],    # end   date with buffer (no matter if avail or not)
                                  label=ilabel,
                                  language=ilanguage,
                                  )
            if not(silent): print("\n\nPlotted: \n  ",file_interpolated['png'])
            result['png-ts'].append(file_interpolated['png'])

        # --------------------
        # Plot data (map with sum of precipitation for time steps identified)
        # --------------------
        languages = ['en_CA','fr_CA']
        for ilanguage in languages:

            # find earliest identified start-date and latest identified end-date
            if highlight_dates_idx != [[]]: # make sure event is found for at least one location
                min_tidx = np.min([ np.min(hh) for hh in highlight_dates_idx if len(hh) > 0 ])
                max_tidx = np.max([ np.max(hh) for hh in highlight_dates_idx if len(hh) > 0 ])

                var          = data["var"][min_tidx:max_tidx+1]
                lat          = data["lat"]
                lon          = data["lon"]
                dates_occur  = data["time"][min_tidx:max_tidx+1]  # only first timestep
                png          = False
                gif          = False
                legend       = True
                pngsum       = True  # PNG with sum of precip
                cities       = True
                basefilename = str(Path(tmpdir+'/analyse_occurrence_'+str(ifeature)+'/occurrence_'+ilanguage+'_'+str(ifeature)+'_'+product.replace(":","-")))
                basefilename = str(Path(tmpdir+'/analyse_occurrence_'+str(ifeature)+'/occurrence_'+ilanguage+'_'+feature['properties']['uuid']+'_'+product.replace(":","-")))
                overwrite    = False

                plots_data = plot_data(var=var,lat=lat,lon=lon,date=dates_occur,
                                           png=png,
                                           pngsum=pngsum,
                                           gif=gif,
                                           legend=legend,
                                           cities=cities,
                                           bbox=bbox,
                                           basefilename=basefilename,
                                           overwrite=overwrite,
                                           language=ilanguage,
                                           locations=locations,
                                           silent=True)
            else:
                plots_data = {'png':[],'gif':[],'legend':[]}

            if not(silent): print("\n\nPlotted (spatial data as PNG, GIF, and/or legend): \n  ")
            if not(silent):
                if len(plots_data['png']) > 0:
                    print("   Number of PNGs    plotted: {:4d}, e.g., {}".format(len(plots_data['png']),plots_data['png'][0]))
                else:
                    print("   Number of PNGs    plotted: {:4d}".format(len(plots_data['png'])))
                if len(plots_data['gif']) > 0:
                    print("   Number of GIFs    plotted: {:4d}, e.g., {}".format(len(plots_data['gif']),plots_data['gif'][0]))
                else:
                    print("   Number of GIFs    plotted: {:4d}".format(len(plots_data['gif'])))
                if len(plots_data['legend']) > 0:
                    print("   Number of legends plotted: {:4d}, e.g., {}".format(len(plots_data['legend']),plots_data['legend'][0]))
                else:
                    print("   Number of legends plotted: {:4d}".format(len(plots_data['legend'])))

            result['png'].append(plots_data['png'])
            result['gif'].append(plots_data['gif'])
            result['legend'].append(plots_data['legend'])



        # --------------------
        # Save data that we will need to augment the HFE database
        # --------------------
        if not(silent): print("\n\nSave results as JSON:")

        ilocation = 0

        # find if there are gaps in the highlighted interval
        if len(highlight_dates_idx[ilocation]) > 0:
            dates_all = determine_dates(
                start_date=dates[np.sort(highlight_dates_idx[ilocation])][0],
                end_date=dates[np.sort(highlight_dates_idx[ilocation])][-1],
                product=product,
                dates_buffer=[0.0,0.0],
                silent=True)
            missing_dates_n = len(dates_all) - len(dates[highlight_dates_idx[ilocation]])
            available_dates_n = len(dates[highlight_dates_idx[ilocation]])
            missing_dates_perc = missing_dates_n*100./len(dates_all)
            available_dates_perc = available_dates_n*100./len(dates_all)
            available_timesteps_precip_mm = { str(dates[idx]):round(interpolated_data['var'][idx,ilocation],2) for idx in highlight_dates_idx[ilocation] }
        else:
            missing_dates_n = 0
            available_dates_n = 0
            missing_dates_perc = 0.0
            available_dates_perc = 0.0
            available_timesteps_precip_mm = {}

        results_dict = {}
        results_dict['accumulated_mm']                   = round(sum_prec[ilocation],2)
        results_dict['start_date_w_buffer']              = str(date[0])
        results_dict['end_date_w_buffer']                = str(date[-1])
        results_dict['missing_timesteps_n']              = missing_dates_n
        results_dict['missing_timesteps_%']              = round(missing_dates_perc,2)
        results_dict['available_timesteps_n']            = available_dates_n
        results_dict['available_timesteps_%']            = round(available_dates_perc,2)
        results_dict['available_timesteps_precip_mm/dt'] = available_timesteps_precip_mm

        results_dict['data_source'] = product

        jsonfile = str(Path(tmpdir+'/analyse_occurrence_'+str(ifeature)+'/occurrence_'+str(ifeature)+'_'+product.replace(":","-")+'.json'))
        jsonfile = str(Path(tmpdir+'/analyse_occurrence_'+str(ifeature)+'/occurrence_'+feature['properties']['uuid']+'_'+product.replace(":","-")+'.json'))
        json_dump = json.dumps(results_dict)
        ff = open(jsonfile, "w")
        ff.write(json_dump)
        ff.close()

        if not(silent): print("Saved: \n  ",jsonfile)
        result['json'].append([jsonfile])


        # # --------------------
        # # Plot data
        # # --------------------
        # if not(silent): print("\n\nPlot data:")

        # var          = data["var"]
        # lat          = data["lat"]
        # lon          = data["lon"]
        # date         = data["time"]
        # png          = True
        # gif          = True
        # legend       = True
        # cities       = True
        # basefilename = str(Path(tmpdir+'/analyse_occurrence_'+str(ifeature)+'/plot_occurrence_'+str(ifeature)))
        # overwrite    = False

        # plots_data = plot_data(var=var,lat=lat,lon=lon,date=date,
        #                            png=png,
        #                            gif=gif,
        #                            legend=legend,
        #                            cities=cities,
        #                            bbox=bbox,
        #                            basefilename=basefilename,
        #                            overwrite=overwrite,
        #                            silent=True)

        # if not(silent): print("   Number of PNGs    plotted: {}".format(len(plots_data['png'])))
        # if not(silent): print("   Number of GIFs    plotted: {}".format(len(plots_data['gif'])))
        # if not(silent): print("   Number of Legends plotted: {}".format(len(plots_data['legend'])))



    # # all NC files from CaSPAr
    # #    total: 9507
    # files_needed_from_caspar = np.unique([item for sublist in [ list(all_files_nc[ff]) for ff in all_files_nc ] for item in sublist ])
    # nfiles_needed_from_caspar = len(files_needed_from_caspar)

    return result




if __name__ == '__main__':
    # import doctest
    # doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)


    ifeatures     = ['873,1092']
    tmpdir        = ['/tmp/']
    bbox_buffer   = ['0.5']
    dates_buffer  = ['5.0,0.0']
    silent        = False
    parser        = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                            description='''Analyse single-point feature(s) (occurrence) from HFE database.''')
    parser.add_argument('-i', '--ifeatures', action='store',
                        default=ifeatures, dest='ifeatures', metavar='ifeatures', nargs=1,
                        help='Index of feature(s) to analyse. Comma separted if several. Numbering starts with 0. (default: "873,1092")')
    parser.add_argument('-t', '--tmpdir', action='store',
                        default=tmpdir, dest='tmpdir', metavar='tmpdir', nargs=1,
                        help='Name of temporary directory to store outputs, i.e., PNG and JSON file (default: /tmp/).')
    parser.add_argument('-b', '--bbox_buffer', action='store',
                        default=bbox_buffer, dest='bbox_buffer', metavar='bbox_buffer', nargs=1,
                        help='Buffer of bounding box around feature. Given in [degrees] (default: "0.5").')
    parser.add_argument('-d', '--dates_buffer', action='store',
                        default=dates_buffer, dest='dates_buffer', metavar='dates_buffer', nargs=1,
                        help='Buffer of time period around feature start and end date. Given in [days] (default: "5.0,0.0").')
    parser.add_argument('-s', '--silent', action='store_true', default=silent, dest="silent",
                        help="If set nothing will be printed to terminal. Default: false.")

    args          = parser.parse_args()
    ifeatures     = [ int(ii) for ii in args.ifeatures[0].split(',') ]
    tmpdir        = args.tmpdir[0]
    bbox_buffer   = float(args.bbox_buffer[0])
    dates_buffer  = [ float(ii) for ii in args.dates_buffer[0].split(',') ]
    silent        = args.silent

    del parser, args

    files_produced = analyse_occurrence(ifeatures=ifeatures,tmpdir=tmpdir,bbox_buffer=bbox_buffer,dates_buffer=dates_buffer,silent=silent)

    print("\n\nAll files produced = ",files_produced)


    # for example, run for all Geomet features:
    # python analyse_occurrence.py --ifeatures "873, 1092, 1098, 1139, 1140, 1174, 1200, 1211, 1247, 1248, 1249, 1270, 1277, 1278, 1293, 1298, 1434, 1446, 1457, 1465, 1512, 1558, 1564, 1646, 1647, 1648, 1651, 1658, 1759, 1794" --bbox_buffer 0.5 --dates_buffer 5.0,0.0 --tmpdir "/project/6070465/julemai/nrcan-hfe/data/output/"
