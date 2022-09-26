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




__all__ = ['analyse_event']

def analyse_event(ifeatures=None,tmpdir='/tmp/',bbox_buffer=0.5,dates_buffer=[5.0,5.0],silent=True):
    """
        Analyses multi-point feature (event) of HFE database.


        Definition
        ----------
        analyse_event(ifeatures=None,tmpdir='/tmp/',bbox_buffer=0.5,dates_buffer=[5.0,5.0])


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
                                        Default: [3.0,1.0]

        silent          Boolean         If set to True, nothing will be printed to terminal.
                                        Default: True


        Output          Format          Description
        -----           -----           -----------
        {json}          dict            Dictionary of files that are produced. Has keys 'png' and 'json'.


        Description
        -----------
        The function summarizes the workflow to analyse multi-point events (events) for the HFE database, i.e.,
        "historical_flood_event.json". The user needs to specify a feature index (or several). The other arguments are optional and
        set to default values.


        Restrictions
        ------------
        None.


        Examples
        --------

        >>> # Analyse events

        >>> files_event = analyse_event(ifeatures=[316],tmpdir='/tmp/',bbox_buffer=0.5,dates_buffer=[5.0,5.0],silent=True)
        >>> print("files_event['png'] = "+str(files_event['png'][0][0]))
        files_event['png'] = /tmp/analyse_event_316/event_dbae8959-f2e0-4bdc-a5d6-84007dec140c_rdpa-10km-6f_2020112900.png
        >>> print("files_event['png'] = "+str(files_event['png'][0][-1]))
        files_event['png'] = /tmp/analyse_event_316/event_dbae8959-f2e0-4bdc-a5d6-84007dec140c_rdpa-10km-6f_2020120718.png
        >>> print("files_event['legend'] = "+str(files_event['legend'][0][0]))
        files_event['legend'] = /tmp/analyse_event_316/event_dbae8959-f2e0-4bdc-a5d6-84007dec140c_rdpa-10km-6f_legend.png
        >>> print("files_event['gif'] = "+str(files_event['gif'][0][0]))
        files_event['gif'] = /tmp/analyse_event_316/event_dbae8959-f2e0-4bdc-a5d6-84007dec140c_rdpa-10km-6f.gif
        >>> print("files_event['png-ts'] = "+str(files_event['png-ts'][0][0]))
        files_event['png-ts'] = /tmp/analyse_event_316/event_dbae8959-f2e0-4bdc-a5d6-84007dec140c_rdpa-10km-6f.png



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
        raise ValueError("analyse_event: List of feature index(es), i.e., ifeatures, must be specified.")

    # initialize return
    result = {}
    result['png-ts'] = []
    result['png'] = []
    result['gif'] = []
    result['legend'] = []
    result['json'] = []

    # # all features
    # ifeatures = range(363)

    # # all GEOMET
    # ifeatures = [0, 1, 62, 168, 173, 178, 179, 192, 202, 210, 215, 240, 241, 245, 247, 277, 283, 294, 316, 339]

    # --------------------
    # Load HFE database
    # --------------------
    if not(silent): print("\n\nReading HFE database")

    # ignore warnings produced
    warnings.filterwarnings("ignore", category=UserWarning, message='read_hfe_json')

    filename        = str(Path(dir_path).parent)+'/data/hfe/historical_flood_event.json'
    filtering       = True
    polygon         = None
    return_filtered = False

    data_hfe = read_hfe_json(filename=filename,filtering=filtering,polygon=polygon,return_filtered=return_filtered,silent=True)
    nfeatures = len(data_hfe['data']['features'])
    if not(silent): print("   Number of flood occurrences/events found: {}".format(nfeatures))

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

        if not(feature['properties']['end_date'] is None):
            length_event = (end_date-start_date).days+(end_date-start_date).seconds/60./60./24.
            if (length_event > 90.):
                print("Event will NOT be analysed because it is TOO LONG:")
                print(">>> Length event {} (idx={}): {} [days]".format(
                    feature['properties']['event_id'],
                    ifeature,
                    length_event))
                return result

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
            raise ValueError("analyse_event: start date is '{}' which is before 1980-01-01 which is when precip data will become available".format(
                feature['properties']['start_date']))
        elif (start_date - datetime.datetime(2018,1,1)).days < 0:
            product  = 'RDRS_v2.1'
            # if not(silent): print(">>>> CASPAR: ",ifeature)
        else:
            product  = 'rdpa:10km:6f'
            # if not(silent): print(">>>> GEOMET: ",ifeature)

        if feature['properties']['end_date'] is None:
            # if there is no end-date make end buffer a bit larger
            date = determine_dates(feature=feature,product=product,dates_buffer=list(np.array(dates_buffer)+[0.0,3.0]),silent=True)
        else:
            date = determine_dates(feature=feature,product=product,dates_buffer=dates_buffer,silent=True)
        if not(silent): print("   date : [ {}, {}, ..., {}, {} ] (in total {} time steps)".format(date[0], date[1],date[-2],date[-1],len(date)))

        # --------------------
        # Request data
        # --------------------
        if not(silent): print("\n\nRequest data:")

        if (product == 'rdpa:10km:6f') or (product == 'rdpa:10km:24f'):

            # ignore warnings produced
            warnings.filterwarnings("ignore", category=UserWarning, message="request_geomet_grib2: File '.*' already exists. Will not be downloaded again.")

            filename  = str(Path(tmpdir+'/analyse_event_'+str(ifeature)+'/geomet'))
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

            raise ValueError("analyse_event: product not known")


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

            raise ValueError("analyse_event: product not known")



        # --------------------
        # Interpolate data at all locations of event
        # --------------------
        if not(silent): print("\n\nInterpolate data:")

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
        if not(silent): print("   Sum of precipitation [mm] at all {} locations over the time period evaluated: {}".format(nlocations,np.round(sum_prec,2)))

        # # plot interpolated data w/o highlighted time steps
        # pngfile = str(Path(tmpdir+'/analyse_event_'+str(ifeature)+'/event_'+str(ifeature)+'_'+product.replace(":","-")+'_no-highlight.png'))
        # pngfile = str(Path(tmpdir+'/analyse_event_'+str(ifeature)+'/event_'+feature['properties']['event_id']+'_'+product.replace(":","-")+'_no-highlight.png'))
        # file_interpolated = plot_interpolated(locations=locations,
        #                       dates=dates,
        #                       data=interpolated_data,
        #                       pngfile=pngfile,
        #                       start_date=start_date,
        #                       end_date=end_date,
        #                       start_date_buffer=date[0],   # start date with buffer (no matter if avail or not)
        #                       end_date_buffer=date[-1],    # end   date with buffer (no matter if avail or not)
        #                       )
        # if not(silent): print("\n\nPlotted: \n  ",file_interpolated['png'])
        # result['png'].append(file_interpolated['png'])


        # --------------------
        # Find cluster of large precipitation for event (within start to end incl. buffer)
        # --------------------
        if not(silent): print("\n\nIdentify large precipitation cluster:")

        highlight_dates_idx = identify_precipitation_event(feature=feature,product=product,dates=dates,data=interpolated_data,length_window_d=2,min_prec_window=3.0,min_prec=0.001,silent=True)

        sum_prec = [ np.sum(interpolated_data['var'][highlight_dates_idx[ilocation],ilocation]) for ilocation in range(nlocations) ]
        if not(silent): print("   Sum of precipitation [mm] at all {} locations over the time period identified: {}".format(nlocations,np.round(sum_prec,2)))

        # plot interpolated data w/ highlighted time steps
        pngfile = str(Path(tmpdir+'/analyse_event_'+str(ifeature)+'/event_'+str(ifeature)+'_'+product.replace(":","-")+'.png'))
        pngfile = str(Path(tmpdir+'/analyse_event_'+str(ifeature)+'/event_'+feature['properties']['event_id']+'_'+product.replace(":","-")+'.png'))
        file_interpolated = plot_interpolated(locations=locations,
                              dates=dates,
                              data=interpolated_data,
                              highlight_dates_idx=highlight_dates_idx,
                              pngfile=pngfile,
                              start_date=start_date,
                              end_date=end_date,
                              start_date_buffer=date[0],   # start date with buffer (no matter if avail or not)
                              end_date_buffer=date[-1],    # end   date with buffer (no matter if avail or not)
                              #label="event_id = '{}'\nevent precip. considered = {:.2f} mm".format(feature['properties']['event_id'],sum_prec[0]),
                              label="event_id = '{}'".format(feature['properties']['event_id']),
                              )
        if not(silent): print("\n\nPlotted (interpolated time series): \n  ",file_interpolated['png'])
        result['png-ts'].append(file_interpolated['png'])



        # --------------------
        # Plot data
        # --------------------

        # find earliest identified start-date and latest identified end-date
        min_tidx = np.min([ np.min(hh) for hh in highlight_dates_idx if len(hh) > 0 ])
        max_tidx = np.max([ np.max(hh) for hh in highlight_dates_idx if len(hh) > 0 ])

        var          = data["var"][min_tidx:max_tidx+1]
        lat          = data["lat"]
        lon          = data["lon"]
        dates_event  = data["time"][min_tidx:max_tidx+1]
        png          = True
        gif          = True
        legend       = True
        cities       = True
        basefilename = str(Path(tmpdir+'/analyse_event_'+str(ifeature)+'/event_'+str(ifeature)+'_'+product.replace(":","-")))
        basefilename = str(Path(tmpdir+'/analyse_event_'+str(ifeature)+'/event_'+feature['properties']['event_id']+'_'+product.replace(":","-")))
        overwrite    = False

        plots_data = plot_data(var=var,lat=lat,lon=lon,date=dates_event,
                                   png=png,
                                   gif=gif,
                                   legend=legend,
                                   cities=cities,
                                   bbox=bbox,
                                   basefilename=basefilename,
                                   overwrite=overwrite,
                                   silent=True)

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

        missing_dates_n = []
        available_dates_n = []
        missing_dates_perc = []
        available_dates_perc = []
        available_timesteps_precip_mm = []
        for ilocation in range(nlocations):

            # find if there are gaps in the highlighted interval
            if len(highlight_dates_idx[ilocation]) > 0:
                dates_all = determine_dates(
                    start_date=dates[np.sort(highlight_dates_idx[ilocation])][0],
                    end_date=dates[np.sort(highlight_dates_idx[ilocation])][-1],
                    product=product,
                    dates_buffer=[0.0,0.0],
                    silent=True)
                missing_dates_n.append( len(dates_all) - len(dates[highlight_dates_idx[ilocation]]) )
                available_dates_n.append( len(dates[highlight_dates_idx[ilocation]]) )
                missing_dates_perc.append( missing_dates_n[-1]*100./len(dates_all) )
                available_dates_perc.append( available_dates_n[-1]*100./len(dates_all) )
                available_timesteps_precip_mm.append( { str(dates[idx]): round(interpolated_data['var'][idx,ilocation],2) for idx in highlight_dates_idx[ilocation] } )
            else:
                missing_dates_n.append( 0 )
                available_dates_n.append( 0 )
                missing_dates_perc.append( 0.0 )
                available_dates_perc.append( 0.0 )
                available_timesteps_precip_mm.append( {} )

        results_dict = {}
        results_dict['accumulated_mm']                   = list(np.round(sum_prec,2))
        results_dict['start_date_w_buffer']              = str(date[0])
        results_dict['end_date_w_buffer']                = str(date[-1])
        results_dict['missing_timesteps_n']              = list(missing_dates_n)
        results_dict['missing_timesteps_%']              = list(np.round(missing_dates_perc,2))
        results_dict['available_timesteps_n']            = list(available_dates_n)
        results_dict['available_timesteps_%']            = list(np.round(available_dates_perc,2))
        results_dict['available_timesteps_precip_mm/dt'] = available_timesteps_precip_mm

        results_dict['data_source'] = product

        jsonfile = str(Path(tmpdir+'/analyse_event_'+str(ifeature)+'/event_'+str(ifeature)+'_'+product.replace(":","-")+'.json'))
        jsonfile = str(Path(tmpdir+'/analyse_event_'+str(ifeature)+'/event_'+feature['properties']['event_id']+'_'+product.replace(":","-")+'.json'))
        json_dump = json.dumps(results_dict)
        ff = open(jsonfile, "w")
        ff.write(json_dump)
        ff.close()

        if not(silent): print("Saved: \n  ",jsonfile)
        result['json'].append([jsonfile])

    return result






if __name__ == '__main__':
    # import doctest
    # doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)


    ifeatures     = ['316']
    tmpdir        = ['/tmp/']
    bbox_buffer   = ['0.5']
    dates_buffer  = ['5.0,5.0']
    silent        = False
    parser        = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                            description='''Analyse multi-point feature(s) (event) from HFE database.''')
    parser.add_argument('-i', '--ifeatures', action='store',
                        default=ifeatures, dest='ifeatures', metavar='ifeatures', nargs=1,
                        help='Index of feature(s) to analyse. Comma separted if several. Numbering starts with 0. (default: "873,1092")')
    parser.add_argument('-t', '--tmpdir', action='store',
                        default=tmpdir, dest='tmpdir', metavar='tmpdir', nargs=1,
                        help='Name of temporary directory to store outputs, i.e., PNG, GIF, and legend (default: /tmp/).')
    parser.add_argument('-b', '--bbox_buffer', action='store',
                        default=bbox_buffer, dest='bbox_buffer', metavar='bbox_buffer', nargs=1,
                        help='Buffer of bounding box around feature. Given in [degrees] (default: "0.5").')
    parser.add_argument('-d', '--dates_buffer', action='store',
                        default=dates_buffer, dest='dates_buffer', metavar='dates_buffer', nargs=1,
                        help='Buffer of time period around feature start and end date. Given in [days] (default: "5.0,5.0").')
    parser.add_argument('-s', '--silent', action='store_true', default=silent, dest="silent",
                        help="If set nothing will be printed to terminal. Default: false.")

    args          = parser.parse_args()
    ifeatures     = [ int(ii) for ii in args.ifeatures[0].split(',') ]
    tmpdir        = args.tmpdir[0]
    bbox_buffer   = float(args.bbox_buffer[0])
    dates_buffer  = [ float(ii) for ii in args.dates_buffer[0].split(',') ]
    silent        = args.silent

    del parser, args

    files_produced = analyse_event(ifeatures=ifeatures,tmpdir=tmpdir,bbox_buffer=bbox_buffer,dates_buffer=dates_buffer,silent=silent)

    print("\n\nAll files produced = ",files_produced)


    # for example, run for all Geomet features:
    # python analyse_event.py --ifeatures "0, 1, 62, 168, 173, 178, 179, 192, 202, 210, 215, 240, 241, 245, 247, 277, 283, 294, 316, 339" --bbox_buffer 0.5 --dates_buffer 5.0,5.0 --tmpdir "/project/6070465/julemai/nrcan-hfe/data/output/"
