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

from fx_determine_dates import determine_dates


__all__ = ['identify_precipitation_event']


def identify_precipitation_event(feature=None,product=None,dates=None,data=None,length_window_d=2,min_prec_window=3.0,min_prec=0.001,silent=True):
    """
        Identified precipitation event in a precipitation time series provided for a given feature.


        Definition
        ----------
        def identify_precipitation_event(feature=None,product=None,dates=None,data=None,length_window_d=2,min_prec_window=3.0,min_prec=0.001,silent=True)


        Input             Format           Description
        -----             -----            -----------

        feature           dict             Dictionary containing a feature that has properties of start and end date,e.g.,
                                                {    "type" : "Feature",
                                                     "geometry" : { ... },
                                                     "properties" : {
                                                        "start_date" : "2004-11-04",
				                                        "end_date" : "2004-11-09",
                                                        ...
                                                     }
                                                }
                                           In cases where "end_date" is None, it will be set to "start_date".
                                           The time of the start date is assumed to be 00:00:00 (UTC). The time of end
                                           date is assumed to be 23:59:59 (UTC).
                                           If "feature" is not provided, "start_date" and "end_date" needs to be specified.
                                           Default: None

        product           string           Specifying the product that will be used to request precipitation data for this
                                           event. This input is required since the different products have different time
                                           steps where they are available, e.g., "RDRS_v2.1" is available every hour while
                                           "rdpa:10km:24f" is only available 12:00 pm each day and "rdpa:10km:6f" is
                                           available 06:00 am/pm and 12:00 am/pm every day.
                                           Must be one of the following:
                                           [ "RDRS_v2.1", "rdpa:10km:24f", "rdpa:10km:6f" ]
                                           Default: None

        dates             list(datetime)   List of datetime objects specifying the time the variable ("data") is valid
                                           for.
                                           Default: None

        data              dict             Interpolated data as returned from "interpolate_data()" which is a dictionary
                                           containing keys "var", "lat", and "lon".
                                           Default: None

        length_window_d   float            Length of rolling window precipitation is analysed. Given in [days].
                                           Default: 2.0 [days]

        min_prec_window   float            Minimum precipitation accumulated over rolling window to be considered
                                           precipitation of event. Given in [mm] within length_window_d.
                                           Default: 3.0 [mm] (over 2.0 days)

        min_prec          float            Minimal amount of precipitation to be considered at one time step to be non-zero.
                                           These timesteps will be added at beginning and end until first value is smaller
                                           than this. Given in [mm/h].
                                           Default: 0.001 [mm/h] (at one time step)

        silent            Boolean          If set to True, nothing will be printed to terminal.
                                           Default: True

        Output            Format           Description
        -----             -----            -----------
        [idx1,            list             List of indexes (of "dates") that should be considered the major precipitation event.
        idx2, ..]                          If the list of indexes contains first (0) or last (-1) index of dates, one should
                                           consider to increase time window where data are provided to avoid start or end of
                                           precipitation event are not cut off.


        Description
        -----------
        It is assumed that a precipitation time series is provided ("data") at given dates ("dates") for a feature ("feature")
        as it is stored in the HFE database. The main event will then be identified starting at the start date given in the
        "feature". A rolling window of a certain length ("length_window_d") will be analysed. If the accumulated precipitation of the
        current window is large enough (larger than "min_prec_window"), all time steps of the window will be added and the window moved
        by one timestep (earlier). This is continued until no new time step is added prior to the start date (because the cumulative
        precipitation is too small, i.e., event start seems to have been found). Then only single time steps (prior to first current
        index) are checked and added until the first time step smaller than a threshold ("min_prec") is detected. This guarantees
        that the entire onset of the event is added. The same approach is then repeated starting the rolling window with the current
        last detected index (current end of event). After no more indexes are added by using the window, single time steps larger
        than a threshold ("min_prec") are added.


        Restrictions
        ------------
        No warning issue if provided time series (data) is too short and precipitation event might be cut off. Check if index 0 or
        n-1 of dates is contained and widen date range before running this function again.


        Examples
        --------

        >>> # --------------------------------------
        >>> # Read HFE database and select one feature
        >>> # --------------------------------------

        >>> from b3_read_hfe_json import read_hfe_json

        >>> filename        = '../data/hfe/historical_flood.json'
        >>> filtering       = True
        >>> polygon         = None
        >>> return_filtered = False
        >>> silent          = True

        >>> data_hfe = read_hfe_json(filename=filename,filtering=filtering,polygon=polygon,return_filtered=return_filtered,silent=silent)

        >>> ifeature = 1345    # uuid = e4440f9a-2c04-4735-840c-d2ba065d3e30    # was 1277 in "HFE_data_v1"
        >>> feature = data_hfe['data']['features'][ifeature]

        >>> # --------------------------------------
        >>> # Determine dates and bbox to request
        >>> # --------------------------------------

        >>> from ex_determine_bbox import determine_bbox
        >>> from fx_determine_dates import determine_dates

        >>> product = 'rdpa:10km:6f'
        >>> crs = 'EPSG:4326'
        >>> bbox_buffer  = 0.5
        >>> dates_buffer = [5.0,5.0]

        >>> bbox = determine_bbox(feature=feature,bbox_buffer=bbox_buffer,silent=True)
        >>> dates = determine_dates(feature=feature,product=product,dates_buffer=dates_buffer,silent=True)

        >>> # --------------------------------------
        >>> # Request and read data from Geomet
        >>> # --------------------------------------

        >>> from a1_request_geomet_grib2 import request_geomet_grib2
        >>> from b1_read_geomet_grib2 import read_geomet_grib2

        >>> # request data
        >>> filename = '/tmp/pytest_rdpa_10km_6f_feature_'+str(ifeature)
        >>> files_geomet = request_geomet_grib2(product=product,date=dates,bbox=bbox,crs=crs,filename=filename,silent=True)

        >>> # read data
        >>> data_geomet = read_geomet_grib2(filenames=files_geomet,silent=True)

        >>> # --------------------------------------
        >>> # Interpolate data from Geomet
        >>> # --------------------------------------

        >>> from dx_interpolate_data import interpolate_data, plot_interpolated

        >>> var       = data_geomet["var"]
        >>> lat       = data_geomet["lat"]
        >>> lon       = data_geomet["lon"]
        >>> dates     = data_geomet["time"]   # does not contain time steps with missing data anymore

        >>> locations = {'lon':np.array([feature['geometry']['coordinates'][0]]),'lat':np.array([feature['geometry']['coordinates'][1]])}
        >>> interpolate_geomet = interpolate_data(var=var,lat=lat,lon=lon,locations=locations,bbox=bbox,post_process=True,silent=True)

        >>> # --------------------------------------
        >>> # Identify event
        >>> # --------------------------------------
        >>> idx_event = identify_precipitation_event(feature=feature,product=product,dates=dates,data=interpolate_geomet,length_window_d=2,min_prec_window=3.0,min_prec=0.001,silent=True)
        >>> print("Indexes of dates identified to be considered for thie event: {} (len = {})".format(idx_event,len(idx_event[0])))
        Indexes of dates identified to be considered for thie event: [[15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28]] (len = 14)

        >>> ilocation = 0
        >>> nlocations = 1
        >>> sum_prec = [ np.sum(interpolate_geomet['var'][idx_event[ilocation],ilocation]) for ilocation in range(nlocations) ]
        >>> print("Sum of precipitation [mm] at all {} locations over the time period identified: {}".format(nlocations,sum_prec))
        Sum of precipitation [mm] at all 1 locations over the time period identified: [121.89938715300597]

        >>> # --------------------------------------
        >>> # Plot interpolated data with identified event highlighted
        >>> # --------------------------------------
        >>> pngfile = '/tmp/pytest_rdpa_10km_6f_interpolated_at_stations_occurrence_'+str(ifeature)+'_identified-timesteps_'+product+'.png'
        >>> label = 'no label'
        >>> file_interpolated = plot_interpolated(locations=locations,dates=dates,data=interpolate_geomet,highlight_dates_idx=idx_event,pngfile=pngfile,start_date_buffer=dates[0],end_date_buffer=dates[-1],label=label)
        >>> print("Plot created: {}".format(file_interpolated['png']))
        Plot created: ['/tmp/pytest_rdpa_10km_6f_interpolated_at_stations_occurrence_1345_identified-timesteps_rdpa:10km:6f.png']



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

    # checking inputs
    if (feature is None):
        raise ValueError("identify_precipitation_event: feature needs to be specified")
    if (product is None):
        raise ValueError("identify_precipitation_event: product needs to be specified")
    if (dates is None):
        raise ValueError("identify_precipitation_event: dates needs to be specified")
    if (data is None):
        raise ValueError("identify_precipitation_event: dataneeds to be specified")


    # initialize return
    result = []



    var       = data["var"]
    lat       = data["lat"]
    lon       = data["lon"]
    # dates     = data["time"]

    if feature['geometry']['type'] == 'Point':
        locations = {'lon':np.array([feature['geometry']['coordinates'][0]]),'lat':np.array([feature['geometry']['coordinates'][1]])}
    else:
        locations = {'lon':np.array(feature['geometry']['coordinates'])[:,0],'lat':np.array(feature['geometry']['coordinates'])[:,1]}
    nlocations = len(locations['lon'])


    idate_wbuffer = dates
    idate_wobuffer = determine_dates(feature=feature,product=product,dates_buffer=[0.0,0.0],silent=True)
    # idx of dates in date range w/o buffer
    idx_dates = [ list(idate_wbuffer).index(idate) for idate in idate_wobuffer if idate in idate_wbuffer ]

    highlight_dates_idx = []
    for ilocation in range(nlocations):

        idx_dates_loc = idx_dates

        ivar = var[:,ilocation]
        ilat = lat[ilocation]
        ilon = lon[ilocation]

        prec = np.sum(ivar[idx_dates_loc])  # initialize precipitation to be only values in exactly specified range
        idates = idate_wobuffer             # initialize dates         to be only values in exactly specified range

        length_window_d = 2       # in [days]
        min_prec_window = 3.0     # in [mm] within time window
        min_prec = 0.001            # in [mm/h] (minimal amount of precipitation to be condidered and added at beginning and end)

        dt = ((idate_wobuffer[1]-idate_wobuffer[0]).seconds)
        length_window_n = int(length_window_d*24*60*60 / dt)

        # (A1) check if EARLIER cluster of time steps lead to significant addition in precipitation amount
        idx_start = max(0,min(idx_dates_loc)-length_window_n)
        stop_search = False

        while (idx_start > 0) and not(stop_search):

            idx_search = list(range(idx_start,idx_start+length_window_n,1))
            add_prec = np.sum(ivar[idx_search])
            if (add_prec > min_prec_window):

                if not(silent): print("   --> A1: adding cluster before: ",add_prec,"(date = ",np.array(idate_wbuffer)[idx_search],")")

                # update variables (add all tested time steps)
                idx_dates_loc = list(np.unique( list(idx_search) + idx_dates_loc ))
                prec = np.sum(ivar[idx_dates_loc])
                idates = list(np.unique( list(np.array(idate_wbuffer)[idx_search]) + idates ))
                idx_start = max(0,idx_search[-1]-length_window_n)

            else:

                # no new prior event found --> stop search
                stop_search = True

        # (A2) add single EARLIER time steps in case they are non-zero as well
        idx_start = max(0,min(idx_dates_loc)-1)
        stop_search = False

        while (idx_start > 0) and not(stop_search):

            idx_search = idx_start
            if ivar[idx_search]/(dt/60/60.) > min_prec:

                add_prec = ivar[idx_search]

                if not(silent): print("   --> A2: adding single before: ",add_prec,"(date = ",idate_wbuffer[idx_search],")")

                # update variables (add single tested time steps)
                idx_dates_loc = list(np.unique( [ idx_search ] + idx_dates_loc ))
                prec = np.sum(ivar[idx_dates_loc])
                idates = list(np.unique( [np.array(idate_wbuffer)[idx_search]] + idates ))
                idx_start = max(0,idx_search-1)

            else:

                # no new non-negative value found --> stop search
                stop_search = True

        # (A3) add one single EARLIER time steps more such that step-wise plot looks ok
        if (idx_start > 0):

            add_prec = ivar[idx_search]

            if not(silent): print("   --> A3: adding single before more: ",add_prec,"(date = ",np.array(idate_wbuffer)[idx_search],")")

            # update variables (add single tested time steps)
            idx_dates_loc = list(np.unique( [ idx_search ] + idx_dates_loc ))
            prec = np.sum(ivar[idx_dates_loc])
            idates = list(np.unique( [np.array(idate_wbuffer)[idx_search]] + idates ))
            # idx_start = max(0,idx_search-1)



        # (B1) check if LATER cluster of time steps lead to significant addition in precipitation amount
        idx_end = min(len(idate_wbuffer),max(idx_dates_loc)+1+length_window_n)
        stop_search = False

        while (idx_end < len(idate_wbuffer)) and not(stop_search):

            idx_search = list(range(idx_end-length_window_n,idx_end,1))
            add_prec = np.sum(ivar[idx_search])
            if (add_prec > min_prec_window):

                if not(silent): print("   --> B1: adding cluster after: ",add_prec,"(date = ",idate_wbuffer[idx_search],")")

                # update variables (add all tested time steps)
                idx_dates_loc = list(np.unique( idx_dates_loc + list(idx_search) ))
                prec = np.sum(ivar[idx_dates_loc])
                idates = list(np.unique( idates + list(np.array(idate_wbuffer)[idx_search]) ))
                idx_end = min(len(idate_wbuffer),idx_search[0]+1+length_window_n)

            else:

                # no new prior event found --> stop search
                stop_search = True

        # (B2) add single LATER time steps in case they are non-zero as well
        idx_end = min(len(idate_wbuffer),max(idx_dates_loc)+1)
        stop_search = False

        while (idx_end < len(idate_wbuffer)) and not(stop_search):

            idx_search = idx_end
            if ivar[idx_search]/(dt/60/60.) > min_prec:

                add_prec = ivar[idx_search]

                if not(silent): print("   --> B2: adding single after: ",add_prec,"(date = ",idate_wbuffer[idx_search],")")

                 # update variables (add single tested time steps)
                idx_dates_loc = list(np.unique( idx_dates_loc + [ idx_search ] ))
                prec = np.sum(ivar[idx_dates_loc])
                idates = list(np.unique( idates + [np.array(idate_wbuffer)[idx_search]] ))
                idx_end = min(len(idate_wbuffer),idx_search+1)

            else:

                # no new non-negative value found --> stop search
                stop_search = True

        # (C) Cleanup; go once again through first and last time steps and remove all that might be almost zero
        #     This happens if A2 and B2 are actually not performed but
        while (len(idx_dates_loc) > 0) and (ivar[idx_dates_loc[0]]/(dt/60/60.) <= min_prec):
            remove_idx = idx_dates_loc.pop(0)
            if not(silent): print("   --> C1: remove_idx because close-to-zero = ",remove_idx)

        while (len(idx_dates_loc) > 0) and (ivar[idx_dates_loc[-1]]/(dt/60/60.) <= min_prec):
            remove_idx = idx_dates_loc.pop(-1)
            if not(silent): print("   --> C2: remove_idx because close-to-zero = ",remove_idx)


        highlight_dates_idx.append(idx_dates_loc)

    result = highlight_dates_idx

    # ----------------------------------------------------
    # Done.
    # ----------------------------------------------------
    return result



if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)


    # # --------------------------------------
    # # Read HFE database and select one feature
    # # --------------------------------------

    # from b3_read_hfe_json import read_hfe_json

    # filename        = '../data/hfe/historical_flood.json'
    # filtering       = True
    # polygon         = None
    # return_filtered = False
    # silent          = True

    # data_hfe = read_hfe_json(filename=filename,filtering=filtering,polygon=polygon,return_filtered=return_filtered,silent=silent)

    # ifeature = 1277
    # feature = data_hfe['data']['features'][ifeature]

    # # --------------------------------------
    # # Determine dates and bbox to request
    # # --------------------------------------

    # from ex_determine_bbox import determine_bbox
    # from fx_determine_dates import determine_dates

    # product = 'rdpa:10km:6f'
    # crs = 'EPSG:4326'
    # bbox_buffer  = 0.5
    # dates_buffer = [5.0,5.0]

    # bbox = determine_bbox(feature=feature,bbox_buffer=bbox_buffer,silent=True)
    # dates = determine_dates(feature=feature,product=product,dates_buffer=dates_buffer,silent=True)

    # # --------------------------------------
    # # Request and read data from Geomet
    # # --------------------------------------

    # from a1_request_geomet_grib2 import request_geomet_grib2
    # from b1_read_geomet_grib2 import read_geomet_grib2

    # # request data
    # filename = '/tmp/pytest_rdpa_10km_6f_feature_'+str(ifeature)
    # files_geomet = request_geomet_grib2(product=product,date=dates,bbox=bbox,crs=crs,filename=filename,silent=True)

    # # read data
    # data_geomet = read_geomet_grib2(filenames=files_geomet,silent=True)

    # # --------------------------------------
    # # Interpolate data from Geomet
    # # --------------------------------------

    # from dx_interpolate_data import interpolate_data, plot_interpolated

    # var       = data_geomet["var"]
    # lat       = data_geomet["lat"]
    # lon       = data_geomet["lon"]
    # dates     = data_geomet["time"]   # does not contain time steps with missing data anymore

    # locations = {'lon':np.array([feature['geometry']['coordinates'][0]]),'lat':np.array([feature['geometry']['coordinates'][1]])}
    # interpolate_geomet = interpolate_data(var=var,lat=lat,lon=lon,locations=locations,bbox=bbox,post_process=True,silent=True)

    # # --------------------------------------
    # # Identify event
    # # --------------------------------------
    # idx_event = identify_precipitation_event(feature=feature,product=product,dates=dates,data=interpolate_geomet,length_window_d=2,min_prec_window=3.0,min_prec=0.001,silent=True)
    # print("Indexes of dates identified to be considered for thie event: {} (len = {})".format(idx_event,len(idx_event[0])))

    # ilocation = 0
    # nlocations = 1
    # sum_prec = [ np.sum(interpolate_geomet['var'][idx_event[ilocation],ilocation]) for ilocation in range(nlocations) ]
    # print("Sum of precipitation [mm] at all {} locations over the time period identified: {}".format(nlocations,sum_prec))

    # # --------------------------------------
    # # Plot interpolated data with identified event highlighted
    # # --------------------------------------
    # pngfile = '/tmp/pytest_rdpa_10km_6f_interpolated_at_stations_occurrence_'+str(ifeature)+'_identified-timesteps_'+product+'.png'
    # file_interpolated = plot_interpolated(locations=locations,
    #                       dates=dates,
    #                       data=interpolate_geomet,
    #                       highlight_dates_idx=idx_event,
    #                       pngfile=pngfile,
    #                       start_date_buffer=dates[0],   # start date with buffer (no matter if avail or not)
    #                       end_date_buffer=dates[-1],    # end   date with buffer (no matter if avail or not)
    #                       label="uuid = '{}'\nevent precip. considered = {:.2f} mm".format(feature['properties']['uuid'],sum_prec[0]),
    #                       )
    # print("Plot created: {}".format(file_interpolated['png']))
