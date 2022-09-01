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


__all__ = ['determine_dates']


def determine_dates(start_date=None,end_date=None,feature=None,product=None,dates_buffer=[3.0,1.0],silent=True):
    """
        Derives a list of time steps during a specified time period of a flood occurrence.


        Definition
        ----------
        def determine_dates(lat=None,lon=None,bbox_buffer=0.5,silent=True)


        Input           Format         Description
        -----           -----          -----------
        start_date      datetime       Start date of a flood occurrence as datetime object. The time of end
                                       date is assumed to be 00:00:00 (UTC).
                                       If "start_date" is not specified, "feature" needs to be provided.
                                       Default: None

        end_date        datetime       End date of a flood occurrence as datetime object. The time of end
                                       date is assumed to be 23:59:59 (UTC).
                                       If "end_date" is not specified, it will be set to start_date.
                                       Default: None

        feature         dict           Dictionary containing a feature that has properties of start and end date,e.g.,
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

        product         string         Specifying the product that will be used to request precipitation data for this
                                       event. This input is required since the different products have different time
                                       steps where they are available, e.g., "RDRS_v2.1" is available every hour while
                                       "rdpa:10km:24f" is only available 12:00 pm each day and "rdpa:10km:6f" is
                                       available 06:00 am/pm and 12:00 am/pm every day.
                                       Must be one of the following:
                                       [ "RDRS_v2.1", "rdpa:10km:24f", "rdpa:10km:6f" ]
                                       Default: None

        dates_buffer   [float,float]   Buffer around period [start_date,end_date] specified (in days). The earliest timesteps
                                       that will be returned will be "start_date-dates_buffer[0]" while the latest
                                       timestep returned will be "end_date+dates_buffer[1]".
                                       Values are in [days].
                                       Negative values of dates_buffer are not allowed.
                                       Default: [3.0,1.0]

        silent          Boolean        If set to True, nothing will be printed to terminal.
                                       Default: True

        Output          Format         Description
        -----           -----          -----------
        [tstep1,        list           List of datetime objects specifying which timesteps should be requested for the
         tstep2, ... ]                 flood occurence starting at "start_date" and ending at "end_date".


        Description
        -----------
        Derives the timesteps that should be requested from the "product" to evaluate the flood that occurred between "start_date"
        and "end_date" (either specified as properties of "feature" or as separate arguments). If "end_date" is None, it will be set to
        "start_date". The time of any start date is assumed to be 00:00:00 (UTC) while the time of any end date is assumed to be
        23:59:59 (UTC).


        Restrictions
        ------------
        Only individual features can be handled; no list of features.


        Examples
        --------

        >>> period = [datetime.datetime(2000,5,28,3,0),datetime.datetime(2000,5,29,20,59)]
        >>> feature_w_end  = {"type" : "Feature", "geometry" : {"type" : "Point","coordinates" : [ -73.5, 45.25 ]}, "properties" : { "start_date" : "2000-05-28", "end_date" : "2000-05-29" }}
        >>> feature_wo_end = {"type" : "Feature", "geometry" : {"type" : "Point","coordinates" : [ -73.5, 45.25 ]}, "properties" : { "start_date" : "2000-05-28", "end_date" : None }}

        >>> # --------------------------------------
        >>> # Derive time steps for RDRS_v2.1 w/ start_date and end-date specified; w/o buffer
        >>> # --------------------------------------
        >>> tsteps = determine_dates(start_date=period[0],end_date=period[1],feature=None,product='RDRS_v2.1',dates_buffer=[0.0,0.0],silent=True)
        >>> print("tsteps[0] = ",tsteps[0:1])
        tsteps[0] = [datetime.datetime(2000, 5, 28, 3, 0)]
        >>> print("tsteps[10] = ",tsteps[10:11])
        tsteps[10] = [datetime.datetime(2000, 5, 28, 13, 0)]
        >>> print("tsteps[-1] = ",tsteps[-1:])
        tsteps[-1] = [datetime.datetime(2000, 5, 29, 20, 0)]
        >>> print("number of timesteps = ",len(tsteps))
        number of timesteps = 42

        >>> # --------------------------------------
        >>> # Derive time steps for RDRS_v2.1 w/ start_date and end-date specified; w/ buffer
        >>> # --------------------------------------
        >>> tsteps = determine_dates(start_date=period[0],end_date=period[1],feature=None,product='RDRS_v2.1',dates_buffer=[3.0,1.0],silent=True)
        >>> print("tsteps[0] = ",tsteps[0:1])
        tsteps[0] = [datetime.datetime(2000, 5, 25, 3, 0)]
        >>> print("tsteps[10] = ",tsteps[10:11])
        tsteps[10] = [datetime.datetime(2000, 5, 25, 13, 0)]
        >>> print("tsteps[-1] = ",tsteps[-1:])
        tsteps[-1] = [datetime.datetime(2000, 5, 30, 20, 0)]
        >>> print("number of timesteps = ",len(tsteps))
        number of timesteps = 138

        >>> # --------------------------------------
        >>> # Derive time steps for RDRS_v2.1 w/ start_date; w/ buffer
        >>> # --------------------------------------
        >>> tsteps = determine_dates(start_date=period[0],end_date=None,feature=None,product='RDRS_v2.1',dates_buffer=[3.0,1.0],silent=True)
        >>> print("tsteps[0] = ",tsteps[0:1])
        tsteps[0] = [datetime.datetime(2000, 5, 25, 3, 0)]
        >>> print("tsteps[10] = ",tsteps[10:11])
        tsteps[10] = [datetime.datetime(2000, 5, 25, 13, 0)]
        >>> print("tsteps[-1] = ",tsteps[-1:])
        tsteps[-1] = [datetime.datetime(2000, 5, 30, 3, 0)]
        >>> print("number of timesteps = ",len(tsteps))
        number of timesteps = 121

        >>> # --------------------------------------
        >>> # Derive time steps for RDRS_v2.1 w/ end-date specified in feature; w/ buffer
        >>> # --------------------------------------
        >>> tsteps = determine_dates(start_date=None,end_date=None,feature=feature_w_end,product='RDRS_v2.1',dates_buffer=[3.0,1.0],silent=True)
        >>> print("tsteps[0] = ",tsteps[0:1])
        tsteps[0] = [datetime.datetime(2000, 5, 25, 0, 0)]
        >>> print("tsteps[10] = ",tsteps[10:11])
        tsteps[10] = [datetime.datetime(2000, 5, 25, 10, 0)]
        >>> print("tsteps[-1] = ",tsteps[-1:])
        tsteps[-1] = [datetime.datetime(2000, 5, 30, 23, 0)]
        >>> print("number of timesteps = ",len(tsteps))
        number of timesteps = 144

        >>> # --------------------------------------
        >>> # Derive time steps for RDRS_v2.1 w/o end-date specified in feature; w/ buffer
        >>> # --------------------------------------
        >>> tsteps = determine_dates(start_date=None,end_date=None,feature=feature_wo_end,product='RDRS_v2.1',dates_buffer=[3.0,1.0],silent=True)
        >>> print("tsteps[0] = ",tsteps[0:1])
        tsteps[0] = [datetime.datetime(2000, 5, 25, 0, 0)]
        >>> print("tsteps[10] = ",tsteps[10:11])
        tsteps[10] = [datetime.datetime(2000, 5, 25, 10, 0)]
        >>> print("tsteps[-1] = ",tsteps[-1:])
        tsteps[-1] = [datetime.datetime(2000, 5, 30, 0, 0)]
        >>> print("number of timesteps = ",len(tsteps))
        number of timesteps = 121

        >>> # --------------------------------------
        >>> # Derive time steps for RDPA-24h w/ start_date and end-date specified; w/o buffer
        >>> # --------------------------------------
        >>> tsteps = determine_dates(start_date=period[0],end_date=period[1],feature=None,product='rdpa:10km:24f',dates_buffer=[0.0,0.0],silent=True)
        >>> print("tsteps = ",tsteps)
        tsteps =  [datetime.datetime(2000, 5, 28, 12, 0), datetime.datetime(2000, 5, 29, 12, 0)]

        >>> # --------------------------------------
        >>> # Derive time steps for RDPA-24h w/ start_date and end-date specified; w/ buffer
        >>> # --------------------------------------
        >>> tsteps = determine_dates(start_date=period[0],end_date=period[1],feature=None,product='rdpa:10km:24f',dates_buffer=[3.0,1.0],silent=True)
        >>> print("tsteps = ",tsteps)
        tsteps =  [datetime.datetime(2000, 5, 25, 12, 0), datetime.datetime(2000, 5, 26, 12, 0), datetime.datetime(2000, 5, 27, 12, 0), datetime.datetime(2000, 5, 28, 12, 0), datetime.datetime(2000, 5, 29, 12, 0), datetime.datetime(2000, 5, 30, 12, 0)]

        >>> # --------------------------------------
        >>> # Derive time steps for RDPA-6h w/ start_date and end-date specified; w/o buffer
        >>> # --------------------------------------
        >>> tsteps = determine_dates(start_date=period[0],end_date=period[1],feature=None,product='rdpa:10km:6f',dates_buffer=[0.0,0.0],silent=True)
        >>> print("tsteps = ",tsteps)
        tsteps = [datetime.datetime(2000, 5, 28, 6, 0), datetime.datetime(2000, 5, 28, 12, 0), datetime.datetime(2000, 5, 28, 18, 0), datetime.datetime(2000, 5, 29, 0, 0), datetime.datetime(2000, 5, 29, 6, 0), datetime.datetime(2000, 5, 29, 12, 0), datetime.datetime(2000, 5, 29, 18, 0)]

        >>> # --------------------------------------
        >>> # Derive time steps for RDPA-6h w/ start_date and end-date specified; w/ buffer
        >>> # --------------------------------------
        >>> tsteps = determine_dates(start_date=period[0],end_date=period[1],feature=None,product='rdpa:10km:6f',dates_buffer=[3.0,1.0],silent=True)
        >>> print("tsteps = ",tsteps)
        tsteps = [datetime.datetime(2000, 5, 25, 6, 0), datetime.datetime(2000, 5, 25, 12, 0), datetime.datetime(2000, 5, 25, 18, 0), datetime.datetime(2000, 5, 26, 0, 0), datetime.datetime(2000, 5, 26, 6, 0), datetime.datetime(2000, 5, 26, 12, 0), datetime.datetime(2000, 5, 26, 18, 0), datetime.datetime(2000, 5, 27, 0, 0), datetime.datetime(2000, 5, 27, 6, 0), datetime.datetime(2000, 5, 27, 12, 0), datetime.datetime(2000, 5, 27, 18, 0), datetime.datetime(2000, 5, 28, 0, 0), datetime.datetime(2000, 5, 28, 6, 0), datetime.datetime(2000, 5, 28, 12, 0), datetime.datetime(2000, 5, 28, 18, 0), datetime.datetime(2000, 5, 29, 0, 0), datetime.datetime(2000, 5, 29, 6, 0), datetime.datetime(2000, 5, 29, 12, 0), datetime.datetime(2000, 5, 29, 18, 0), datetime.datetime(2000, 5, 30, 0, 0), datetime.datetime(2000, 5, 30, 6, 0), datetime.datetime(2000, 5, 30, 12, 0), datetime.datetime(2000, 5, 30, 18, 0)]



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
        Written,  Juliane Mai, August 2022
    """

    # checking inputs
    dates_buffer = np.array(dates_buffer)

    if ((start_date is None) and (feature is None)):
        raise ValueError("determine_dates: either start_date or feature needs to be specified")
    if (not(start_date is None) and not(feature is None)):
        raise ValueError("determine_dates: either start_date or feature needs to be specified")
    if ((start_date is None) and not(end_date is None)):
        raise ValueError("determine_dates: if end_date is specified, start_date also needs to be specified")
    if (np.any(dates_buffer < 0.0)):
        raise ValueError("determine_dates: all values in dates_buffer needs to be non-negative")
    if (product is None):
        raise ValueError("determine_dates: product needs to be specified")
    if ( not(product in [ "RDRS_v2.1", "rdpa:10km:24f", "rdpa:10km:6f" ]) ):
        raise ValueError("determine_dates: product needs to be one of the following: [ 'RDRS_v2.1', 'rdpa:10km:24f', 'rdpa:10km:6f' ]")


    # initialize return
    result = []

    # set the start_date and end_date properly
    if not(start_date is None):
        if end_date is None:
            end_date = start_date + datetime.timedelta(days=1)

    if not(feature is None):
        start_date = feature['properties']['start_date']   # YYYY-MM-DD (format checked while "read_hfe_json()")
        start_date = datetime.datetime(int(start_date[0:4]),int(start_date[5:7]),int(start_date[8:10]),0,0)

        end_date = feature['properties']['end_date']   # YYYY-MM-DD (format checked while "read_hfe_json()") but might be "None"
        if end_date is None:
            end_date = start_date + datetime.timedelta(days=1)
        else:
            end_date = datetime.datetime(int(end_date[0:4]),int(end_date[5:7]),int(end_date[8:10]),23,59)

        if not(silent): print("start_date = ",start_date)
        if not(silent): print("end_date = ",end_date)

    # apply buffer
    start_date = start_date - datetime.timedelta(days=dates_buffer[0])
    end_date   =   end_date + datetime.timedelta(days=dates_buffer[1])

    # adjust start_date to earliest time step available in chosen product after current start_date
    if product == 'RDRS_v2.1':
        tstep = 1   # [hours]
        tini  = 0   # first issue per day [hours]
    elif product == 'rdpa:10km:24f':
        tstep = 24   # [hours]
        tini  = 12   # first issue per day [hours]
    elif product == 'rdpa:10km:6f':
        tstep = 6   # [hours]
        tini  = 0   # first issue per day [hours]
    else:
        raise ValueError("determine_dates: product needs to be one of the following: [ 'RDRS_v2.1', 'rdpa:10km:24f', 'rdpa:10km:6f' ]")

    if (start_date.hour+start_date.minute/60.)%tstep != 0: # not already at a valid time step
        start_date = start_date+datetime.timedelta(hours=1.0*tstep-(start_date.hour+start_date.minute/60.)%tstep)-datetime.timedelta(hours=tini)

    # first time step not to include anymore (easier for creating a list)
    end_date = end_date+datetime.timedelta(hours=1.0*tstep-(end_date.hour+end_date.minute/60.)%tstep)+datetime.timedelta(hours=tini)


    # timedelta only stores days, seconds and microseconds
    timediff_hour = (end_date-start_date).days*24. + (end_date-start_date).seconds/60./60. + (end_date-start_date).microseconds/1000000./60./60.
    tsteps = [ start_date + datetime.timedelta(hours=itime) for itime in np.arange(0.0,timediff_hour,tstep*1.0) ]

    if not(silent): result = [start_date,end_date]   # end_date is first time step NOT to include anymore
    result = tsteps

    # ----------------------------------------------------
    # Done.
    # ----------------------------------------------------
    return result



if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # period = [datetime.datetime(2000,5,28,3,0),datetime.datetime(2000,5,29,20,59)]
    # feature_w_end  = {"type" : "Feature", "geometry" : {"type" : "Point","coordinates" : [ -73.5, 45.25 ]}, "properties" : { "start_date" : "2000-05-28", "end_date" : "2000-05-29" }}
    # feature_wo_end = {"type" : "Feature", "geometry" : {"type" : "Point","coordinates" : [ -73.5, 45.25 ]}, "properties" : { "start_date" : "2000-05-28", "end_date" : None }}

    # # --------------------------------------
    # # Derive time steps for RDRS_v2.1 w/ start_date and end-date specified; w/o buffer
    # # --------------------------------------
    # tsteps = determine_dates(start_date=period[0],end_date=period[1],feature=None,product='RDRS_v2.1',dates_buffer=[0.0,0.0],silent=True)
    # print("tsteps = ",tsteps)
