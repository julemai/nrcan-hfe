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
sys.path.append(dir_path+'/lib')


import numpy as np
import json as json
import warnings
import datetime as datetime
from in_poly import in_poly     # in lib/



__all__ = ['read_hfe_json']

def read_hfe_json(filename=None,filtering=False,polygon=None,return_filtered=False,silent=True):
    """
        Read HFE database file in JSON format and checks that all relevant attributes are present.

        Definition
        ----------
        def read_hfe_json(filename=None,filtering=False,polygon=None,return_filtered=False,silent=True)


        Input            Format          Description
        -----            -----           -----------
        filename         string          Name JSON file to read.
                                         Default: None

        filtering        Boolean         If True, all features will be filtered such that:
                                         1. required keys are available
                                         2. start_date and end_date are YYYY-MM-DD
                                            (not just MM or YYYY, etc.)
                                         3. start_date is 1980-01-01 or later
                                         4. If end_date given it must be after start_date
                                         5. points are within "polygon" if not None
                                         Default: False

        polygon          list            List of lat/lon pairs to check whether points in JSON
                                         features are located within this polygon.
                                         Ignored if filtering=False.
                                         Default: None

        return_filtered  dict            Dictionary returning list of filtered features (uuid)
                                         and the reason why they were filtered.
                                         Ignored if filtering=False.
                                         Default: False

        silent           Boolean         If set to True, nothing will be printed to terminal.
                                         Default: True


        Output           Format          Description
        -----            -----           -----------
        {data,filtered}  dict            Dictionary containing content of JSON file (data) and a
                                         list of filtered features ("uuid" if occurence JSON file,
                                         "event_id" if event JSON file). The list will be empty if
                                         "return_filtered" is False. The data are the structure
                                         the input JSON looked like.


        Description
        -----------
        Read HFE database file in JSON format.

        Most importantly the content can be checked that:
        1. required keys are available
        2. start_date and end_date are YYYY-MM-DD
           (not just MM or YYYY, etc.)
        3. start_date is 1980-01-01 or later
        4. if end_date given it must be after start_date
        5. points are within "polygon" if not None

        These checks will make sure that any processing of data for those events and occurrences will
        run smoothly hereafter.


        Restrictions
        ------------
        List of required_keys should be updated whenever another key attribute is required in later processing.


        Examples
        --------

        >>> # Polygon that was used to request data from CaSPAr
        >>> polygon = [[-30,30], [-130,30], [-130,80], [-30,80]]

        >>> # ------------------------
        >>> # Read flood events (MULTIPOINTS)
        >>> # ------------------------

        >>> filename='../data/hfe/historical_flood_event.json'

        >>> # no filtering; just reading
        >>> data_hfe = read_hfe_json(filename=filename,filtering=False,polygon=None,silent=True)
        >>> print("Number of flood occurrences/events found: {}".format(len(data_hfe['data']['features'])))
        Number of flood occurrences/events found: 487

        >>> # filtering w/o checking if points are in polygon
        >>> data_hfe = read_hfe_json(filename=filename,filtering=True,polygon=None,silent=True)
        >>> print("Number of flood occurrences/events found: {}".format(len(data_hfe['data']['features'])))
        Number of flood occurrences/events found: 363

        >>> # filtering w/  checking if points are in polygon
        >>> data_hfe = read_hfe_json(filename=filename,filtering=True,polygon=polygon,silent=True)
        >>> print("Number of flood occurrences/events found: {}".format(len(data_hfe['data']['features'])))
        Number of flood occurrences/events found: 361

        >>> # ------------------------
        >>> # Read flood occurences (POINTS)
        >>> # ------------------------

        >>> filename='../data/hfe/historical_flood.json'

        >>> # no filtering; just reading
        >>> data_hfe = read_hfe_json(filename=filename,filtering=False,polygon=None,silent=True)
        >>> print("Number of flood occurrences/events found: {}".format(len(data_hfe['data']['features'])))
        Number of flood occurrences/events found: 1904

        >>> # filtering w/o checking if points are in polygon
        >>> data_hfe = read_hfe_json(filename=filename,filtering=True,polygon=None,silent=True)
        >>> print("Number of flood occurrences/events found: {}".format(len(data_hfe['data']['features'])))
        Number of flood occurrences/events found: 1854

        >>> # filtering w/  checking if points are in polygon
        >>> data_hfe = read_hfe_json(filename=filename,filtering=True,polygon=polygon,silent=True)
        >>> print("Number of flood occurrences/events found: {}".format(len(data_hfe['data']['features'])))
        Number of flood occurrences/events found: 1852


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
    if filename is None:
        raise ValueError("read_hfe_json: filename needs to be specified")

    if not(filtering):
        # overwrite some settings because filtering will not be applied anyway
        polygon = None
        return_filtered = False

    # initialize return
    result = {}
    result_filtered = {}

    # read data
    with open(filename) as json_file:
        data = json.load(json_file)

    # filter
    if filtering:

        required_keys = ['type', 'geometry', 'properties']
        required_keys_geom = ['type', 'coordinates']

        file_is_for = None
        if np.any([ 'uuid' in dd['properties'].keys() for dd in data['features'] ]):
            # at lest one feature has 'uuid' --> is file of occurrences
            file_is_for = 'occurrences'
            required_keys_prop = ['event_id', 'start_date', 'end_date', 'season', 'flood_cause','uuid', 'rainfall_mm']
        else:
            # no feature has 'uuid' --> is file with events
            file_is_for = 'events'
            required_keys_prop = ['event_id', 'start_date', 'end_date', 'season', 'flood_cause']


        # check 1. required keys are available
        filter_features = []
        for iifeature,ifeature in enumerate(data['features']):

            all_there = [ kk in ifeature.keys() for kk in required_keys ]
            if not(all_there):
                raise ValueError('read_hfe_json: File {} does not contain following keys: {}'.format(filename,required_keys))

            all_there = [ kk in ifeature['geometry'].keys() for kk in required_keys_geom ]
            if not(all_there):
                raise ValueError('read_hfe_json: File {} does not contain following keys for "geometry": {}'.format(filename,required_keys_geom))

            all_there = [ kk in ifeature['properties'].keys() for kk in required_keys_prop ]
            if not(all_there):
                filter_features.append(iifeature) # cant use uuid or event_id because this might be the one actually missing

        result_filtered['check_keys'] = filter_features

        # check 2. start_date and end_date are YYYY-MM-DD (not just MM or YYYY, etc.)
        filter_features = []
        for iifeature,ifeature in enumerate(data['features']):

            tofilter = False

            start_date = ifeature['properties']['start_date']
            if len(start_date.strip().split(" ")) > 1:
                tofilter = True
            if len(start_date.strip().split("-")) != 3:
                tofilter = True
            try:
                start_datetime = datetime.datetime(int(start_date[0:4]),int(start_date[5:7]),int(start_date[8:10]))
            except:
                tofilter = True

            if tofilter:

                if file_is_for == 'occurrences':
                    warnings.warn("read_hfe_json: start_date of feature[uuid={}] is '{}' which is not conform to expected YYYY-MM-DD.".format(ifeature['properties']['uuid'],start_date))
                    filter_features.append(ifeature['properties']['uuid'])
                elif file_is_for == 'events':
                    warnings.warn("read_hfe_json: start_date of feature[eventid={}] is '{}' which is not conform to expected YYYY-MM-DD.".format(ifeature['properties']['event_id'],start_date))
                    filter_features.append(ifeature['properties']['event_id'])
                else:
                    raise ValueError('read_hfe_json: What kind of file is this?')

                # do not even bother checking end_date
                # also makes sure "filter_feature" only appears once
                continue

            end_date   = ifeature['properties']['end_date']

            if (end_date != "") and not(end_date is None):

                if len(end_date.strip().split(" ")) > 1:
                    tofilter = True
                if len(end_date.strip().split("-")) != 3:
                    tofilter = True
                try:
                    end_datetime = datetime.datetime(int(end_date[0:4]),int(end_date[5:7]),int(end_date[8:10]))
                except:
                    tofilter = True

                if tofilter:
                    if file_is_for == 'occurrences':
                        warnings.warn("read_hfe_json: end_date of feature[uuid={}] is '{}' which is not conform to expected YYYY-MM-DD.".format(ifeature['properties']['uuid'],end_date))
                        filter_features.append(ifeature['properties']['uuid'])
                    elif file_is_for == 'events':
                        warnings.warn("read_hfe_json: end_date of feature[eventid={}] is '{}' which is not conform to expected YYYY-MM-DD.".format(ifeature['properties']['event_id'],end_date))
                        filter_features.append(ifeature['properties']['event_id'])
                    else:
                        raise ValueError('read_hfe_json: What kind of file is this?')

        result_filtered['check_date_format'] = filter_features

        # check 3. start_date is 1980-01-01 or later (because that's when RDRS-v2.1 becomes available; for events before no data could be retrieved anyway)
        filter_features = []
        for iifeature,ifeature in enumerate(data['features']):

            already_filtered = False
            if file_is_for == 'occurrences':
                if (ifeature['properties']['uuid'] in result_filtered['check_date_format']):
                    already_filtered = True
            elif file_is_for == 'events':
                if (ifeature['properties']['event_id'] in result_filtered['check_date_format']):
                    already_filtered = True
            else:
                raise ValueError('read_hfe_json: What kind of file is this?')

            if not(already_filtered):

                start_date = ifeature['properties']['start_date']
                start_datetime = datetime.datetime(int(start_date[0:4]),int(start_date[5:7]),int(start_date[8:10]))

                if (start_datetime - datetime.datetime(1980,1,1)).days < 0:
                    if file_is_for == 'occurrences':
                        warnings.warn("read_hfe_json: start_date of feature[uuid={}] is '{}' which is before 1980-01-01 which is when precip data will become available.".format(ifeature['properties']['uuid'],start_date))
                        filter_features.append(ifeature['properties']['uuid'])
                    elif file_is_for == 'events':
                        warnings.warn("read_hfe_json: start_date of feature[eventid={}] is '{}' which is before 1980-01-01 which is when precip data will become available.".format(ifeature['properties']['event_id'],start_date))
                        filter_features.append(ifeature['properties']['event_id'])
                    else:
                        raise ValueError('read_hfe_json: What kind of file is this?')

        result_filtered['check_too_early'] = filter_features

        # check 4. if end_date given, it must be after start_date
        filter_features = []
        for iifeature,ifeature in enumerate(data['features']):

            already_filtered = False
            if file_is_for == 'occurrences':
                if (ifeature['properties']['uuid'] in result_filtered['check_date_format']):
                    already_filtered = True
                if (ifeature['properties']['uuid'] in result_filtered['check_too_early']):
                    already_filtered = True
            elif file_is_for == 'events':
                if (ifeature['properties']['event_id'] in result_filtered['check_date_format']):
                    already_filtered = True
                if (ifeature['properties']['event_id'] in result_filtered['check_too_early']):
                    already_filtered = True
            else:
                raise ValueError('read_hfe_json: What kind of file is this?')

            if not(already_filtered):

                tofilter = False
                start_date = ifeature['properties']['start_date']
                end_date   = ifeature['properties']['end_date']
                if (end_date != "") and not(end_date is None):

                    start_datetime = datetime.datetime(int(start_date[0:4]),int(start_date[5:7]),int(start_date[8:10]))
                    end_datetime = datetime.datetime(int(end_date[0:4]),int(end_date[5:7]),int(end_date[8:10]))

                    if start_datetime > end_datetime:
                        tofilter = True

                if tofilter:
                    if file_is_for == 'occurrences':
                        warnings.warn("read_hfe_json: end_date of feature[uuid={}] is '{}' and this is before start_date '{}'.".format(ifeature['properties']['uuid'],end_date,start_date))
                        filter_features.append(ifeature['properties']['uuid'])
                    elif file_is_for == 'events':
                        warnings.warn("read_hfe_json: end_date of feature[eventid={}] is '{}' and this is before start_date '{}'.".format(ifeature['properties']['event_id'],end_date,start_date))
                        filter_features.append(ifeature['properties']['event_id'])
                    else:
                        raise ValueError('read_hfe_json: What kind of file is this?')

        result_filtered['check_enddate_too_early'] = filter_features

        # check 5. points are within "polygon" if not None
        if not(polygon is None):

            lon_poly = np.array(polygon)[:,0]
            lat_poly = np.array(polygon)[:,1]

            filter_features = []
            for iifeature,ifeature in enumerate(data['features']):

                tofilter = False
                coords = ifeature['geometry']['coordinates']
                geomtype = ifeature['geometry']['type']
                if geomtype == "Point":
                    coords = [ coords ]
                for coord in coords:
                    # in_poly returns ::  1 if point inside polygon
                    #                    -1 if point outside polygon
                    #                     0 if point on vertex/edge of polygon
                    inpoly = in_poly( coord, lon_poly, lat_poly )
                    if inpoly < 1:
                        tofilter = True
                        break  # stop checking other ccordinates

                if tofilter:
                    if file_is_for == 'occurrences':
                        warnings.warn("read_hfe_json: coordinate of feature[uuid={}] is {} is outside of provided polygon.".format(ifeature['properties']['uuid'],coord))
                        filter_features.append(ifeature['properties']['uuid'])
                    elif file_is_for == 'events':
                        warnings.warn("read_hfe_json: coordinate(s) of feature[eventid={}] is {} which is outside of provided polygon.".format(ifeature['properties']['event_id'],coord))
                        filter_features.append(ifeature['properties']['event_id'])
                    else:
                        raise ValueError('read_hfe_json: What kind of file is this?')

            result_filtered['check_outside_polygon'] = filter_features
        else:
            result_filtered['check_outside_polygon'] = []

    if filtering:

        # get all features to filter in one list and make sure each only appears once
        all_filtered = []
        all_filtered += result_filtered['check_keys']
        all_filtered += result_filtered['check_date_format']
        all_filtered += result_filtered['check_too_early']
        all_filtered += result_filtered['check_enddate_too_early']
        all_filtered += result_filtered['check_outside_polygon']
        all_filtered = np.unique(all_filtered)

        # add everything from data read into tmp structure
        tmp = {}
        tmp['type'] = data['type']
        tmp['name'] = data['name']
        tmp['features'] = []
        for iifeature,ifeature in enumerate(data['features']):

            if file_is_for == 'occurrences':
                if not(ifeature['properties']['uuid'] in all_filtered):
                    tmp['features'].append( data['features'][iifeature] )
            elif file_is_for == 'events':
                if not(ifeature['properties']['event_id'] in all_filtered):
                    tmp['features'].append( data['features'][iifeature] )
            else:
                raise ValueError('read_hfe_json: What kind of file is this?')

        # overwrite read data
        data = tmp

    result['data'] = data
    result['filtered'] = result_filtered

    # return dictionary
    return result


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)



    # filename='../data/hfe/historical_flood.json'
    # filename='../data/hfe/historical_flood_event.json'
    # filtering=True
    # polygon=[[-61.43692000000001, 42.107902], [-69.679413, 43.733895], [-75.802917, 41.948767], [-81.803513, 40.336600000000004], [-85.34317, 40.905729], [-86.717148, 43.926583], [-90.00412, 46.520603], [-98.191681, 46.843286], [-109.592743, 46.635294], [-121.171646, 46.26819], [-127.60894799999998, 47.161241], [-133.259354, 50.607646], [-136.426163, 54.110943], [-135.518814, 56.883876], [-119.99999999999999, 59.99999999999999], [-110, 59.99999999999999], [-90, 52], [-70, 52], [-52.67807, 53.008173], [-50.820694, 50.59152], [-49.156265, 47.78086600000001], [-52.49267600000001, 44.657909], [-57.761306999999995, 45.09049], [-61.43692000000001, 42.107902]]
    # return_filtered=False
    # silent=True

    # data_hfe = read_hfe_json(filename=filename,filtering=filtering,polygon=polygon,silent=silent,return_filtered=True)

    # print("Number of flood occurrences/events found: {}".format(len(data_hfe['data']['features'])))
