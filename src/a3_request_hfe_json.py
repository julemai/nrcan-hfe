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

import warnings
from pathlib import Path

__all__ = ['request_hfe_json']

def request_hfe_json(filename=None, jsonfilebase=None,silent=True):
    """
        Request JSON dump of HFE database.

        This is a dummy function (at the moment).


        Definition
        ----------
        request_hfe_json(filename=None, jsonfilebase=None, silent=True)


        Input           Format          Description
        -----           -----           -----------
        filename        string          Name of HFE database
                                        (currently not used)

        jsonfilebase    string          Basename (no file extension) of the JSON files to
                                        create, e.g.,
                                           jsonfilebase = '../data/hfe/historical_flood'
                                        Filenames will be as follows:
                                           File 1: <jsonfilebase>.json
                                           File 2: <jsonfilebase>_event.json

        silent          Boolean        If set to True, nothing will be printed to terminal.
                                       Default: True


        Output          Format          Description
        -----           -----           -----------
        {json}          dict            Dictionary of files that are produced


        Description
        -----------
        The JSON version of the HFE database used here was provided by Philippe Aussant
        (philippe.aussant@NRCan-RNCan.gc.ca) and has been placed under:
        1. data/hfe/historical_flood.json         (points)
        2. data/hfe/historical_flood_event.json   (multipoints)

        To be consistent and assure the codes using this JSON file will work also with
        later versions, it is recommended to populate this function with the code Philippe
        used to obtain the JSON files. The function "read_hfe_json()" will check later
        that all required attributes are found in the JSON files.

        This function currently only checks that the two JSON files exist at the indicated
        location.


        Restrictions
        ------------
        Dummy function is currently just checking if files exist but should be the proper
        code of how they were produced.


        Examples
        --------

        >>> # Request data

        >>> jsonfilebase = '../data/hfe/historical_flood'
        >>> file_hfe = request_hfe_json(filename=None, jsonfilebase=jsonfilebase,silent=True)
        >>> print('file_hfe = '+str(file_hfe))
        file_hfe = {'json': ['../data/hfe/historical_flood.json', '../data/hfe/historical_flood.json']}


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

    if filename is None:
        warnings.warn("request_hfe_json: HFE data in JSON format are currently not produced. The function only checks if files exist at indicated location.")

    if jsonfilebase is None:
        raise ValueError('a3_request_hfe_json: Basename (w/o file extension) of JSON files to create needs to be specified.')


    # initialize return
    filenames = {'json':[]}

    # for now, just checking if files exist
    # --> Philippe, you might want to place your code right here
    filename_occurrence = Path(jsonfilebase+'.json')
    if filename_occurrence.exists():
        if not(silent): print('Checked file {}.'.format(filename_occurrence))
        filenames['json'].append(str(filename_occurrence))

    filename_event = Path(jsonfilebase+'.json')
    if filename_event.exists():
        if not(silent): print('Checked file {}.'.format(filename_event))
        filenames['json'].append(str(filename_event))


    return filenames



if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)


    # jsonfilebase = '../data/hfe/historical_flood'
    # file_hfe = request_hfe_json(filename=None, jsonfilebase=jsonfilebase, silent=True)
    # print('file_hfe = '+str(file_hfe))
