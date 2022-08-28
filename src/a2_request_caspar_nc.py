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

__all__ = ['request_caspar_nc']

def request_caspar_nc():
    """
        Request NetCDF from CaSPAr.

        This is a dummy function. Requests of data from CaSPAr have to be made manually.
        A description is provided here:
        https://github.com/julemai/nrcan-hfe/blob/main/data/caspar/README.md

        Or directly under CaSPAr:
        https://github.com/julemai/CaSPAr/wiki/How-to-get-started-and-download-your-first-data


        Definition
        ----------
        def request_caspar_nc()


        Input           Format          Description
        -----           -----           -----------
        none            none            none


        Output          Format          Description
        -----           -----           -----------
        none            none            none


        Description
        -----------
        A dummy function to provide some information on how to request data from CaSPAr.


        Restrictions
        ------------
        Requests need to be made manually.


        Examples
        --------

        Request data

        >>> file_caspar = request_caspar_nc()
        >>> print('file_caspar = '+str(file_caspar))
        file_caspar = None


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

    warnings.warn("CaSPAr data need to be requested manually. Find information under 'https://github.com/julemai/nrcan-hfe/blob/main/data/caspar/README.md'.")

    filenames = None

    return filenames



if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)


    # date = datetime.datetime(2022,8,24,12,0)
    # bbox = {"lat":{"min":45.0,"max":46.0},"lon":{"min":-74.0,"max":-73.0}}
    # filename = '/tmp/tmf/test'
    # file_geomet = request_geomet_grib2(product='rdpa:10km:6f',date=date,bbox=bbox,crs='EPSG:4326',filename=filename)
