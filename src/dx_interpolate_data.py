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
import scipy.interpolate


__all__ = ['interpolate_data']


def interpolate_data(var=None,lat=None,lon=None,locations=None,bbox=None,return_tmp=False,post_process=False,silent=True):
    """
        Plot data retrieved from either Geomet or CaSPAr (data will be all formatted the same by now).

        Definition
        ----------
        def interpolate_data(var=None,lat=None,lon=None,locations=None, bbox=None,silent=True)


        Input           Format         Description
        -----           -----          -----------
        var             2D array or    Variable to plot; needs to have same shape as "lat" and "lon" (2D)
                        3D array       but can be a list of those variables (3D) in case variable is available
                                       for several time steps. First dimension is then considered the time dimension.
                                       In this case "date" needs to be a list of dates (length of list equals
                                       size of first dimension equals number of days provided)
                                       Default: None

        lat             2D array       Latitudes of each grid cell. Shape needs to be the same as for "var".
                                       In case "var" is 3-dimensional, the shape of "lat" is the same as the
                                       second and third dimension of "var" since first dimension of "var" is time.
                                       Default: None

        lon             2D array       Longitudes of each grid cell. Shape needs to be the same as for "var".
                                       In case "var" is 3-dimensional, the shape of "lat" is the same as the
                                       second and third dimension of "var" since first dimension of "var" is time.
                                       Default: None

        locations       dict           Dictionary providing attributes "lat" as list of latitudes and "lon" as list
                                       of "longitudes" for locations where raster shall be interpolated. Lists of
                                       "lat" and "lon" are assumed to be of the same length. If only one location
                                       is requested, instead of lists a scalar can be provided

        bbox            dict           Dictionary specifying the bounding box, i.e.,
                                         {"lat":{"min":min_lat,"max":max_lat},
                                          "lon":{"min":min_lon,"max":max_lon}}
                                       If provided, the location will be checked to fall into the bounding box.
                                       Error generated if not; even if data might be available in vicinity of the
                                       bounding box. In any case the code will check that data are available for
                                       locations. The bounding box check is only an additional safety check.
                                       Default: None

        return_tmp      Boolean        If set to True a "tmp" attribute will be returned containing a dictionary
                                       tmp = {"lat":tmp_lat, "lon":tmp_lon, "var":tmp_var} with information of the
                                       nine cells used for interpolation of each station. Hence the shape of "lat_tmp"
                                       is (nlocations,3,3), "lon_tmp" is (nlocations,3,3), "var_tmp" is
                                       (nlocations,ntime,3,3).
                                       Default: False

        post_process    Boolean        If set to True, all interpolated values will be checked to be non-negative.
                                       If negative they will be set to zero; unless they are smaller than -0.1 which
                                       will raise an error. The error is raise because it is likely that the interpolation
                                       is not working properly when returning large negative values (for precipitation).
                                       Default: False

        silent          Boolean        If set to True, nothing will be printed to terminal.
                                       Default: True

        Output          Format         Description
        -----           -----          -----------
        {var,lat,lon}   dict           Dictionary {"var":var,"lat":lat,"lon":lon}
        or                             containing interpolated variable "var" at locations requested (2D in any case). For
        tmp}                           consistency, also "lat" and "lon" at requested points are returned but
                                       they are the same as the input "lat" and "lon" variables. "tmp" will only be returned
                                       if "return_tmp" is set to True.


        Description
        -----------
        Interpolates variable (var) in space at locations (locations) using bilinear interpolation. The base of this
        function is "scipy.interpolate.interp2d". All coordinates are transformed internally to meter (instead of
        degree north and east) to assure that the weights for the interpolation are based on a true geographical distance.


        Restrictions
        ------------
        None.


        Examples
        --------

        >>> # --------------------------------------
        >>> # Request and read data from Geomet
        >>> # --------------------------------------

        >>> from a1_request_geomet_grib2 import request_geomet_grib2
        >>> from b1_read_geomet_grib2 import read_geomet_grib2

        >>> product = 'rdpa:10km:6f'
        >>> crs = 'EPSG:4326'
        >>> bbox = {"lat":{"min":45.0,"max":46.0},"lon":{"min":-74.0,"max":-73.0}}

        >>> # four dates 6h apart from each other
        >>> dates = [ datetime.datetime(2022,8,24,12,0) + datetime.timedelta(hours=24*ii) for ii in range(2) ]

        >>> # request data
        >>> filename = '/tmp/pytest_rdpa_10km_6f'
        >>> files_geomet = request_geomet_grib2(product=product,date=dates,bbox=bbox,crs=crs,filename=filename,silent=True)

        >>> # read data
        >>> data_geomet = read_geomet_grib2(filenames=files_geomet,silent=True)

        >>> # --------------------------------------
        >>> # Interpolate data from Geomet
        >>> # --------------------------------------

        >>> # 4 locations ; w/ post-process (negative values set to zero)
        >>> locations = {"lat":[45.25,45.5,45.75,45.75], "lon":[-73.5,-73.5,-73.5,-73.75]}
        >>> interpolate_geomet = interpolate_data(var=data_geomet["var"],lat=data_geomet["lat"],lon=data_geomet["lon"],locations=locations,bbox=bbox,post_process=True,silent=True)
        >>> print("Interpolated data location #1: (", interpolate_geomet['lat'][0],",",interpolate_geomet['lon'][0],") --> time[0:2] = ",interpolate_geomet['var'][:,0])
        Interpolated data location #1: ( 45.25 , -73.5 ) --> time[0:2] =  [0.06768641 2.84033923]
        >>> print("Interpolated data location #2: (", interpolate_geomet['lat'][1],",",interpolate_geomet['lon'][1],") --> time[0:2] = ",interpolate_geomet['var'][:,1])
        Interpolated data location #2: ( 45.5 , -73.5 ) --> time[0:2] =  [0.03125    0.13948481]
        >>> print("Interpolated data location #3: (", interpolate_geomet['lat'][2],",",interpolate_geomet['lon'][2],") --> time[0:2] = ",interpolate_geomet['var'][:,2])
        Interpolated data location #3: ( 45.75 , -73.5 ) --> time[0:2] =  [0.00104299 0.        ]
        >>> print("Interpolated data location #4: (", interpolate_geomet['lat'][3],",",interpolate_geomet['lon'][3],") --> time[0:2] = ",interpolate_geomet['var'][:,3])
        Interpolated data location #4: ( 45.75 , -73.75 ) --> time[0:2] =  [0.00281244 0.        ]

        >>> # 1 location (as list) ; w/ post-process (negative values set to zero)
        >>> locations = {"lat":[45.25], "lon":[-73.5]}
        >>> interpolate_geomet = interpolate_data(var=data_geomet["var"],lat=data_geomet["lat"],lon=data_geomet["lon"],locations=locations,bbox=bbox,post_process=True,silent=True)
        >>> print("Interpolated data location #1: (", interpolate_geomet['lat'][0],",",interpolate_geomet['lon'][0],") --> time[0:2] = ",interpolate_geomet['var'][:,0])
        Interpolated data location #1: ( 45.25 , -73.5 ) --> time[0:2] =  [0.06768641 2.84033923]

        >>> # 1 location (as scalar) ; w/ post-process (negative values set to zero)
        >>> locations = {"lat":45.25, "lon":-73.5}
        >>> interpolate_geomet = interpolate_data(var=data_geomet["var"],lat=data_geomet["lat"],lon=data_geomet["lon"],locations=locations,bbox=bbox,post_process=True,silent=True)
        >>> print("Interpolated data location #1: (", interpolate_geomet['lat'][0],",",interpolate_geomet['lon'][0],") --> time[0:2] = ",interpolate_geomet['var'][:,0])
        Interpolated data location #1: ( 45.25 , -73.5 ) --> time[0:2] =  [0.06768641 2.84033923]

        >>> # --------------------------------------
        >>> # Request and read data from CaSPAr
        >>> # --------------------------------------

        >>> from a2_request_caspar_nc import request_caspar_nc
        >>> from b2_read_caspar_nc import read_caspar_nc

        >>> product='RDRS_v2.1'
        >>> variable='RDRS_v2.1_A_PR0_SFC'
        >>> dates=[ datetime.datetime(2018,8,8,15,0) + datetime.timedelta(hours=ii) for ii in range(10) ]
        >>> bbox={"lat":{"min":45.0,"max":46.0},"lon":{"min":-74.0,"max":-73.0}}
        >>> foldername='test-data/'
        >>> lintransform={'a':1000.0,'b':0.0}  # convert from m/h to mm/h
        >>> silent=False

        >>> # request data
        >>> files_caspar = request_caspar_nc(product=product,variable=variable,date=dates,foldername=foldername,silent=True)

        >>> # read data
        >>> data_caspar = read_caspar_nc( variable=variable, filenames=files_caspar, bbox=bbox, lintransform=lintransform, silent=True)

        >>> # --------------------------------------
        >>> # Interpolate data from CaSPAr
        >>> # --------------------------------------

        >>> # 4 locations; w/o post-process (could lead to negative values)
        >>> locations = {"lat":[45.25,45.5,45.75,45.75], "lon":[-73.5,-73.5,-73.5,-73.75]}
        >>> interpolate_caspar = interpolate_data(var=data_caspar["var"],lat=data_caspar["lat"],lon=data_caspar["lon"],locations=locations,bbox=bbox,post_process=False,silent=True)
        >>> print("Interpolated data location #1: (", interpolate_caspar['lat'][0],",",interpolate_caspar['lon'][0],") --> time[0:2] = ",interpolate_caspar['var'][:4,0])
        Interpolated data location #1: ( 45.25 , -73.5 ) --> time[0:2] =  [0. 0. 0. 0.]
        >>> print("Interpolated data location #2: (", interpolate_caspar['lat'][1],",",interpolate_caspar['lon'][1],") --> time[0:2] = ",interpolate_caspar['var'][:4,1])
        Interpolated data location #2: ( 45.5 , -73.5 ) --> time[0:2] =  [0.03943589 0.25160196 0.32552247 0.0252854 ]
        >>> print("Interpolated data location #3: (", interpolate_caspar['lat'][2],",",interpolate_caspar['lon'][2],") --> time[0:2] = ",interpolate_caspar['var'][:4,2])
        Interpolated data location #3: ( 45.75 , -73.5 ) --> time[0:2] =  [0.00069719 0.03622241 0.20670044 0.00618395]
        >>> print("Interpolated data location #4: (", interpolate_caspar['lat'][3],",",interpolate_caspar['lon'][3],") --> time[0:2] = ",interpolate_caspar['var'][:4,3])
        Interpolated data location #4: ( 45.75 , -73.75 ) --> time[0:2] =  [9.51688978e-05 1.66070234e-02 1.09431941e-01 5.25377496e-02]


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
    if var is None:
        raise ValueError("interpolate_data: variable needs to be specified")
    if lat is None:
        raise ValueError("interpolate_data: latitude needs to be specified")
    if lon is None:
        raise ValueError("interpolate_data: longitude needs to be specified")
    if locations is None:
        raise ValueError("interpolate_data: locations needs to be specified")

    location_was_scalar = False
    if type(locations["lat"]) == float:
        # scalar specified --> make sure it's a list
        locations["lon"] = [ locations["lon"] ]
        locations["lat"] = [ locations["lat"] ]
        location_was_scalar = True

    if len(locations["lat"]) != len(locations["lon"]):
        raise ValueError("interpolate_data: Lat (n={}) and lon (n={}) of locations are not of the same shape.".format(len(locations["lat"]),len(locations["lon"])))
    else:
        nlocations = len(locations["lat"])

    if len(np.shape(var)) == 3: # several time steps for variable provided
        # make sure everything is an np array
        var  = np.array( var )

        # dimensions
        ntime = np.shape(var)[0]
        nlat = np.shape(var)[1]
        nlon = np.shape(var)[2]

        if list(np.shape(lat)) != [nlat,nlon]:
            raise ValueError("interpolate_data: Shape of 'lat' ({}) not matching shape of 'var' ({}).".format(np.shape(lat),np.shape(var)[1:]))
        if list(np.shape(lon)) != [nlat,nlon]:
            raise ValueError("interpolate_data: Shape of 'lon' ({}) not matching shape of 'var' ({}).".format(np.shape(lon),np.shape(var)[1:]))

    elif len(np.shape(var)) == 2: # only one time step provided
        # if only one date provided make a dummy dimension is added where needed
        var  = np.array( [var] )

        # dimensions
        ntime = 1
        nlat = np.shape(var)[1]
        nlon = np.shape(var)[2]

        if list(np.shape(lat)) != [nlat,nlon]:
            raise ValueError("interpolate_data: Shape of 'lat' ({}) not matching shape of 'var' ({}).".format(np.shape(lat),np.shape(var)[1:]))
        if list(np.shape(lon)) != [nlat,nlon]:
            raise ValueError("interpolate_data: Shape of 'lon' ({}) not matching shape of 'var' ({}).".format(np.shape(lon),np.shape(var)[1:]))

    else:
        raise ValueError("interpolate_data: Unknown shape of variable. Needs to be 2D or 3D.")

    # check if locations are in bounding box
    if not(bbox is None):

        for ilocation in range(nlocations):

            if (        (locations["lat"][ilocation] < bbox["lat"]["min"]) or
                        (locations["lat"][ilocation] > bbox["lat"]["max"]) or
                        (locations["lon"][ilocation] < bbox["lon"]["min"]) or
                        (locations["lon"][ilocation] > bbox["lon"]["max"])):
                raise ValueError("interpolate_data: location (lat={},lon={}) is not within the bounding box {}. The data might still contain the grid cell. Try to run without using 'bbox' argument in 'interplolate_data'. The function will crash if the location is also not available in the dataset.".format(locations["lat"][ilocation],locations["lon"][ilocation],bbox))



    if not(silent): print("nlat=",nlat)
    if not(silent): print("nlon=",nlon)

    # initialize return
    result = {}
    interpol_all = []
    if return_tmp:
        results_tmp = {'lat':[],'lon':[],'var':[]}

    if not(silent): print("nlocations = ",nlocations)

    for ilocation in range(nlocations):

        [ii, jj] = get_index_of_latlon_location(lat, lon, locations["lat"][ilocation], locations["lon"][ilocation])

        if ii < 1 or ii > nlat-1:
            raise ValueError("interpolate_data: Location too close to edge. (lat, ii={})".format(ii))
        if jj < 1 or jj > nlon-1:
            raise ValueError("interpolate_data: Location too close to edge. (lon, jj={})".format(jj))

        neighbors_lon = np.array([ [ lon[iii,jjj] for jjj in range(jj-1,jj+2) ] for iii in range(ii-1,ii+2) ])
        neighbors_lat = np.array([ [ lat[iii,jjj] for jjj in range(jj-1,jj+2) ] for iii in range(ii-1,ii+2) ])
        neighbors_var = np.array([ [ [ var[itime,iii,jjj] for jjj in range(jj-1,jj+2) ] for iii in range(ii-1,ii+2) ] for itime in range(ntime) ])

        # save variables that are used for interpolation
        if return_tmp:
            results_tmp['lat'].append(neighbors_lat)
            results_tmp['lon'].append(neighbors_lon)
            results_tmp['var'].append(neighbors_var)

        xx = neighbors_lon.flatten()
        yy = neighbors_lat.flatten()

        if not(silent): print("lat_loc = ",locations["lat"][ilocation]," ---> idx = [row=",ii,",col=",jj,"]")
        if not(silent): print("lon_loc = ",locations["lon"][ilocation]," ---> idx = [row=",ii,",col=",jj,"]")
        if not(silent): print("neighbor_lat = ",neighbors_lat)
        if not(silent): print("neighbor_lon = ",neighbors_lon)

        interpol = []
        for itime in range(ntime):

            xy = np.asarray([(ii,jj) for ii, jj in zip(xx,yy)])
            zz = neighbors_var[itime].flatten()
            ff = scipy.interpolate.LinearNDInterpolator(xy, zz)
            izz = ff(locations["lon"][ilocation], locations["lat"][ilocation])

            if post_process:
                # check if value is "too negative"
                if izz < -0.1:
                    raise ValueError("interpolate_data: very negative value obtained for interpolated data for location #{} (lat={},lon={}): value = {}. It seems that the interpolation is not working properly as it is based on 9 values that are all non-negative and therefore shouldn't produce large negative interpolated values.".format(ilocation, locations["lat"][ilocation], locations["lon"][ilocation], izz))

                # set negative value to zero
                if izz < 0.0:
                    izz = 0.0

            interpol.append( izz ) # interpol.append( np.diagonal(ii) )

        interpol_all.append(interpol)

    interpol_all = np.transpose(np.array(interpol_all))   # shape (ntime,nloc)
    if location_was_scalar:
        interpol_all = interpol_all[:,:]

    result['var'] = interpol_all
    result['lat'] = locations["lat"]
    result['lon'] = locations["lon"]
    if return_tmp:
        results_tmp['lat'] = np.array(results_tmp['lat'])
        results_tmp['lon'] = np.array(results_tmp['lon'])
        results_tmp['var'] = np.array(results_tmp['var'])
        result['tmp'] = results_tmp

    # ----------------------------------------------------
    # Done.
    # ----------------------------------------------------
    return result


def get_index_of_latlon_location(lat_2d, lon_2d, lat_location, lon_location):

    # determines the index of the location (lat_location, lon_location) in the 2D map (lat_2d, lon_2d)
    #
    # return: [row_idx, col_idx]
    #
    # restriction: all longitudes are in range [-180, 180]
    #
    dist   = (np.abs( lat_2d - lat_location ))**2 + np.abs( lon_2d - lon_location )**2
    minpos = np.transpose(np.array(np.where(dist == np.min(dist))))[0]

    return minpos



if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # # --------------------------------------
    # # Request and read data from Geomet
    # # --------------------------------------

    # from a1_request_geomet_grib2 import request_geomet_grib2
    # from b1_read_geomet_grib2 import read_geomet_grib2

    # product = 'rdpa:10km:24f'
    # crs = 'EPSG:4326'
    # bbox = {"lat":{"min":45.0,"max":46.0},"lon":{"min":-74.0,"max":-73.0}}

    # # four dates 6h apart from each other
    # dates = [ datetime.datetime(2022,8,24,12,0) + datetime.timedelta(hours=24*ii) for ii in range(2) ]

    # # request data
    # filename = '/tmp/test'
    # files_geomet = request_geomet_grib2(product=product,date=dates,bbox=bbox,crs=crs,filename=filename,silent=True)

    # # read data
    # data_geomet = read_geomet_grib2(files_geomet,silent=True)

    # # --------------------------------------
    # # Interpolate data from Geomet
    # # --------------------------------------

    # # 4 locations
    # locations = {"lat":[45.25,45.5,45.75,45.75], "lon":[-73.5,-73.5,-73.5,-73.75]}
    # interpolate_geomet = interpolate_data(var=data_geomet["var"],lat=data_geomet["lat"],lon=data_geomet["lon"],locations=locations,bbox=bbox,silent=True)
    # print("Interpolated data location #1: (", interpolate_geomet['lat'][0],",",interpolate_geomet['lon'][0],") --> time[0:2] = ",interpolate_geomet['var'][:,0])
    # print("Interpolated data location #2: (", interpolate_geomet['lat'][1],",",interpolate_geomet['lon'][1],") --> time[0:2] = ",interpolate_geomet['var'][:,1])
    # print("Interpolated data location #3: (", interpolate_geomet['lat'][2],",",interpolate_geomet['lon'][2],") --> time[0:2] = ",interpolate_geomet['var'][:,2])
    # print("Interpolated data location #4: (", interpolate_geomet['lat'][3],",",interpolate_geomet['lon'][3],") --> time[0:2] = ",interpolate_geomet['var'][:,3])

    # # 1 location (as list)
    # locations = {"lat":[45.25], "lon":[-73.5]}
    # interpolate_geomet = interpolate_data(var=data_geomet["var"],lat=data_geomet["lat"],lon=data_geomet["lon"],locations=locations,bbox=bbox,silent=True)
    # print("Interpolated data location #1: (", interpolate_geomet['lat'][0],",",interpolate_geomet['lon'][0],") --> time[0:2] = ",interpolate_geomet['var'][:,0])

    # # 1 location (as scalar)
    # locations = {"lat":45.25, "lon":-73.5}
    # interpolate_geomet = interpolate_data(var=data_geomet["var"],lat=data_geomet["lat"],lon=data_geomet["lon"],locations=locations,bbox=bbox,silent=True)
    # print("Interpolated data location #1: (", interpolate_geomet['lat'][0],",",interpolate_geomet['lon'][0],") --> time[0:2] = ",interpolate_geomet['var'][:,0])

    # # --------------------------------------
    # # Request and read data from CaSPAr
    # # --------------------------------------

    # from b2_read_caspar_nc import read_caspar_nc
    # product='RDRS_v2.1'
    # variable='RDRS_v2.1_A_PR0_SFC'
    # dates=[ datetime.datetime(2018,8,8,15,0) + datetime.timedelta(hours=ii) for ii in range(10) ]
    # bbox={"lat":{"min":45.0,"max":46.0},"lon":{"min":-74.0,"max":-73.0}}
    # foldername='test-data/'
    # lintransform={'a':1000.0,'b':0.0}  # convert from m/h to mm/h
    # silent=False

    # data_caspar = read_caspar_nc( product=product, variable=variable, date=dates, bbox=bbox, foldername=foldername, lintransform=lintransform, silent=True)

    # # --------------------------------------
    # # Interpolate data from CaSPAr
    # # --------------------------------------

    # # 4 locations
    # locations = {"lat":[45.25,45.5,45.75,45.75], "lon":[-73.5,-73.5,-73.5,-73.75]}
    # interpolate_caspar = interpolate_data(var=data_caspar["var"],lat=data_caspar["lat"],lon=data_caspar["lon"],locations=locations,bbox=bbox,silent=True)
    # print("Interpolated data location #1: (", interpolate_caspar['lat'][0],",",interpolate_caspar['lon'][0],") --> time[0:2] = ",interpolate_caspar['var'][:4,0])
    # print("Interpolated data location #2: (", interpolate_caspar['lat'][1],",",interpolate_caspar['lon'][1],") --> time[0:2] = ",interpolate_caspar['var'][:4,1])
    # print("Interpolated data location #3: (", interpolate_caspar['lat'][2],",",interpolate_caspar['lon'][2],") --> time[0:2] = ",interpolate_caspar['var'][:4,2])
    # print("Interpolated data location #4: (", interpolate_caspar['lat'][3],",",interpolate_caspar['lon'][3],") --> time[0:2] = ",interpolate_caspar['var'][:4,3])

    # # 1 location (as list)
    # locations = {"lat":[45.25], "lon":[-73.5]}
    # interpolate_caspar = interpolate_data(var=data_caspar["var"],lat=data_caspar["lat"],lon=data_caspar["lon"],locations=locations,bbox=bbox,silent=True)
    # print("Interpolated data location #1: (", interpolate_caspar['lat'][0],",",interpolate_caspar['lon'][0],") --> time[0:2] = ",interpolate_caspar['var'][:4,0])

    # # 1 location (as scalar)
    # locations = {"lat":45.25, "lon":-73.5}
    # interpolate_caspar = interpolate_data(var=data_caspar["var"],lat=data_caspar["lat"],lon=data_caspar["lon"],locations=locations,bbox=bbox,return_tmp=True,silent=True)
    # print("Interpolated data location #1: (", interpolate_caspar['lat'][0],",",interpolate_caspar['lon'][0],") --> time[0:2] = ",interpolate_caspar['var'][:4,0])
