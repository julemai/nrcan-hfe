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
from pathlib import Path


__all__ = ['interpolate_data', 'plot_interpolated']


def interpolate_data(var=None,lat=None,lon=None,locations=None,bbox=None,return_tmp=False,post_process=False,silent=True):
    """
        Interpolate gridded data to specific location(s).

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
        {var,lat,lon,                  consistency, also "lat" and "lon" at requested points are returned but
        tmp}                           they are the same as the input "lat" and "lon" variables. "tmp" will only be returned
                                       if "return_tmp" is set to True.


        Description
        -----------
        Interpolates variable (var) in space at locations (locations) using bilinear interpolation. The base of this
        function is "scipy.interpolate.LinearNDInterpolator". The grid cell containing the exact location as well as the 8
        neighboring cells are used for interpolation.


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



def plot_interpolated(locations=None,dates=None,data=None,start_date_buffer=None,end_date_buffer=None,pngfile=None,highlight_dates_idx=None,start_date=None,end_date=None,label=None,silent=True):

    """
        Plot interpolated data.

        Definition
        ----------
        def plot_interpolated(locations=None,dates=None,data=None,start_date_buffer=None,end_date_buffer=None,pngfile=None,
                              highlight_dates_idx=None,start_date=None,end_date=None,label=None,silent=True)


        Input           Format         Description
        -----           -----          -----------
        locations       dict           Dictionary providing attributes "lat" as list of latitudes and "lon" as list
                                       of "longitudes" for locations where raster shall be interpolated. Lists of
                                       "lat" and "lon" are assumed to be of the same length. If only one location
                                       is requested, instead of lists a scalar can be provided
                                       Default: None

        dates           list(datetime) List of datetime objects specifying the time the variable ("data") is valid
                                       for.
                                       Default: None

        data            dict           Interpolated data as returned from "interpolate_data()" which is a dictionary
                                       containing keys "var", "lat", and "lon".
                                       Default: None

        start_date_buffer  datetime    Start date of event (e.g., date that is specified in HFE database) including
                                       buffer considered to be used to draw vertical dashed line.
                                       Default: None

        end_date_buffer    datetime    End date of event (e.g., date that is specified in HFE database) including
                                       buffer considered to be used to draw vertical dashed line.
                                       Default: None

        highlight_dates_idx  list      List of indexes of dates ("date") that should be highlighted in plot (e.g.,
                                       because those are the time steps that were considered the major precipitation
                                       event).
                                       Default: None

        pngfile         string         Name of PNG file to create
                                       Default: None

        start_date      datetime       Start date of event (e.g., date that is specified in HFE database) to be used
                                       to draw vertical dotted line.
                                       Default: None

        end_date        datetime       End date of event (e.g., date that is specified in HFE database) to be used
                                       to draw vertical dotted line.
                                       Default: None

        label           string         Label to be added below the figure to contain, for example, station name or
                                       additional custom information.
                                       Default: None

        silent          Boolean        If set to True, nothing will be printed to terminal.
                                       Default: True

        Output          Format         Description
        -----           -----          -----------
        filenames       dict           Dictionary containing files produced. Key is "png".


        Description
        -----------
        Plots interpolated data. Can be interpolated data at several locations. All will be plotted in one figure.
        Output format is PNG.


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
        >>> from dx_interpolate_data import interpolate_data

        >>> product = 'rdpa:10km:6f'
        >>> crs = 'EPSG:4326'
        >>> bbox = {"lat":{"min":45.0,"max":46.0},"lon":{"min":-74.0,"max":-73.0}}

        >>> # four dates 6h apart from each other
        >>> dates = [ datetime.datetime(2022,8,24,12,0) + datetime.timedelta(hours=6*ii) for ii in range(4) ]

        >>> # request data
        >>> filename = '/tmp/pytest_rdpa_10km_6f'
        >>> files_geomet = request_geomet_grib2(product=product,date=dates,bbox=bbox,crs=crs,filename=filename,silent=True)

        >>> # read data
        >>> data_geomet = read_geomet_grib2(filenames=files_geomet,silent=True)

        >>> # interpolate data
        >>> locations = {"lat":[45.25,45.5,45.75,45.75], "lon":[-73.5,-73.5,-73.5,-73.75]}
        >>> interpolate_geomet = interpolate_data(var=data_geomet["var"],lat=data_geomet["lat"],lon=data_geomet["lon"],locations=locations,bbox=bbox,post_process=True,silent=True)

        >>> # --------------------------------------
        >>> # Plot interpolated data
        >>> # --------------------------------------

        >>> filenames = plot_interpolated(locations=locations, dates=dates, data=interpolate_geomet, start_date_buffer=dates[0], end_date_buffer=dates[-1], pngfile='/tmp/test-interpolate-data.png',silent=True)
        >>> print("PNG files interpolated data: ", filenames['png'])
        PNG files interpolated data: ['/tmp/test-interpolate-data.png']


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
    if (locations is None):
        raise ValueError("plot_interpolated: location needs to be specified")
    else:
        nlocations = np.shape(locations['lon'])[0]
    if (dates is None):
        raise ValueError("plot_interpolated: dates needs to be specified")
    if (data is None):
        raise ValueError("plot_interpolated: data needs to be specified")
    if not(highlight_dates_idx is None):
        if len(highlight_dates_idx) != nlocations:
            raise ValueError("plot_interpolated: if values at specific dates should be highlighted in plot, the indexes of those dates need to specified for each location individually")
    if (pngfile is None):
        raise ValueError("plot_interpolated: pngfile needs to be specified")
    if (start_date_buffer is None):
        raise ValueError("plot_interpolated: start_date_buffer needs to be specified")
    if (end_date_buffer is None):
        raise ValueError("plot_interpolated: end_date_buffer needs to be specified")

    # initialize return
    result = {}
    result["png"] = []

    # make sure inputs are of right type
    dates = np.array(dates)

    # make sure missing time steps are masked
    timedelta = (
        ((dates[1]-dates[0]).days)*24+
        ((dates[1]-dates[0]).seconds)/60/60)

    ntime = int(((end_date_buffer-start_date_buffer).days*24+(end_date_buffer-start_date_buffer).seconds/60/60)/timedelta)+1
    dates_all = np.array([start_date_buffer + datetime.timedelta(hours=timedelta*ii) for ii in range(ntime)])

    var = np.ones([ntime,nlocations]) * -9999.
    for iitime,itime in enumerate(dates):
        tidx = list(dates_all).index(itime)
        var[tidx] = data['var'][iitime,:]
    data_var_ma = np.ma.array(var,mask=(var==-9999.))
    data_all = {'lat':data['lat'], 'lon':data['lon'],'var':data_var_ma}

    if not(highlight_dates_idx is None):
        highlight_dates_idx_all = [ [ list(dates_all).index(dates[iidx]) for iidx in iidx_loc ] for iidx_loc in highlight_dates_idx ]

    # make sure folder to store file exists; otherwise create
    Path(pngfile).parent.mkdir(parents=True, exist_ok=True)

    import matplotlib as mpl
    import matplotlib.pyplot as plt


    # -------------------------------------------------------------------------
    # create plot
    # -------------------------------------------------------------------------

    # plot settings
    textsize    = 10          # standard text size
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

    llxbbox     = 0.5         # x-anchor legend bounding box
    llybbox     = 1.0        # y-anchor legend bounding box
    llrspace    = 0.25        # spacing between rows in legend
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
    ifig += 1
    iplot = 0
    if not(silent): print('     Plot - Fig ', ifig, ' ::  ',pngfile)
    fig = plt.figure(ifig)

    # -------------------------------------------------------------------------
    # Create line plots on figure
    # -------------------------------------------------------------------------

    sub    = fig.add_axes( [0.0,0.0,1.0,0.25] )  # [left, bottom, width, height]

    if not(label is None):
        sub.text( 0.5, -0.25, label,
                    ha = 'center', va = 'top',
                    fontweight='bold',
                    transform=sub.transAxes,
                    fontsize=textsize+2 )

    # -------------------------------------------------------------------------
    # Regular step-wise plot (for entire time series)
    # -------------------------------------------------------------------------
    ilabel = [ "Loc #"+str(iloc+1)+"  ("+str(locations["lat"][iloc])+","+str(locations["lon"][iloc])+")" for iloc in range(len(locations["lat"])) ]
    sub.step(dates_all,data_all['var']/timedelta,color='0.4',where='pre',linewidth=1.0*lwidth,linestyle='-',zorder=100) #label=ilabel,

    xmin, xmax = sub.get_xlim()
    ymin, ymax = sub.get_ylim()

    # -------------------------------------------------------------------------
    # Plot some vertical lines to visualize time period requested (incl. buffer) and start and end date of event
    # -------------------------------------------------------------------------
    sub.plot([start_date_buffer,start_date_buffer],[ymin,ymax],label='Start date incl. buffer',color='0.8',linewidth=1.0*lwidth,linestyle='--',zorder=50)

    if not(start_date is None):
        sub.plot([start_date,start_date],[ymin,ymax],label='Event start date',color='0.8',linewidth=1.0*lwidth,linestyle='dotted',zorder=50)
    if not(end_date is None):
        sub.plot([end_date,end_date],[ymin,ymax],label='Event end date',color='0.8',linewidth=1.0*lwidth,linestyle='dotted',zorder=50)

    sub.plot([end_date_buffer,end_date_buffer],[ymin,ymax],label='End date incl. buffer',color='0.8',linewidth=1.0*lwidth,linestyle='--',zorder=50)

    # -------------------------------------------------------------------------
    # Highlight some time steps if requested
    # -------------------------------------------------------------------------
    if not(highlight_dates_idx is None):

        for ilocation in range(nlocations):
            if ilocation == 0:
                ilabel = "Event precip. considered"
            else:
                ilabel = ""

            if len(highlight_dates_idx_all[ilocation]) > 0:
                # set all missing values to 0.0 such that fill-between looks appropriate (i.e., fill under step-plot)
                tmp = data_all['var'][np.min(highlight_dates_idx_all[ilocation]):np.max(highlight_dates_idx_all[ilocation])+1,ilocation]/timedelta
                tmp[tmp.mask] = 0.0

                sub.fill_between(    dates_all[np.min(highlight_dates_idx_all[ilocation]):np.max(highlight_dates_idx_all[ilocation])+1],
                                     tmp,
                                     color=ocean_color, step="pre", interpolate=True,alpha=0.4,linewidth=1.0*lwidth,linestyle='-',zorder=50, label=ilabel)

    # -------------------------------------------------------------------------
    # Format date axis
    # -------------------------------------------------------------------------
    import matplotlib.dates as mdates
    from matplotlib.ticker import FormatStrFormatter
    myFmt = mdates.DateFormatter("%d %h\n%Y\n%H:%M")
    sub.xaxis.set_major_formatter(myFmt)
    sub.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    # -------------------------------------------------------------------------
    # Label on axes
    # -------------------------------------------------------------------------
    plt.setp(sub,  ylabel="Precipitation Rate [mm h$^{-1}$]")   # [mm/6h]

    # -------------------------------------------------------------------------
    # plot legend
    # -------------------------------------------------------------------------
    ll = sub.legend(frameon=frameon, ncol=5,
                            labelspacing=llrspace, handletextpad=llhtextpad, handlelength=llhlength,
                            loc='lower center', bbox_to_anchor=(llxbbox,llybbox), scatterpoints=1, numpoints=1,
                            fontsize = 'small')

    # -------------------------------------------------------------------------
    # Done.
    # -------------------------------------------------------------------------

    fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
    plt.close(fig)

    # add file to return
    result["png"].append(pngfile)


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





    # # --------------------------------------
    # # Plotting
    # # --------------------------------------

    # from a1_request_geomet_grib2 import request_geomet_grib2
    # from b1_read_geomet_grib2 import read_geomet_grib2
    # from dx_interpolate_data import interpolate_data, plot_interpolated

    # product = 'rdpa:10km:6f'
    # crs = 'EPSG:4326'
    # bbox = {"lat":{"min":45.0,"max":46.0},"lon":{"min":-74.0,"max":-73.0}}

    # # four dates 6h apart from each other
    # dates = [ datetime.datetime(2022,8,24,12,0) + datetime.timedelta(hours=24*ii) for ii in range(2) ]

    # # request data
    # filename = '/tmp/pytest_rdpa_10km_6f'
    # files_geomet = request_geomet_grib2(product=product,date=dates,bbox=bbox,crs=crs,filename=filename,silent=True)

    # # read data
    # data_geomet = read_geomet_grib2(filenames=files_geomet,silent=True)

    # # interpolate data
    # locations = {"lat":[45.25,45.5,45.75,45.75], "lon":[-73.5,-73.5,-73.5,-73.75]}
    # interpolate_geomet = interpolate_data(var=data_geomet["var"],lat=data_geomet["lat"],lon=data_geomet["lon"],locations=locations,bbox=bbox,post_process=True,silent=True)

    # filenames = plot_interpolated(locations=locations, dates=dates, data=interpolate_geomet, start_date_buffer=dates[0], end_date_buffer=dates[-1], pngfile='/tmp/test-interpolate-data.png',silent=True)
    # print("PNG files interpolated data: ", filenames['png'])
