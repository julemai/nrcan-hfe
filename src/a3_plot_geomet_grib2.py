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
import json as json
from pathlib import Path

# to create PNGs and legend
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# to create GIF
from PIL import Image


__all__ = ['plot_geomet_grib2']

def plot_geomet_grib2(var=None,lat=None,lon=None,date=None,png=True,gif=False,legend=False,cities=True,basefilename='/tmp/test',silent=True):
    """
        Read GRIB2 retrieved from Geomet .

        Definition
        ----------
        def read_geomet_grib2(filename)


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

        date            datetime or    Datetime object specifying the time the variable is valid for. If "var"
                        list(datetime) is provided for several time steps, a list of datetime objects needs to be provided.
                                       Default: None

        png             Boolean        If true, a PNG file for each time step will be created. If false,
                                       no PNG file will be created. Value will be ignored in case "gif" is true as the
                                       PNGs are required to build the GIF.
                                       Default: True

        gif             Boolean        If true, an animated GIF file is created. If false, no GIF file will be created.
                                       Creating a GIF makes most sense when variable "var" is given for multiple time steps
                                       (i.e., var is #D and date is a list of datetime objects). The time steps will be sorted
                                       before creating the GIF.
                                       Default: False

        legend          Boolean        If set to true, an additional PNG will be created containing the legend of the plots.
                                       The PNGs/GIF will not contain any legend themselves. The legend is the same as used by
                                       Geomet and returned when requesting, for example, "https://geo.weather.gc.ca/geomet? \
                                       SERVICE=WMS&VERSION=1.3.0&REQUEST=GetLegendGraphic&LAYERS=RDPA.6F_PR&STYLES=RDPA-WXO& \
                                       CRS=EPSG:4326&BBOX=45,-74,46,-73&WIDTH=400&HEIGHT=400&FORMAT=image/png".
                                       Default: False

        cities          Boolean        If set to true, the map will contain cities that are located in the domain the data are
                                       provided. The list of cities is provided through "lib/canadian-cities.csv". The setting
                                       of true/false will have no effect if there is no city located in the lat/lon domain.
                                       Default: True

        basefilename    string         Basename of files created. The basename will be appended by:
                                       - PNG files: datetime in format YYYYMMDDHH
                                         (i.e, default = "/tmp/test_<YYYYMMDDHH>.png")
                                       - Legend file: "_legend"
                                         (i.e, default = "/tmp/test_legend.png")
                                       - GIF file: Nothing
                                         (i.e, default = "/tmp/test.gif")
                                       Default: '/tmp/test'
        silent          Boolean        If set to True, nothing will be printed to terminal.

        Output          Format         Description
        -----           -----          -----------
        filenames       dict           Dictionary of all files created, e.g.:
                                          {"png": [file1, file2,...],
                                           "gif": [],
                                           "legend": [file3]}


        Description
        -----------
        Plots variable(s) on a map with the possibility to include nearby cities ("cities=True").
        The plots can be PNG files (png=True) or a GIF (gif=True).
        The latter is most useful when the variable is provided for several time steps.
        A legend can be plotted separately ("legend=True")


        Restrictions
        ------------
        None.

        Examples
        --------

        Request and read data


        >>> from a1_request_geomet_grib2 import request_geomet_grib2
        >>> from a2_read_geomet_grib2 import read_geomet_grib2

        >>> product = 'rdpa:10km:24f'
        >>> crs = 'EPSG:4326'
        >>> bbox = {"lat":{"min":45.0,"max":46.0},"lon":{"min":-74.0,"max":-73.0}}

        >>> # four dates 6h apart from each other
        >>> dates = [ datetime.datetime(2022,8,24,12,0) + datetime.timedelta(hours=24*ii) for ii in range(2) ]

        >>> # request data
        >>> filename = '/tmp/test'
        >>> files_geomet = request_geomet_grib2(product=product,date=dates,bbox=bbox,crs=crs,filename=filename,silent=True)

        >>> # read data
        >>> data_geomet = read_geomet_grib2(files_geomet,silent=True)

        >>> # plot (2 time steps)
        >>> plot_geomet = plot_geomet_grib2(var=data_geomet["var"],lat=data_geomet["lat"],lon=data_geomet["lon"],date=dates,png=True,gif=True,legend=True,cities=True,basefilename='/tmp/test_4dates',silent=True)
        >>> print("All PNG files created: ", plot_geomet['png'])
        All PNG files created:  ['/tmp/test_4dates_2022082412.png', '/tmp/test_4dates_2022082512.png']

        >>> # plot (1 time step)
        >>> plot_geomet = plot_geomet_grib2(var=data_geomet["var"][0],lat=data_geomet["lat"],lon=data_geomet["lon"],date=dates[0],png=True,gif=False,legend=False,cities=True,basefilename='/tmp/test_1date',silent=True)
        >>> print("All PNG files created: ", plot_geomet['png'])
        All PNG files created:  ['/tmp/test_1date_2022082412.png']


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
        raise ValueError("plot_geomet_grib2: variable needs to be specified")
    if lat is None:
        raise ValueError("plot_geomet_grib2: latitude needs to be specified")
    if lon is None:
        raise ValueError("plot_geomet_grib2: longitude needs to be specified")
    if date is None:
        raise ValueError("plot_geomet_grib2: date(s) needs to be specified")

    if len(np.shape(var)) == 3: # several time steps for variable provided
        # make sure everything is an np array
        var  = np.array( var )
        date = np.array( date )

        # dimensions
        ntime = np.shape(var)[0]
        nlat = np.shape(var)[1]
        nlon = np.shape(var)[2]

        if list(np.shape(lat)) != [nlat,nlon]:
            raise ValueError("Shape of 'lat' ({}) not matching shape of 'var' ({}).".format(np.shape(lat),np.shape(var)[1:]))
        if list(np.shape(lon)) != [nlat,nlon]:
            raise ValueError("Shape of 'lon' ({}) not matching shape of 'var' ({}).".format(np.shape(lon),np.shape(var)[1:]))
        if list(np.shape(date)) != [ntime]:
            raise ValueError("Shape of 'date' ({}) not matching time steps in 'var' ({}).".format(np.shape(date),np.shape(var)[0]))

    elif len(np.shape(var)) == 2: # only one time step provided
        # if only one date provided make a dummy dimension is added where needed
        var  = np.array( [var] )
        date = np.array( [date] )

        # dimensions
        ntime = np.shape(var)[0]
        nlat = np.shape(var)[1]
        nlon = np.shape(var)[2]

        if list(np.shape(lat)) != [nlat,nlon]:
            raise ValueError("Shape of 'lat' ({}) not matching shape of 'var' ({}).".format(np.shape(lat),np.shape(var)[1:]))
        if list(np.shape(lon)) != [nlat,nlon]:
            raise ValueError("Shape of 'lon' ({}) not matching shape of 'var' ({}).".format(np.shape(lon),np.shape(var)[1:]))
        if list(np.shape(date)) != [ntime]:
            raise ValueError("Shape of 'date' ({}) not matching time steps in 'var' ({}).".format(np.shape(date),0))

    else:
        raise ValueError("Unknown shape of variable. Needs to be 2D or 3D.")

    # overwrite png in case gif is set to True
    if gif:
        png = True

    # initialize return
    result = {}
    result["png"] = []
    result["gif"] = []
    result["legend"] = []

    # plot settings
    textsize    = 12          # standard text size
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

    # ----------------------------------------------------
    # Read and prepare data if PNG/GIF files will be created
    # ----------------------------------------------------
    if png or gif:

        # boundaries between lats and lons
        lonh = np.empty((nlat+1,nlon+1), dtype=float)
        lath = np.empty((nlat+1,nlon+1), dtype=float)

        tmp = [ [ (lat[ii+1,jj+1]-lat[ii,jj])/2 for jj in range(nlon-1) ] + [ (lat[ii+1,nlon-1]-lat[ii,nlon-2])/2 ] for ii in range(nlat-1) ]
        dlat = np.array(tmp + [ tmp[-1] ])

        tmp = [ [ (lon[ii+1,jj+1]-lon[ii,jj])/2 for jj in range(nlon-1) ] + [ (lon[ii+1,nlon-1]-lon[ii,nlon-2])/2 ] for ii in range(nlat-1) ]
        dlon = np.array(tmp + [ tmp[-1] ])

        lonh[0:nlat,0:nlon] = lon - dlon
        lath[0:nlat,0:nlon] = lat - dlat

        # make lat and lon one column and row wider such that all
        lonh[nlat,0:nlon] = lonh[nlat-1,0:nlon] + (lonh[nlat-1,0:nlon] - lonh[nlat-2,0:nlon])
        lath[nlat,0:nlon] = lath[nlat-1,0:nlon] + (lath[nlat-1,0:nlon] - lath[nlat-2,0:nlon])
        lonh[0:nlat,nlon] = lonh[0:nlat,nlon-1] + (lonh[0:nlat,nlon-1] - lonh[0:nlat,nlon-2])
        lath[0:nlat,nlon] = lath[0:nlat,nlon-1] + (lath[0:nlat,nlon-1] - lath[0:nlat,nlon-2])
        lonh[nlat,nlon]   = lonh[nlat-1,nlon-1] + (lonh[nlat-1,nlon-1] - lonh[nlat-2,nlon-2])
        lath[nlat,nlon]   = lath[nlat-1,nlon-1] + (lath[nlat-1,nlon-1] - lath[nlat-2,nlon-2])

        # bounding box
        llcrnrlon_raw =  np.min(lonh)
        urcrnrlon_raw =  np.max(lonh)
        llcrnrlat_raw =  np.min(lath)
        urcrnrlat_raw =  np.max(lath)
        llcrnrlon_buf =  np.min(lonh) - (np.max(lonh) - np.min(lonh))*.05
        urcrnrlon_buf =  np.max(lonh) + (np.max(lonh) - np.min(lonh))*.05
        llcrnrlat_buf =  np.min(lath) - (np.max(lath) - np.min(lath))*.05
        urcrnrlat_buf =  np.max(lath) + (np.max(lath) - np.min(lath))*.05

        # load cities if needed
        if cities:

            if not( Path('lib/canadian-cities.json').is_file() ):
                dict_cities = {}
                file_cities = 'lib/canadian-cities.csv'
                ff = open(file_cities, "r")
                lines = ff.readlines()
                ff.close()

                for ll in lines[1:]:
                    tmp = ll.strip().split(';')
                    tmp2 = tmp[0].split(',')
                    city = ','.join(tmp2[0:-2])
                    province = tmp2[-2]
                    country = tmp2[-1]
                    dict_cities[city] = {"province": province, "country":country, "lat":float(tmp[1]),"lon":float(tmp[2])}

                json_dump = json.dumps(dict_cities)
                ff = open('lib/canadian-cities.json', "w")
                ff.write(json_dump)
                ff.close()


            with open('lib/canadian-cities.json', 'r') as ff:
                dict_cities = json.load(ff)

            # find all cities in bounding box
            city_inregion = {}
            for cc in dict_cities:
                if (    (dict_cities[cc]["lat"] > llcrnrlat_raw) and
                        (dict_cities[cc]["lat"] < urcrnrlat_raw) and
                        (dict_cities[cc]["lon"] > llcrnrlon_raw) and
                        (dict_cities[cc]["lon"] < urcrnrlon_raw)):
                    city_inregion[cc] = dict_cities[cc]

            if not(silent): print("Found {} cities in region.".format(len(city_inregion)))


    # ----------------------------------------------------
    # Create PNG files
    # ----------------------------------------------------
    if png:

        if not(silent): print('Create PNG ')

        for iidate,idate in enumerate(date):

            timestep_title = idate.strftime('%d %h %Y %H:%M:%S')+' UTC'   # '02 Oct 2017 18:00:00 UTC'
            timestep_filename = idate.strftime("%Y%m%d%H")                # '2017100218'
            pngfile = "{0}_{1}.png".format(basefilename,timestep_filename)

            # -------------------------------------------------------------------------
            # Create figure object
            # -------------------------------------------------------------------------
            ifig += 1
            iplot = 0
            if not(silent): print('     Plot - Fig ', ifig, ' ::  ',pngfile)
            fig = plt.figure(ifig)

            # -------------------------------------------------------------------------
            # Create map on figure
            # -------------------------------------------------------------------------
            iplot += 1
            sub    = fig.add_axes( [0.0,0.0,1.0,1.0] )  # [left, bottom, width, height]

            # Basemap
            lat_1     =  (llcrnrlat_buf+urcrnrlat_buf)/2  # first  "equator"
            lat_2     =  (llcrnrlat_buf+urcrnrlat_buf)/2  # second "equator"
            lat_0     =  (llcrnrlat_buf+urcrnrlat_buf)/2  # center of the map
            lon_0     =  (llcrnrlon_buf+urcrnrlon_buf)/2  # center of the map
            bmap = Basemap(projection='merc',area_thresh=2000.,
                    llcrnrlon=llcrnrlon_buf, urcrnrlon=urcrnrlon_buf, llcrnrlat=llcrnrlat_buf, urcrnrlat=urcrnrlat_buf,
                    lat_1=lat_1, lat_2=lat_2, lat_0=lat_0, lon_0=lon_0,
                    resolution='i') # Lambert conformal

            # plot coastlines
            bmap.drawcoastlines(linewidth=0.3)

            #bmap.drawmapboundary(color='black', fill_color=ocean_color, linewidth=0.3)
            #bmap.drawcountries(color='black', linewidth=0.3)
            bmap.fillcontinents(color='white', lake_color=ocean_color)

            # latitudes - draw 3 parallels, labels = [left, right, top, bottom]
            nparallels = 3.
            bmap.drawparallels(np.arange(np.round(np.min(lath),1),np.round(np.max(lath),1)+0.1,np.round((np.max(lath)+0.1-np.min(lath))/nparallels,1)),
                                   labels=[1,0,0,0], dashes=[1,1], linewidth=0.25, color='0.5')

            # longitudes - draw 3 meridians, labels = [left, right, top, bottom]
            nmeridians = 3.
            bmap.drawmeridians(np.arange(np.round(np.min(lonh),1),np.round(np.max(lonh),1)+0.1,np.round((np.max(lonh)+0.1-np.min(lonh))/nmeridians,1)),
                                   labels=[0,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5')

            # geo-referenced
            # xx, yy = bmap(lon,lat)
            xxh, yyh = bmap(lonh,lath)
            zz = var[iidate,:,:]
            variable_plot = bmap.pcolor(xxh, yyh, zz, cmap=cmap, norm=norm, zorder = 100, alpha=0.7)

            # plot cities
            if cities:

                labelshiftx = 0.0 #(urcrnrlon_raw - llcrnrlon_raw)*.01
                labelshifty = (urcrnrlat_raw - llcrnrlat_raw)*.005
                for icity in city_inregion:
                    cityx, cityy           = bmap(city_inregion[icity]["lon"],city_inregion[icity]["lat"])
                    xshift, yshift = bmap(city_inregion[icity]["lon"]+labelshiftx,city_inregion[icity]["lat"]-labelshifty)
                    sub.plot(cityx, cityy, 'ok',      markersize=3, zorder = 200)
                    sub.text(xshift, yshift, icity , fontsize=textsize-1, zorder = 200, horizontalalignment="center", verticalalignment="top")

            # set title as time step
            sub.set_title(timestep_title,fontsize=textsize+2)

            fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
            plt.close(fig)

            # add file to return
            result["png"].append(pngfile)

    # ----------------------------------------------------
    # Create legend
    # ----------------------------------------------------
    if legend:

        if not(silent): print('Create Legend ')

        legendfile = "{0}_legend.png".format(basefilename)

        # -------------------------------------------------------------------------
        # Create figure object
        # -------------------------------------------------------------------------
        ifig += 1
        iplot = 0
        if not(silent): print('     Plot - Fig ', ifig, ' ::  ',legendfile)
        fig = plt.figure(ifig)

        # -------------------------------------------------------------------------
        # Create map on figure
        # -------------------------------------------------------------------------
        iplot += 1
        # sub    = fig.add_axes( [0.0,0.0,1.0,1.0] )  # [left, bottom, width, height]

        f = lambda m,c: plt.plot([],[],marker=m, color=c, ls="none")[0]
        handles = [f("s", colors[icol]) for icol in range(len(colors))]
        labels = [ "{0:.2f} - {1:.2f}".format(bounds[ibb],bounds[ibb+1]) for ibb,bb in enumerate(bounds[:-1]) ]

        # no frame (frameon), wider handle (numpoints)
        legend = plt.legend(handles, labels, loc=3, framealpha=1, frameon=False, numpoints=4)

        fig  = legend.figure
        fig.canvas.draw()

        # no axis
        plt.axis('off')

        # limiting the bounding box to the legend area,
        # slighly extended to ensure the surrounding rounded corner box of
        # is not cropped. Transparency is enabled, so it is not an issue.
        bbox  = legend.get_window_extent().padded(2)
        bbox = bbox.transformed(fig.dpi_scale_trans.inverted())
        fig.savefig(legendfile, transparent=transparent, bbox_inches=bbox)

        plt.close(fig)

        # add file to return
        result["legend"].append(legendfile)

    # ----------------------------------------------------
    # Create GIF files (assumes that PNGs are available)
    # ----------------------------------------------------
    if gif:

        if not(silent): print('Create GIF ')
        giffile = "{0}.gif".format(basefilename)
        ifig += 1

        if not(silent): print('     Plot - Fig ', ifig, ' ::  ',pngfile)

        # Create frames
        frames = []
        pngs = np.sort(result["png"])   # makes sure files are sorted
        for png in pngs:
            frame = Image.open(png)
            frames.append(frame)

        # Save into GIF file that loops forever
        frames[0].save(giffile, format='GIF',
                       append_images=frames[1:],
                       save_all=True,
                       duration=300, loop=0)

        # add file to return
        result["gif"].append(pngfile)


    # ----------------------------------------------------
    # Done.
    # ----------------------------------------------------
    return result



if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # from a1_request_geomet_grib2 import request_geomet_grib2
    # from a2_read_geomet_grib2 import read_geomet_grib2

    # product = 'rdpa:10km:24f'
    # crs = 'EPSG:4326'
    # bbox = {"lat":{"min":45.0,"max":46.0},"lon":{"min":-74.0,"max":-73.0}}

    # # four dates 6h apart from each other
    # dates = [ datetime.datetime(2022,8,24,12,0) + datetime.timedelta(hours=24*ii) for ii in range(2) ]

    # # request data
    # filename = '/tmp/test'
    # files_geomet = request_geomet_grib2(product=product,date=dates,bbox=bbox,crs=crs,filename=filename)

    # # read data
    # data_geomet = read_geomet_grib2(files_geomet)

    # # plot (4 time steps)
    # plot_geomet = plot_geomet_grib2(    var=data_geomet["var"],
    #                       lat=data_geomet["lat"],
    #                       lon=data_geomet["lon"],
    #                       date=dates,
    #                       png=True,
    #                       gif=True,
    #                       legend=True,
    #                       cities=True,
    #                       basefilename='/tmp/test_4dates',
    #                       silent=True)
    # print("All files created: ", plot_geomet)

    # # plot (1 time step)
    # plot_geomet = plot_geomet_grib2(    var=data_geomet["var"][0],
    #                       lat=data_geomet["lat"],
    #                       lon=data_geomet["lon"],
    #                       date=dates[0],
    #                       png=True,
    #                       gif=False,
    #                       legend=False,
    #                       cities=True,
    #                       basefilename='/tmp/test_1date',
    #                       silent=True)
    # print("All files created: ", plot_geomet)
