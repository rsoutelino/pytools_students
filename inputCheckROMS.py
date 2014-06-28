#!/usr/bin/env python
import cdms2, glob, os # cdms2 comes from cdat-lite (pip install cdat-lite)
import numpy as np
from pymo.models.roms import ROMS
import datetime as dt
from mpl_toolkits.basemap import Basemap
from pymo.core.lib import dt2ncep, ncep2dt
import matplotlib.pyplot as plt
import glob

from variable_ranges_ROMS_dicts import *

# TO-DO list ################################################################
# - make a better function to plot the boundary conditions
#############################################################################


BASEDIR = "/home/rsoutelino/myroms/phd_run/"


def init_basemap(lon, lat, proj):
    m = Basemap(projection=proj, llcrnrlat=lat.min(), 
                urcrnrlat=lat.max(), llcrnrlon=lon.min(),
                urcrnrlon=lon.max(), lat_1=lat.min(),
                lon_0=0, resolution='i')
    return m


def mask_fields(var):
    return np.ma.masked_where(var > 1.0e+10, var)


class CheckROMSInputFile(object):
    """
    Loads a ROMS implementation in a particular date and time
    Input arguments:

    implementation : implementation name (cabofrio, phd, etc)
    filename       : input ROMS file to be plotted (fullpath, please)
    tstart         : python datetime object of the first time stamp of the run
    """
    def __init__(self, implementation, filename, tstart):
        self.filename = filename
        self.implementation = implementation
        self.startstr = "%s%02d%02d" %(tstart.year, tstart.month, tstart.day)
        self.tstart = tstart
        self.taste_file()
        self.load_ROMS_grid()
        self.load_file()

    def taste_file(self):
        if "clm" in self.filename:
            self.datatype = "clm"
            self.rangeDict = OCEAN
        elif "ini" in self.filename:
            self.datatype = "clim"
            self.rangeDict = OCEAN
        elif "frc" in self.filename:
            self.datatype = "meteo"
            self.rangeDict = METEO
        elif "tide" in self.filename:
            self.datatype = "tide"
            self.rangeDict = TIDE
        elif "bry" in self.filename:
            self.datatype = "bry"
            self.rangeDict = OCEAN
        else:
            print "ERROR: data file not recognized"


    def load_ROMS_grid(self):
        gridpath = os.path.join(BASEDIR, self.implementation, "whatever/is/your/standard/path")
        gridfile = glob.glob(gridpath + "/*.nc")[0]
        self.grid = cdms2.open(gridfile)
        self.xg = self.grid['lon_rho'][:]
        self.yg = self.grid['lat_rho'][:]
        self.m = init_basemap(self.xg, self.yg, 'cyl')
        self.mxg, self.myg = self.m(self.xg, self.yg)

    def load_file(self):
        self.data = cdms2.open(self.filename)
        print "\n    You loaded %s, here is a summary: \n" %self.filename
        varlist = self.data.listvariables()
        print varlist
        print "\n\n\n"
            
    def plot_variable(self, var, mask=False, level=-1):
        """
        var    :  name of the netcdf variable [str]  
        level  :  only if it's a 3D var
        mask   :  True if you want to mask the array 
        """
        self.var = self.data[var]
        self.varname = var
        # getting time dimension
        if self.datatype == "tide":
            tidefile = "tide%s.nc" %(self.startstr)
            tidepath = os.path.join(BASEDIR, self.implementation, "in", tidefile)
            try:
                tide = cdms2.open(tidepath)
            except cdms2.error.CDMSError:
                print "\n DID YOU CHOOSE THE RIGHT FORECAST START DAY? \n"

            cons = tide['cons'][:]
            self.cons = []
            for k in range(cons.shape[0]):
                consk = cons[k]
                string = ""
                for s in consk:
                    string += s
                self.cons.append(string)
        else:
            time = self.var.getAxis(0)
            self.time = self.time_handler(time)

        if mask:
            self.var[:] = mask_fields(self.var[:])
            print "\nWARNING: variable was masked\n"


        shape = self.var.shape
        if len(shape) > 3:
            self.levels = self.var.getAxis(1)[:]
            print "\nThis is a 3D var, here are the levels:\n"
            print self.levels 
            print "\nPlotting surface as default, but pick a level if you wish:\n"
            self.plot_3Dvars(mask, level)
        elif len(shape) > 2:
            self.plot_2Dvars(mask) 
    
    def plot_2Dvars(self, mask):
        vmin, vmax = self.get_min_max(mask)
        if np.abs(vmin) == np.abs(vmax):
            cmap = plt.cm.RdBu
        else:
            cmap = plt.cm.jet

        nfields = self.var.shape[0]
        if nfields > 20:
            skip = 5
            if nfields > 100:
                skip = 10
        else:
            skip = 1

        for t in range(0, nfields, skip):
            plt.figure(facecolor='w')
            if self.datatype == "meteo":
                lon, lat = self.var.getAxis(2), self.var.getAxis(1)
                lon, lat = np.meshgrid(lon, lat)
                mlon, mlat = self.m(lon, lat)
                self.m.pcolormesh(mlon, mlat, self.var[t,...], 
                                  vmin=vmin, vmax=vmax, cmap=cmap)
            else:
                self.m.pcolormesh(self.mxg, self.myg, self.var[t,...], 
                                  vmin=vmin, vmax=vmax, cmap=cmap)
            self.m.fillcontinents()
            if self.datatype == "tide":
                plt.title( "%s: %s --> %s" %(self.source, self.varname, self.cons[t]) )
            else:
                print self.time[t]
                plt.title( "%s: %s --> %s" %(self.source, self.varname, self.time[t]) )
            plt.colorbar()
            plt.show()

    def plot_3Dvars(self, mask, level):
        if level == -1:
            lev = "Surface"
        elif level == 0:
            lev = "Bottom"
        vmin, vmax = self.get_min_max(mask)
        if np.abs(vmin) == np.abs(vmax):
            cmap = plt.cm.RdBu
        else:
            cmap = plt.cm.jet

        nfields = self.var.shape[0]
        if nfields > 20:
            skip = 5
            if nfields > 100:
                skip = 10
        else:
            skip = 1

        print skip

        for t in range(0, nfields, skip):
            plt.figure(facecolor='w')
            self.m.pcolormesh(self.mxg, self.myg, self.var[t, level, ...], 
                              vmin=vmin, vmax=vmax, cmap=cmap)
            self.m.fillcontinents()

            print self.time[t], lev
            plt.title( "%s: %s --> %s, %s" %(self.source, self.varname, self.time[t], lev))
            plt.colorbar()
            plt.show()

    def time_handler(self, time):
        timelist = []
        for k in range(time[:].size):
            timestamp = dt2ncep(self.tstart) + time[k]
            timelist.append( ncep2dt(timestamp) )
        return timelist

    def get_min_max(self, mask):
        if mask:
            vmin = None
            vmax = None
        else:
            vmin = self.rangeDict[self.varname][0] 
            vmax = self.rangeDict[self.varname][1] 
        return vmin, vmax



#############################################################################


clm = CheckROMSInputFile('nsea', '/home/rsoutelino/myroms/cf_run/cf_frc.nc', 
                           dt.datetime(2014, 6, 17))

clm.plot_variable('sustr', level=-1)




