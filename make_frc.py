#!/usr/bin/env python
# -*- coding: utf-8 -*-


# IMPORTING MODULES #################################################
print ' \n' + '==> ' + '  IMPORTING MODULES ...\n' + ' ' 

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import delaunay
from matplotlib.mlab import griddata
from mpl_toolkits.basemap import Basemap
import datetime as dt
import netCDF4
import scipy.io as sp

# classes and functions to the computings
from roms_setup import run_setup, zlev, ztosigma 


# SCRIPT START ######################################################

# Basic Settings:

frc_filenamestr = '_cfsr_frc.nc'
frc_filetypestr = 'ROMS Forcing file'


# READING PREVIOUSLY BUILT RELEVANT FILES: ###########################
# metadata ascii file
# grid netcdf file

print ' \n' + '==> ' + '  READING ASCII METADATA FILE ...\n' + ' '
Run = run_setup('piloto.setup')

print ' \n' + '==> ' + '  READING GRID NETCDF FILE ...\n' + ' '
grdfile  = netCDF4.Dataset(Run.datadir + Run.run_name + '_grd.nc')

#######################################################################
##                                                                   ##
##					!!!!	BULK FLUXES		!!!!                     ##
##                                                                   ##
#######################################################################

# some computings regarding BULK FLUXES:

# frc_time = np.arange(0, 90)
# lon = grdfile.variables['lon_rho'][:]
# Mp, Lp = lon.shape

# # WIND

# wnd_ms = 5.5             		   ## 10 knot wind = 5.5 ms-1

# u = wnd_ms * np.cos( np.pi/4 )     ## decomposing wnd_ms vector to zonal(u) component

# # wind ramping from 0 to 'u' in the first 2 days, mantaining constant for
# # the rest 88 days --> 90 days wind period
# uwnd = np.concatenate( ( np.linspace(0,u,4), u*np.ones(86)  ) )

# uwnda = np.dot( np.ones((Lp, Mp, 1)), uwnd.reshape(1,frc_time.size) )  #creating wind array
# UWND  = uwnda.transpose()
# VWND  = UWND.copy()

# ATMOSPHERIC PRESSURE
# pair = 1015

# AIR TEMPERATURE
# tair = 25

# AIR HUMIDITY
# qair = 50

#######################################################################
##                                                                   ##
##			 	!!!!	SMS (UVstress) FLUXES		!!!!             ##
##                                                                   ##
#######################################################################


## some computings regarding SMS FLUXES:

# sms_time = np.arange(0, 90)
# lon = grdfile.variables['lon_rho'][:]
# Mp, Lp = lon.shape

# # U, V STRESS

# wnd_str = -0.1             ## 0.1Pa = ~6.0 ms-1

# # wind ramping from 0 to 'u' in the first 2 days, mantaining constant for
# # the rest 86 days --> 90 days wind period
# ustr = np.concatenate( ( np.linspace(0,wnd_str,4), wnd_str*np.ones(86)  ) )
# vstr = np.concatenate( ( np.linspace(0,wnd_str,4), wnd_str*np.ones(86)  ) )



# ustra = np.dot( np.ones((Lp-1, Mp, 1)), ustr.reshape(1,sms_time.size) )  #creating wind array
# USTR  = ustra.transpose()

# vstra = np.dot( np.ones((Lp, Mp-1, 1)), vstr.reshape(1,sms_time.size) )  #creating wind array
# VSTR  = vstra.transpose()


#######################################################################
##                                                                   ##
##		!!!!	CFSR --> SMS (UVstress) FLUXES		!!!!             ##
##                                                                   ##
#######################################################################
gridfile = netCDF4.Dataset('init/piloto_grd.nc')
lon      = grdfile.variables['lon_rho'][:]
Mp, Lp   = lon.shape



cfsr = netCDF4.Dataset('init/wind_jan2010.nc')

u = cfsr.variables['u_wind'][:]
v = cfsr.variables['v_wind'][:]

sms_time = np.arange(0, 90)

## some computings regarding the conversio from wind to surface momentum stress:

Cd   = 0.0022   ## Wind-drag coefficient
Pair = 1.275    ## Air density

ustra = Pair * Cd * u[11,...]**2    # Parametrization of shear stress (sustr) as a function of wind speed (u)
vstra = Pair * Cd * v[11,...]**2    # v[11,...] -> S wind/western domain - NE wind/eastern domain 

ustr = ustra*np.ones([ len(sms_time), ustra.shape[0], ustra.shape[1] ])
vstr = vstra*np.ones([ len(sms_time), ustra.shape[0] ,ustra.shape[1] ])


# WRITING THE NETCDF FILE ####################################################
# Based on "sms_uvstress.cdl" NETCDF sample structure


print ' \n' + '==> ' + '  WRITING NETCDF Frc FILE ...\n' + ' '

ncfile = netCDF4.Dataset(Run.datadir + Run.run_name + frc_filenamestr, mode='w',
    clobber='true', format='NETCDF3_CLASSIC')

# creating DIMENSIONS
ncfile.createDimension('xi_u',  size= ustr.shape[1])
ncfile.createDimension('eta_u', size= ustr.shape[0])
ncfile.createDimension('xi_v',  size= ustr.shape[1])
ncfile.createDimension('eta_v', size= ustr.shape[0])
ncfile.createDimension('sms_time', size= len(sms_time))


# creating GLOBAL ATTRIBUTES
setattr(ncfile, 'type', frc_filetypestr)
setattr(ncfile, 'title', Run.frc_info)
setattr(ncfile, 'out_file', Run.run_name + frc_filenamestr)
setattr(ncfile, 'grd_file', Run.run_name + '_grd.nc')
setattr(ncfile,'history',np.str( dt.datetime.now() ))

# creating VARIABLES, ATTRIBUTES and ASSIGNING VALUES

# ---------------------------------------------------------------------------
ncfile.createVariable('sms_time', 'd', dimensions=('sms_time'))
setattr(ncfile.variables['sms_time'], 'long_name', 'surface momentum stress time')
setattr(ncfile.variables['sms_time'], 'units', 'time-units since 0001-01-01 00:00:00')
setattr(ncfile.variables['sms_time'], 'calendar', 'gregorian')
#setattr(ncfile.variables['sms_time'], 'cycle_length', 30.)
ncfile.variables['sms_time'][:] = sms_time

# ---------------------------------------------------------------------------
ncfile.createVariable('sustr', 'd', dimensions=('sms_time', 'eta_u', 'xi_u'))
setattr(ncfile.variables['sustr'], 'long_name', 'surface u-momentum stress')
setattr(ncfile.variables['sustr'], 'units', 'Newton meter-2')
setattr(ncfile.variables['sustr'], 'time', 'sms_time')
#setattr(ncfile.variables['sustr'], '', '')
ncfile.variables['sustr'][:] = ustr

# ---------------------------------------------------------------------------
ncfile.createVariable('svstr', 'd', dimensions=('sms_time', 'eta_v', 'xi_v'))
setattr(ncfile.variables['svstr'], 'long_name', 'surface v-momentum stress')
setattr(ncfile.variables['svstr'], 'units', 'Newton meter-2')
setattr(ncfile.variables['svstr'], 'time', 'sms_time')
#setattr(ncfile.variables['svstr'], '', '')
ncfile.variables['svstr'][:] = vstr

if

# # ---------------------------------------------------------------------------
# ncfile.createVariable('frc_time', 'd', dimensions=('frc_time'))
# setattr(ncfile.variables['frc_time'], 'long_name', 'atmospheric forcing time')
# setattr(ncfile.variables['frc_time'], 'units', 'time-units since 0001-01-01 00:00:00')
# setattr(ncfile.variables['frc_time'], 'calendar', 'gregorian')
# #setattr(ncfile.variables['frc_time'], 'cycle_length', 30.)
# ncfile.variables['frc_time'][:] = frc_time

# # ---------------------------------------------------------------------------
# ncfile.createVariable('Uwind', 'd', dimensions=('frc_time', 'eta_rho', 'xi_rho'))
# setattr(ncfile.variables['Uwind'], 'long_name', 'surface u-wind component')
# setattr(ncfile.variables['Uwind'], 'units', 'meter second-1')
# setattr(ncfile.variables['Uwind'], 'time', 'frc_time')
# setattr(ncfile.variables['Uwind'], 'coordinates', 'lon_rho lat_rho s_rho ocean_time')
# ncfile.variables['Uwind'][:] = UWND*0

# # ---------------------------------------------------------------------------
# ncfile.createVariable('Vwind', 'd', dimensions=('frc_time', 'eta_rho', 'xi_rho'))
# setattr(ncfile.variables['Vwind'], 'long_name', 'surface v-wind component')
# setattr(ncfile.variables['Vwind'], 'units', 'meter second-1')
# setattr(ncfile.variables['Vwind'], 'time', 'frc_time')
# setattr(ncfile.variables['Vwind'], 'coordinates', 'lon_rho lat_rho s_rho ocean_time')
# ncfile.variables['Vwind'][:] = VWND*0

# # ---------------------------------------------------------------------------
# ncfile.createVariable('Pair', 'd', dimensions=('frc_time', 'eta_rho', 'xi_rho'))
# setattr(ncfile.variables['Pair'], 'long_name', 'surface air pressure')
# setattr(ncfile.variables['Pair'], 'units', 'millibar')
# setattr(ncfile.variables['Pair'], 'time', 'frc_time')
# ncfile.variables['Pair'][:] = UWND*0 + pair

# # ---------------------------------------------------------------------------
# ncfile.createVariable('Tair', 'd', dimensions=('frc_time', 'eta_rho', 'xi_rho'))
# setattr(ncfile.variables['Tair'], 'long_name', 'surface air temperature')
# setattr(ncfile.variables['Tair'], 'units', 'Celsius')
# setattr(ncfile.variables['Tair'], 'time', 'frc_time')
# ncfile.variables['Tair'][:] = UWND*0 + tair

# # ---------------------------------------------------------------------------
# ncfile.createVariable('Qair', 'd', dimensions=('frc_time', 'eta_rho', 'xi_rho'))
# setattr(ncfile.variables['Qair'], 'long_name', 'surface air relative humidity')
# setattr(ncfile.variables['Qair'], 'units', 'percentage')
# setattr(ncfile.variables['Qair'], 'time', 'frc_time')
# ncfile.variables['Qair'][:] = UWND*0

# # ---------------------------------------------------------------------------
# ncfile.createVariable('rain', 'd', dimensions=('frc_time', 'eta_rho', 'xi_rho'))
# setattr(ncfile.variables['rain'], 'long_name', 'rain fall rate')
# setattr(ncfile.variables['rain'], 'units', 'kilogram meter-2 second-1')
# setattr(ncfile.variables['rain'], 'time', 'frc_time')
# ncfile.variables['rain'][:] = UWND*0

# # ---------------------------------------------------------------------------
# ncfile.createVariable('swrad', 'd', dimensions=('frc_time', 'eta_rho', 'xi_rho'))
# setattr(ncfile.variables['swrad'], 'long_name', 'solar shortwave radiation')
# setattr(ncfile.variables['swrad'], 'units', 'Watts meter-2')
# setattr(ncfile.variables['swrad'], 'positive_value', 'downward flux, heating')
# setattr(ncfile.variables['swrad'], 'negative_value', 'upward flux, cooling')
# setattr(ncfile.variables['swrad'], 'time', 'frc_time')
# ncfile.variables['swrad'][:] = UWND*0

# # ---------------------------------------------------------------------------
# ncfile.createVariable('lwrad', 'd', dimensions=('frc_time', 'eta_rho', 'xi_rho'))
# setattr(ncfile.variables['lwrad'], 'long_name', 'net longwave radiation flux')
# setattr(ncfile.variables['lwrad'], 'units', 'Watts meter-2')
# setattr(ncfile.variables['lwrad'], 'positive_value', 'downward flux, heating')
# setattr(ncfile.variables['lwrad'], 'negative_value', 'upward flux, cooling')
# setattr(ncfile.variables['lwrad'], 'time', 'frc_time')
# ncfile.variables['lwrad'][:] = UWND*0

ncfile.sync()

print ' \n' + '==> ' + '  #############################  ...\n' + ' '
print ' \n' + '==> ' + '  Frc FILE SUCCESSFULLY CREATED  ...\n' + ' '
print ' \n' + '==> ' + '  #############################  ...\n' + ' '

