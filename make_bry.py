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

bry_filenamestr = '_bry.nc'
bry_filetypestr = 'ROMS Boundary Forcing file'


# READING PREVIOUSLY BUILT RELEVANT FILES: ###########################
# metadata ascii file
# grid netcdf file

print ' \n' + '==> ' + '  READING ASCII METADATA FILE ...\n' + ' '
Run = run_setup('piloto.setup')

print ' \n' + '==> ' + '  READING GRID NETCDF FILE ...\n' + ' '
grdfile  = netCDF4.Dataset(Run.datadir + Run.run_name + '_grd.nc')


# WRITING THE NETCDF FILE ####################################################
# Based on "bry_unlimit.cdl" NETCDF sample structure

# some computings regarding boundary limits:

print ' \n' + '==> ' + '  WRITING NETCDF Bry FILE ...\n' + ' '

ncfile = netCDF4.Dataset(run.datadir + run.run_name + bry_filenamestr, mode='w',
    clobber='true', format='NETCDF3_CLASSIC')

bry_time = np.arange(1, 31)

# creating DIMENSIONS
ncfile.createDimension('xi_rho', size=Lp)
ncfile.createDimension('xi_u', size=L)
ncfile.createDimension('xi_v', size=Lp)
ncfile.createDimension('eta_rho', size=Mp)
ncfile.createDimension('eta_u', size=Mp)
ncfile.createDimension('eta_v', size=M)
ncfile.createDimension('s_rho', size=N)
ncfile.createDimension('s_w', size=Np)
ncfile.createDimension('one', size=1)
ncfile.createDimension('bry_time', size=bry_time.size)
#ncfile.createDimension('zeta_time', size= )
#ncfile.createDimension('v2d_time',  size= )
#ncfile.createDimension('v3d_time',  size= )
#ncfile.createDimension('temp_time', size= )
#ncfile.createDimension('salt_time', size= )
#ncfile.createDimension('tracer', size=2)
#ncfile.createDimension('time', size=1)


# creating GLOBAL ATTRIBUTES
setattr(ncfile, 'type', bry_filetypestr)
setattr(ncfile, 'title', run.bry_info)
setattr(ncfile, 'out_file', run.run_name + bry_filenamestr)
setattr(ncfile, 'grd_file', run.run_name + '_grd.nc')
now = dt.datetime.now()
setattr(ncfile,'history',np.str(now))


# creating VARIABLES, ATTRIBUTES and ASSIGNING VALUES

# ---------------------------------------------------------------------------
ncfile.createVariable('spherical', 'c')
setattr(ncfile.variables['spherical'], 'long_name', 'grid type logical switch')
setattr(ncfile.variables['spherical'], 'flag_values', '0, 1')
setattr(ncfile.variables['spherical'], 'flag_meanings', 'cartesian spherical')
ncfile.variables['spherical'][:]  = spherical

# ---------------------------------------------------------------------------
ncfile.createVariable('Vtransform', 'd', dimensions=('one'))
setattr(ncfile.variables['Vtransform'], 'long_name',
    'vertical terrain-following transformation equation')
ncfile.variables['Vtransform'][:]  = run.vtransform

# ---------------------------------------------------------------------------
ncfile.createVariable('Vstretching', 'd', dimensions=('one'))
setattr(ncfile.variables['Vstretching'], 'long_name',
    'vertical terrain-following stretching function')
ncfile.variables['Vstretching'][:]  = run.vstretching

# ---------------------------------------------------------------------------
ncfile.createVariable('theta_s', 'd', dimensions=('one'))
setattr(ncfile.variables['theta_s'], 'long_name',
    'S-coordinate surface control parameter')
ncfile.variables['theta_s'][:] = run.theta_s

# ---------------------------------------------------------------------------
ncfile.createVariable('theta_b', 'd', dimensions=('one'))
setattr(ncfile.variables['theta_b'], 'long_name',
    'S-coordinate bottom control parameter')
ncfile.variables['theta_b'][:] = run.theta_b

# ---------------------------------------------------------------------------
ncfile.createVariable('Tcline', 'd', dimensions=('one'))
setattr(ncfile.variables['Tcline'], 'long_name',
    'S-coordinate surface/bottom layer width')
setattr(ncfile.variables['Tcline'], 'units', 'meter')
ncfile.variables['Tcline'][:] = run.tcline

# ---------------------------------------------------------------------------
ncfile.createVariable('hc', 'd', dimensions=('one'))
setattr(ncfile.variables['hc'],'long_name',
    'S-coordinate parameter, critical depth')
setattr(ncfile.variables['hc'], 'units', 'meter')
ncfile.variables['hc'][:] = run.hc

# ---------------------------------------------------------------------------
ncfile.createVariable('s_rho', 'd', dimensions=('s_rho'))
setattr(ncfile.variables['s_rho'], 'long_name', 'S-coordinate at RHO-points')
setattr(ncfile.variables['s_rho'], 'valid_min', -1.0)
setattr(ncfile.variables['s_rho'], 'valid_max', 0.0)
setattr(ncfile.variables['s_rho'], 'positive', 'up')
setattr(ncfile.variables['s_rho'], 'standard_name', 'ocean_s_coordinate_g1')
setattr(ncfile.variables['s_rho'], 'formula_terms', 
    's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc')
ncfile.variables['s_rho'][:] = sc

# ---------------------------------------------------------------------------
ncfile.createVariable('s_w', 'd', dimensions=('s_w'))
setattr(ncfile.variables['s_w'], 'long_name', 'S-coordinate at W-points')
setattr(ncfile.variables['s_w'], 'valid_min', -1.0)
setattr(ncfile.variables['s_w'], 'valid_max', 0.0)
setattr(ncfile.variables['s_w'], 'positive', 'up')
setattr(ncfile.variables['s_w'], 'standard_name', 'ocean_s_coordinate_g1')
setattr(ncfile.variables['s_w'], 'formula_terms', 
    's: s_rho C: Cs_w eta: zeta depth: h depth_c: hc')
ncfile.variables['s_w'][:] = scw


# ---------------------------------------------------------------------------
ncfile.createVariable('Cs_r', 'd', dimensions=('s_rho'))
setattr(ncfile.variables['Cs_r'], 'long_name',
    'S-coordinate stretching curves at RHO-points')
setattr(ncfile.variables['Cs_r'], 'valid_min', -1.0)
setattr(ncfile.variables['Cs_r'], 'valid_max', 0.0)
ncfile.variables['Cs_r'][:] = Cs

# ---------------------------------------------------------------------------
ncfile.createVariable('Cs_w', 'd', dimensions=('s_w'))
setattr(ncfile.variables['Cs_w'], 'long_name',
    'S-coordinate stretching curves at W-points')
setattr(ncfile.variables['Cs_w'], 'valid_min', -1.0)
setattr(ncfile.variables['Cs_w'], 'valid_max', 0.0)
ncfile.variables['Cs_w'][:] = Csw

# ---------------------------------------------------------------------------
ncfile.createVariable('h', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['h'], 'long_name', 'bathymetry at RHO-points')
setattr(ncfile.variables['h'], 'units', 'meter')
setattr(ncfile.variables['h'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['h'][:] = h2

# ---------------------------------------------------------------------------
ncfile.createVariable('lon_rho', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['lon_rho'], 'long_name', 'longitude of RHO-points')
setattr(ncfile.variables['lon_rho'], 'units', 'degree_east')
setattr(ncfile.variables['lon_rho'], 'standard_name', 'longitude')
ncfile.variables['lon_rho'][:] = grdfile.variables['lon_rho'][:]

# ---------------------------------------------------------------------------
ncfile.createVariable('lat_rho', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['lat_rho'], 'long_name', 'latitude of RHO-points')
setattr(ncfile.variables['lat_rho'], 'units', 'degree_north')
setattr(ncfile.variables['lat_rho'], 'standard_name', 'latitude')
ncfile.variables['lat_rho'][:] = grdfile.variables['lat_rho'][:]

# ---------------------------------------------------------------------------
ncfile.createVariable('lon_u', 'd', dimensions=('eta_u', 'xi_u'))
setattr(ncfile.variables['lon_u'], 'long_name', 'longitude of U-points')
setattr(ncfile.variables['lon_u'], 'units', 'degree_east')
setattr(ncfile.variables['lon_u'], 'standard_name', 'longitude')
ncfile.variables['lon_u'][:] = grdfile.variables['lon_u'][:]

# ---------------------------------------------------------------------------
ncfile.createVariable('lat_u', 'd', dimensions=('eta_u', 'xi_u'))
setattr(ncfile.variables['lat_u'], 'long_name', 'latitude of U-points')
setattr(ncfile.variables['lat_u'], 'units', 'degree_north')
setattr(ncfile.variables['lat_u'], 'standard_name', 'latitude')
ncfile.variables['lat_u'][:] = grdfile.variables['lat_u'][:]

# ---------------------------------------------------------------------------
ncfile.createVariable('lon_v', 'd', dimensions=('eta_v', 'xi_v'))
setattr(ncfile.variables['lon_v'], 'long_name', 'longitude of V-points')
setattr(ncfile.variables['lon_v'], 'units', 'degree_east')
setattr(ncfile.variables['lon_v'], 'standard_name', 'lonitude')
ncfile.variables['lon_v'][:] = grdfile.variables['lon_v'][:]

# ---------------------------------------------------------------------------
ncfile.createVariable('lat_v', 'd', dimensions=('eta_v', 'xi_v'))
setattr(ncfile.variables['lat_v'], 'long_name', 'latitude of V-points')
setattr(ncfile.variables['lat_v'], 'units', 'degree_north')
setattr(ncfile.variables['lat_v'], 'standard_name', 'latitude')
ncfile.variables['lat_v'][:] = grdfile.variables['lat_v'][:]


# ---------------------------------------------------------------------------
ncfile.createVariable('ocean_time', 'd', dimensions=('bry_time'))
setattr(ncfile.variables['ocean_time'], 'long_name', 'open boundary conditions time')
setattr(ncfile.variables['ocean_time'], 'units', 'time-units since 0001-01-01 00:00:00')
setattr(ncfile.variables['ocean_time'], 'calendar', 'gregorian')
setattr(ncfile.variables['ocean_time'], 'cycle_length', 30.)
ncfile.variables['ocean_time'][:] = bry_time


# ---------------------------------------------------------------------------
ncfile.createVariable('temp_west', 'd', dimensions=('bry_time', 's_rho', 'eta_rho'))
setattr(ncfile.variables['temp_west'], 'long_name', 
    'potential temperature western boundary condition ')
setattr(ncfile.variables['temp_west'], 'units', 'celcius')
setattr(ncfile.variables['temp_west'], 'time', 'ocean_time')
setattr(ncfile.variables['temp_west'], 'coordinates', 'lon_rho lat_rho s_rho ocean_time')
ncfile.variables['temp_west'][:] = TEMP[:,:,0].reshape(1, N, Mp).repeat(bry_time.size, axis=0)

# ---------------------------------------------------------------------------
ncfile.createVariable('temp_east', 'd', dimensions=('bry_time', 's_rho', 'eta_rho'))
setattr(ncfile.variables['temp_east'], 'long_name', 
    'potential temperature eastern boundary condition')
setattr(ncfile.variables['temp_east'], 'units', 'celcius')
setattr(ncfile.variables['temp_east'], 'time', 'ocean_time')
setattr(ncfile.variables['temp_east'], 'coordinates', 'lon_rho lat_rho s_rho ocean_time')
ncfile.variables['temp_east'][:] = TEMP[:,:,-1].reshape(1, N, Mp).repeat(bry_time.size, axis=0)

# ---------------------------------------------------------------------------
ncfile.createVariable('temp_south', 'd', dimensions=('bry_time', 's_rho', 'xi_rho'))
setattr(ncfile.variables['temp_south'], 'long_name', 
    'potential temperature southern boundary condition')
setattr(ncfile.variables['temp_south'], 'units', 'celcius')
setattr(ncfile.variables['temp_south'], 'time', 'ocean_time')
setattr(ncfile.variables['temp_south'], 'coordinates', 'lon_rho lat_rho s_rho ocean_time')
ncfile.variables['temp_south'][:] = TEMP[:,0,:].reshape(1, N, Lp).repeat(bry_time.size, axis=0)

# ---------------------------------------------------------------------------
ncfile.createVariable('salt_west', 'd', dimensions=('bry_time', 's_rho', 'eta_rho'))
setattr(ncfile.variables['salt_west'], 'long_name', 'salinity western boundary condition')
setattr(ncfile.variables['salt_west'], 'time', 'ocean_time')
setattr(ncfile.variables['salt_west'], 'coordinates', 'lon_rho lat_rho s_rho ocean_time')
ncfile.variables['salt_west'][:] = SALT[:,:,0].reshape(1, N, Mp).repeat(bry_time.size, axis=0)

# ---------------------------------------------------------------------------
ncfile.createVariable('salt_east', 'd', dimensions=('bry_time', 's_rho', 'eta_rho'))
setattr(ncfile.variables['salt_east'], 'long_name', 'salinity eastern boundary condition')
setattr(ncfile.variables['salt_east'], 'time', 'ocean_time')
setattr(ncfile.variables['salt_east'], 'coordinates', 'lon_rho lat_rho s_rho ocean_time')
ncfile.variables['salt_east'][:] = SALT[:,:,-1].reshape(1, N, Mp).repeat(bry_time.size, axis=0)

# ---------------------------------------------------------------------------
ncfile.createVariable('salt_south', 'd', dimensions=('bry_time', 's_rho', 'xi_rho'))
setattr(ncfile.variables['salt_south'], 'long_name', 'salinity southern boundary condition')
setattr(ncfile.variables['salt_south'], 'time', 'ocean_time')
setattr(ncfile.variables['salt_south'], 'coordinates', 'lon_rho lat_rho s_rho ocean_time')
ncfile.variables['salt_south'][:] = SALT[:,0,:].reshape(1, N, Lp).repeat(bry_time.size, axis=0)


ncfile.sync()

print ' \n' + '==> ' + '  #############################  ...\n' + ' '
print ' \n' + '==> ' + '  Bry FILE SUCCESSFULLY CREATED  ...\n' + ' '
print ' \n' + '==> ' + '  #############################  ...\n' + ' '

