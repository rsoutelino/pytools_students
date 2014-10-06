#!/usr/bin/env python
######################################################################
#
#  Build a ROMS grid file
#
#  Further Information:  
#  http://www.brest.ird.fr/Roms_tools/
#  
#  This file is part of ROMSTOOLS
#
#  ROMSTOOLS is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published
#  by the Free Software Foundation; either version 2 of the License,
#  or (at your option) any later version.
#
#  ROMSTOOLS is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#  MA  02111-1307  USA
#
#  Copyright (c) 2002-2006 by Pierrick Penven 
#  e-mail:Pierrick.Penven@ird.fr  
#
#  Contributions of P. Marchesiello (IRD) and X. Capet (UCLA)
#
#  Updated    Aug-2006 by Pierrick Penven
#  Updated    24-Oct-2006 by Pierrick Penven (mask correction)
#
#  Translated to Python by Rafael Soutelino: rsoutelino@gmail.com 
#  Last Modification: Aug, 2010
################################################################
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import delaunay
from matplotlib.mlab import griddata, find
from mpl_toolkits.basemap import Basemap
from scipy.io import loadmat
import datetime as dt
import netCDF4 as nc

# classes and functions to the computings
from roms.romslab import RunSetup, rho2uvp, get_metrics, spheric_dist
from roms.romslab import get_angle, add_topo, process_mask, uvp_mask, smoothgrid
from roms.romslab import rotfilter, rfact, hanning_smoother
from roms.romslab import hanning_smoother_coef2d, FX, FY
import roms.romslab as rl

from roms.core import get_zlev, stretching

parser = argparse.ArgumentParser(description='Runs a ROMS simulation.')
parser.add_argument('imp', metavar='imp', type=str,
                   help='ROMS implementation name: Ex: [nsea, makas]')
parser.add_argument('imp_type', metavar='type', type=str,
                   help='Simulation type: [forecast, hindcast]')
args = parser.parse_args()

# sponge/nudging layer defaults
SPONGE_SIZE   = 6
SPONGE_FACTOR = 4 
INNER_NUDGING = 1/60.                        
OUTER_NUDGING = 1/5.                        

# SCRIPT START ######################################################

# Basic Settings:
filetypestr = 'ROMS Grid file - %s %s' %(args.imp, args.imp_type)
rootdir = "/home/rsoutelino/metocean/github/roms/static"

# READING PREVIOUSLY BUILT RELEVANT FILES: ###########################
# metadata ascii file
# OA-created netcdf initial T, S file 
# grid netcdf file
print ' \n' + '==> ' + '  READING ASCII METADATA FILE ...\n' + ' '
roms = RunSetup( os.path.join(rootdir, "domains_%s.yaml" %args.imp_type), args.imp )

dl = roms.res
lonr  = np.arange(roms.lonmin, roms.lonmax + dl, dl)

i = 0; 
latr = np.array([roms.latmin])
while latr[i] <= roms.latmax:
    i    = i + 1
    tmp  = latr[i-1] + dl * np.cos( latr[i-1]*np.pi/180 )
    latr = np.hstack([latr, tmp])


Lonr, Latr        = np.meshgrid(lonr, latr)
Lonu, Lonv, Lonp  = rho2uvp(Lonr)
Latu, Latv, Latp  = rho2uvp(Latr)

M, L = Latp.shape

print ' \n' + '==> ' + '  COMPUTING METRICS  ...\n' + ' '


print ' \n' + '==> ' + '  LLm = ' + np.str(L-1) + ' ...\n' + ' '
print ' \n' + '==> ' + '  MMm = ' + np.str(M-1) + ' ...\n' + ' '

# !!!!!!!!!!!!!!!!!!!!!
### CODE SOMETHING HERE TO WRITE THIS INFORMATION IN THE METADATA FILE
# !!!!!!!!!!!!!!!!!!!!!

pm, pn, dndx, dmde = get_metrics(Latu, Lonu, Latv, Lonv)
xr = 0*pm
yr = xr

for i in np.arange(0, L):
    xr[:, i+1] = xr[:, i] + 2 / ( pm[:, i+1] + pm[:, i] )

for j in np.arange(0, M):
    yr[j+1, :] = yr[j, :] + 2 / ( pn[j+1, :] + pn[j, :] )

xu, xv, xp = rho2uvp(xr)
yu, yv, yp = rho2uvp(yr)

dx    = 1 / pm
dy    = 1 / pn
dxmax = np.max( dx/1000 )
dxmin = np.min( dx/1000 )
dymax = np.max( dy/1000 )
dymin = np.min( dy/1000 )

angle = get_angle(Latu, Lonu)
angle = np.rad2deg(angle)

f0 = 4 * np.pi * np.sin( np.pi * Latr/180 ) / ( 24*3600 )


print ' \n' + '==> ' + '  ADDING TOPOGRAPHY ...\n' + ' '

h = add_topo(Lonr, Latr, pm, pn, roms.bathy)
hraw = h.copy()

print ' \n' + '==> ' + '  COMPUTING THE MASK ...\n' + ' '

maskr = h*0
maskr[ np.where(h > 0) ] = 1 
maskr = process_mask(maskr)
[masku, maskv, maskp] = uvp_mask(maskr) 

print ' \n' + '==> ' + '  FILTERING THE TOPOGRAPHY ...\n' + ' '
h[h > roms.hmax] = roms.hmax 
h = smoothgrid(h, maskr, roms.hmin, roms.hmaxc,
             roms.slope, roms.npass, roms.nfinal)

stop
# computing rx1 factor based on N, Ts, Tb, Hc choices: should be 3 < rx1 < 7
sc = ( np.arange(1, roms.N + 1) - roms.N - 0.5 ) / roms.N
sigma = stretching(sc, roms.Vstret, roms.theta_s, roms.theta_b)
zlev = get_zlev(h, sigma, roms.hc, sc, Vtransform=roms.Vtrans)
rx1 = rl.rx1(zlev, maskr)

if rx1.max() < 3 or rx1.max() > 7:
    print "\n\n\n WARNING: rx1 = %s values violates thresholds !" %rx1.max()
    m = Basemap(resolution='i', llcrnrlon=Lonr.min(), urcrnrlon=Lonr.max(),
                                llcrnrlat=Latr.min(), urcrnrlat=Latr.max())
    m.pcolormesh(Lonr, Latr, rx1, vmin=3, vmax=7)
    m.drawcoastlines()
    plt.colorbar()
    plt.show()


# making sponge layer configuration
fac = np.linspace(1, SPONGE_FACTOR, SPONGE_SIZE)
ss = SPONGE_SIZE
visc = 0*h + 1

visc[:,:ss] = fac[::-1].reshape(1, ss).repeat(visc.shape[0], axis=0) # east
visc[:,-ss:] = fac.reshape(1, ss).repeat(visc.shape[0], axis=0) # west
visc[:ss, :] = fac[::-1].reshape(ss, 1).repeat(visc.shape[1], axis=1) # south
visc[-ss:, :] = fac.reshape(ss, 1).repeat(visc.shape[1], axis=1) # north

diff = visc 

# making nudging layers
nud2d = 0*h
# Initializing coefficients with zeros.  Recall that division by zero is not
# defined and we will give "Inf".  Therefore, we just need to set
# only the values in the desired areas and the rest can be zero
# (for no nud2dging) because the nud2dging in ROMS is:
#
#      F(...,new) = F(...,new) +
#                   dt * F_nud2dgcoef * (Fclm - F(...,new))
fac = np.linspace(INNER_NUDGING, OUTER_NUDGING, SPONGE_SIZE)
nud2d[:,:ss] = fac[::-1].reshape(1, ss).repeat(nud2d.shape[0], axis=0) # east
nud2d[:,-ss:] = fac.reshape(1, ss).repeat(nud2d.shape[0], axis=0) # west
nud2d[:ss, :] = fac[::-1].reshape(ss, 1).repeat(nud2d.shape[1], axis=1) # south
nud2d[-ss:, :] = fac.reshape(ss, 1).repeat(nud2d.shape[1], axis=1) # north

nud3d = nud2d.reshape(1, nud2d.shape[0], nud2d.shape[1]).repeat(roms.N, axis=0)


####################################################################
####################################################################

print ' \n' + '==> ' + '  WRITING NETCDF GRID FILE ...\n' + ' '

today = dt.date.today()
Lp = L + 1
Mp = M + 1

spherical = 'T'

if not os.path.exists(os.path.dirname(roms.grdfile)):
    os.makedirs(os.path.dirname(roms.grdfile))

ncfile = nc.Dataset(roms.grdfile, mode='w',
    clobber='true', format='NETCDF3_CLASSIC')

# creating DIMENSIONS
ncfile.createDimension('xi_u', size=L)
ncfile.createDimension('eta_u', size=Mp)
ncfile.createDimension('xi_v', size=Lp)
ncfile.createDimension('eta_v', size=M)
ncfile.createDimension('xi_rho', size=Lp)
ncfile.createDimension('eta_rho', size=Mp)
ncfile.createDimension('xi_psi', size=L)
ncfile.createDimension('eta_psi', size=M)
ncfile.createDimension('one', size=1)
ncfile.createDimension('two', size=2)
ncfile.createDimension('four', size=4)
ncfile.createDimension('bath', size=1)


# creating GLOBAL ATTRIBUTES
setattr(ncfile, 'title', roms.name)
setattr(ncfile, 'date', str(today))
setattr(ncfile, 'type', filetypestr)

# creating VARIABLES, ATTRIBUTES and ASSIGNING VALUES

# ---------------------------------------------------------------------------
ncfile.createVariable('xl', 'd', dimensions=('one'))
setattr(ncfile.variables['xl'], 'long_name', 'domain length in XI-direction')
setattr(ncfile.variables['xl'], 'units', 'meter')
ncfile.variables['xl'][:]  = xr.max()

# ---------------------------------------------------------------------------
ncfile.createVariable('el', 'd', dimensions=('one'))
setattr(ncfile.variables['el'], 'long_name', 'domain length in ETA-direction')
setattr(ncfile.variables['el'], 'units', 'meter')
ncfile.variables['el'][:]  = yr.max()

# ---------------------------------------------------------------------------
ncfile.createVariable('depthmin', 'd', dimensions=('one'))
setattr(ncfile.variables['depthmin'], 'long_name', 'Shallow bathymetry clipping depth')
setattr(ncfile.variables['depthmin'], 'units', 'meter')
ncfile.variables['depthmin'][:]  = h.min()

# ---------------------------------------------------------------------------
ncfile.createVariable('depthmax', 'd', dimensions=('one'))
setattr(ncfile.variables['depthmax'], 'long_name', 'Deep bathymetry clipping depth')
setattr(ncfile.variables['depthmax'], 'units', 'meter')
ncfile.variables['depthmax'][:]  = h.max()

# ---------------------------------------------------------------------------
ncfile.createVariable('spherical', 'c')
setattr(ncfile.variables['spherical'], 'long_name', 'Grid type logical switch')
setattr(ncfile.variables['spherical'], 'option_T', 'spherical')
setattr(ncfile.variables['spherical'], 'option_F', 'cartesian')
ncfile.variables['spherical'][:]  = spherical

# ---------------------------------------------------------------------------
ncfile.createVariable('angle', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['angle'], 'long_name', 'angle between XI-axis and EAST')
setattr(ncfile.variables['angle'], 'units', 'degree')
# setattr(ncfile.variables['angle'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['angle'][:]  = angle

# ---------------------------------------------------------------------------
ncfile.createVariable('h', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['h'], 'long_name', 'Final bathymetry at RHO-points')
setattr(ncfile.variables['h'], 'units', 'meter')
# setattr(ncfile.variables['h'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['h'][:]  = h

# ---------------------------------------------------------------------------
ncfile.createVariable('hraw', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['hraw'], 'long_name', 'Working bathymetry at RHO-points')
setattr(ncfile.variables['hraw'], 'units', 'meter')
# setattr(ncfile.variables['hraw'], 'coordinates', 'lon_rho lat_rho bath')
ncfile.variables['hraw'][:]  = hraw

# ---------------------------------------------------------------------------
ncfile.createVariable('f', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['f'], 'long_name', 'Coriolis parameter at RHO-points')
setattr(ncfile.variables['f'], 'units', 'second-1')
# setattr(ncfile.variables['f'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['f'][:]  = f0

# ---------------------------------------------------------------------------
ncfile.createVariable('pm', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['pm'], 'long_name', 'Curvilinear coordinate metric in XI')
setattr(ncfile.variables['pm'], 'units', 'meter-1')
# setattr(ncfile.variables['pm'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['pm'][:]  = pm

# ---------------------------------------------------------------------------
ncfile.createVariable('pn', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['pn'], 'long_name', 'Curvilinear coordinate metric in ETA')
setattr(ncfile.variables['pn'], 'units', 'meter-1')
# setattr(ncfile.variables['pn'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['pn'][:]  = pn

# ---------------------------------------------------------------------------
ncfile.createVariable('dndx', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['dndx'], 'long_name', 
    'XI derivative of inverse metric factor pn')
setattr(ncfile.variables['dndx'], 'units', 'meter')
# setattr(ncfile.variables['dndx'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['dndx'][:]  = dndx

# ---------------------------------------------------------------------------
ncfile.createVariable('dmde', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['dmde'], 'long_name', 
    'ETA derivative of inverse metric factor pm')
setattr(ncfile.variables['dmde'], 'units', 'meter')
# setattr(ncfile.variables['dmde'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['dmde'][:]  = dmde

# ---------------------------------------------------------------------------
ncfile.createVariable('x_rho', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['x_rho'], 'long_name', 'x location of RHO-points')
setattr(ncfile.variables['x_rho'], 'units', 'meter')
ncfile.variables['x_rho'][:]  = xr

# ---------------------------------------------------------------------------
ncfile.createVariable('y_rho', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['y_rho'], 'long_name', 'y location of RHO-points')
setattr(ncfile.variables['y_rho'], 'units', 'meter')
ncfile.variables['y_rho'][:]  = yr

# ---------------------------------------------------------------------------
ncfile.createVariable('x_psi', 'd', dimensions=('eta_psi', 'xi_psi'))
setattr(ncfile.variables['x_psi'], 'long_name', 'x location of PSI-points')
setattr(ncfile.variables['x_psi'], 'units', 'meter')
ncfile.variables['x_psi'][:]  = xp

# ---------------------------------------------------------------------------
ncfile.createVariable('y_psi', 'd', dimensions=('eta_psi', 'xi_psi'))
setattr(ncfile.variables['y_psi'], 'long_name', 'y location of PSI-points')
setattr(ncfile.variables['y_psi'], 'units', 'meter')
ncfile.variables['y_psi'][:]  = yp

# ---------------------------------------------------------------------------
ncfile.createVariable('x_u', 'd', dimensions=('eta_u', 'xi_u'))
setattr(ncfile.variables['x_u'], 'long_name', 'x location of U-points')
setattr(ncfile.variables['x_u'], 'units', 'meter')
ncfile.variables['x_u'][:]  = xu

# ---------------------------------------------------------------------------
ncfile.createVariable('y_u', 'd', dimensions=('eta_u', 'xi_u'))
setattr(ncfile.variables['y_u'], 'long_name', 'y location of U-points')
setattr(ncfile.variables['y_u'], 'units', 'meter')
ncfile.variables['y_u'][:]  = yu

# ---------------------------------------------------------------------------
ncfile.createVariable('x_v', 'd', dimensions=('eta_v', 'xi_v'))
setattr(ncfile.variables['x_v'], 'long_name', 'x location of V-points')
setattr(ncfile.variables['x_v'], 'units', 'meter')
ncfile.variables['x_v'][:]  = xv

# ---------------------------------------------------------------------------
ncfile.createVariable('y_v', 'd', dimensions=('eta_v', 'xi_v'))
setattr(ncfile.variables['y_v'], 'long_name', 'y location of V-points')
setattr(ncfile.variables['y_v'], 'units', 'meter')
ncfile.variables['y_v'][:]  = yv

# ---------------------------------------------------------------------------
ncfile.createVariable('lon_rho', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['lon_rho'], 'long_name', 'longitude of RHO-points')
setattr(ncfile.variables['lon_rho'], 'units', 'degree east')
ncfile.variables['lon_rho'][:]  = Lonr

# ---------------------------------------------------------------------------
ncfile.createVariable('lat_rho', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['lat_rho'], 'long_name', 'latitude of RHO-points')
setattr(ncfile.variables['lat_rho'], 'units', 'degree north')
ncfile.variables['lat_rho'][:]  = Latr

# ---------------------------------------------------------------------------
ncfile.createVariable('lon_psi', 'd', dimensions=('eta_psi', 'xi_psi'))
setattr(ncfile.variables['lon_psi'], 'long_name', 'longitude of PSI-points')
setattr(ncfile.variables['lon_psi'], 'units', 'degree east')
ncfile.variables['lon_psi'][:]  = Lonp

# ---------------------------------------------------------------------------
ncfile.createVariable('lat_psi', 'd', dimensions=('eta_psi', 'xi_psi'))
setattr(ncfile.variables['lat_psi'], 'long_name', 'latitude of PSI-points')
setattr(ncfile.variables['lat_psi'], 'units', 'degree north')
ncfile.variables['lat_psi'][:]  = Latp

# ---------------------------------------------------------------------------
ncfile.createVariable('lon_u', 'd', dimensions=('eta_u', 'xi_u'))
setattr(ncfile.variables['lon_u'], 'long_name', 'longitude of U-points')
setattr(ncfile.variables['lon_u'], 'units', 'degree east')
ncfile.variables['lon_u'][:]  = Lonu

# ---------------------------------------------------------------------------
ncfile.createVariable('lat_u', 'd', dimensions=('eta_u', 'xi_u'))
setattr(ncfile.variables['lat_u'], 'long_name', 'latitude of U-points')
setattr(ncfile.variables['lat_u'], 'units', 'degree north')
ncfile.variables['lat_u'][:]  = Latu

# ---------------------------------------------------------------------------
ncfile.createVariable('lon_v', 'd', dimensions=('eta_v', 'xi_v'))
setattr(ncfile.variables['lon_v'], 'long_name', 'longitude of V-points')
setattr(ncfile.variables['lon_v'], 'units', 'degree east')
ncfile.variables['lon_v'][:]  = Lonv

# ---------------------------------------------------------------------------
ncfile.createVariable('lat_v', 'd', dimensions=('eta_v', 'xi_v'))
setattr(ncfile.variables['lat_v'], 'long_name', 'latitude of V-points')
setattr(ncfile.variables['lat_v'], 'units', 'degree north')
ncfile.variables['lat_v'][:]  = Latv

# ---------------------------------------------------------------------------
ncfile.createVariable('mask_rho', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['mask_rho'], 'long_name', 'mask on RHO-points')
setattr(ncfile.variables['mask_rho'], 'flag_values', '0, 1')
setattr(ncfile.variables['mask_rho'], 'flag_meanings', 'land, water')
# setattr(ncfile.variables['mask_rho'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['mask_rho'][:]  = maskr

# ---------------------------------------------------------------------------
ncfile.createVariable('mask_u', 'd', dimensions=('eta_u', 'xi_u'))
setattr(ncfile.variables['mask_u'], 'long_name', 'mask on U-points')
setattr(ncfile.variables['mask_u'], 'flag_values', '0, 1')
setattr(ncfile.variables['mask_u'], 'flag_meanings', 'land, water')
# setattr(ncfile.variables['mask_u'], 'coordinates', 'lon_u lat_u')
ncfile.variables['mask_u'][:]  = masku

# ---------------------------------------------------------------------------
ncfile.createVariable('mask_v', 'd', dimensions=('eta_v', 'xi_v'))
setattr(ncfile.variables['mask_v'], 'long_name', 'mask on V-points')
setattr(ncfile.variables['mask_v'], 'flag_values', '0, 1')
setattr(ncfile.variables['mask_v'], 'flag_meanings', 'land, water')
# setattr(ncfile.variables['mask_v'], 'coordinates', 'lon_v lat_v')
ncfile.variables['mask_v'][:]  = maskv

# ---------------------------------------------------------------------------
ncfile.createVariable('mask_psi', 'd', dimensions=('eta_psi', 'xi_psi'))
setattr(ncfile.variables['mask_psi'], 'long_name', 'mask on PSI-points')
setattr(ncfile.variables['mask_psi'], 'flag_values', '0, 1')
setattr(ncfile.variables['mask_psi'], 'flag_meanings', 'land, water')
# setattr(ncfile.variables['mask_psi'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['mask_psi'][:]  = maskp

# ---------------------------------------------------------------------------
ncfile.createVariable('visc_factor', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['visc_factor'], 'long_name', 'horizontal viscosity sponge factor')
setattr(ncfile.variables['visc_factor'], 'valid_min', 0.)
setattr(ncfile.variables['visc_factor'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['visc_factor'][:]  = visc

# ---------------------------------------------------------------------------
ncfile.createVariable('diff_factor', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['diff_factor'], 'long_name', 'horizontal diffusivity sponge factor')
setattr(ncfile.variables['diff_factor'], 'valid_min', 0.)
setattr(ncfile.variables['diff_factor'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['diff_factor'][:]  = diff

ncfile.sync()

print ' \n' + '==> ' + '  ############################################  ...\n' + ' '
print ' \n' + '==> ' + '        GRID FILE SUCCESSFULLY CREATED          ...\n' + ' '
print ' \n' + '==> ' + '  ############################################  ...\n' + ' '


print ' \n' + '==> ' + '  WRITING NUDGING COEFICIENT FILE ...\n' + ' '

roms.nudgfile = roms.grdfile.replace("grd", "nudg")
filetypestr = filetypestr.replace("Grid", "Nudging Coefficients")
today = dt.date.today()
Lp = L + 1
Mp = M + 1

spherical = 'T'

if not os.path.exists(os.path.dirname(roms.nudgfile)):
    os.makedirs(os.path.dirname(roms.nudgfile))

ncfile = nc.Dataset(roms.nudgfile, mode='w',
    clobber='true', format='NETCDF3_CLASSIC')

# creating DIMENSIONS

ncfile.createDimension('xi_rho', size=Lp)
ncfile.createDimension('eta_rho', size=Mp)
ncfile.createDimension('s_rho', size=roms.N)



# creating GLOBAL ATTRIBUTES
setattr(ncfile, 'title', roms.name)
setattr(ncfile, 'date', str(today))
setattr(ncfile, 'type', filetypestr)


# creating VARIABLES, ATTRIBUTES and ASSIGNING VALUES

# ---------------------------------------------------------------------------
ncfile.createVariable('M2_NudgeCoef', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['M2_NudgeCoef'], 'long_name', '2D momentum inverse nudging coefficients')
setattr(ncfile.variables['M2_NudgeCoef'], 'units', 'day-1')
setattr(ncfile.variables['M2_NudgeCoef'], 'coordinates', 'xi_rho eta_rho')
ncfile.variables['M2_NudgeCoef'][:]  = nud2d

# ---------------------------------------------------------------------------
ncfile.createVariable('M3_NudgeCoef', 'd', dimensions=('s_rho', 'eta_rho', 'xi_rho'))
setattr(ncfile.variables['M3_NudgeCoef'], 'long_name', '3D momentum inverse nudging coefficients')
setattr(ncfile.variables['M3_NudgeCoef'], 'units', 'day-1')
setattr(ncfile.variables['M3_NudgeCoef'], 'coordinates', 'xi_rho eta_rho s_rho')
ncfile.variables['M3_NudgeCoef'][:]  = nud3d

# ---------------------------------------------------------------------------
ncfile.createVariable('tracer_NudgeCoef', 'd', dimensions=('s_rho', 'eta_rho', 'xi_rho'))
setattr(ncfile.variables['tracer_NudgeCoef'], 'long_name', 'generic tracer inverse nudging coefficients')
setattr(ncfile.variables['tracer_NudgeCoef'], 'units', 'day-1')
setattr(ncfile.variables['tracer_NudgeCoef'], 'coordinates', 'xi_rho eta_rho s_rho')
ncfile.variables['tracer_NudgeCoef'][:]  = nud3d

# ---------------------------------------------------------------------------
ncfile.createVariable('temp_NudgeCoef', 'd', dimensions=('s_rho', 'eta_rho', 'xi_rho'))
setattr(ncfile.variables['temp_NudgeCoef'], 'long_name', 'temp inverse nudging coefficients')
setattr(ncfile.variables['temp_NudgeCoef'], 'units', 'day-1')
setattr(ncfile.variables['temp_NudgeCoef'], 'coordinates', 'xi_rho eta_rho s_rho')
ncfile.variables['temp_NudgeCoef'][:]  = nud3d

# ---------------------------------------------------------------------------
ncfile.createVariable('salt_NudgeCoef', 'd', dimensions=('s_rho', 'eta_rho', 'xi_rho'))
setattr(ncfile.variables['salt_NudgeCoef'], 'long_name', 'salt inverse nudging coefficients')
setattr(ncfile.variables['salt_NudgeCoef'], 'units', 'day-1')
setattr(ncfile.variables['salt_NudgeCoef'], 'coordinates', 'xi_rho eta_rho s_rho')
ncfile.variables['salt_NudgeCoef'][:]  = nud3d


ncfile.sync()

print ' \n' + '==> ' + '  ############################################  ...\n' + ' '
print ' \n' + '==> ' + '        NUDGING FILE SUCCESSFULLY CREATED          ...\n' + ' '
print ' \n' + '==> ' + '  ############################################  ...\n' + ' '












































































