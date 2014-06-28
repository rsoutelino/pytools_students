#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
import glob
import cdms2
import mpl_util
from matplotlib.path import Path

# NW Australia
xmin, xmax, ymin, ymax = (-46, -38, -28, -21)

grid1 = ( [xmin+2, xmax-2, xmax-2, xmin+2, xmin+2], 
	      [ymax-2, ymax-2, ymin+2, ymin+2, ymax-2] )
# grid2 = ( [119, 138, 138, 119, 119], [-5, -5, +21, -21, -5] )

m = Basemap(projection='cyl', llcrnrlat=ymin, urcrnrlat=ymax,
    llcrnrlon=xmin, urcrnrlon=xmax, lat_1=ymin,
    lon_0=xmin, resolution='l')

levels = [0, 500, 1000, 1500, 2000, 2500, 3000, 4000, 4500, 5000, 5500, 6000]

h = cdms2.open("cf_grd.nc")['h'].getValue()
lon = cdms2.open("cf_grd.nc")['lon_rho'].getValue()
lat = cdms2.open("cf_grd.nc")['lat_rho'].getValue()
mlon, mlat = m(lon, lat)


plt.figure(facecolor="w")
m.fillcontinents()
m.contourf(mlon, mlat, h, levels,
             cmap=mpl_util.LevelColormap(levels, cmap=plt.cm.Blues), 
             extend='both')
m.plot(grid1[0], grid1[1], 'k', linewidth=2)
# m.plot(grid2[0], grid2[1], 'r', linewidth=2)
plt.colorbar()
# m.contour(mlon, mlat, h, [1000, 2000], colors='0.5')

plt.show()
