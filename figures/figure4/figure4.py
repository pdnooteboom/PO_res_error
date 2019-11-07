#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 15:16:58 2019

@author: nooteboom
"""

import numpy as np
import matplotlib.pylab as plt
from netCDF4 import Dataset
import cartopy.crs as ccrs
import seaborn as sns
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

sns.set_style("darkgrid")
sns.set_context("paper")

projection = ccrs.PlateCarree()
size = 17
cmap='rainbow'

lat= -47
lon = 42

lat2= -46
lon2 = 39


width = 55
height = 30
width2 = width 
height2 = height 

# Manage the location on the earth and extents of the subfigures
r1 = -19
u1 = 0
exte = [-360+lon-width/2.+r1, -360+lon+width/2.+r1, lat-height/2.+u1, lat+height/2.+u1]
r2 = -13
u2 = 2
exte2 = np.array([-360+lon2-width2/2.+r2, -360+lon2+width2/2.+r2, lat2-height2/2.+u2, lat2+height2/2.+u2])

xloc = np.arange(-180, 180, 20)
yloc = np.arange(-90,90,15)

ddeg = 1
sp = 6
dd = 10
cs = 3.0

fs = 18

csea = 'whitesmoke'
cland = 'k'

ch = 'red'
cl = 'royalblue'
cd = 'y'

opac = 0.5
#%% First the highres
dirReadhigh = '/Volumes/HardDisk/POP/output/highres/timeseries/'
nc_hr = Dataset(dirReadhigh + 'timeseries_per_location_ddeg%d_sp%d_dd%d_tempres5.nc'%(ddeg,sp,dd))

Lons = nc_hr['Lons'][:]
Lats = nc_hr['Lats'][:]
vLons = nc_hr['vLons'][:]
vLats = nc_hr['vLats'][:]

hlo = find_nearest(Lons,lon)
hla = find_nearest(Lats,lat)
hlo2 = find_nearest(Lons,lon2)
hla2 = find_nearest(Lats,lat2)

idx = np.where(np.logical_and(vLons==hlo, vLats==hla))
idx2 = np.where(np.logical_and(vLons==hlo2, vLats==hla2))
assert len(idx[0]==1)

hxs = nc_hr['lon'][idx]
hys = nc_hr['lat'][idx]
hstemp = nc_hr['temp'][idx]

ml = len(hstemp[0])

hxs2 = nc_hr['lon'][idx2]
hys2 = nc_hr['lat'][idx2]
hstemp2 = nc_hr['temp'][idx2]

ml2 = len(hstemp2[0])
#%% Then the highres monthly
dirReadlow = '/Volumes/HardDisk/POP/output/highres/timeseries/'
nc_lr = Dataset(dirReadlow + 'timeseries_per_location_ddeg1_sp6_dd10_tempresmonmean.nc')

Lons = nc_lr['Lons'][:]
Lats = nc_lr['Lats'][:]
vLons = nc_lr['vLons'][:]
vLats = nc_lr['vLats'][:]

llo = find_nearest(Lons,lon)
lla = find_nearest(Lats,lat)
llo2 = find_nearest(Lons,lon2)
lla2 = find_nearest(Lats,lat2)
idx2 = np.where(np.logical_and(vLons==llo2, vLats==lla2))
idx = np.where(np.logical_and(vLons==llo, vLats==lla))

assert len(idx2[0]==1)
lxs2 = nc_lr['lon'][idx2]
lys2 = nc_lr['lat'][idx2]
lstemp2 = nc_lr['temp'][idx2][:,:ml2]
lxs = nc_lr['lon'][idx]
lys = nc_lr['lat'][idx]
lstemp = nc_lr['temp'][idx][:,:ml]

#%% The plot
fig = plt.figure(figsize=(15,7))

ax0 = plt.subplot(1,2,2, projection=projection)
plt.title('(b)', fontsize=fs)
ax0.add_feature(cartopy.feature.LAND, color=cland)
ax0.add_feature(cartopy.feature.OCEAN, color=csea, zorder=0)
g = ax0.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')
g.xlabels_top = False
g.ylabels_left = False
g.xformatter = LONGITUDE_FORMATTER
g.yformatter = LATITUDE_FORMATTER
g.xlabel_style = {'fontsize': fs-2}
g.ylabel_style = {'fontsize': fs-2}
g.xlocator = mticker.FixedLocator(xloc)
g.ylocator = mticker.FixedLocator(yloc)
ax0.set_extent(exte, ccrs.PlateCarree())

plt.scatter(hxs, hys, c=ch, s=size, label='1', alpha=opac)
plt.scatter(lxs, lys, c=cl, s=size, label='2', alpha=opac)
plt.scatter(llo+0.5, lla+0.5 ,c='k', marker='P', s=150)
plt.legend(loc=4, fontsize=fs-4, title='configuration')

ax0 = plt.subplot(1,2,1, projection=projection)
plt.title('(a)', fontsize=fs)
ax0.add_feature(cartopy.feature.LAND, color=cland)
ax0.add_feature(cartopy.feature.OCEAN, color=csea, zorder=0)
g = ax0.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')
g.xlabels_top = False
g.ylabels_right = False
g.xformatter = LONGITUDE_FORMATTER
g.yformatter = LATITUDE_FORMATTER
g.xlocator = mticker.FixedLocator(xloc)
g.ylocator = mticker.FixedLocator(yloc)
g.xlabel_style = {'fontsize': fs}
g.ylabel_style = {'fontsize': fs}
ax0.set_extent(exte2, ccrs.PlateCarree())

plt.scatter(lxs2, lys2, c=cl, s=size, label='1', alpha=opac)
plt.scatter(hxs2, hys2, c=ch, s=size, label='0.1', alpha=opac)
plt.scatter(llo2+0.5, lla2+0.5 ,c='k', marker='P', s=150)#'+'

#%%
plt.savefig('figure4.pdf',bbox_inches='tight',pad_inches=0)

plt.show()
