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
size = 11
cmap='rainbow'

lat= -49.5
lon = 240.5
#lat2= 40
#lon2 = 220
lat2= 40.5
lon2 = 219.5

#lat3= -40
#lon3 = 25
#lat3= -37
#lon3 = 19
lat3= -44.5
lon3 = 20.5

#lat4= -55
#lon4 = 150
#lat4= -55
#lon4 = 140
#lat4= -47
#lon4 = 130

lon4 = 142.5
lat4 = -52.5

width = 50
hei = 15
width2 = 8
hei2 = 6
exte = [-360+(lon-42), -360+lon+6, lat-hei, lat+hei]
exte2 = [-360+(lon2-width2), -360+lon2+10, lat2-4, lat2+hei2]
exte3 = [(lon3-26), lon3+38, lat3-hei, lat3+25]
exte4 = [(lon4-55)-4, lon4+20, lat4-15, lat4+21]


width = 80/1.5
height = 40/1.5
width2 = 80
height2 = 40

r1 = -15
u1 = 0
exte = [-360+lon-width/2.+r1, -360+lon+width/2.+r1, lat-height/2.+u1, lat+height/2.+u1]
r2 = -0
u2 = 2
exte2 = [-360+lon2-width/2+r2, -360+lon2+width/2+r2, lat2-height/2+u2, lat2+height/2+u2]
r1 = 0
u1 = 7
exte3 = [-360+lon3-width2/2.+r1, -360+lon3+width2/2.+r1, lat3-height2/2.+u1, lat3+height2/2.+u1]
r2 = -22
u2 = 2
exte4 = [-360+lon4-width2/2.+r2, -360+lon4+width2/2.+r2, lat4-height2/2.+u2, lat4+height2/2.+u2]

ddeg = 1
sp = 6
dd = 10
cs = 2.0

fs = 16

vmin= -2; vmax = 28;

csea = 'whitesmoke'
cland = 'k'

ch = 'red'
cl = 'royalblue'
cd = 'y'

opac = 0.5

xloc = np.arange(-180, 180, 20)
yloc = np.arange(-90,90,20)
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
hlo4 = find_nearest(Lons,lon4)
hla4 = find_nearest(Lats,lat4)
hlo3 = find_nearest(Lons,lon3)
hla3 = find_nearest(Lats,lat3)

idx = np.where(np.logical_and(vLons==hlo, vLats==hla))
idx2 = np.where(np.logical_and(vLons==hlo2, vLats==hla2))
idx3 = np.where(np.logical_and(vLons==hlo3, vLats==hla3))
idx4 = np.where(np.logical_and(vLons==hlo4, vLats==hla4))
assert len(idx[0]==1)

hxs = nc_hr['lon'][idx]
hys = nc_hr['lat'][idx]
hstemp = nc_hr['temp'][idx]

hxs2 = nc_hr['lon'][idx2]
hys2 = nc_hr['lat'][idx2]
hstemp2 = nc_hr['temp'][idx2]

hxs3 = nc_hr['lon'][idx3]
hys3 = nc_hr['lat'][idx3]
hstemp3 = nc_hr['temp'][idx3]

hxs4 = nc_hr['lon'][idx4]
hys4 = nc_hr['lat'][idx4]
hstemp4 = nc_hr['temp'][idx4]

print('# particles in distribution: ',len(hstemp), len(hstemp2),len(hstemp3),len(hstemp4))
#%% Then the lowres with cs

dirReadlow = '/Volumes/HardDisk/POP/output/lowres/timeseries/'
nc_lr = Dataset(dirReadlow + 'timeseries_per_location_smagorinksi_wn_gm_Cs%.1f_ddeg%d_sp%d_dd%d.nc'%(cs,ddeg,sp,dd))

Lons = nc_lr['Lons'][:]
Lats = nc_lr['Lats'][:]
vLons = nc_lr['vLons'][:]
vLats = nc_lr['vLats'][:]

llo = find_nearest(Lons,lon)
lla = find_nearest(Lats,lat)
llo2 = find_nearest(Lons,lon2)
lla2 = find_nearest(Lats,lat2)
llo3 = find_nearest(Lons,lon3)
lla3 = find_nearest(Lats,lat3)
llo4 = find_nearest(Lons,lon4)
lla4 = find_nearest(Lats,lat4)
idx2 = np.where(np.logical_and(vLons==llo2, vLats==lla2))
idx = np.where(np.logical_and(vLons==llo, vLats==lla))
idx3 = np.where(np.logical_and(vLons==llo3, vLats==lla3))
idx4 = np.where(np.logical_and(vLons==llo4, vLats==lla4))

assert len(idx[0]==1)
lxs = nc_lr['lon'][idx]
lys = nc_lr['lat'][idx]
lstemp = nc_lr['temp'][idx]
lxs2 = nc_lr['lon'][idx2]
lys2 = nc_lr['lat'][idx2]
lstemp2 = nc_lr['temp'][idx2]
lxs3 = nc_lr['lon'][idx3]
lys3 = nc_lr['lat'][idx3]
lstemp3 = nc_lr['temp'][idx3]
lxs4 = nc_lr['lon'][idx4]
lys4 = nc_lr['lat'][idx4]
lstemp4 = nc_lr['temp'][idx4]
print('# particles in distribution: ',np.sum(0<lstemp), np.sum(0<lstemp2),np.sum(0<lstemp3),np.sum(0<lstemp4))
#%% Then the lowres with cs==0
cs = 0.0
dirReadlow = '/Volumes/HardDisk/POP/output/lowres/timeseries/'
nc_lr = Dataset(dirReadlow + 'timeseries_per_location_smagorinksi_wn_gm_Cs%.1f_ddeg%d_sp%d_dd%d.nc'%(cs,ddeg,sp,dd))

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
assert len(idx[0]==1)
lxs0 = nc_lr['lon'][idx]
lys0 = nc_lr['lat'][idx]
lstemp0 = nc_lr['temp'][idx]
lxs20 = nc_lr['lon'][idx2]
lys20 = nc_lr['lat'][idx2]
lstemp20 = nc_lr['temp'][idx2]

llo3 = find_nearest(Lons,lon3)
lla3 = find_nearest(Lats,lat3)
llo4 = find_nearest(Lons,lon4)
lla4 = find_nearest(Lats,lat4)
idx3 = np.where(np.logical_and(vLons==llo3, vLats==lla3))
idx4 = np.where(np.logical_and(vLons==llo4, vLats==lla4))
assert len(idx[0]==1)
lxs30 = nc_lr['lon'][idx3]
lys30 = nc_lr['lat'][idx3]
lstemp30 = nc_lr['temp'][idx3]
lxs40 = nc_lr['lon'][idx4]
lys40 = nc_lr['lat'][idx4]
lstemp40 = nc_lr['temp'][idx4]
print('# particles in distribution: ',np.sum(0<lstemp0), np.sum(0<lstemp20),np.sum(0<lstemp30),np.sum(0<lstemp40))
#%% The plot
fig = plt.figure(figsize=(10,6))

ax0 = plt.subplot(2,2,2, projection=projection)
plt.title('(b)', fontsize=fs)
ax0.add_feature(cartopy.feature.LAND, color=cland)
ax0.add_feature(cartopy.feature.OCEAN, color=csea, zorder=0)
g = ax0.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')
g.xlabels_top = False
g.ylabels_left = False
g.xformatter = LONGITUDE_FORMATTER
g.yformatter = LATITUDE_FORMATTER
g.xlabel_style = {'fontsize': fs-3}
g.ylabel_style = {'fontsize': fs-3}
g.xlocator = mticker.FixedLocator(xloc)
g.ylocator = mticker.FixedLocator(yloc)
ax0.set_extent(exte, ccrs.PlateCarree())

plt.scatter(hxs, hys, c=ch, s=size, label='0.1', alpha=opac)
plt.scatter(lxs, lys, c=cl, s=size, label='1', alpha=opac)
plt.scatter(lxs0, lys0, c=cd, s=size, label='1, cs=0', alpha=opac)
plt.scatter(llo+0.5, lla +0.5,c='k', marker='P', s=150)

print('# particles in distribution: ',np.sum(0<hxs[0]), np.sum(0<lxs[0]), np.sum(0<lxs0[0]))

ax0 = plt.subplot(2,2,1, projection=projection)
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
g.xlabel_style = {'fontsize': fs-3}
g.ylabel_style = {'fontsize': fs-3}
ax0.set_extent(exte2, ccrs.PlateCarree())

plt.scatter(hxs2, hys2, c=ch, s=size, label='0.1', alpha=opac)
plt.scatter(lxs2, lys2, c=cl, s=size, label='1', alpha=opac)
plt.scatter(lxs20, lys20, c=cd, s=size, label='1, cs=0', alpha=opac)
plt.scatter(llo2+0.5, lla2+0.5 ,c='k', marker='P', s=150)#'+'

print('# particles in distribution: ',np.sum(0<hxs2[0]), np.sum(0<lxs2[0]), np.sum(0<lxs20[0]))
#%%
plt.subplots_adjust(hspace=0.35)
ax0 = plt.subplot(2,2,3, projection=projection)
plt.title('(c)', fontsize=fs)
ax0.add_feature(cartopy.feature.LAND, color=cland)
ax0.add_feature(cartopy.feature.OCEAN, color=csea, zorder=0)
g = ax0.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')
g.xlabels_top = False
g.ylabels_right = False
g.xformatter = LONGITUDE_FORMATTER
g.yformatter = LATITUDE_FORMATTER
g.xlabel_style = {'fontsize': fs-3}
g.ylabel_style = {'fontsize': fs-3}
g.xlocator = mticker.FixedLocator(xloc)
g.ylocator = mticker.FixedLocator(yloc)
ax0.set_extent(exte3, ccrs.PlateCarree())

plt.scatter(hxs3, hys3, c=ch, s=size, label='$R_{0.1}$', alpha=opac)
plt.scatter(lxs3, lys3, c=cl, s=size, label='$R_{1md}$', alpha=opac)
plt.scatter(lxs30, lys30, c=cd, s=size, label='$R_{1m}$', alpha=opac)
plt.scatter(llo3+0.5, lla3+0.5 ,c='k', marker='P', s=150)#'+'

print('# particles in distribution: ',np.sum(0<hxs3[0]), np.sum(0<lxs3[0]), np.sum(0<lxs30[0]))
      
plt.legend(loc=4, prop={'size': 10}, title='Configuration')

ax0 = plt.subplot(2,2,4, projection=projection)
plt.title('(d)', fontsize=fs)
ax0.add_feature(cartopy.feature.LAND, color=cland)
ax0.add_feature(cartopy.feature.OCEAN, color=csea, zorder=0)
g = ax0.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')
g.xlabels_top = False
g.ylabels_left = False
g.xformatter = LONGITUDE_FORMATTER
g.yformatter = LATITUDE_FORMATTER
g.xlocator = mticker.FixedLocator(np.arange(exte4[0], exte4[1]+20, 20))
g.xlabel_style = {'fontsize': fs-3}
g.ylabel_style = {'fontsize': fs-3}
g.xlocator = mticker.FixedLocator(xloc)
g.ylocator = mticker.FixedLocator(yloc)
ax0.set_extent(exte4, ccrs.PlateCarree())

plt.scatter(hxs4, hys4, c=ch, s=size, label='0.1', alpha=opac)
plt.scatter(lxs4, lys4, c=cl, s=size, label='1', alpha=opac)
plt.scatter(lxs40, lys40, c=cd, s=size, label='1, cs=0', alpha=opac)
plt.scatter(llo4+0.5, lla4+0.5 ,c='k', marker='P', s=150)#'+'
print('# particles in distribution: ',np.sum(0<hxs4[0]), np.sum(0<lxs4[0]), np.sum(0<lxs40[0]))
#%%

plt.savefig('figure6.pdf',bbox_inches='tight',pad_inches=0)

plt.show()
