#!/usr/bin/env python3
import numpy as np
import matplotlib.pylab as plt
from netCDF4 import Dataset
import matplotlib
import cartopy.crs as ccrs
import seaborn as sns
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import math
from numba import jit
import cartopy.mpl.ticker as cticker
#from matplotlib.lines import Line2D

sns.set_style("whitegrid")
sns.set_context("paper")
fs = 18
font = {'size'   : fs}
matplotlib.rc('font', **font)

#variables
sp = 6
dd = 10
projection = ccrs.PlateCarree(180)
exte = [1, 360, -74, 71]
exte2 = [-179, 181, -74, 71]
Cs = 2.0
ddeg = 1


dist =True
surfa = True#False
cmap2 = 'seismic'#'magma'# For the average travel distance
cmap3 = 'BrBG_r'# 'summer'#'For the surface area
vssurf = [0,3]# [0,17]#
vsdist = [0,5]#[0,27]#
vssurf = [-2.5,2.5]# [0,17]#
vsdist = [-5,5]#[0,27]#
#vssurf = [0,17]#
#vsdist = [0,27]#

sns.set_style("darkgrid")
sns.set_context("paper")
si = 141
lw = 2 # linewidth

sp1 = 6
sp2 = 25

#%%
@jit(nopython=True)
def find_nearest_index(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

@jit(nopython=True)
def distance(origin, destination):
    lat1, lon1 = origin
    lat2, lon2 = destination
    radius = 6371.1 # km

    dlat = math.radians(lat2-lat1)
    dlon = math.radians(lon2-lon1)
    a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1)) \
        * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    d = radius * c

    return d

@jit(nopython=True)
def avgdistf(lat, lon, avgdist, tsl, vLons, vLats, ml):
    for i in range(len(vLons)):
        dist = 0.
        n = 0.
        for j in range(ml):
           if(lat[i,j]<1000):
               n += 1.
               dist += distance([vLats[i], vLons[i]], [lat[i,j], lon[i,j]])
        if(n>0):
            avgdist[i] = dist / n
        else:
            avgdist[i] = np.nan
    return avgdist

@jit(nopython=True)
def uniqueloc(lats, lons, la, lo):
    bo = True
    if(len(lons)>0):
        for i in range(len(lons)):
            if(lons[i]==lo):
                if(lats[i]==la):
                    bo = False
    return bo

@jit(nopython=True)
def surfacef(lat, lon, surface, tsl, vLons, vLats, Lons, Lats, ml):
    res = 1. # resolution of the binning
    for i in range(len(vLons)):
        surf = 0.
        lons = []
        lats = []
        for j in range(ml):
            if(lat[i,j]<1000):
                lo = find_nearest_index(Lons, lon[i,j])
                la = find_nearest_index(Lats, lat[i,j])
                if(uniqueloc(lats, lons, la, lo)):
                    lats.append(la)
                    lons.append(lo)
        for j in range(len(lons)):
            surf += distance([lats[j]-res/2.,lons[j]-res/2.],[lats[j]-res/2.,lons[j]+res/2.]) * distance([lats[j]+res/2.,lons[j]-res/2.],[lats[j]-res/2.,lons[j]-res/2.])        
        surface[i] = surf         
    return surface

def calc_fields(name = '', ml=131):
    ncf = Dataset(name)
    Lons = ncf['Lons'][:]
    Lats = ncf['Lats'][:]
    vLons = ncf['vLons'][:]
    vLats = ncf['vLats'][:]
    assert len(Lons)*len(Lats)==len(vLons)    
    lat = ncf['lat'][:]
    lon = ncf['lon'][:]
    tsl = ncf['temp'][:].shape[1]

    avgdist = np.zeros(len(vLons))
    surface = np.zeros(len(vLons))
    avgdist = avgdistf(lat, lon, avgdist, tsl, vLons, vLats, ml)
    surface = surfacef(lat, lon, surface, tsl, vLons, vLats, Lons, Lats, ml)
    avgdist, surface = avgdist.reshape(len(Lats), len(Lons)), surface.reshape(len(Lats), len(Lons))
    avgdist = np.concatenate((avgdist[:,180:-2], avgdist[:,:180]),axis=1);
    surface = np.concatenate((surface[:,180:-2],surface[:,:180]),axis=1);
    Lons = np.append(Lons[180:-2]-360,Lons[:180])
    
    return avgdist,surface, Lons, Lats

#%%
maxminlat = -76
minmaxlat = 71
minmaxlon = 359 - 180
maxminlon = 0 - 180  
    
#avgd, surf, Lons, Lats = calc_fields(name = '/Volumes/HardDisk/POP/output/highres/timeseries/timeseries_per_location_ddeg1_sp6_dd10_tempresmonmean.nc')
#avgd, surf = np.flip(avgd,0), np.flip(surf,0)
#land = np.full(avgd.shape, np.nan); land[surf==0] = 1;
#surf[surf==0] = np.nan
#surf_temp = np.nanmean(surf) / 10**5.

##%%
#avgd, surf, Lons, Lats = calc_fields(name = '/Volumes/HardDisk/POP/output/highres/timeseries/timeseries_per_location_ddeg1_sp25_dd10_tempresmonmean.nc')
#avgd, surf = np.flip(avgd,0), np.flip(surf,0)
#land = np.full(avgd.shape, np.nan); land[surf==0] = 1;
#surf[surf==0] = np.nan
#surf_temp2 = np.nanmean(surf) / 10**5.
#%%


#%% start figure
#%
avgd, surf, Lons, Lats = calc_fields(name = '/Volumes/HardDisk/POP/output/highres/timeseries/timeseries_per_location_ddeg%d_sp%d_dd%d_tempres5_ds2.nc'%(ddeg,sp,dd))
idxlat = np.logical_and(Lats>=maxminlat, Lats<=minmaxlat)
idxlon = np.logical_and(Lons>=maxminlon, Lons<=minmaxlon)
assert (idxlon==True).all()
highres_surf = surf[idxlat].copy()
highres_avgd = avgd[idxlat].copy()
    
#avgd, surf = np.flip(avgd,0), np.flip(surf,0)

avgd, surf, Lons, Lats = calc_fields(name = '/Volumes/HardDisk/POP/output/highres/timeseries/timeseries_per_location_ddeg1_sp%d_dd10_tempresmonmean.nc'%(sp))
idxlat = np.logical_and(Lats>=maxminlat, Lats<=minmaxlat)
idxlon = np.logical_and(Lons>=maxminlon, Lons<=minmaxlon)
assert (idxlon==True).all()
surf = -surf[idxlat] + highres_surf
avgd = -avgd[idxlat] + highres_avgd

avgd, surf = np.flip(avgd,0), np.flip(surf,0)
#highres_avgd, highres_surf = np.flip(highres_avgd,0), np.flip(highres_surf,0)
#%%
fig = plt.figure(figsize=(18,6))
#% subplot (a)
ax = plt.subplot(1,2,1, projection=projection)
plt.title('(a)', fontsize=fs)

g = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')
g.xlabels_top = False
g.ylabels_right = False
g.ylabels_left = True
g.xlabels_bottom = False
g.xlabel_style = {'fontsize': fs}
g.ylabel_style = {'fontsize': fs}
g.xformatter = LONGITUDE_FORMATTER
g.yformatter = LATITUDE_FORMATTER
g.xlocator = mticker.FixedLocator([-180,-90, -0, 90, 180])
g.ylocator = mticker.FixedLocator([-75,-50,-25, 0, 25, 50, 75, 100])
ax.set_extent(exte, ccrs.PlateCarree())

ax.set_xticks([0., 90., 180., 270., 360.], crs=ccrs.PlateCarree())
ax.set_xticklabels([0., 90., 180., 270., 360.], fontsize=fs)
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.grid(linewidth=2, color='black', alpha=0., linestyle='--')

im = plt.imshow(surf/10**5., vmin=vssurf[0], vmax=vssurf[1], extent = exte2, transform=ccrs.PlateCarree(), 
           cmap=cmap3, zorder = 0)
    
land = np.full(avgd.shape, np.nan);
land[np.isnan(avgd)] = 1
plt.imshow(land, vmin=0, vmax=1.9, extent = exte2, transform=ccrs.PlateCarree(), cmap='binary', zorder = 0)

#% subplot (b)
ax = plt.subplot(1,2,2, projection=projection)
plt.title('(b)', fontsize=fs)

g = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')
g.xlabels_top = False
g.ylabels_right = False
g.ylabels_left = False
g.xlabels_bottom = False
g.xlabel_style = {'fontsize': fs}
g.ylabel_style = {'fontsize': fs}
g.xformatter = LONGITUDE_FORMATTER
g.yformatter = LATITUDE_FORMATTER
g.xlocator = mticker.FixedLocator([-180,-90, -0, 90, 180])
g.ylocator = mticker.FixedLocator([-75,-50,-25, 0, 25, 50, 75, 100])
ax.set_extent(exte, ccrs.PlateCarree())

ax.set_xticks([0., 90., 180., 270., 360.], crs=ccrs.PlateCarree())
ax.set_xticklabels([0., 90., 180., 270., 360.], fontsize=fs)
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.grid(linewidth=2, color='black', alpha=0., linestyle='--')

im2 = plt.imshow(avgd/100., vmin=vsdist[0], vmax=vsdist[1], extent = exte2, transform=ccrs.PlateCarree(), cmap=cmap2, zorder = 0)
    
plt.imshow(land, vmin=0, vmax=1.9, extent = exte2, transform=ccrs.PlateCarree(), cmap='binary', zorder = 0)

#% final
#fig.subplots_adjust(bottom=0.1)

cbar_ax = fig.add_axes([0.13, 0.05, 0.35, 0.07])
cbar_ax.set_visible(False)
cbar = fig.colorbar(im, ax=cbar_ax, orientation = 'horizontal', fraction = 1.2, aspect=17)
cbar.ax.xaxis.set_label_position('bottom')
cbar.ax.set_xlabel('$10^5$ km$^2$', fontsize=fs)
cbar.ax.tick_params(labelsize=fs) 
#cbar.set_ticklabels([1,2,3,4]) 

cbar_ax = fig.add_axes([0.55, 0.05, 0.35, 0.07])
cbar_ax.set_visible(False)
cbar = fig.colorbar(im2, ax=cbar_ax, orientation = 'horizontal', fraction = 1.2, aspect=17)
cbar.ax.xaxis.set_label_position('bottom')
cbar.ax.set_xlabel('$10^2$ km', fontsize=fs)
cbar.ax.tick_params(labelsize=fs) 


plt.savefig('figure5.pdf',bbox_inches='tight',pad_inches=0)
plt.show()