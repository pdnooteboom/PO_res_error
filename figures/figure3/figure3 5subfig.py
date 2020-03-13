#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 15:41:33 2019

@author: nooteboom
"""

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
from matplotlib.lines import Line2D

sns.set_style("whitegrid")
sns.set_context("paper")
fs = 17
font = {'size'   : fs}
matplotlib.rc('font', **font)

#variables
sp = 6
dd = 10
projection = ccrs.PlateCarree(180)
exte = [1, 360, -75, 72]
exte2 = [-179, 181, -75, 72]
Cs = 2.0
ddeg = 1
cmap2 = 'coolwarm' # For the surface area
cmap3 = 'hot'# For the average travel distance
vssurf = [0,17]

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
    
avgd, surf, Lons, Lats = calc_fields(name = '/Volumes/HardDisk/POP/output/highres/timeseries/timeseries_per_location_ddeg1_sp6_dd10_tempresmonmean.nc')
idxlat = np.logical_and(Lats>=maxminlat, Lats<=minmaxlat)
idxlon = np.logical_and(Lons>=maxminlon, Lons<=minmaxlon)
assert (idxlon==True).all()
surf =  np.flip(surf[idxlat],0)
land = np.full(surf.shape, np.nan); land[surf==0] = 1;
surf[surf==0] = np.nan
surf_temp = np.nanmean(surf) / 10**5.
#%%
avgd, surf, Lons, Lats = calc_fields(name = '/Volumes/HardDisk/POP/output/highres/timeseries/timeseries_per_location_ddeg1_sp25_dd10_tempresmonmean.nc')
idxlat = np.logical_and(Lats>=maxminlat, Lats<=minmaxlat)
idxlon = np.logical_and(Lons>=maxminlon, Lons<=minmaxlon)
assert (idxlon==True).all()
surf = np.flip(surf[idxlat],0)
surf[surf==0] = np.nan
surf_temp2 = np.nanmean(surf) / 10**5.

avgd50, surf50hr, Lons, Lats = calc_fields(name = '/Volumes/HardDisk/POP/output/highres/timeseries/timeseries_per_location_ddeg1_sp25_dd10_tempres5.nc')
surf50hr[surf50hr==0] = np.nan
surf50mean = np.nanmean(surf50hr)  / 10**5.

avgd, surf, Lons, Lats = calc_fields(name = '/Volumes/HardDisk/POP/output/highres/timeseries/timeseries_per_location_ddeg%d_sp%d_dd%d_tempres5_ds2.nc'%(ddeg,sp,dd))
idxlat = np.logical_and(Lats>=maxminlat, Lats<=minmaxlat)
idxlon = np.logical_and(Lons>=maxminlon, Lons<=minmaxlon)
assert (idxlon==True).all()
surf_highres = np.flip(surf[idxlat],0)
highres_surf = surf_highres.copy()
highres_surf[highres_surf==0] = np.nan
highres_surf = np.nanmean(highres_surf) / 10**5.
print('highres_surf: ',highres_surf)

avgd, surf, Lons, Lats = calc_fields(name = '/Volumes/HardDisk/POP/output/lowres/timeseries/timeseries_per_location_smagorinksi_wn_Cs%.1f_ddeg%d_sp%d_dd%d.nc'%(0.0,ddeg,sp,dd))
idxlat = np.logical_and(Lats>=maxminlat, Lats<=minmaxlat)
idxlon = np.logical_and(Lons>=maxminlon, Lons<=minmaxlon)
assert (idxlon==True).all()
surf_lr = np.flip(surf[idxlat],0)

avgd, surf, Lons, Lats = calc_fields(name = '/Volumes/HardDisk/POP/output/lowres/timeseries/timeseries_per_location_smagorinksi_wn_Cs%.1f_ddeg%d_sp%d_dd%d.nc'%(Cs,ddeg,sp,dd))
idxlat = np.logical_and(Lats>=maxminlat, Lats<=minmaxlat)
idxlon = np.logical_and(Lons>=maxminlon, Lons<=minmaxlon)
assert (idxlon==True).all()
surf_lr2 = np.flip(surf[idxlat],0)

sns.set_style("darkgrid")
sns.set_context("paper")
fs = 14 # fontsize
si = 141
lw = 2 # linewidth

sp1 = 6
sp2 = 25

color1 = 'k'
color2 = 'red'
color3 = 'k'
#% Load the data
CS =  np.array([0., 0.25, 0.5, 1.0, 2.0, 5.0])
CS50 = np.array([0., 0.25, 0.5, 1.0, 2.0, 5.0])
cs = Cs

sur  = np.zeros(len(CS))
sur50  = np.zeros(len(CS50))
surgm  = np.zeros(len(CS))
sur50gm  = np.zeros(len(CS50))
for j in range(len(CS)):
    if(CS[j]!=0.25):
        avgd, surf, Lons, Lats = calc_fields(name = '/Volumes/HardDisk/POP/output/lowres/timeseries/timeseries_per_location_smagorinksi_wn_Cs%.1f_ddeg1_sp%d_dd10.nc'%(CS[j],sp1))
    else:
        avgd, surf, Lons, Lats = calc_fields(name = '/Volumes/HardDisk/POP/output/lowres/timeseries/timeseries_per_location_smagorinksi_wn_Cs%.2f_ddeg1_sp%d_dd10.nc'%(CS[j],sp1))    
    idxlat = np.logical_and(Lats>=maxminlat, Lats<=minmaxlat)
    idxlon = np.logical_and(Lons>=maxminlon, Lons<=minmaxlon)
    assert (idxlon==True).all()
    if(CS[j]==cs):
        surf_cs =  np.flip(surf[idxlat],0)
    surf = surf[idxlat]
    surf[surf==0] = np.nan
    sur[j] = np.nanmean(surf) / 10**5.
    
    if(CS[j]!=0.25):
        avgd, surf, Lons, Lats = calc_fields(name = '/Volumes/HardDisk/POP/output/lowres/timeseries/timeseries_per_location_smagorinksi_wn_gm_Cs%.1f_ddeg1_sp%d_dd10.nc'%(CS[j],sp1))
    else:
        avgd, surf, Lons, Lats = calc_fields(name = '/Volumes/HardDisk/POP/output/lowres/timeseries/timeseries_per_location_smagorinksi_wn_gm_Cs%.2f_ddeg1_sp%d_dd10.nc'%(CS[j],sp1)) 
    idxlat = np.logical_and(Lats>=maxminlat, Lats<=minmaxlat)
    idxlon = np.logical_and(Lons>=maxminlon, Lons<=minmaxlon)
    assert (idxlon==True).all()
    if(CS[j]==cs):
        surf_gm =  np.flip(surf[idxlat],0)
    surf = surf[idxlat]   
    surf[surf==0] = np.nan
    surgm[j] = np.nanmean(surf) / 10**5.


for j in range(len(CS50)):
    if(CS50[j] != 0.25):
        avgd, surf, Lons, Lats = calc_fields(name = '/Volumes/HardDisk/POP/output/lowres/timeseries/timeseries_per_location_smagorinksi_wn_Cs%.1f_ddeg1_sp%d_dd10.nc'%(CS50[j],sp2))          
    else:   
        avgd, surf, Lons, Lats = calc_fields(name = '/Volumes/HardDisk/POP/output/lowres/timeseries/timeseries_per_location_smagorinksi_wn_Cs%.2f_ddeg1_sp%d_dd10.nc'%(CS50[j],sp2))    
    idxlat = np.logical_and(Lats>=maxminlat, Lats<=minmaxlat)
    idxlon = np.logical_and(Lons>=maxminlon, Lons<=minmaxlon)
    assert (idxlon==True).all()
    surf = surf[idxlat] 
    surf[surf==0] = np.nan
    sur50[j] = np.nanmean(surf) / 10**5.

    if(CS[j]!=0.25):
        avgd, surf, Lons, Lats = calc_fields(name = '/Volumes/HardDisk/POP/output/lowres/timeseries/timeseries_per_location_smagorinksi_wn_gm_Cs%.1f_ddeg1_sp%d_dd10.nc'%(CS50[j],sp2))
    else:
        avgd, surf, Lons, Lats = calc_fields(name = '/Volumes/HardDisk/POP/output/lowres/timeseries/timeseries_per_location_smagorinksi_wn_gm_Cs%.2f_ddeg1_sp%d_dd10.nc'%(CS50[j],sp2))  
    idxlat = np.logical_and(Lats>=maxminlat, Lats<=minmaxlat)
    idxlon = np.logical_and(Lons>=maxminlon, Lons<=minmaxlon)
    assert (idxlon==True).all()
    surf = surf[idxlat]  
    surf[surf==0] = np.nan
    sur50gm[j] = np.nanmean(surf) / 10**5.

plt.plot(CS, sur)
plt.plot(CS50, sur50)
plt.plot(CS, surgm, '--')
plt.plot(CS50, sur50gm, '--')
plt.scatter([0], [surf_temp])
plt.show()
#%% start figure
fig = plt.figure(figsize=(19,15))
grid = plt.GridSpec(3, 24, wspace=0., hspace=0.4)
#% subplot (a)
ax = plt.subplot(grid[0, :12], projection=projection)#plt.subplot(2,2,1, projection=projection)
plt.title('(a) $R_{0.1}$', fontsize=fs)

g = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')
g.xlocator = mticker.FixedLocator([-180,-90, -0, 90, 180])
g.xlabels_top = False
g.ylabels_right = False
g.xlabels_bottom = False
g.xlabel_style = {'fontsize': fs}
g.ylabel_style = {'fontsize': fs}
g.xformatter = LONGITUDE_FORMATTER
g.yformatter = LATITUDE_FORMATTER
g.ylocator = mticker.FixedLocator([-75,-50,-25, 0, 25, 50, 75, 100])
ax.set_extent(exte, ccrs.PlateCarree())

plt.imshow(surf_highres/10.**5, vmin=vssurf[0], vmax=vssurf[1], extent = exte2, transform=ccrs.PlateCarree(), 
           cmap=cmap3, zorder = 0)
#land = np.full(avgd.shape, np.nan); land[surf==0] = 1;
plt.imshow(land, vmin=0, vmax=1.6, extent = exte2, transform=ccrs.PlateCarree(), cmap='binary', zorder = 0)

#% subplot (b)
ax = plt.subplot(grid[0, 12:], projection=projection)#plt.subplot(2,2,2, projection=projection)
plt.title('(b) $R_{1m}$', fontsize=fs)

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

#ax.set_xticks([0., 90., 180., 270., 360.], crs=ccrs.PlateCarree())
#ax.set_xticklabels([0., 90., 180., 270., 360.], fontsize=fs)
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.grid(linewidth=2, color='black', alpha=0., linestyle='--')

plt.imshow(surf_lr/10**5., vmin=vssurf[0], vmax=vssurf[1], extent = exte2, transform=ccrs.PlateCarree(), 
           cmap=cmap3, zorder = 0)
#land = np.full(avgd.shape, np.nan); land[surf==0] = 1;
plt.imshow(land, vmin=0, vmax=1.6, extent = exte2, transform=ccrs.PlateCarree(), cmap='binary', zorder = 0)

#% subplot (c)
ax = plt.subplot(grid[1, :12], projection=projection)#plt.subplot(2,2,3, projection=projection)
plt.title('(c) $R_{1md}$, $c_s$=%.1f'%(Cs), fontsize=fs)

g = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')
g.xlabels_top = False
g.ylabels_right = False
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

im2 = plt.imshow(surf_lr2/10**5., vmin=vssurf[0], vmax=vssurf[1], extent = exte2, transform=ccrs.PlateCarree(), 
           cmap=cmap3, zorder = 0)
#land = np.full(avgd.shape, np.nan); land[surf==0] = 1;
plt.imshow(land, vmin=0, vmax=1.6, extent = exte2, transform=ccrs.PlateCarree(), cmap='binary', zorder = 0)

#% subplot (d)

ax = plt.subplot(grid[1, 12:], projection=projection)#plt.subplot(2,2,3, projection=projection)
plt.title('(d) $R_{1mdb}$, $c_s$=%.1f'%(Cs), fontsize=fs)

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

im2 = plt.imshow(surf_gm/10**5., vmin=vssurf[0], vmax=vssurf[1], extent = exte2, transform=ccrs.PlateCarree(), 
           cmap=cmap3, zorder = 0)
#land = np.full(avgd.shape, np.nan); land[surf==0] = 1;
plt.imshow(land, vmin=0, vmax=1.6, extent = exte2, transform=ccrs.PlateCarree(), cmap='binary', zorder = 0)
    
#%
dsWD = [4,7,4,7] # the line dash of the first configuration
dsWD2 = [4,7,4,7] # the line dash of the bolus configuration

ax = plt.subplot(grid[2, 12:-1])
plt.title('(e)', fontsize=fs)

plt.xlabel('$c_s$', fontsize=fs)

a0 = sns.lineplot(x=CS, y=np.full(len(CS),highres_surf), linewidth=lw,
                   color=color1, zorder=1)
a0.lines[0].set_dashes(dsWD)
a0 = sns.lineplot(x=CS, y=np.full(len(CS),surf_temp), linewidth=lw,
                   color=color1, zorder=1)
a0.lines[1].set_linestyle(":")
a0 = sns.lineplot(x=CS, y=np.full(len(CS),surf50mean), linewidth=lw,
                   color=color2, zorder=1)
a0.lines[2].set_dashes(dsWD)
a0 = sns.lineplot(x=CS, y=np.full(len(CS),surf_temp2), linewidth=lw,
                   color=color2, zorder=1)
a0.lines[3].set_linestyle(":")

sns.lineplot(x=CS, y=sur, color=color1, linewidth=lw, zorder=10)
sns.scatterplot(x=CS, y=sur, color=color1, s=si, zorder=11)

sns.lineplot(x=CS,y=surgm, linewidth = lw, ax=ax, color=color1, 
             zorder=9)
sns.scatterplot(x=CS,y=surgm, ax=ax, color=color1, s=si, zorder=12,
                   legend=False, marker="^")

sns.lineplot(x=CS50, y=sur50, color=color2, linewidth=lw, zorder=10)
sns.scatterplot(x=CS50, y=sur50, color=color2, s=si, zorder=11)

sns.lineplot(x=CS50,y=sur50gm, linewidth = lw, ax=ax, color=color2, 
             zorder=9)
sns.scatterplot(x=CS50,y=sur50gm, ax=ax, color=color2, s=si, zorder=12,
                   legend=False, marker="^")

for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(fs)
ax.set_ylim(0,7.5)

ax.set_ylabel('surface area (10$^5$ km$^2$)', fontsize=fs, color = color1)


for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(fs) 

lw = 2
colo = 'k'
legend_el = [Line2D([0], [0], dashes=dsWD, color=colo, lw=lw, label='$R_{0.1}$'), 
             Line2D([0], [0], linestyle=':', color=colo, lw=lw, label='$R_{0.1m}$'), 
             Line2D([0], [0], linestyle='-', marker='o', markersize=8, color=colo, lw=lw, label='$W_d(R_{0.1}$, $R_{1m}$/ $R_{1md}$)'), 
             Line2D([0], [0], linestyle='-', marker='^', markersize=8, color=colo, lw=lw, label='$W_d(R_{0.1}$, $R_{1mb}$/ $R_{1mdb}$)')]
#first_legend = plt.legend(handles=legend_el, title='Configuration',loc=4, fontsize=fs, bbox_to_anchor=(0., .01, 1., .022))
first_legend = ax.legend(handles=legend_el, title='Configuration', fontsize=fs, loc='center right', bbox_to_anchor=(-0.1, 0.2))


ax2 = plt.gca().add_artist(first_legend)

legend_el = [Line2D([0], [0], linestyle='solid', color=color1, lw=lw, label='$w_f=6$'), 
             Line2D([0], [0], linestyle='solid', color=color2, lw=lw, label='$w_f=25$')]
#plt.legend(handles=legend_el, title='Sinking speed (m/day)',loc=4, fontsize=fs, bbox_to_anchor=(0., .52, 1., .102))
ax.legend(handles=legend_el, title='Sinking speed (m/day)', fontsize=fs, loc='center right', bbox_to_anchor=(-0.1, 0.65))


#% final
#fig.subplots_adjust(bottom=0.17)
#cbar_ax = fig.add_axes([0.11, 0.05, 0.35, 0.07])
#cbar_ax.set_visible(False)
#cbar = fig.colorbar(im2, ax=cbar_ax, orientation = 'horizontal', fraction = 1.2)
#cbar.ax.xaxis.set_label_position('bottom')
#cbar.ax.set_xlabel('$10^5$ km$^2$', fontsize=fs)
#cbar.ax.tick_params(labelsize=fs) 
#cbar.set_ticklabels([1,2,3,4]) 

#fig.subplots_adjust(bottom=0.17)
cbar_ax = fig.add_axes([0.135, 0.285, 0.35, 0.07])
cbar_ax.set_visible(False)
cbar = fig.colorbar(im2, ax=cbar_ax, orientation = 'horizontal', fraction = 1.2,
                    aspect=18)
cbar.ax.xaxis.set_label_position('bottom')
cbar.ax.set_xlabel('$10^5$ km$^2$', fontsize=fs)
cbar.ax.tick_params(labelsize=fs)

plt.savefig('figure3_withandwithoutbolus.pdf',bbox_inches='tight',pad_inches=0)
plt.show()