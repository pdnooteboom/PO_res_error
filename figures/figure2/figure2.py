#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 10:14:33 2019

@author: nooteboom
"""
import matplotlib.pylab as plt
import numpy as np
from netCDF4 import Dataset
import cartopy.crs as ccrs
import seaborn as sns
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import matplotlib.colors as colors
import cartopy.mpl.ticker as cticker
from matplotlib.lines import Line2D


#Some parameters for the figure
projection = ccrs.PlateCarree(180)
exte = [1, 360, -75, 72]#[-179, 181, -76, 86]#[1,359,-75,81]#
exte2 = exte
extemm = [-179, 181, -76, 86]
exte3 = [-179, 181, -76, 71]
cmap = 'viridis'
  
sns.set_style("darkgrid")
sns.set_context("paper")
fs = 17 # fontsize
si = 140
lw = 2 # linewidth

sp1 = 6
sp2 = 25

color1 = 'k'
color2 = 'red'
color3 = 'dodgerblue'
#%% Load the data
CS =  np.array([0., 0.25, 0.5, 1.0, 2.0, 5.0])
CS50 = np.array([0., 0.25, 0.5, 1.0, 2.0, 5.0])
cs = 2.0

dirRead = 'Wd_output/' 
lsq = np.zeros(len(CS))
lsq50 = np.zeros(len(CS50))
lsqgm = np.zeros(len(CS))
lsq50gm = np.zeros(len(CS50))

maxminlat = -1000000
minmaxlat = 10000000
maxminlon = -1000000
minmaxlon = 10000000
for j in range(len(CS)):
    if(CS[j]!=0.25):
        ncf = Dataset(dirRead + 'OTs_sp%d_cs%.1f.nc'%(sp1,CS[j]))
    else:
        ncf = Dataset(dirRead + 'OTs_sp%d_cs%.2f.nc'%(sp1,CS[j]))
    maxminlat = max(maxminlat, np.min(ncf['Lats']))
    minmaxlat = min(minmaxlat, np.max(ncf['Lats']))
    maxminlon = max(maxminlon, np.min(ncf['Lons']))
    minmaxlon = min(minmaxlon, np.max(ncf['Lons']))
    ncf.close()
    
    if(CS[j]!=0.25):
        ncf = Dataset(dirRead + 'OTs_sp%d_cs%.1fgm.nc'%(sp1,CS[j]))
    else:
        ncf = Dataset(dirRead + 'OTs_sp%d_cs%.2fgm.nc'%(sp1,CS[j]))
    maxminlat = max(maxminlat, np.min(ncf['Lats']))
    minmaxlat = min(minmaxlat, np.max(ncf['Lats']))
    maxminlon = max(maxminlon, np.min(ncf['Lons']))
    minmaxlon = min(minmaxlon, np.max(ncf['Lons']))
    ncf.close()
for j in range(len(CS50)):
    if(CS[j]!=0.25):
        ncf = Dataset(dirRead + 'OTs_sp%d_cs%.1f.nc'%(sp2,CS50[j]))
    else:
        ncf = Dataset(dirRead + 'OTs_sp%d_cs%.2f.nc'%(sp2,CS50[j]))
    maxminlat = max(maxminlat, np.min(ncf['Lats']))
    minmaxlat = min(minmaxlat, np.max(ncf['Lats']))
    maxminlon = max(maxminlon, np.min(ncf['Lons']))
    minmaxlon = min(minmaxlon, np.max(ncf['Lons']))
    ncf.close()
    if(CS[j]!=0.25):
        ncf = Dataset(dirRead + 'OTs_sp%d_cs%.1fgm.nc'%(sp2,CS[j]))
    else:
        ncf = Dataset(dirRead + 'OTs_sp%d_cs%.2fgm.nc'%(sp2,CS[j]))
    maxminlat = max(maxminlat, np.min(ncf['Lats']))
    minmaxlat = min(minmaxlat, np.max(ncf['Lats']))
    maxminlon = max(maxminlon, np.min(ncf['Lons']))
    minmaxlon = min(minmaxlon, np.max(ncf['Lons']))
    ncf.close()
ncf = Dataset(dirRead + 'OTs_nultest.nc')
maxminlat = max(maxminlat, np.min(ncf['Lats']))
minmaxlat = min(minmaxlat, np.max(ncf['Lats']))
maxminlon = max(maxminlon, np.min(ncf['Lons']))
minmaxlon = min(minmaxlon, np.max(ncf['Lons']))
ncf.close()
ncf = Dataset(dirRead + 'OTs_monmean_sp25.nc')
maxminlat = max(maxminlat, np.min(ncf['Lats']))
minmaxlat = min(minmaxlat, np.max(ncf['Lats']))
maxminlon = max(maxminlon, np.min(ncf['Lons']))
minmaxlon = min(minmaxlon, np.max(ncf['Lons']))
ncf = Dataset(dirRead + 'OTs_monmean_sp6.nc')
maxminlat = max(maxminlat, np.min(ncf['Lats']))
minmaxlat = min(minmaxlat, np.max(ncf['Lats']))
maxminlon = max(maxminlon, np.min(ncf['Lons']))
minmaxlon = min(minmaxlon, np.max(ncf['Lons']))
minmaxlon -= 1

for j in range(len(CS)):
    if(CS[j]!=0.25):
        ncf = Dataset(dirRead + 'OTs_sp%d_cs%.1f.nc'%(sp1,CS[j]))
    else:
        ncf = Dataset(dirRead + 'OTs_sp%d_cs%.2f.nc'%(sp1,CS[j]))
    idxlat = np.logical_and(ncf['Lats']>=maxminlat, ncf['Lats']<=minmaxlat)
    idxlon = np.logical_and(ncf['Lons']>=maxminlon, ncf['Lons']<=minmaxlon)
    if(CS[j]==cs):
        distance = np.flip(ncf['Wasserstein distance'][idxlat, idxlon],0)
    if(CS[j]==0):
        distance_cs0 = np.flip(ncf['Wasserstein distance'][idxlat, idxlon],0)
#        distance_cs0 = ncf['Wasserstein distance'][idxlat, idxlon]
    lsq[j] = np.nanmean(ncf['Wasserstein distance'][idxlat, idxlon])
    ncf.close()
    
    if(CS[j]!=0.25):
        ncf = Dataset(dirRead + 'OTs_sp%d_cs%.1fgm.nc'%(sp1,CS[j]))
    else:
        ncf = Dataset(dirRead + 'OTs_sp%d_cs%.2fgm.nc'%(sp1,CS[j]))
    idxlat = np.logical_and(ncf['Lats']>=maxminlat, ncf['Lats']<=minmaxlat)
    idxlon = np.logical_and(ncf['Lons']>=maxminlon, ncf['Lons']<=minmaxlon)
    if(CS[j]==cs):
        distancegm = np.flip(ncf['Wasserstein distance'][idxlat, idxlon],0)
    lsqgm[j] = np.nanmean(ncf['Wasserstein distance'][idxlat, idxlon])
    ncf.close()

for j in range(len(CS50)):
    if(CS[j]!=0.25):
        ncf = Dataset(dirRead + 'OTs_sp%d_cs%.1f.nc'%(sp2,CS50[j]))
    else:
        ncf = Dataset(dirRead + 'OTs_sp%d_cs%.2f.nc'%(sp2,CS50[j]))
    idxlat = np.logical_and(ncf['Lats']>=maxminlat, ncf['Lats']<=minmaxlat)
    idxlon = np.logical_and(ncf['Lons']>=maxminlon, ncf['Lons']<=minmaxlon)
    lsq50[j] = np.nanmean(ncf['Wasserstein distance'][idxlat, idxlon])
    ncf.close()
    if(CS[j]!=0.25):
        ncf = Dataset(dirRead + 'OTs_sp%d_cs%.1fgm.nc'%(sp2,CS[j]))
    else:
        ncf = Dataset(dirRead + 'OTs_sp%d_cs%.2fgm.nc'%(sp2,CS[j]))
    idxlat = np.logical_and(ncf['Lats']>=maxminlat, ncf['Lats']<=minmaxlat)
    idxlon = np.logical_and(ncf['Lons']>=maxminlon, ncf['Lons']<=minmaxlon)
    lsq50gm[j] = np.nanmean(ncf['Wasserstein distance'][idxlat, idxlon])
    ncf.close()

ncf = Dataset(dirRead + 'OTs_monmean_sp25.nc')
idxlat = np.logical_and(ncf['Lats']>=maxminlat, ncf['Lats']<=minmaxlat)
idxlon = np.logical_and(ncf['Lons']>=maxminlon, ncf['Lons']<=minmaxlon)
distance_mm = np.flip(ncf['Wasserstein distance'][idxlat,idxlon],0)
mean_distance_mm2 = np.nanmean(distance_mm)
ncf = Dataset(dirRead + 'OTs_monmean_sp6.nc')
idxlat = np.logical_and(ncf['Lats']>=maxminlat, ncf['Lats']<=minmaxlat)
idxlon = np.logical_and(ncf['Lons']>=maxminlon, ncf['Lons']<=minmaxlon)
distance_mm = np.flip(ncf['Wasserstein distance'][idxlat, idxlon],0)
#distance_mm = ncf['Wasserstein distance'][idxlat, idxlon]
mean_distance_mm = np.nanmean(distance_mm)

ncf = Dataset(dirRead + 'OTs_nultest.nc')
idxlat = np.logical_and(ncf['Lats']>=maxminlat, ncf['Lats']<=minmaxlat)
idxlon = np.logical_and(ncf['Lons']>=maxminlon, ncf['Lons']<=minmaxlon)
lsqnull = np.nanmean(ncf['Wasserstein distance'][idxlat, idxlon])
#distance_mm = np.flip(ncf['Wasserstein distance'][idxlat, idxlon],0)
ncf.close()

extenew = [maxminlon, minmaxlon, maxminlat, minmaxlat]
Lats = np.arange(maxminlat, minmaxlat+1)
Lons = np.arange(maxminlon, minmaxlon+1)

#%%
fig = plt.figure(figsize=(19,9))
grid = plt.GridSpec(2, 24, wspace=0., hspace=0.35)
#%% First subplot for the highres monthly Wd
ax0 = plt.subplot(grid[0, 12:], projection=projection)#plt.subplot(2,2,2, projection=projection)
for tick in ax0.xaxis.get_major_ticks():
                tick.label.set_fontsize(fs) 
plt.title('(b) Eddying vs. non-eddying ($W_d(R_{0.1},R_{1m})$)', fontsize=fs)
#ax0.add_feature(cartopy.feature.LAND, color='gray')
g = ax0.gridlines( draw_labels=True,crs=ccrs.PlateCarree(central_longitude=180),
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')
g.xlabels_top = False
g.ylabels_right = False
g.ylabels_left = False
g.xlabels_bottom = False
g.xformatter = LONGITUDE_FORMATTER
g.yformatter = LATITUDE_FORMATTER
g.xlocator = mticker.FixedLocator([0, 90, 180, -90])
g.ylocator = mticker.FixedLocator([-75,-50,-25, 0, 25, 50, 75, 100])

g.ylabel_style = {'fontsize': fs}
ax0.set_extent(exte, ccrs.PlateCarree())

#distance_cs0 = distance_cs0[:,:-1]
#plt.pcolormesh(Lons,Lats, distance_cs0, cmap=cmap, zorder = 0,
#                norm=colors.LogNorm(vmin=10000, vmax=5*10**6))
#land = np.full(distance_cs0.shape, np.nan); land[np.isnan(distance_cs0)] = 1;
#plt.pcolormesh(Lons,Lats,land, vmin=0, vmax=2.3, cmap='binary', zorder = 0)
im = plt.imshow(distance_cs0, extent = exte2,
                transform=ccrs.PlateCarree(), cmap=cmap, zorder = 0,
                norm=colors.LogNorm(vmin=10000, vmax=5*10**6))
land = np.full(distance_cs0.shape, np.nan); land[np.isnan(distance_cs0)] = 1;
plt.imshow(land, vmin=0, vmax=2.3, extent = exte, transform=ccrs.PlateCarree(), cmap='binary', zorder = 0)

ax0.set_xticks([0., 90., 180., 270., 360.], crs=ccrs.PlateCarree())
ax0.set_xticklabels([0., 90., 180., 270., 360.])

lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax0.xaxis.set_major_formatter(lon_formatter)
ax0.yaxis.set_major_formatter(lat_formatter)
ax0.grid(linewidth=2, color='black', alpha=0., linestyle='--')
#%% Second subplot for the lowres monthly Wd
ax0 = plt.subplot(grid[0, :12], projection=projection)#plt.subplot(2,2,1, projection=projection)
for tick in ax0.xaxis.get_major_ticks():
                tick.label.set_fontsize(fs)

plt.title('(a) 5-daily vs. monthly model output ($W_d(R_{0.1},R_{0.1m})$)', fontsize=fs)
ax0.add_feature(cartopy.feature.LAND, color='gray')
g = ax0.gridlines(crs=ccrs.PlateCarree(central_longitude=180), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')
g.xlabels_top = False
g.ylabels_right = False
g.xlabels_bottom = False
g.xformatter = LONGITUDE_FORMATTER
g.yformatter = LATITUDE_FORMATTER
#g.xlocator = mticker.FixedLocator([-180,-90, -0, 90, 180])
g.ylocator = mticker.FixedLocator([-50,-25, 0, 25, 50, 75, 100])

g.xlabel_style = {'fontsize': fs}
g.ylabel_style = {'fontsize': fs}
ax0.set_extent(exte, ccrs.PlateCarree())

#distance_mm = np.concatenate((distance_mm[:,180:-2], distance_mm[:,:180]),axis=1);
#distance_mm = np.concatenate((distance_mm, np.zeros((147,2))),axis=1);
#Lons_mm = np.append(np.append(Lons[180:-2], Lons[:180]),np.array([360,361]))
im = plt.imshow(distance_mm, extent = exte2,
                transform=ccrs.PlateCarree(), cmap=cmap, zorder = 0,
                norm=colors.LogNorm(vmin=10000, vmax=5*10**6))
land = np.full(distance_mm.shape, np.nan); land[np.isnan(distance_mm)] = 1;
plt.imshow(land, vmin=0, vmax=2.3, extent = exte2, transform=ccrs.PlateCarree(), cmap='binary', zorder = 0)

#plt.pcolormesh(Lons,Lats, distance_mm, cmap=cmap, zorder = 0,
#                norm=colors.LogNorm(vmin=10000, vmax=5*10**6))
#land = np.full(distance_mm.shape, np.nan); land[np.isnan(distance_mm)] = 1;
#plt.pcolormesh(Lons,Lats, land, vmin=0, vmax=2.3, cmap='binary', zorder = 0)

#ax0.set_xticks([0., 90., 180., 270., 360.], crs=ccrs.PlateCarree())
#ax0.set_xticklabels([0., 90., 180., 270., 360.])

lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax0.xaxis.set_major_formatter(lon_formatter)
ax0.yaxis.set_major_formatter(lat_formatter)
ax0.grid(linewidth=2, color='black', alpha=0., linestyle='--')
#%%
# The fourth subplot, which shows the optimization of Cs
def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

dsWD = [4,7,4,7] # the line dash of the first configuration

ax = plt.subplot(grid[1, 13:])#plt.subplot(2,4,7)
plt.title('(d)', fontsize=fs)

plt.xlabel('$c_s$                   ', fontsize=fs)
plt.ylabel('Mean $W_d$ (10$^{5}$ km)', fontsize=fs, color=color1)

# the wasserstein axis
#ax.set_yticklabels(labels=[-2,-1,0,1,2,3,4,5],fontsize=fs, color='k')

#a0 = sns.lineplot(x=CS,y=lsq/10**5, linewidth = lw, ax=ax, color=color1, 
#             zorder=9)

a0 = sns.lineplot(x=CS,y=np.full(len(lsq), mean_distance_mm/10**5), 
                  linewidth = lw, ax=ax, color=color1, zorder=10)
a0.lines[0].set_linestyle(":")
sns.lineplot(x=CS,y=lsqnull/10**5, linewidth = lw, ax=ax, color=color1, 
             zorder=10)
a0.lines[1].set_dashes(dsWD)
a0 = sns.lineplot(x=CS,y=np.full(len(lsq), mean_distance_mm2/10**5), 
                  linewidth = lw, ax=ax, color=color2, zorder=10)
a0.lines[2].set_linestyle(":")
#sns.scatterplot(x=CS,y=lsq/10**5, ax=ax, color=color1, s=si+20, zorder=9,
#                   legend=False)
sns.lineplot(x=CS,y=lsqgm/10**5, linewidth = lw, ax=ax, color=color1, 
             zorder=9)
sns.scatterplot(x=CS,y=lsqgm/10**5, ax=ax, color=color1, s=si, zorder=10,
                   legend=False)#, marker="^")

#sns.lineplot(x=CS50,y=lsq50/10**5, linewidth = lw, ax=ax, color=color2, 
#             zorder=9)
#sns.scatterplot(x=CS50,y=lsq50/10**5, ax=ax, color=color2, s=si+20, zorder=10,
#                   legend=False)#marker="^")
print('Wds: ', lsq/10**5, lsqnull/10**5, lsq50/10**5)
sns.lineplot(x=CS50,y=lsq50gm/10**5, linewidth = lw, ax=ax, color=color2, 
             zorder=9)
sns.scatterplot(x=CS50,y=lsq50gm/10**5, ax=ax, color=color2, s=si, zorder=10,
                   legend=False)#, marker="^")


for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(fs) 
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(fs) 

lw = 2
colo = 'k'
legend_el = [Line2D([0], [0], dashes=dsWD, color=colo, lw=lw, label='$W_d(R_{0.1}$, $R_{0.1})$'), 
             Line2D([0], [0], linestyle=':', color=colo, lw=lw, label='$W_d(R_{0.1}$, $R_{0.1m})$'), 
             Line2D([0], [0], linestyle='-', marker='o', markersize=8+1, color=colo, lw=lw, label='$W_d(R_{0.1}$, $R_{1m}$/ $R_{1md}$)'), 
#             Line2D([0], [0], linestyle='-', marker='^', markersize=8, color=colo, lw=lw, label='$W_d(R_{0.1}$, $R_{1mb}$/ $R_{1mdb}$)')
             ]
#first_legend = ax.legend(handles=legend_el,loc=4, fontsize=fs, bbox_to_anchor=(0., .15, 1., .102))#, title='Configuration'
first_legend = ax.legend(handles=legend_el, fontsize=fs, loc='upper right', bbox_to_anchor=(1, -0.1))

ax2 = plt.gca().add_artist(first_legend)

legend_el = [Line2D([0], [0], linestyle='solid', color=color1, lw=lw, label='$w_f$ = 6 m day$^{-1}$'), 
             Line2D([0], [0], linestyle='solid', color=color2, lw=lw, label='$w_f$ = 25 m day$^{-1}$')]
#plt.legend(handles=legend_el, title='Sinking speed (m/day)',loc=4, fontsize=fs, bbox_to_anchor=(0., .65, 1., .102))
ax.legend(handles=legend_el, fontsize=fs, loc='upper right', bbox_to_anchor=(0.3, -0.1))

# %%third subplot
ax0 = plt.subplot(grid[1, :12], projection=projection)#plt.subplot(2,2,3, projection=projection)
for tick in ax0.xaxis.get_major_ticks():
                tick.label.set_fontsize(fs)
plt.title('(c) Eddying vs. non-eddying with difussion ($W_d(R_{0.1},R_{1md})$)', fontsize=fs)
g = ax0.gridlines(crs=ccrs.PlateCarree(central_longitude=180), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')
g.xlabels_top = False
g.ylabels_right = False
g.xlabels_bottom = False
g.xformatter = LONGITUDE_FORMATTER
g.yformatter = LATITUDE_FORMATTER
g.xlocator = mticker.FixedLocator([-180,-90, -0, 90, 180])
g.ylocator = mticker.FixedLocator([-75,-50,-25, 0, 25, 50, 75, 100])

g.xlabel_style = {'fontsize': fs}
g.ylabel_style = {'fontsize': fs}
ax0.set_extent(exte, ccrs.PlateCarree())

#distance = distance[:,:-1]
im = plt.imshow(distance, extent = exte2,
                transform=ccrs.PlateCarree(), cmap=cmap, zorder = 0,
                norm=colors.LogNorm(vmin=10000, vmax=5*10**6))
land = np.full(distance.shape, np.nan); land[np.isnan(distance)] = 1;
plt.imshow(land, vmin=0, vmax=2.3, extent = exte2, transform=ccrs.PlateCarree(), cmap='binary', zorder = 0)

ax0.set_xticks([0., 90., 180., 270., 360.], crs=ccrs.PlateCarree())
ax0.set_xticklabels([0., 90., 180., 270., 360.])

lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax0.xaxis.set_major_formatter(lon_formatter)
ax0.yaxis.set_major_formatter(lat_formatter)
ax0.grid(linewidth=2, color='black', alpha=0., linestyle='--')

#%%
fig.subplots_adjust(bottom=0.17)
cbar_ax = fig.add_axes([0.135, 0., 0.35, 0.07])
cbar_ax.set_visible(False)
cbar = fig.colorbar(im, ax=cbar_ax, orientation = 'horizontal', fraction = 1.2,
                    aspect=15)
cbar.ax.xaxis.set_label_position('bottom')
cbar.set_label(label='W$_d$ (km)', fontsize=fs)
cbar.ax.tick_params(labelsize=fs)

plt.savefig('figure2.pdf',bbox_inches='tight',pad_inches=0)
plt.show()