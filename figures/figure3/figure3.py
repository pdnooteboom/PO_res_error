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
#exte = [1, 360, -74, 81]
exte = [1, 360, -75, 72]
exte2 = exte
#exte2 = [-179, 181, -74, 81]
Cs = 5.0
ddeg = 1
cmap2 = 'coolwarm' # For the surface area
cmap3 = 'hot'# For the average travel distance
vssurf = [0,17]
#%%
def distance_nojit(origin, destination):
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


def mean_latitudeweigthed(avgd, lats):
#    result = np.nanmean(avgd)
    weights = np.zeros(avgd.shape)
    for i in range(weights.shape[0]):
        for j in range(weights.shape[1]):
            if(avgd[i,j]>0):
                weights[i,j] = distance_nojit([lats[i], -0.5], [lats[i], 0.5])
    avgd[np.isnan(avgd)] = 0
    result = np.average(avgd, weights=weights)
    return result
#%%
ncf = Dataset('avgd_surf.nc')
Lats = ncf['Lats'][:]
CS = ncf['CS'][:]
idcs = np.where(CS==Cs)[0][0]

surf =  ncf['surftemp_hr'][:]
surf_temp = mean_latitudeweigthed(surf, Lats)/ 10**5.
#surf_temp = np.nanmean(surf) / 10**5.

surf =  ncf['surftemp_hr25'][:]
surf_temp2 = mean_latitudeweigthed(surf, Lats)/ 10**5.
#surf_temp2 = np.nanmean(surf) / 10**5.

surf =  ncf['surf_hr25'][:]
surf50mean = mean_latitudeweigthed(surf, Lats)/ 10**5.
#surf50mean = np.nanmean(surf)  / 10**5.

surf_highres =  ncf['surf_hr'][:]
highres_surf = mean_latitudeweigthed(surf_highres, Lats) / 10**5.
#highres_surf = np.nanmean(surf_highres) / 10**5.

surf_lr =  np.flip(ncf['surf_lr'][0],0)

surf_cs =  np.flip(ncf['surf_lr'][idcs],0)

#sur = np.nanmean(np.nanmean(ncf['surf_lr'][:],axis=1),axis=1)  / 10**5.
#sur50 = np.nanmean(np.nanmean(ncf['surf_lr25'][:],axis=1),axis=1)  / 10**5.
#surgm = np.nanmean(np.nanmean(ncf['surf_gm_lr'][:],axis=1),axis=1)  / 10**5.
#sur50gm = np.nanmean(np.nanmean(ncf['surf_gm_lr25'][:],axis=1),axis=1)  / 10**5.

sur50 = np.zeros(len(CS))
sur = np.zeros(len(CS))
sur50gm = np.zeros(len(CS))
surgm = np.zeros(len(CS))

for i in range(sur.shape[0]):
    sur[i] = mean_latitudeweigthed(ncf['surf_lr'][i], Lats) / 10**5.
    sur50[i] = mean_latitudeweigthed(ncf['surf_lr25'][i], Lats)/ 10**5.
    surgm[i] = mean_latitudeweigthed(ncf['surf_gm_lr'][i], Lats) / 10**5.
    sur50gm[i] = mean_latitudeweigthed(ncf['surf_gm_lr25'][i], Lats) / 10**5.


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
#%% start figure
fig = plt.figure(figsize=(19,9))
grid = plt.GridSpec(2, 24, wspace=0., hspace=0.2)
#% subplot (a)
ax = plt.subplot(grid[0, :12], projection=projection)#plt.subplot(2,2,1, projection=projection)
plt.title('(a) $R_{0.1}$, eddying', fontsize=fs)

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

plt.imshow(np.flip(surf_highres,0)/10.**5, vmin=vssurf[0], vmax=vssurf[1], extent = exte2, transform=ccrs.PlateCarree(), 
           cmap=cmap3, zorder = 0)
land = np.full(surf_highres.shape, np.nan); land[0==surf_highres] = 1;
plt.imshow(np.flip(land,0), vmin=0, vmax=1.6, extent = exte2, transform=ccrs.PlateCarree(), cmap='binary', zorder = 0)


#%subplot (b)
ax = plt.subplot(grid[0, 12:], projection=projection)#plt.subplot(2,2,2, projection=projection)
plt.title('(b) $R_{1m}$, non-eddying', fontsize=fs)

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

plt.imshow(surf_lr/10**5., vmin=vssurf[0], vmax=vssurf[1], extent = exte2, transform=ccrs.PlateCarree(), 
           cmap=cmap3, zorder = 0)
land = np.full(surf_lr.shape, np.nan); land[np.isnan(surf_lr)] = 1;
plt.imshow(land, vmin=0, vmax=1.6, extent = exte2, transform=ccrs.PlateCarree(), cmap='binary', zorder = 0)

#% subplot (c)
ax = plt.subplot(grid[1, :12], projection=projection)#plt.subplot(2,2,3, projection=projection)
plt.title('(c) $R_{1md}$, non-eddying with diffusion', fontsize=fs)

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

im2 = plt.imshow(surf_cs/10**5., vmin=vssurf[0], vmax=vssurf[1], extent = exte2, transform=ccrs.PlateCarree(), 
           cmap=cmap3, zorder = 0)
land = np.full(surf_cs.shape, np.nan); land[np.isnan(surf_cs)] = 1;
plt.imshow(land, vmin=0, vmax=1.6, extent = exte2, transform=ccrs.PlateCarree(), cmap='binary', zorder = 0)
    
#%
dsWD = [4,7,4,7] # the line dash of the first configuration
dsWD2 = [4,7,4,7] # the line dash of the bolus configuration

ax = plt.subplot(grid[1, 13:])#plt.subplot(2,2,3, projection=projection)
plt.title('(d)', fontsize=fs)

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
sns.scatterplot(x=CS,y=surgm, ax=ax, color=color1, s=si, zorder=10,
                   legend=False, marker="^")

sns.lineplot(x=CS, y=sur50, color=color2, linewidth=lw, zorder=10)
sns.scatterplot(x=CS, y=sur50, color=color2, s=si, zorder=11)

sns.lineplot(x=CS,y=sur50gm, linewidth = lw, ax=ax, color=color2, 
             zorder=9)
sns.scatterplot(x=CS,y=sur50gm, ax=ax, color=color2, s=si, zorder=10,
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
             Line2D([0], [0], linestyle='-', marker='o', markersize=8+1, color=colo, lw=lw, label='$R_{1m}$/ $R_{1md}$'), 
#             Line2D([0], [0], linestyle='-', marker='^', markersize=8, color=colo, lw=lw, label='$W_d(R_{0.1}$, $R_{1mb}$/ $R_{1mdb}$)')
             ]
#first_legend = ax.legend(handles=legend_el,loc=4, fontsize=fs, bbox_to_anchor=(0., .15, 1., .102))#, title='Configuration'
first_legend = ax.legend(handles=legend_el, fontsize=fs, loc='upper right', bbox_to_anchor=(1, -0.1))

ax2 = plt.gca().add_artist(first_legend)

legend_el = [Line2D([0], [0], linestyle='solid', color=color1, lw=lw, label='$w_f$ = 6 m day$^{-1}$'), 
             Line2D([0], [0], linestyle='solid', color=color2, lw=lw, label='$w_f$ = 25 m day$^{-1}$')]
#plt.legend(handles=legend_el, title='Sinking speed (m/day)',loc=4, fontsize=fs, bbox_to_anchor=(0., .65, 1., .102))
ax.legend(handles=legend_el,  fontsize=fs, loc='upper right', bbox_to_anchor=(0.3, -0.1))

#% final
#fig.subplots_adjust(bottom=0.17)
#cbar_ax = fig.add_axes([0.11, 0.05, 0.35, 0.07])
#cbar_ax.set_visible(False)
#cbar = fig.colorbar(im2, ax=cbar_ax, orientation = 'horizontal', fraction = 1.2)
#cbar.ax.xaxis.set_label_position('bottom')
#cbar.ax.set_xlabel('$10^5$ km$^2$', fontsize=fs)
#cbar.ax.tick_params(labelsize=fs) 
#cbar.set_ticklabels([1,2,3,4]) 

fig.subplots_adjust(bottom=0.17)
cbar_ax = fig.add_axes([0.135, 0., 0.35, 0.07])
cbar_ax.set_visible(False)
cbar = fig.colorbar(im2, ax=cbar_ax, orientation = 'horizontal', fraction = 1.2,
                    aspect=18)
cbar.ax.xaxis.set_label_position('bottom')
cbar.ax.set_xlabel('$10^5$ km$^2$', fontsize=fs)
cbar.ax.tick_params(labelsize=fs)

plt.savefig('figure3.pdf',bbox_inches='tight',pad_inches=0)
plt.show()