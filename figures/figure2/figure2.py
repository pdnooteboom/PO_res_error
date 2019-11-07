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
exte = [1, 360, -75, 81]
exte2 = exte
extemm = [-179, 181, -76, 86]
cmap = 'viridis'
  
sns.set_style("dark")
sns.set_context("paper")
fs = 14 # fontsize
si = 140
lw = 2 # linewidth

sp1 = 6
sp2 = 25

color1 = 'k'
color2 = 'red'
color3 = 'dodgerblue'
#%% Load the data
CS =  np.array([0.0,0.5, 2.0 ,  3.0,   6.0, 15.0])
CS50 = np.array([0.0, 1.5, 3.0, 4.5])
cs = 3.0

highres_avgd = 5535.848817 / 1000.
highres_surf =  5.550313492752732        

dirRead = '/Volumes/HardDisk/POP/output/OT/' 
lsq = np.zeros(len(CS))
lsq50 = np.zeros(len(CS50))
avd = np.zeros(len(CS))
sur  = np.zeros(len(CS))

for j in range(len(CS)):
    ncf = Dataset(dirRead + 'OTs_sp%d_cs%.1f.nc'%(sp1,CS[j]))
    if(CS[j]==cs):
        distance = np.flip(ncf['Wasserstein distance'][:],0)
    if(CS[j]==0):
        distance_cs0 = np.flip(ncf['Wasserstein distance'][:],0)
    
    lsq[j] = np.nanmean(ncf['Wasserstein distance'][:])
    avd[j] = np.nanmean(ncf['Average horizontal travel distance'][:])  /100.
    surf = ncf['Surface area'][:]
    sur[sur==0] = np.nan
    sur[j] = np.nanmean(surf) / 10**4
    ncf.close()

for j in range(len(CS50)):
    ncf = Dataset(dirRead + 'OTs_sp%d_cs%.1f.nc'%(sp2,CS50[j]))
    lsq50[j] = np.nanmean(ncf['Wasserstein distance'][:])
    ncf.close()

ncf = Dataset(dirRead + 'OTs_nultest.nc')
lsqnull = np.nanmean(ncf['Wasserstein distance'][:])
ncf.close()
#%%
fig = plt.figure(figsize=(18,8))
#%% First subplot for the highres monthly Wd
ax0 = plt.subplot(2,2,2, projection=projection)
for tick in ax0.xaxis.get_major_ticks():
                tick.label.set_fontsize(fs) 
plt.title('(b) 1$^{\circ}$, monthly, cs=0', fontsize=fs)
ax0.add_feature(cartopy.feature.LAND, color='gray')
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

im = plt.imshow(distance_cs0, extent = exte2,
                transform=ccrs.PlateCarree(), cmap=cmap, zorder = 0,
                norm=colors.LogNorm(vmin=10000, vmax=5*10**6))
land = np.full(distance_cs0.shape, np.nan); land[np.isnan(distance_cs0)] = 1;
plt.imshow(land, vmin=0, vmax=2.3, extent = exte2, transform=ccrs.PlateCarree(), cmap='binary', zorder = 0)

ax0.set_xticks([0., 90., 180., 270., 360.], crs=ccrs.PlateCarree())
ax0.set_xticklabels([0., 90., 180., 270., 360.])

lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax0.xaxis.set_major_formatter(lon_formatter)
ax0.yaxis.set_major_formatter(lat_formatter)
ax0.grid(linewidth=2, color='black', alpha=0., linestyle='--')

#%% Second subplot for the lowres monthly Wd
ncf_mm = Dataset(dirRead + 'OTs_monmean4.nc')
distance_mm = np.flip(ncf_mm['Wasserstein distance'][:],0)   
mean_distance_mm = np.nanmean(distance_mm)

ax0 = plt.subplot(2,2,1, projection=projection)
for tick in ax0.xaxis.get_major_ticks():
                tick.label.set_fontsize(fs) 

plt.title('(a) 0.1$^{\circ}$, monthly, cs=0', fontsize=fs)
ax0.add_feature(cartopy.feature.LAND, color='gray')
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
ax0.set_extent(extemm, ccrs.PlateCarree())

distance_mm = np.concatenate((distance_mm[:,180:-2], distance_mm[:,:180]),axis=1);

im = plt.imshow(distance_mm, extent = extemm,
                transform=ccrs.PlateCarree(), cmap=cmap, zorder = 0,
                norm=colors.LogNorm(vmin=10000, vmax=5*10**6))
land = np.full(distance_mm.shape, np.nan); land[np.isnan(distance_mm)] = 1;
plt.imshow(land, vmin=0, vmax=2.3, extent = extemm, transform=ccrs.PlateCarree(), cmap='binary', zorder = 0)

ax0.set_xticks([0., 90., 180., 270., 360.], crs=ccrs.PlateCarree())
ax0.set_xticklabels([0., 90., 180., 270., 360.])

lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax0.xaxis.set_major_formatter(lon_formatter)
ax0.yaxis.set_major_formatter(lat_formatter)
ax0.grid(linewidth=2, color='black', alpha=0., linestyle='--')

#%%
# The third subplot, which is the global Wd with optimized Cs
# The fourth subplot, which shows the optimization of Cs
def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)


fig.subplots_adjust(right=0.85)

dsWD = [4,7,4,7] # the line dash of the first configuration

ax3 = plt.subplot(2,2,4)
ax2 = ax3.twinx()
ax = ax3.twinx()
plt.title('(d)', fontsize=fs)
for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(fs) 
for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(fs) 

plt.xlabel('$c_s$', fontsize=fs)
plt.ylabel('Mean $W_d$ (10$^{5}$ km)', fontsize=fs, color=color1)
ax.set_yticks([0,1.,2.,3.,4,5,6])


# the travel distance axis
ax2.spines["right"].set_position(("axes", 1.1))
a0 = sns.lineplot(x=CS/3., y=np.full(len(CS),highres_avgd), linewidth=lw,
                  ax=ax2, color=color2, zorder=1)
a0.lines[0].set_dashes(dsWD)

a0 = sns.lineplot(x=CS/3., y=np.full(len(CS),np.nanmean(ncf_mm['Average horizontal travel distance'][:]/1000.)), linewidth=lw,
                  ax=ax2, color=color2, zorder=1)
a0.lines[1].set_linestyle(":")


sns.lineplot(x=CS/3., y=avd, ax=ax2, color=color2, linewidth=lw, zorder=10)
sns.scatterplot(x=CS/3., y=avd, ax=ax2, color=color2, s=si, zorder=11)
make_patch_spines_invisible(ax2)
ax2.spines["right"].set_visible(True)
ax2.set_ylabel('average travel distance (10$^2$ km)', fontsize=fs, color = color2)
for tick in ax2.yaxis.get_major_ticks():
                tick.label.set_fontsize(fs)
ax2.set_yticklabels([4,4.5,5,5.5,6], fontsize=fs, color=color2)
ax2.set_ylim(4.,6.4)


print('average distances:  ', highres_avgd,'   ', np.nanmean(ncf_mm['Average horizontal travel distance'][:]/1000.), '   ', avd)



# the surface area axis
ax3.set_ylabel('surface area (10$^5$ km$^2$)', fontsize=fs, color = color3)
ax3.set_yticklabels(labels=[1,2,3,4,5,6],fontsize=fs, color=color3)#

a0 = sns.lineplot(x=CS/3., y=np.full(len(CS),highres_surf), linewidth=lw, ax=ax3, 
                  color=color3, label='1', legend=False, zorder=10)
a0.lines[0].set_dashes(dsWD)
a0 = sns.lineplot(x=CS/3., y=np.full(len(CS),4.845305778217904), ax=ax3, linewidth=lw,
                  color=color3, label='2', 
                  legend=False, zorder=5)
a0.lines[1].set_linestyle(":")


sns.scatterplot(x=CS/3., y=sur, ax=ax3, color=color3, s=si, zorder=11)
sns.lineplot(x=CS/3., y=sur, ax=ax3, color=color3, linewidth=lw, markers = True,
             label='3/4', zorder=10)

print('surface areas:  ', highres_surf,'   ', np.nanmean(ncf_mm['Surface area'][:]/10**5), '   ', sur)

lw = 2
colo = 'k'
legend_el = [Line2D([0], [0], dashes=dsWD, color=colo, lw=lw, label='1'), 
             Line2D([0], [0], linestyle=':', color=colo, lw=lw, label='2'), 
             Line2D([0], [0], color=colo, lw=lw, label='3/4')]
a0.legend(handles=legend_el, title='Configuration',loc=4, fontsize=fs, bbox_to_anchor=(0., .1, 1., .102))



# the wasserstein axis
ax.set_yticklabels(labels=[0,1,2,3,4,5],fontsize=fs, color='k')

sns.lineplot(x=CS/3.,y=lsq/10**5, linewidth = lw, ax=ax, color=color1, 
             zorder=10)

a0 = sns.lineplot(x=CS/3.,y=np.full(len(lsq), mean_distance_mm/10**5), 
                  linewidth = lw, ax=ax, color=color1, zorder=10)
a0.lines[1].set_linestyle(":")
sns.lineplot(x=CS/3.,y=lsqnull/10**5, linewidth = lw, ax=ax, color=color1, 
             zorder=10)
a0.lines[2].set_dashes(dsWD)
sns.scatterplot(x=CS/3.,y=lsq/10**5, ax=ax, color=color1, s=si, zorder=11,
                   legend=False)
sns.lineplot(x=CS50/3.,y=lsq50/10**5, linewidth = lw, ax=ax, color=color1, 
             zorder=11)
sns.scatterplot(x=CS50/3.,y=lsq50/10**5, ax=ax, color=color1, s=si, zorder=10,
                   legend=False, marker="^")
print('Wds: ', lsq/10**5, lsqnull/10**5, lsq50/10**5)

# %%third subplot
ax0 = plt.subplot(2,2,3, projection=projection)
for tick in ax0.xaxis.get_major_ticks():
                tick.label.set_fontsize(fs)
plt.title('(c) 1$^{\circ}$, monthly, $c_s$=%.1f'%(cs/3.), fontsize=fs)
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

fig.subplots_adjust(bottom=0.17)
cbar_ax = fig.add_axes([0.11, 0.0, 0.35, 0.07])
cbar_ax.set_visible(False)
cbar = fig.colorbar(im, ax=cbar_ax, orientation = 'horizontal', fraction = 1.2)
cbar.ax.xaxis.set_label_position('bottom')
cbar.set_label(label='W$_d$ (km)', fontsize=fs)
cbar.ax.tick_params(labelsize=fs)

plt.savefig('figure2.pdf',bbox_inches='tight',pad_inches=0)
plt.show()