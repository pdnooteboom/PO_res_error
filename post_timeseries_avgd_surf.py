#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 16:14:35 2020

@author: nooteboom
"""

import numpy as np
from netCDF4 import Dataset
import math
from numba import jit

# variables
sp = 6
dd = 10
ddeg = 1

sp1 = 6
sp2 = 25

CS =  np.array([0., 0.25, 0.5, 1.0, 2.0, 5.0])
CS25 = np.array([0., 0.25, 0.5, 1.0, 2.0, 5.0])
#%% functions

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
maxminlat = -76
minmaxlat = 71
minmaxlon = 359 - 180
maxminlon = 0 - 180

readhr = '/projects/0/palaeo-parcels/POP/POPres/0.1degree/particlefiles/'
readlr = '/projects/0/palaeo-parcels/POP/POPres/particlefiles/'

import matplotlib.pylab as plt
def plotim(field):
    plt.imshow(field)
    plt.colorbar()
    plt.show()

avgd, surf, Lons, Lats = calc_fields(name = readhr + 'sp6_dd10/timeseries_per_location_ddeg1_sp6_dd10_tempresmonmean.nc')
idxlat = np.logical_and(Lats>=maxminlat, Lats<=minmaxlat)
idxlon = np.logical_and(Lons>=maxminlon, Lons<=minmaxlon)
assert (idxlon==True).all()
Lonsf = Lons[idxlon]; Latsf = Lats[idxlat];
surf_temp6 =  surf[idxlat]
surf_temp6[surf_temp6==0] = np.nan
avgd_temp6= avgd[idxlat]

avgd, surf, Lons, Lats = calc_fields(name = readhr + 'sp25_dd10/timeseries_per_location_ddeg1_sp25_dd10_tempresmonmean.nc')
idxlat = np.logical_and(Lats>=maxminlat, Lats<=minmaxlat)
idxlon = np.logical_and(Lons>=maxminlon, Lons<=minmaxlon)
assert (idxlon==True).all()
assert (Lons==Lonsf).all()
surf_temp25 = surf[idxlat]
surf_temp25[surf_temp25==0] = np.nan
avgd_temp25 = avgd[idxlat]

avgd, surf, Lons, Lats = calc_fields(name = readhr + 'sp25_dd10/timeseries_per_location_ddeg1_sp25_dd10_tempres5.nc')
idxlat = np.logical_and(Lats>=maxminlat, Lats<=minmaxlat)
idxlon = np.logical_and(Lons>=maxminlon, Lons<=minmaxlon)
assert (idxlon==True).all()
assert (Lons==Lonsf).all()
surf_25hr = surf[idxlat]
surf_25hr[surf_25hr==0] = np.nan
avgd_25hr = avgd[idxlat]

avgd, surf, Lons, Lats = calc_fields(name = readhr + 'sp6_dd10/timeseries_per_location_ddeg%d_sp%d_dd%d_tempres5_ds2.nc'%(ddeg,sp1,dd))
idxlat = np.logical_and(Lats>=maxminlat, Lats<=minmaxlat)
idxlon = np.logical_and(Lons>=maxminlon, Lons<=minmaxlon)
assert (idxlon==True).all()
assert (Lons==Lonsf).all()
surf_6hr = surf[idxlat]
surf_6hr[surf_6hr==0] = np.nan
avgd_6hr = avgd[idxlat]

sur  = np.zeros((len(CS), len(Latsf), len(Lonsf)))
sur25  = np.zeros((len(CS), len(Latsf), len(Lonsf)))
surgm  = np.zeros((len(CS), len(Latsf), len(Lonsf)))
sur25gm  = np.zeros((len(CS), len(Latsf), len(Lonsf)))
av  = np.zeros((len(CS), len(Latsf), len(Lonsf)))
av25  = np.zeros((len(CS), len(Latsf), len(Lonsf)))
avgm  = np.zeros((len(CS), len(Latsf), len(Lonsf)))
av25gm  = np.zeros((len(CS), len(Latsf), len(Lonsf)))

for j in range(len(CS)):
    print(j/len(CS25))
    if(CS[j]!=0.25):
        avgd, surf, Lons, Lats = calc_fields(name = readlr + 'sp6_dd10/timeseries_per_location_smagorinksi_wn_Cs%.1f_ddeg1_sp%d_dd10.nc'%(CS[j],sp1))
    else:
        avgd, surf, Lons, Lats = calc_fields(name = readlr + 'sp6_dd10/timeseries_per_location_smagorinksi_wn_Cs%.2f_ddeg1_sp%d_dd10.nc'%(CS[j],sp1))    
    idxlat = np.logical_and(Lats>=maxminlat, Lats<=minmaxlat)
    idxlon = np.logical_and(Lons>=maxminlon, Lons<=minmaxlon)
    assert (idxlon==True).all()
    assert (Lons==Lonsf).all()
    surf = surf[idxlat]
    surf[surf==0] = np.nan
    sur[j] = surf
    av[j] = avgd[idxlat]
#    if(j==0):
#        plotim(sur[j])
#        plotim(av[j])

    
    if(CS[j]!=0.25):
        avgd, surf, Lons, Lats = calc_fields(name = readlr + 'sp6_dd10/timeseries_per_location_smagorinksi_wn_gm_Cs%.1f_ddeg1_sp%d_dd10.nc'%(CS[j],sp1))
    else:
        avgd, surf, Lons, Lats = calc_fields(name = readlr + 'sp6_dd10/timeseries_per_location_smagorinksi_wn_gm_Cs%.2f_ddeg1_sp%d_dd10.nc'%(CS[j],sp1)) 
    idxlat = np.logical_and(Lats>=maxminlat, Lats<=minmaxlat)
    idxlon = np.logical_and(Lons>=maxminlon, Lons<=minmaxlon)
    assert (idxlon==True).all()
    assert (Lons==Lonsf).all()
    surf = surf[idxlat]
    surf[surf==0] = np.nan
    surgm[j] = surf
    avgm[j] = avgd[idxlat]


for j in range(len(CS25)):
    print(j/len(CS25))
    if(CS25[j] != 0.25):
        avgd, surf, Lons, Lats = calc_fields(name = readlr + 'sp25_dd10/timeseries_per_location_smagorinksi_wn_Cs%.1f_ddeg1_sp%d_dd10.nc'%(CS25[j],sp2))          
    else:   
        avgd, surf, Lons, Lats = calc_fields(name = readlr + 'sp25_dd10/timeseries_per_location_smagorinksi_wn_Cs%.2f_ddeg1_sp%d_dd10.nc'%(CS25[j],sp2))    
    idxlat = np.logical_and(Lats>=maxminlat, Lats<=minmaxlat)
    idxlon = np.logical_and(Lons>=maxminlon, Lons<=minmaxlon)
    assert (idxlon==True).all()
    assert (Lons==Lonsf).all()
    surf = surf[idxlat]
    surf[surf==0] = np.nan
    sur25[j] = surf
    av25[j] = avgd[idxlat]

    if(CS[j]!=0.25):
        avgd, surf, Lons, Lats = calc_fields(name = readlr + 'sp25_dd10/timeseries_per_location_smagorinksi_wn_gm_Cs%.1f_ddeg1_sp%d_dd10.nc'%(CS25[j],sp2))
    else:
        avgd, surf, Lons, Lats = calc_fields(name = readlr + 'sp25_dd10/timeseries_per_location_smagorinksi_wn_gm_Cs%.2f_ddeg1_sp%d_dd10.nc'%(CS25[j],sp2))  
    idxlat = np.logical_and(Lats>=maxminlat, Lats<=minmaxlat)
    idxlon = np.logical_and(Lons>=maxminlon, Lons<=minmaxlon)
    assert (idxlon==True).all()
    assert (Lons==Lonsf).all()
    surf = surf[idxlat]
    surf[surf==0] = np.nan
    sur25gm[j] = surf
    av25gm[j] = avgd[idxlat]

#%% Write netCDF file

ncf = Dataset('avgd_surf.nc', 'w')

ncf.createDimension('lon', len(Lonsf))
ncf.createDimension('lat', len(Latsf))
ncf.createDimension('cs', len(CS))

Lonss = ncf.createVariable('Lons', np.float32,('lon',))
Latss = ncf.createVariable('Lats', np.float32,('lat',))
CSs = ncf.createVariable('CS', np.float32,('cs',))

Lonss[:] = Lonsf
Latss[:] = Latsf
CSs[:] = CS

shr = ncf.createVariable('surf_hr',np.float32,('lat', 'lon',))
ahr = ncf.createVariable('avgd_hr',np.float32,('lat', 'lon',))
shr25 = ncf.createVariable('surf_hr25',np.float32,('lat', 'lon',))
ahr25 = ncf.createVariable('avgd_hr25', np.float32,('lat', 'lon',))
shrt = ncf.createVariable('surftemp_hr', np.float32,('lat', 'lon',))
ahrt = ncf.createVariable('avgdtemp_hr', np.float32,('lat', 'lon',))
shrt25 = ncf.createVariable('surftemp_hr25', np.float32, ('lat', 'lon',))
ahrt25 = ncf.createVariable('avgdtemp_hr25', np.float32,('lat', 'lon',))

slr = ncf.createVariable('surf_lr', np.float32, ('cs','lat', 'lon',))
slr25 = ncf.createVariable('surf_lr25', np.float32,('cs','lat', 'lon',))
slrgm = ncf.createVariable('surf_gm_lr',np.float32, ('cs','lat', 'lon',))
slr25gm =ncf.createVariable('surf_gm_lr25', np.float32, ('cs','lat', 'lon',))

alr = ncf.createVariable('avgd_lr', np.float32, ('cs','lat', 'lon',))
alr25 = ncf.createVariable('avgd_lr25', np.float32, ('cs','lat', 'lon',))
alrgm = ncf.createVariable('avgd_gm_lr', np.float32,('cs','lat', 'lon',))
alr25gm =ncf.createVariable('avgd_gm_lr25', np.float32,('cs','lat', 'lon',))

shrt[:] = surf_temp6
ahrt[:] = avgd_temp6
shrt25[:] = surf_temp25
ahrt25[:] = avgd_temp25
shr[:] = surf_6hr
ahr[:] = avgd_6hr
shr25[:] = surf_25hr
ahr25[:] = avgd_25hr

slr[:] = sur
slr25[:] = sur25
slrgm[:] = surgm
slr25gm[:] = sur25gm

alr[:] = av
alr25[:] = av25
alrgm[:] = avgm
alr25gm[:] = av25gm

ncf.close()