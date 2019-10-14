# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 16:27:01 2018

@author: nooteboom
"""
import numpy as np
from numba import jit
from netCDF4 import Dataset

def find_down(array,value):
    if(value>=array[0]):
        array = [n-value for n in array]
        idx = np.array([n for n in array if n<0]).argmax()
    else:
        idx = np.nan
    return idx  

#%% To open a file
res = 1
sp = 6.
dd = 10. 
tempres = 'monmean'#5

ddeg = 1

dirRead = '/projects/0/palaeo-parcels/POP/POPres/0.1degree/particlefiles/'

dirWrite = '/projects/0/palaeo-parcels/POP/POPres/0.1degree/particlefiles/sp%d_dd%d/'%(int(sp),int(dd))

adv = True
nadv = False#True

#%%
pfile = Dataset(dirRead + 'sp%d_dd%d/'%(int(sp),int(dd)) + 'concatenated_sp%d_dd%d_res%d_tempres'%(int(sp),int(dd),res)+str(tempres)+'.nc')

if(nadv):
    pfilefix = Dataset(dirRead + 'surface/' + 'concatenatedsurface_dd10_res1.nc')

#%%
minlon = max(0, min(pfile['lon0'][:]))
maxlon = max(pfile['lon0'][:])+1
minlat = min(pfile['lat0'][:])
maxlat = max(pfile['lat0'][:])+1
#%% Now loop over all regions to define the T matrix
    
#Define the grid at the bottom:
Lons = np.arange(minlon,maxlon+ddeg,ddeg)-0.5
Lats = np.arange(minlat,maxlat+ddeg,ddeg)-0.5
Lonss, Latss = np.meshgrid(Lons, Lats)

vLons = Lonss.reshape(Lonss.size) 
vLats = Latss.reshape(Latss.size)   

@jit
def makets_ex():
    ts = []
    for i in range(len(vLons)): ts.append([]);
    return ts
    
ts = makets_ex()

def construct_ts(tstemp, tssalin, tslon, tslat, tslon0, tslat0, tsage, lats0, lons0,  vLons, vLats, Lons, Lats, temp, salin, lonadv, latadv, ageadv):

    """Return a list with first dimension len(VLons) for all the bottom grid boxes
    and second dimension the temperatures at the surface related to the bottom box"""  

    le = len(lons0)

    tslens = np.zeros(len(tstemp))
    
    maxlents = 0    
    #Look in which temperature region the particle surfaces and grid box region
    # released
    for i in range(le):
        if(lons0[i]>0):        
            lo = find_down(Lons,lons0[i]); lon = Lons[lo];
            la = find_down(Lats,lats0[i]); lat = Lats[la];
       
            j = np.where(np.logical_and(vLons==lon,vLats==lat))[0][0]  
        
            tssalin[j] = np.append(tssalin[j], salin[i])
            tstemp[j] = np.append(tstemp[j], temp[i])   
            tslon[j] = np.append(tslon[j], lonadv[i])
            tslat[j] = np.append(tslat[j], latadv[i])  
            tslon0[j] = np.append(tslon0[j], lons0[i])
            tslat0[j] = np.append(tslat0[j], lats0[i])
            tsage[j] = np.append(tsage[j], ageadv[i])   
            tslens[j] += 1

    for j in range(len(vLons)):
        maxlents = max(maxlents,len(tstemp[j]))                    
                    
    return tstemp, tssalin, tslon, tslat, tslon0, tslat0, tsage, maxlents, tslens

def construct_tsfix(tsfixtemp, tsfixsalin, fixlon, fixlat,  vLons, vLats, Lons, Lats, fixtemp, fixsalin):
    """Return a list with first dimension len(VLons) for all the bottom grid boxes
    and second dimension the temperatures at the surface related to the bottom box"""  
    sle = len(fixlon)
    
    maxlentsfix = 0
    
    for i in range(sle):
        lo = find_down(Lons,fixlon[i]);lon = Lons[lo];
        la = find_down(Lats,fixlat[i]); lat = Lats[la];
        
        j = np.where(np.logical_and(vLons==lon,vLats==lat))[0][0]

        tsfixtemp[j] = np.append(tsfixtemp[j],fixtemp[i,:])     
        tsfixsalin[j] = np.append(tsfixsalin[j],fixsalin[i,:])
        
    for i in range(len(vLons)):
        maxlentsfix = max(maxlentsfix,len(tsfixtemp[i]))  
                    
    return tsfixtemp, tsfixsalin, maxlentsfix


lat0 = np.array(pfile['lat0'][:])
lon0 = np.array(pfile['lon0'][:])
lon0[np.where(lon0<0)] += 360
latadv = np.array(pfile['lat'][:])
lonadv = np.array(pfile['lon'][:])
lonadv[lonadv<0] += 360
ageadv = np.array(pfile['age'][:])
if(nadv):
    lonfix = pfilefix['lon'][:,0]
    lonfix[lonfix<0] += 360
    latfix = pfilefix['lat'][:,0]
temp = np.array(pfile['temp'][:])
salin = np.array(pfile['salin'][:])
if(nadv):
    salinfix = pfilefix['salin'][:,:]
    tempfix = pfilefix['temp'][:,:]

tstemp = np.empty((len(vLons),), dtype=object)  
for i in range(len(tstemp)):tstemp[i] = [];
tssalin = np.empty((len(vLons),), dtype=object)  
for i in range(len(tssalin)):tssalin[i] = [];
tslon = np.empty((len(vLons),), dtype=object)  
for i in range(len(tslon)):tslon[i] = [];
tslat = np.empty((len(vLons),), dtype=object)  
for i in range(len(tslat)):tslat[i] = [];
tslon0 = np.empty((len(vLons),), dtype=object)
for i in range(len(tslon0)):tslon0[i] = [];
tslat0 = np.empty((len(vLons),), dtype=object)
for i in range(len(tslon0)):tslon0[i] = [];
tsage = np.empty((len(vLons),), dtype=object)  
for i in range(len(tsage)):tsage[i] = [];

tstemp, tssalin, tslon, tslat, tslon0, tslat0, tsage, maxlents, tslens =  construct_ts(tstemp, tssalin, tslon, tslat, tslon0, tslat0, tsage, lat0, lon0, vLons, vLats, Lons, Lats, temp, salin, lonadv, latadv, ageadv)

ts = np.array(ts)
if(nadv):   
    tsfixtemp = np.empty((len(vLons),), dtype=object)  
    for i in range(len(tsfixtemp)):tsfixtemp[i] = [];
    tsfixsalin = np.empty((len(vLons),), dtype=object)  
    for i in range(len(tsfixsalin)):tsfixsalin[i] = [];
    tsfixtemp, tsfixsalin, maxlentsfix =  construct_tsfix(tsfixtemp, tsfixsalin, lonfix,latfix,  vLons, vLats, Lons, Lats, tempfix, salinfix)

maxlents = 0

for i in range(len(tstemp)):

    if(len(tstemp[i])>0): 
        maxlents = max(maxlents,len(tstemp[i]))
  
dataset = Dataset(dirWrite  + 'timeseries_per_location_ddeg%d_sp%d_dd%d'%(ddeg, int(sp),int(dd))+'_tempres'+str(tempres)+'.nc','w',format='NETCDF4_CLASSIC')

traj = dataset.createDimension('tslen', maxlents)
if(nadv):
    trajfix = dataset.createDimension('tslenfix', maxlentsfix)
vlatlon = dataset.createDimension('vLons', len(vLons))
lons = dataset.createDimension('Lons', len(Lons))
lats = dataset.createDimension('Lats', len(Lats))

vLon = dataset.createVariable('vLons', np.float64, ('vLons',))
vLat = dataset.createVariable('vLats', np.float64, ('vLons',))
Lon = dataset.createVariable('Lons', np.float64, ('Lons',))
Lat = dataset.createVariable('Lats', np.float64, ('Lats',))

temps = dataset.createVariable('temp', np.float64, ('vLons','tslen',))
salins = dataset.createVariable('salin', np.float64, ('vLons','tslen',))
lat = dataset.createVariable('lat', np.float64, ('vLons','tslen',))
lon = dataset.createVariable('lon', np.float64, ('vLons','tslen',))
ages = dataset.createVariable('age', np.float64, ('vLons','tslen',))

tslenss = dataset.createVariable('tslens', np.float64, ('vLons',))
tslenss[:] = tslens
if(nadv):
    fixtemps = dataset.createVariable('fixtemp', np.float64, ('vLons','tslenfix',))
    fixsalins = dataset.createVariable('fixsalin', np.float64, ('vLons','tslenfix',))

vLon[:] = vLons
vLat[:] = vLats
Lon[:] = Lons
Lat[:] = Lats

for j in range(len(vLons)):
    le = len(tstemp[j])
    if(le>0):
        temps[j,:le] = np.array(tstemp[j])
        salins[j,:le] = np.array(tssalin[j])
        lon[j,:le] = np.array(tslon[j])
        lat[j,:le] = np.array(tslat[j])
        ages[j,:le] = np.array(tsage[j])

    if(nadv):
        fle = len(tsfixtemp[j])
        if(fle>0):
            fixtemps[j,:fle] = np.array(tsfixtemp[j])
            fixsalins[j,:fle] = np.array(tsfixsalin[j])

dataset.close()
