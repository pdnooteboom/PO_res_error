# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 15:30:32 2018

@author: nooteboom
"""

from __future__ import division
import numpy as np
from netCDF4 import Dataset
from numba import jit

#%% Some functions
from datetime import date, datetime, timedelta

Y = 2000 # dummy leap year to allow input X-02-29 (leap day)
seasons = [('winter', (date(Y,  1,  1),  date(Y,  3, 20))),
           ('spring', (date(Y,  3, 21),  date(Y,  6, 20))),
           ('summer', (date(Y,  6, 21),  date(Y,  9, 22))),
           ('autumn', (date(Y,  9, 23),  date(Y, 12, 20))),
           ('winter', (date(Y, 12, 21),  date(Y, 12, 31)))]

def get_season(now):
    if isinstance(now, datetime):
        now = now.date()
    now = now.replace(year=Y)
    return next(season for season, (start, end) in seasons
                if start <= now <= end)

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]
    
def find_nearest_index(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx    
    
def find_down(array,value):
    if(value>=array[0]):
        array = [n-value for n in array]
        idx = np.array([n for n in array if n<0]).argmax()
    else:
        idx = np.nan
    return idx     

#%% To open a file
adv = True#False
nadv = False#True


res = 1
sp = 6.
dd = 10. 
Cs = 3.0

ddeg = 1

dirRead = '/projects/0/palaeo-parcels/POP/POPres/particlefiles/'
dirWrite2 = '/projects/0/palaeo-parcels/POP/POPres/particlefiles/surface/'

#%% Getting the average and variance differences 
 
pltse = 'annual' 
if(pltse=='winter'):
    se = 0
elif(pltse=='spring'):
    se = 1
elif(pltse=='summer'):
    se = 2
elif(pltse=='autumn'):
    se = 3
else:
    se = 4

if(adv):
    lons0 = np.array([])
    lats0 = np.array([])
    lons = np.array([])
    lats = np.array([])
    season = np.array([])
    temp = np.array([])
    salin = np.array([])
    time = np.array([])
    age = np.array([])
    zs = np.array([])

ncn = Dataset(dirRead + 'surface/surfacegrid' +'_id'+'_dd'+str(int(dd))+"_res"+str(res) + '.nc')
#Arrays for fixed surface locations
if(nadv):
    fixlon = np.empty([0,ncn['lon'][:].shape[1]])
    fixlat = np.empty([0,ncn['lon'][:].shape[1]])
    fixtemp = np.empty([0,ncn['lon'][:].shape[1]])
    fixsalin = np.empty([0,ncn['lon'][:].shape[1]])
    fixtime = np.empty([0,ncn['lon'][:].shape[1]])
    print 'shapes: ',fixtime.shape, ncn['lon'][:].shape
for posidx in range(1,19):
    if(posidx%10==0):
        print posidx
    if(adv):
        print dirRead + 'sp%d_dd%d/grid_smagorinski_Cs%.1f_id%d'%(int(sp),int(dd),Cs,posidx) +'_dd'+str(int(dd)) +'_sp'+str(int(sp))+"_res"+str(res) + ".nc"
        nc = Dataset(dirRead + 'sp%d_dd%d/grid_smagorinski_Cs%.1f_id%d'%(int(sp),int(dd),Cs,posidx) +'_dd'+str(int(dd)) +'_sp'+str(int(sp))+"_res"+str(res) + ".nc") 
        if(nc['lon0'][0,0]<0):
            lons0 = np.append(lons0,nc['lon0'][:,0]+360)   
        else:
            lons0 = np.append(lons0,nc['lon0'][:,0])         
#        lons0 = np.concatenate((lons0,nc['lon0'][:,0]),axis=0)
        lats0 = np.concatenate((lats0,nc['lat0'][:,0]), axis=0)
        lons = np.concatenate((lons,nc['lon'][:,0]),axis=0)
        lats = np.concatenate((lats,nc['lat'][:,0]), axis=0) 
        temp = np.concatenate((temp,nc['temp'][:,0]), axis=0) 
        salin = np.concatenate((salin,nc['salin'][:,0]), axis=0)
        time = np.concatenate((time,nc['time'][:,0]), axis=0)
        age = np.concatenate((age,nc['age'][:,0]), axis=0)
        zs = np.concatenate((zs,nc['z'][:,0]), axis=0)
    if(nadv):
        ncn = Dataset(dirRead + 'surface/surfacegrid' +'_id'+'_dd'+str(int(dd))+"_res"+str(res) + '.nc')
        fixlon = np.concatenate((fixlon,ncn['lon'][:]),axis=0)
        fixlat = np.concatenate((fixlat,ncn['lat'][:]), axis=0) 
        fixtemp = np.concatenate((fixtemp,ncn['temp'][:]), axis=0)  
        fixsalin = np.concatenate((fixsalin,ncn['salin'][:]), axis=0) 
        fixtime = np.concatenate((fixtime,ncn['time'][:]), axis=0)      
        
if(adv):
    lons = lons%360
    lons0 = lons0%360

if(nadv):
    fixlon = fixlon%360
    
    fixtemp = np.ma.filled(fixtemp, fill_value=np.nan)
    fixsalin = np.ma.filled(fixsalin, fill_value=np.nan)
    fixtime = np.ma.filled(fixtime, fill_value=np.nan)
    fixlon = np.ma.filled(fixlon, fill_value=np.nan)
    fixlat = np.ma.filled(fixlat, fill_value=np.nan) 

#%% Write two netcdf files where all posidx files are concatenated
#First the advected nc file
if(adv):
    dirWrite = '/projects/0/palaeo-parcels/POP/POPres/particlefiles/sp%d_dd%d/'%(int(sp),int(dd))
    
    #np.savez(dirWrite + 'concatenated_sp%d_dd%d_res%d.npy'%(int(sp),int(dd),res), lats0=lats0, lons0=lons0, lons=lons, lats=lats, temp = temp, salin = salin, time = time, fixlon = fixlon, fixlat = fixlat, fixtemp = fixtemp, fixsalin = fixsalin, fixtime = fixtime)   
    
    dataset = Dataset(dirWrite + 'concatenated_smagorinski_Cs%.1f_sp%d_dd%d_res%d.nc'%(Cs,int(sp),int(dd),res),'w',format='NETCDF4_CLASSIC')
    
    traj = dataset.createDimension('traj',lats.shape[0])
    
    times = dataset.createVariable('time', np.float64, ('traj',))
    lat = dataset.createVariable('lat', np.float64, ('traj',))
    lon = dataset.createVariable('lon', np.float64, ('traj',))
    salins = dataset.createVariable('salin', np.float64, ('traj',)) 
    temps = dataset.createVariable('temp', np.float64, ('traj',))
    ages = dataset.createVariable('age', np.float64, ('traj',))
    lat0 = dataset.createVariable('lat0', np.float64, ('traj',))
    lon0 = dataset.createVariable('lon0', np.float64, ('traj',))
    z = dataset.createVariable('z', np.float64, ('traj',))
    
    lat[:] = lats
    lon[:] = lons
    lon0[:] = lons0
    lat0[:] = lats0
    times[:] = time
    salins[:] = salin
    temps[:] = temp
    ages[:] = age
    z[:] = zs
    
    dataset.close()

# Then the fixed surface nc file
if(nadv):
    dataset = Dataset(dirWrite2 + 'concatenatedsurface_dd%d_res%d.nc'%(int(dd),res),'w',format='NETCDF4_CLASSIC')

    fixtemp = fixtemp[:,1:]
    fixlon = fixlon[:,1:]
    fixlat = fixlat[:,1:]
    fixtime = fixtime[:,1:]
    fixsalin = fixsalin[:,1:]
    
    traj = dataset.createDimension('traj',fixlon.shape[0])
    obs = dataset.createDimension('obs',fixlon.shape[1])
    
    fixlons = dataset.createVariable('lon', np.float64, ('traj','obs',))
    fixlats = dataset.createVariable('lat', np.float64, ('traj','obs',))
    fixtimes = dataset.createVariable('time', np.float64, ('traj','obs',))
    fixtemps = dataset.createVariable('temp', np.float64, ('traj','obs',))
    fixsalins = dataset.createVariable('salin', np.float64, ('traj','obs',))
    
    print 'any nan in fixtemp:', np.isnan(fixtemp).any()
    fixlons[:] = fixlon[:]; del fixlon;
    fixlats[:] = fixlat[:]; fixlat
    fixtimes[:] = fixtime[:]; del fixtime
    fixtemps[:] = fixtemp[:];
    fixsalins[:] = fixsalin; del fixsalin
    
    dataset.close()
