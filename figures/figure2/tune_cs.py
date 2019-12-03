#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 11:35:51 2019

@author: nooteboom
"""

import numpy as np
from netCDF4 import Dataset
import calc_metrics2 as cm

sp = 6#50#    # sinking speed
ddeg = 1
dd = 10
CS = [0.0, 3.0, 4.5] # cs values where we calculate all measures from

minmaxlat = 100000.
maxminlat = -100000
allresults = ['monmean', 'null', '0.0','0.25', '0.5', '1.0', '1.5', '2.0', '5.0']

dirWrite = '/Volumes/HardDisk/POP/output/OT/'  #'/projects/0/palaeo-parcels/POP/OT/'
dirReadhigh = '/Volumes/HardDisk/POP/output/highres/timeseries/'   
dirReadlow = '/Volumes/HardDisk/POP/output/lowres/timeseries/'  
nc_target = Dataset(dirReadhigh + 'timeseries_per_location_ddeg%d_sp%d_dd%d_tempres5.nc'%(ddeg,sp,dd))

tslenmin = 10000
hid = 140 # The 'maximum particle cloud size'

nc = Dataset(dirReadhigh + 'timeseries_per_location_ddeg%d_sp%d_dd%d_tempres5.nc'%(ddeg,sp,dd))
minmaxlat = min(minmaxlat,np.max(nc['Lats'][:]))
maxminlat = max(maxminlat, np.min(nc['Lats'][:]))
tslen = nc['tslens'][:]
print('minlat, maxlat  '+str(np.min(nc['Lats'][:])) + '   ' + str(np.max(nc['Lats'][:])))
print(np.sum(tslen>hid)/np.float(len(tslen[tslen>0])))
tslenmin = min(tslenmin, np.min(tslen[tslen>hid]))
for i in allresults:
    print(i)  # Check if all particle clouds are of the 'maximum particle cloud size'
    if('monmean'==i):
        nc = Dataset(dirReadhigh + 'timeseries_per_location_ddeg1_sp6_dd10_tempresmonmean.nc')
        minmaxlat = min(minmaxlat,np.max(nc['Lats'][:]))
        maxminlat = max(maxminlat, np.min(nc['Lats'][:]))
        tslen = nc['tslens'][:]
        print(np.sum(tslen>hid)/np.float(len(tslen[tslen>0])))
        tslenmin = min(tslenmin, np.min(tslen[tslen>hid]))
    elif('null' == i):
        nc = Dataset(dirReadhigh + 'timeseries_per_location_ddeg%d_sp%d_dd%d_tempres5_ds.nc'%(ddeg,sp,dd))
        minmaxlat = min(minmaxlat,np.max(nc['Lats'][:]))
        maxminlat = max(maxminlat, np.min(nc['Lats'][:]))
        tslen = nc['tslens'][:]
        print(np.sum(tslen>hid)/np.float(len(tslen[tslen>0])))
        tslenmin = min(tslenmin, np.min(tslen[tslen>hid]))
    else:
        if(i=='0.0'):
            nc = nc_target = Dataset(dirReadlow + 'timeseries_per_location_smagorinksi_Cs%s_ddeg%d_sp%d_dd%d.nc'%(i,ddeg,sp,dd))
        else:
            nc = nc_target = Dataset(dirReadlow + 'timeseries_per_location_smagorinksi_wn_Cs%s_ddeg%d_sp%d_dd%d.nc'%(i,ddeg,sp,dd))
 #       nc = Dataset(dirReadlow + 'timeseries_per_location_smagorinksi_Cs%s_ddeg%d_sp%d_dd%d.nc'%(i,ddeg,sp,dd))
        minmaxlat = min(minmaxlat,np.max(nc['Lats'][:]))
        maxminlat = max(maxminlat, np.min(nc['Lats'][:]))
        vLons = nc['vLons'][:]
        vLats = nc['vLats'][:]
        tslen = nc['tslens'][:]
        print(np.sum(tslen>hid)/np.float(len(tslen[tslen>0])))
        tslenmin = min(tslenmin, np.min(tslen[tslen>hid]))
    nc.close()
    print(tslenmin)
Lons = np.arange(360)
Lats = np.arange(maxminlat, minmaxlat+1)

nc_target = Dataset(dirReadhigh + 'timeseries_per_location_ddeg%d_sp%d_dd%d_tempres5.nc'%(ddeg,sp,dd))
idxt = np.where(np.logical_and(nc_target['vLats'][:]>=Lats[0],nc_target['vLats'][:]<=Lats[-1]))
lat = nc_target['lat'][:]; lon = nc_target['lon'][:] 
vLons_target =  nc_target['vLons'][:]; vLats_target =  nc_target['vLats'][:]; 

distance = np.full(vLons.shape, np.nan)
lsq = np.array([])

avgd, surf = cm.calc_fields(Lats, Lons, vLats, vLons, name = dirReadhigh + 'timeseries_per_location_ddeg%d_sp%d_dd%d_tempres5.nc'%(ddeg,sp,dd))

print('average travelled distance in the high res simulation: %f'%np.nanmean(avgd))
print('average surface area in the high res simulation: %f'%np.nanmean(surf))

j=0
for i in allresults:
    if('monmean'==i): 
        print('monmean')
        ds = Dataset(dirWrite + 'OTs_monmean4.nc', 'w')
        lats = ds.createDimension('lat',len(Lats))
        latitudes = ds.createVariable('Lats', np.float32, ('lat',))
        latitudes[:] = Lats
        lons = ds.createDimension('lon',len(Lons))
        longitudes = ds.createVariable('Lons', np.float32, ('lon',))
        longitudes[:] = Lons
        #cs = ds.createDimension('cs',len(CS))
        OTS = ds.createVariable('Wasserstein distance', np.float32, ('lat','lon',))
        
        AVGD = ds.createVariable('Average horizontal travel distance', np.float32, ('lat','lon',))
        SURF = ds.createVariable('Surface area', np.float32, ('lat','lon',))
    
        nc_source = Dataset(dirReadhigh + 'timeseries_per_location_ddeg1_sp6_dd10_tempresmonmean.nc')
        lat_s = nc_source['lat'][:]
        lon_s = nc_source['lon'][:]
        vLats_s = nc_source['vLats'][:]
        vLons_s = nc_source['vLons'][:]
    
        for j in range(Lons.shape[0]-1):
            if(j%10==0): print(j/float(len(Lons)));
            for k in range(Lats.shape[0]):
                idx0 = np.where(np.logical_and(vLons%360==Lons[j]%360,vLats==Lats[k]))
                idxt = np.where(np.logical_and(vLons_target%360==Lons[j]%360,vLats_target==Lats[k]))
                idxs = np.where(np.logical_and(vLons_s%360==Lons[j]%360,vLats_s==Lats[k]))
                latu = lat[idxt[0],:].flatten(); latu = latu[np.where(latu>-500)];
                lasu = lat_s[idxs[0],:].flatten(); lasu = lasu[np.where(lasu>-500)];
                lotu = lon[idxt[0],:].flatten(); lotu = lotu[np.where(lotu>-500)];
                losu = lon_s[idxs[0],:].flatten(); losu = losu[np.where(losu>-500)];
                if(len(lasu)>0 and len(latu>0) and len(losu)>0 and len(lotu)>0):
                    if(not (lasu.mask.all() or latu.mask.all())):
                        lasu = lasu[~lasu.mask][:int(tslenmin)];losu = losu[~losu.mask][:int(tslenmin)];
                        latu = latu[~latu.mask][:int(tslenmin)];lotu = lotu[~lotu.mask][:int(tslenmin)];
    
                        yt, xt, ys, xs = cm.transform_latlon_to_km_frame(latu, lotu, lasu, losu)
            
                        t = np.concatenate((xt[:,np.newaxis], yt[:,np.newaxis]),axis=1)
                        s = np.concatenate((xs[:,np.newaxis], ys[:,np.newaxis]),axis=1)
    
                        
                        distance[idx0] = cm.Wdist(t,s)
        OTS[:] = distance.reshape(len(Lats), len(Lons))
        
        avgdt, surft = cm.calc_fields(Lats, Lons, vLats, vLons, name = dirReadhigh + 'timeseries_per_location_ddeg1_sp6_dd10_tempresmonmean.nc',ml=tslenmin)
    
        AVGD[:] = avgdt
        SURF[:] = surft
        ds.close()    
    elif('null'==i): 
        print('nul test')
        ds = Dataset(dirWrite + 'OTs_nultest.nc', 'w')
        lats = ds.createDimension('lat',len(Lats))
        latitudes = ds.createVariable('Lats', np.float32, ('lat',))
        latitudes[:] = Lats
        lons = ds.createDimension('lon',len(Lons))
        longitudes = ds.createVariable('Lons', np.float32, ('lon',))
        longitudes[:] = Lons
        
        OTS = ds.createVariable('Wasserstein distance', np.float32, ('lat','lon',))
        AVGD = ds.createVariable('Average horizontal travel distance', np.float32, ('lat','lon',))
        SURF = ds.createVariable('Surface area', np.float32, ('lat','lon',))
    
        nc_source  = Dataset(dirReadhigh + 'timeseries_per_location_ddeg%d_sp%d_dd%d_tempres5_ds.nc'%(ddeg,sp,dd))
        lat_s = nc_source['lat'][:]
        lon_s = nc_source['lon'][:]
        vLats_s = nc_source['vLats'][:]
        vLons_s = nc_source['vLons'][:]
        
    
        for i in range(Lons.shape[0]-1):
            if(i%10==0): print(i/float(len(Lons)));
            for k in range(Lats.shape[0]):
                idx0 = np.where(np.logical_and(vLons%360==Lons[j]%360,vLats==Lats[k]))
                idxt = np.where(np.logical_and(vLons_target%360==Lons[i]%360,vLats_target==Lats[k]))
                idxs = np.where(np.logical_and(vLons_s%360==Lons[i]%360,vLats_s==Lats[k]))
                latu = lat[idxt[0],:].flatten(); latu = latu[np.where(latu>-500)];
                lasu = lat_s[idxs[0],:].flatten(); lasu = lasu[np.where(lasu>-500)];
                lotu = lon[idxt[0],:].flatten(); lotu = lotu[np.where(lotu>-500)];
                losu = lon_s[idxs[0],:].flatten(); losu = losu[np.where(losu>-500)];
                if(len(lasu)>0 and len(latu>0) and len(losu)>0 and len(lotu)>0):
                    if(not (lasu.mask.all() or latu.mask.all())):
                        lasu = lasu[~lasu.mask][:int(tslenmin)];losu = losu[~losu.mask][:int(tslenmin)];
                        latu = latu[~latu.mask][:int(tslenmin)];lotu = lotu[~lotu.mask][:int(tslenmin)];
    
                        yt, xt, ys, xs = cm.transform_latlon_to_km_frame(latu, lotu, lasu, losu)
            
                        t = np.concatenate((xt[:,np.newaxis], yt[:,np.newaxis]),axis=1)
                        s = np.concatenate((xs[:,np.newaxis], ys[:,np.newaxis]),axis=1)

                        distance[idx0] = cm.Wdist(t,s)
                        
        OTS[:] = distance.reshape(len(Lats), len(Lons))
        
        avgdt, surft = cm.calc_fields(Lats, Lons, vLats, vLons, name = dirReadhigh + 'timeseries_per_location_ddeg%d_sp%d_dd%d_tempres5_ds.nc'%(ddeg,sp,dd), ml=tslenmin)
    
        AVGD[:] = avgdt
        SURF[:] = surft
        ds.close()    
    else:
        print('CS: %s'%i)
        ds = Dataset(dirWrite + 'OTs_sp%d_cs%s.nc'%(sp,i), 'w')
        lats = ds.createDimension('lat',len(Lats))
        latitudes = ds.createVariable('Lats', np.float32, ('lat',))
        latitudes[:] = Lats
        lons = ds.createDimension('lon',len(Lons))
        longitudes = ds.createVariable('Lons', np.float32, ('lon',))
        longitudes[:] = Lons
        OTS = ds.createVariable('Wasserstein distance', np.float32, ('lat','lon',))
        
        AVGD = ds.createVariable('Average horizontal travel distance', np.float32, ('lat','lon',))
        SURF = ds.createVariable('Surface area', np.float32, ('lat','lon',))
        
        cs = np.float(i)
        if(cs==0.):
            nc_source = nc_target = Dataset(dirReadlow + 'timeseries_per_location_smagorinksi_Cs%.1f_ddeg%d_sp%d_dd%d.nc'%(cs,ddeg,sp,dd))
        elif(cs!=0.25):
            nc_source = nc_target = Dataset(dirReadlow + 'timeseries_per_location_smagorinksi_wn_Cs%.1f_ddeg%d_sp%d_dd%d.nc'%(cs,ddeg,sp,dd))
        else:
            nc_source = nc_target = Dataset(dirReadlow + 'timeseries_per_location_smagorinksi_wn_Cs%.2f_ddeg%d_sp%d_dd%d.nc'%(cs,ddeg,sp,dd))
#        nc_source = nc_target = Dataset(dirReadlow + 'timeseries_per_location_smagorinksi_Cs%s_ddeg%d_sp%d_dd%d.nc'%(i,ddeg,sp,dd))
        lat_s = nc_source['lat'][:]
        lon_s = nc_source['lon'][:]
        vLats_s = nc_source['vLats'][:]
        vLons_s = nc_source['vLons'][:]
    
        for j in range(Lons.shape[0]):
            if(j%10==0): print(j/float(len(Lons)));
            for k in range(Lats.shape[0]):
                idx0 = np.where(np.logical_and(vLons%360==Lons[j]%360,vLats==Lats[k]))
                idxt = np.where(np.logical_and(vLons_target%360==Lons[j]%360,vLats_target==Lats[k]))
                idxs = np.where(np.logical_and(vLons_s%360==Lons[j]%360,vLats_s==Lats[k]))
                latu = lat[idxt[0],:].flatten(); latu = latu[np.where(latu>-500)];
                lasu = lat_s[idxs[0],:].flatten(); lasu = lasu[np.where(lasu>-500)];
                lotu = lon[idxt[0],:].flatten(); lotu = lotu[np.where(lotu>-500)];
                losu = lon_s[idxs[0],:].flatten(); losu = losu[np.where(losu>-500)];
                if(len(lasu)>0 and len(latu>0) and len(losu)>0 and len(lotu)>0):
                    if(not (lasu.mask.all() or latu.mask.all())):
                        lasu = lasu[~lasu.mask][:int(tslenmin)];losu = losu[~losu.mask][:int(tslenmin)];
                        latu = latu[~latu.mask][:int(tslenmin)];lotu = lotu[~lotu.mask][:int(tslenmin)];
    
                        yt, xt, ys, xs = cm.transform_latlon_to_km_frame(latu, lotu, lasu, losu)
            
                        t = np.concatenate((xt[:,np.newaxis], yt[:,np.newaxis]),axis=1)
                        s = np.concatenate((xs[:,np.newaxis], ys[:,np.newaxis]),axis=1)
    
                        
                        distance[idx0] = cm.Wdist(t,s)
        OTS[:] = distance.reshape(len(Lats), len(Lons))
        lsq = np.append(lsq, np.nanmean(distance))
        
        avgdt, surft= cm.calc_fields(Lats, Lons, vLats, vLons, name = dirReadlow + 'timeseries_per_location_smagorinksi_Cs%s_ddeg%d_sp%d_dd%d.nc'%(i,ddeg,sp,dd))
    
        AVGD[:] = avgdt; SURF[:] = surft;
        ds.close()
        j+=1
    