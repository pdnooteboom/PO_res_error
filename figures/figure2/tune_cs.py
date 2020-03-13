#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 11:35:51 2019

@author: nooteboom
"""

import numpy as np
from netCDF4 import Dataset
import calc_metrics2 as cm

def plotamount(nc, ty='', lim=130):
    import matplotlib.pylab as plt
    Lats = nc['temp'][:]
    sizedist = np.zeros(Lats.shape[0])
    for i in range(Lats.shape[0]):
        sizedist[i] = np.sum((Lats[i]>-200))
    
    sizedist = sizedist.reshape(len(nc['Lats'][:]), len(nc['Lons'][:]))
    plt.imshow(np.flip(sizedist,0))
    plt.colorbar()
    plt.title('sizes of particle distributions ' + ty)
    plt.show()

sp = 6
ddeg = 1
dd = 10

minmaxlat = 100000.
maxminlat = -100000
allresults= ['monmean','null','0.0','0.25', '0.5', '1.0', '2.0', '5.0', '0.0gm','0.25gm','0.5gm', '1.0gm', '2.0gm', '5.0gm']

dirWrite = '/Volumes/HardDisk/POP/output/OT/'  #'/projects/0/palaeo-parcels/POP/OT/'
dirReadhigh = '/Volumes/HardDisk/POP/output/highres/timeseries/'   
dirReadlow = '/Volumes/HardDisk/POP/output/lowres/timeseries/'  
nc_target = Dataset(dirReadhigh + 'timeseries_per_location_ddeg%d_sp%d_dd%d_tempres5_ds2.nc'%(ddeg,sp,dd))

tslenmin = 10000
hid = 130#

nc = Dataset(dirReadhigh + 'timeseries_per_location_ddeg%d_sp%d_dd%d_tempres5_ds2.nc'%(ddeg,sp,dd))
#plotamount(nc, ty='ref')
minmaxlat = min(minmaxlat,np.max(nc['Lats'][:]))
maxminlat = max(maxminlat, np.min(nc['Lats'][:]))
tslen = nc['tslens'][:]
print('minlat, maxlat  '+str(np.min(nc['Lats'][:])) + '   ' + str(np.max(nc['Lats'][:])))
print(np.sum(tslen>hid)/np.float(len(tslen[tslen>0])))
tslenmin = min(tslenmin, np.min(tslen[tslen>hid]))

for i in allresults:
    print(i)
    if('monmean'==i):
        nc = Dataset(dirReadhigh + 'timeseries_per_location_ddeg1_sp%d_dd10_tempresmonmean.nc'%(sp))
#        plotamount(nc, ty=i)
        minmaxlat = min(minmaxlat,np.max(nc['Lats'][:]))
        maxminlat = max(maxminlat, np.min(nc['Lats'][:]))
        tslen = nc['tslens'][:]
        print('minlat, maxlat  '+str(np.min(nc['Lats'][:])) + '   ' + str(np.max(nc['Lats'][:])))
        print(np.sum(tslen>hid)/np.float(len(tslen[tslen>0])))
        tslenmin = min(tslenmin, np.min(tslen[tslen>hid]))
    elif('null' == i):
        nc = Dataset(dirReadhigh + 'timeseries_per_location_ddeg%d_sp%d_dd%d_tempres5.nc'%(ddeg,sp,dd))
#        plotamount(nc, ty=i)
        minmaxlat = min(minmaxlat,np.max(nc['Lats'][:]))
        maxminlat = max(maxminlat, np.min(nc['Lats'][:]))
        tslen = nc['tslens'][:]
        print(np.sum(tslen>hid)/np.float(len(tslen[tslen>0])))
        tslenmin = min(tslenmin, np.min(tslen[tslen>hid]))
    else:
        if(i[-2:]=='gm'):
            nc = Dataset(dirReadlow + 'timeseries_per_location_smagorinksi_wn_gm_Cs%s_ddeg%d_sp%d_dd%d.nc'%(i[:-2],ddeg,sp,dd))
        else:
            nc = Dataset(dirReadlow + 'timeseries_per_location_smagorinksi_wn_Cs%s_ddeg%d_sp%d_dd%d.nc'%(i,ddeg,sp,dd))
#        plotamount(nc, ty=i)       
        minmaxlat = min(minmaxlat,np.max(nc['Lats'][:]))
        maxminlat = max(maxminlat, np.min(nc['Lats'][:]))
        vLons = nc['vLons'][:]
        vLats = nc['vLats'][:]
        tslen = nc['tslens'][:]
        print(np.sum(tslen>hid)/np.float(len(tslen[tslen>0])))
        tslenmin = min(tslenmin, np.min(tslen[tslen>hid]))
    nc.close()
    print(tslenmin)
Lons = np.arange(361)
Lats = np.arange(maxminlat, minmaxlat+1);
vLons, vLats = np.meshgrid(Lons, Lats);  vLons = vLons.flatten(); vLats = vLats.flatten();

nc_target = Dataset(dirReadhigh + 'timeseries_per_location_ddeg%d_sp%d_dd%d_tempres5_ds2.nc'%(ddeg,sp,dd))
idxt = np.where(np.logical_and(nc_target['vLats'][:]>=Lats[0],nc_target['vLats'][:]<=Lats[-1]))
lat = nc_target['lat'][:]; lon = nc_target['lon'][:] 
vLons_target =  nc_target['vLons'][:]; vLats_target =  nc_target['vLats'][:]; 

lsq = np.array([])

j=0
for i in allresults:
    if('monmean'==i): 
        distance = np.full(len(vLons), np.nan)
        print('monmean')
        ds = Dataset(dirWrite + 'OTs_monmean_sp%d.nc'%(sp), 'w')
        lats = ds.createDimension('lat',len(Lats))
        latitudes = ds.createVariable('Lats', np.float32, ('lat',))
        latitudes[:] = Lats
        lons = ds.createDimension('lon',len(Lons))
        longitudes = ds.createVariable('Lons', np.float32, ('lon',))
        longitudes[:] = Lons
        OTS = ds.createVariable('Wasserstein distance', np.float32, ('lat','lon',))

        nc_source = Dataset(dirReadhigh + 'timeseries_per_location_ddeg1_sp%d_dd10_tempresmonmean.nc'%(sp))
        lat_s = nc_source['lat'][:]
        lon_s = nc_source['lon'][:]
        vLats_s = nc_source['vLats'][:]
        vLons_s = nc_source['vLons'][:]
    
        for j in range(Lons.shape[0]-1):
            if(j%30==0): print(j/float(len(Lons)));
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
    
                        distance[idx0] = cm.Wdist(s,t)
        OTS[:] = distance.reshape(len(Lats), len(Lons))
        
        ds.close()    
    elif('null'==i): 
        distance = np.full(len(vLons), np.nan)
        print('nul test')
        ds = Dataset(dirWrite + 'OTs_nultest.nc', 'w')
        lats = ds.createDimension('lat',len(Lats))
        latitudes = ds.createVariable('Lats', np.float32, ('lat',))
        latitudes[:] = Lats
        lons = ds.createDimension('lon',len(Lons))
        longitudes = ds.createVariable('Lons', np.float32, ('lon',))
        longitudes[:] = Lons
        OTS = ds.createVariable('Wasserstein distance', np.float32, ('lat','lon',))
        
    
        nc_source  = Dataset(dirReadhigh + 'timeseries_per_location_ddeg%d_sp%d_dd%d_tempres5.nc'%(ddeg,sp,dd))
        lat_s = nc_source['lat'][:]
        lon_s = nc_source['lon'][:]
        vLats_s = nc_source['vLats'][:]
        vLons_s = nc_source['vLons'][:]

        for j in range(Lons.shape[0]-1):
            if(j%30==0): print(j/float(len(Lons)));
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

                        distance[idx0] = cm.Wdist(s,t)
                        
        OTS[:] = distance.reshape(len(Lats), len(Lons))
  
        ds.close()   
    else:
        print('CS: %s'%i)
        distance = np.full(len(vLons), np.nan)
        ds = Dataset(dirWrite + 'OTs_sp%d_cs%s.nc'%(sp,i), 'w')
        lats = ds.createDimension('lat',len(Lats))
        latitudes = ds.createVariable('Lats', np.float32, ('lat',))
        latitudes[:] = Lats
        lons = ds.createDimension('lon',len(Lons))
        longitudes = ds.createVariable('Lons', np.float32, ('lon',))
        longitudes[:] = Lons
        OTS = ds.createVariable('Wasserstein distance', np.float32, ('lat','lon',))
 
        if(i[-2:]=='gm'):
            cs = np.float(i[:-2])
        else:
            cs = np.float(i)

        if(i[-2:]=='gm'):
            if(cs!=0.25):
                nc_source = nc_target = Dataset(dirReadlow + 'timeseries_per_location_smagorinksi_wn_gm_Cs%.1f_ddeg%d_sp%d_dd%d.nc'%(cs,ddeg,sp,dd))
            else:
                nc_source = nc_target = Dataset(dirReadlow + 'timeseries_per_location_smagorinksi_wn_gm_Cs%.2f_ddeg%d_sp%d_dd%d.nc'%(cs,ddeg,sp,dd))
        else:
            if(cs!=0.25):
                nc_source = nc_target = Dataset(dirReadlow + 'timeseries_per_location_smagorinksi_wn_Cs%.1f_ddeg%d_sp%d_dd%d.nc'%(cs,ddeg,sp,dd))
            else:
                nc_source = nc_target = Dataset(dirReadlow + 'timeseries_per_location_smagorinksi_wn_Cs%.2f_ddeg%d_sp%d_dd%d.nc'%(cs,ddeg,sp,dd))

        lat_s = nc_source['lat'][:]
        lon_s = nc_source['lon'][:]
        vLats_s = nc_source['vLats'][:]
        vLons_s = nc_source['vLons'][:]
    
        for j in range(Lons.shape[0]):
            if(j%30==0): print(j/float(len(Lons)));
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
    
                        
                        distance[idx0] = cm.Wdist(s,t)
        OTS[:] = distance.reshape(len(Lats), len(Lons))
        lsq = np.append(lsq, np.nanmean(distance))
        
        ds.close()
        j+=1