# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 16:00:14 2019

@author: nooteboom
"""

import numpy as np
from  parcels import Field, FieldSet
import parcels#v200beta
print 'parcels version: ',parcels.__version__
import time
#import parcels

start = time.time()
latsmin = np.arange(-76,82,9)
latsmax = np.arange(-67,90,9)

filename = '/projects/0/palaeo-parcels/POP/POPdata/mesh_0.1degree/grid_coordinates_pop_tx0.1.nc'


bfile = filename#mesh#dirread_pop+'bathymetry_POP_0.1res.nc'

filenames = { 'B' : {'lon': [filename],
                    'lat': [filename],
#                    'depth': filename,
                    'data':[filename]}}

variables = {'B':'BOT_DEP'}

dimensions = {'B': {'lon': 'U_LON_2D', 'lat': 'U_LAT_2D'} } 

bla = 25#59
indices1 = {'lon': range(0,3600),'lat': range(1200, 2400) }
indices2 =  {'lon': range(0,3600),'lat': range(0, 2200) }
fieldsetsbig = FieldSet.from_netcdf(filenames, variables, dimensions, indices = indices1, allow_time_extrapolation=False)
fieldsetsmall = FieldSet.from_netcdf(filenames, variables, dimensions, indices = indices2, allow_time_extrapolation=False)
    #print fieldset.B[0,1,la,0]

#print 'min, max release loc latitudes: ', np.min(), np.max()

print 'latitude and longitude ranges from the indices big:'
print 'longitude: ',np.min(fieldsetsbig.B.lon), np.max(fieldsetsbig.B.lon)
print 'latitude: ',np.min(fieldsetsbig.B.lat), np.max(fieldsetsbig.B.lat)

print 'latitude and longitude ranges from the indices small:'
print 'longitude: ',np.min(fieldsetsmall.B.lon), np.max(fieldsetsmall.B.lon)
print 'latitude: ',np.min(fieldsetsmall.B.lat), np.max(fieldsetsmall.B.lat)

for posidx in range(9,len(latsmin)):
    print 'posidx: ', posidx
    startid = time.time()

#Chaning the lon and lat you must also do this within the kernels    
    minlat = latsmin[posidx];maxlat = latsmax[posidx]
    minlon=-180;maxlon=180;

    res = 1
    dd = 10 #forming depth

    grid = np.mgrid[int(minlon):int(maxlon):res,int(minlat):int(maxlat):res]+0.5
    n=grid[0].size
    lonsp = np.reshape(grid[0],n)
    latsp = np.reshape(grid[1],n)

    print latsp>bla
#Delete the particles on the land
#    filename = '/projects/0/palaeo-parcels/POP/POPdata/mesh_0.1degree/grid_coordinates_pop_tx0.1.nc'
#    [lons,lats]=[np.array([lo for lo, la in zip(lons,lats) if Land[0,lo,la,0]>dd  ]),np.array([la for lo, la in zip(lons,lats) if Land[0,lo,la,0]>dd ])]
    [lons,lats]=[np.array([lo for lo, la in zip(lonsp[latsp<bla],latsp[latsp<bla]) if fieldsetsmall.B[0,0,la, lo]>dd  ]),np.array([la for lo, la in zip(lonsp[latsp<bla],latsp[latsp<bla]) if fieldsetsmall.B[0,0,la, lo]>dd ])]
    [lonsbig,latsbig]=[np.array([lo for lo, la in zip(lonsp[latsp>=bla],latsp[latsp>=bla]) if fieldsetsbig.B[0,0,la, lo]>dd  ]),np.array([la for lo, la in zip(lonsp[latsp>=bla],latsp[latsp>=bla]) if fieldsetsbig.B[0,0,la, lo]>dd ])]
    print latsbig
    lons = np.append(lons,lonsbig); lats = np.append(lats,latsbig)
    #print lats
    print 'number of particles beginning: ', len(lats)
#Remove land particles and mediterranean particles
#[lons,lats]=[np.array([lo for lo, la in zip(lons,lats) if Land[0,lo,la,0]>0.0  ]),np.array([la for lo, la in zip(lons,lats) if Land[0,lo,la,0]>0.0 ])]
    print 'Time for one posidx (in minutes):', ((time.time()-startid)/60.)

    np.save('coor/lons_id%d_dd%d.npy'%(posidx,int(dd)),lons)
    np.save('coor/lats_id%d_dd%d.npy'%(posidx,int(dd)),lats)

print 'Time of the whole execution (in minutes):', ((time.time()-start)/60.)
