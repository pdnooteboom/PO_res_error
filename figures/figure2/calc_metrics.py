#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 10:14:27 2019

@author: nooteboom
"""

import numpy as np
import matplotlib.pylab as plt
import ot
from netCDF4 import Dataset
from numba import jit
import math

#%% general
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
def uniqueloc(lats, lons, la, lo):
    bo = True
    if(len(lons)>0):
        for i in range(len(lons)):
            if(lons[i]==lo):
                if(lats[i]==la):
                    bo = False
    return bo
#%% calculate average lateral distance
@jit(nopython=True)
def avgdistf(lat, lon, avgdist, tsl, vLons, vLats,ml):
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


#%% calculate the surface area
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

def calc_fields(Lats, Lons, vLats, vLons, name = '', ml=140 ):
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
    
    print(avgdist.shape)
    
    return avgdist,surface

#%%

def transform_latlon_to_km_frame(lat,lon, lat1, lon1):
    # lat and lon should be 1D arrays
    assert len(lat.shape)==1, 'should be 1D array'
    conc_lon = np.append(lon,lon1)
    if((conc_lon<100).any() and (conc_lon>260).any()):
        lon[np.where(lon<100)] += 360
        lon1[np.where(lon1<100)] += 360   
    
    lat0, lon0 = lat[0], lon[0]
    
    y=np.zeros(lat.shape); x=np.zeros(lat.shape);
    y1=np.zeros(lat1.shape); x1=np.zeros(lat1.shape);
    for i in range(lat.shape[0]):
        y[i] = 1852 * 60 * (lat[i]-lat0) / 1000.
        x[i] = (1852 * 60 * np.cos((lat[i]+lat0)/2.*np.pi/180.))*(lon[i]-lon0) / 1000.
    for i in range(lat1.shape[0]):
        y1[i] = 1852 * 60 * (lat1[i]-lat0) / 1000.
        x1[i] = (1852 * 60 * np.cos((lat1[i]+lat0)/2.*np.pi/180.))*(lon1[i]-lon0) / 1000.
    return y,x,y1,x1

def test_transform_latlon_to_km_frame():
    
    lathh = np.array([0,1,2]); lonhh=np.array([355,2,340]);
    lat1 = np.array([0,1,2])+2; lon1=np.array([1,2,1])+2;    
    
    plt.scatter(lonhh,lathh, label=1)
    plt.scatter(lon1, lat1, label=2)
    plt.legend()
    plt.show()
    
    y,x,y1,x1 = transform_latlon_to_km_frame(lathh,lonhh, lat1, lon1)
    
    plt.scatter(x,y, label=1)
    plt.scatter(x1, y1, label=2)
    plt.legend()
    plt.show() 

#test_transform_latlon_to_km_frame()

def Wdist(t, s):
#    
    M = ot.dist(s, t)
    Wd = ot.emd2([], [], M)

    return Wd
    
def test_Wdist():
    t = np.random.rand(2,2)
    s = np.random.rand(2,2)*2 + 1
    wd = Wdist(t, s)
    print(wd)

#test_Wdist()