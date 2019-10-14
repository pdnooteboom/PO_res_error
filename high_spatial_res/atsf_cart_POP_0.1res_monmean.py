# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 15:31:22 2017

@author: nooteboom
"""

from parcels import (FieldSet,ParticleSet, JITParticle, AdvectionRK4_3D,
                     ErrorCode, ParticleFile, Variable)
from datetime import timedelta as delta
from datetime import  datetime
import numpy as np
import math
from glob import glob
import sys
from netCDF4 import Dataset

dirread_POP = '/projects/0/samoc/pop/tx0.1/output/run_henk_mixedbc_extravars_nooteboom/'
mesh = '/projects/0/palaeo-parcels/POP/POPdata/mesh_0.1degree/'
dirread_pop_monmean = '/projects/0/palaeo-parcels/POP/POPdata/control_PD_0.1degree_monthlymean/'

sp = 6. #The sinkspeed m/day
dd = 10. #The dwelling depth
res = 1 #resolution in degrees

dirwrite = '/projects/0/palaeo-parcels/POP/POPres/0.1degree/particlefiles/sp%d_dd%d/'%(int(sp),int(dd))

posidx = int(sys.argv[1]) #ID of the file to define latitude and longitude ranges

latsz = np.load('release_loc/coor/lats_id%d_dd%d.npy'%(posidx,int(dd)))
lonsz = np.load( 'release_loc/coor/lons_id%d_dd%d.npy'%(posidx,int(dd)))

latsz = latsz; lonsz= lonsz;#np.array([latsz[0]]); lonsz = np.array([lonsz[0]]);

print 'here0 : ',latsz.shape, lonsz.shape
print 'latsz: ' , np.min(latsz), np.max(latsz)
print 'lonsz: ' , np.min(lonsz), np.max(lonsz)

if(not lonsz.size):
    sys.exit("Only land in the run with this idx")

dep = dd * np.ones(latsz.shape)

times = np.array([datetime(1974, 11, 30) - delta(days=x) for x in range(0,365*5,3)])#int(365*8),3)])
time = np.empty(shape=(0));lons = np.empty(shape=(0));lats = np.empty(shape=(0));
for i in range(len(times)):
    lons = np.append(lons,lonsz)
    lats = np.append(lats, latsz)
    time = np.append(time, np.full(len(lonsz),times[i]))


#%%
def set_fieldset(files_highres, hormesh, sfile):
    ufiles = files_highres 
    bfile = hormesh

    filenames = { 'U': {'lon': hormesh,
                        'lat': hormesh,
                        'depth': sfile,
                        'data':ufiles},
                'V' : {'lon': hormesh,
                        'lat': hormesh,
                        'depth': sfile,
                        'data':ufiles},
                'W' : {'lon': hormesh,
                        'lat': hormesh,
                        'depth': sfile,
                        'data':ufiles},  
                'S' : {'lon': hormesh,
                        'lat': hormesh,
#                        'depth': sfile,
                        'data':ufiles},
                'T' : {'lon': hormesh,
                        'lat': hormesh,
#                        'depth': sfile,
                        'data':ufiles},
                'B' : {'lon': hormesh,
                        'lat': hormesh,
                        'depth': sfile,
                        'data':hormesh}#, 
                }

    variables = {'U': 'UVEL',
                 'V': 'VVEL',
                 'W': 'WVEL',
                 'T': 'TEMP_5m',
                 'S': 'SALT_5m',
                 'B':'BOT_DEP'}

    dimensions = {'U':{'lon': 'U_LON_2D', 'lat': 'U_LAT_2D', 'depth': 'BOTTOM_GRIDCELL','time': 'time'},#
                  'V': {'lon': 'U_LON_2D', 'lat': 'U_LAT_2D', 'depth': 'BOTTOM_GRIDCELL','time': 'time'},#
                    'W': {'lon': 'U_LON_2D', 'lat': 'U_LAT_2D', 'depth': 'BOTTOM_GRIDCELL','time': 'time'},#
                    'T': {'lon': 'U_LON_2D', 'lat': 'U_LAT_2D', 'time':'time'},
                    'S': {'lon': 'U_LON_2D', 'lat': 'U_LAT_2D', 'time':'time'},
                    'B': {'lon': 'U_LON_2D', 'lat': 'U_LAT_2D'} }



    #Determine the latitude indices to be loaded
    gf = Dataset(hormesh)
    latsmin = np.arange(-76,82,9)
    latsmax = np.arange(-67,90,9)
    #Chaning the lon and lat you must also do this within the kernels    
    minlat = latsmin[posidx];maxlat = latsmax[posidx]

    if(posidx==0 or (posidx >=9 and posidx<=14)):
        latrange = 20
    else:
        latrange = 15

    for i_i in range(2400):
        minla = min(gf['U_LAT_2D'][i_i])
        maxla = max(gf['U_LAT_2D'][i_i])
        if(maxla<minlat-latrange):
            minlaidx = i_i
        elif(minlat-latrange<-78.45172):
            minlaidx = 0
        if(posidx>=13):
            maxlaidx = 2400
        elif(minla<maxlat+latrange):
            maxlaidx = i_i
    latind = range(minlaidx, maxlaidx)
    print 'minlat, maxlat:', minlat,maxlat
    print 'latind: ', latind[0], np.array(latind)[-1]
    indices = {'lat': latind}
    gf.close()

#old:
#    latrange = 15
#
#    for i_i in range(2400):
#        minla = min(gf['U_LAT_2D'][i_i])
#        maxla = max(gf['U_LAT_2D'][i_i])
#        if(maxla<minlat-latrange):
#            minlaidx = i_i
#        elif(minlat-latrange<-78.45172):
#            minlaidx = 0
#        if(posidx>=13):
#            maxlaidx = 2400
#        elif(minla<maxlat+latrange):
#            maxlaidx = i_i
#    latind = range(minlaidx, maxlaidx)
#    print 'latind: ', latind[0], np.array(latind)[-1]
#    indices = {'lat': latind}
#    gf.close()


    fieldset = FieldSet.from_pop(filenames, variables, dimensions, indices=indices, allow_time_extrapolation=False)
    
    print 'latitude and longitude ranges from the indices:' 
    print 'longitude: ',fieldset.U.lon[0], fieldset.U.lon[-1]
    print 'latitude: ',np.max(fieldset.U.lat[0]), np.min(fieldset.U.lat[-1])
#    print lonsz[0],latsz[0]
#    print 'depth at bottom: ', fieldset.B[0,lonsz[0], latsz[0],0]

    fieldset.U.vmax = 10    # set max of flow to 10 m/s
    fieldset.V.vmax = 10
    fieldset.W.vmax = 10
    fieldset.T.vmin = -5
    return fieldset
        

def periodicBC(particle, fieldSet, time):
    if particle.lon > 180:
        particle.lon -= 360        
    if particle.lon < -180:
        particle.lon += 360   
        
#Sink Kernel if only atsf is saved:
def Sink(particle, fieldset, time):
    if(particle.depth>fieldset.dwellingdepth):
        particle.depth = particle.depth + fieldset.sinkspeed * particle.dt
    elif(particle.depth<=fieldset.dwellingdepth and particle.depth>1):
        particle.depth = fieldset.surface
        particle.temp = fieldset.T[time+particle.dt, fieldset.surface, particle.lat, particle.lon]
        particle.salin = fieldset.S[time+particle.dt, fieldset.surface, particle.lat, particle.lon]
        particle.delete()        

def Age(particle, fieldset, time):
    particle.age = particle.age + math.fabs(particle.dt)  

def DeleteParticle(particle, fieldset, time):
    particle.delete()

def initials(particle, fieldset, time):
    if particle.age==0.:
        particle.depth = fieldset.B[time, fieldset.surface, particle.lat, particle.lon] - 1
        if(particle.depth  > 5370.):
            particle.age = (particle.depth - 5370.)*fieldset.sinkspeed
            particle.depth = 5370. -1       
        particle.lon0 = particle.lon
        particle.lat0 = particle.lat
        particle.depth0 = particle.depth

def run_corefootprintparticles(dirwrite,outfile,lonss,latss,dep):
    hormesh = mesh + 'grid_coordinates_pop_tx0.1.nc'
    Sdepth = mesh + 'bottom_cell.nc'
    files_highres = sorted(glob(dirread_pop_monmean+'t.t0.1_42l_nccs01.03??????_monthly.nc'))
    fieldset = set_fieldset(files_highres, hormesh, Sdepth)
    fieldset.B.allow_time_extrapolation = True

    fieldset.add_periodic_halo(zonal=True)       
    fieldset.add_constant('dwellingdepth', np.float(dd))
    fieldset.add_constant('sinkspeed', sp/86400.)
    fieldset.add_constant('maxage', 300000.*86400)
    fieldset.add_constant('surface', 5.00622)

    class DinoParticle(JITParticle):
        temp = Variable('temp', dtype=np.float32, initial=np.nan)
        age = Variable('age', dtype=np.float32, initial=0.)
        salin = Variable('salin', dtype=np.float32, initial=np.nan)
        lon0 = Variable('lon0', dtype=np.float32, initial=0.)
        lat0 = Variable('lat0', dtype=np.float32, initial=0.)
        depth0 = Variable('depth0',dtype=np.float32, initial=0.) 

    pset = ParticleSet.from_list(fieldset=fieldset, pclass=DinoParticle, lon=lonss.tolist(), lat=latss.tolist(), 
                       time = time)

    pfile = ParticleFile(dirwrite + outfile, pset, write_ondelete=True)

    kernels = pset.Kernel(initials) + Sink  + pset.Kernel(AdvectionRK4_3D) + Age + periodicBC  

    pset.execute(kernels, runtime=delta(days=365*5), dt=delta(minutes=-15), 
        output_file=pfile, verbose_progress=False, 
        recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})

    print 'Execution finished'

outfile = "grid"+'_id'+str(posidx)+'_dd'+str(int(dd)) +'_sp'+str(int(sp)) +'_tempresmonmean'
run_corefootprintparticles(dirwrite,outfile,lons,lats,dep)

