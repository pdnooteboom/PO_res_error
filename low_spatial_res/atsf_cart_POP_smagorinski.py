# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 15:31:22 2017

@author: nooteboom
"""

from parcels import (FieldSet, ParticleSet, JITParticle, AdvectionRK4_3D,
                     ErrorCode, ParticleFile, Variable)
from datetime import timedelta as delta
from datetime import  datetime
import numpy as np
import math
from glob import glob
import sys
import smagorinski

dirread_top = '/projects/0/topios/hydrodynamic_data/NEMO-MEDUSA/ORCA0083-N006/'
dirread_pop = '/projects/0/palaeo-parcels/POP/POPdata/'

sp = 50. #The sinkspeed m/day
dd = 10. #The dwelling depth
res = 1 #resolution in degrees
Cs = float(sys.argv[1]) #Diffusion parameter

posidx = int(sys.argv[2]) #ID of the file to define latitude and longitude ranges

dirwrite = '/projects/0/palaeo-parcels/POP/POPres/particlefiles/sp%d_dd%d/'%(int(sp),int(dd))

latsz = np.load(dirread_pop + 'lats_dd%d.npy'%(int(dd)))
lonsz = np.load(dirread_pop + 'lons_dd%d.npy'%(int(dd)))

# These three lines must be temporary:
w_05 = np.where(lonsz!=0.5)
lonsz = lonsz[w_05]
latsz = latsz[w_05]

if(posidx==0):
    idx = np.where(latsz<-81)
    ind = {'lat':range(0,30)}
elif(posidx==1):
    idx = np.where(np.logical_and(latsz>=-81, latsz<-72))
    ind = {'lat':range(0,40)}
elif(posidx==2):
    idx = np.where(np.logical_and(latsz>=-72, latsz<-63))
    ind = {'lat':range(0,50)}
elif(posidx==3):
    idx = np.where(np.logical_and(latsz>=-63, latsz<-54))
    ind = {'lat':range(10,60)}
elif(posidx==4):
    idx = np.where(np.logical_and(latsz>=-54, latsz<-45))
    ind = {'lat':range(10,70)}
elif(posidx==5):
    idx = np.where(np.logical_and(latsz>=-45, latsz<-36))
    ind = {'lat':range(20,80)}
elif(posidx==6):
    idx = np.where(np.logical_and(latsz>=-36, latsz<-27))
    ind = {'lat':range(30,90)}
elif(posidx==7):
    idx = np.where(np.logical_and(latsz>=-27, latsz<-18))
    ind = {'lat':range(40,90)}
elif(posidx==8):
    idx = np.where(np.logical_and(latsz>=-18, latsz<-9))
    ind = {'lat':range(50,100)}
elif(posidx==9):
    idx = np.where(np.logical_and(latsz>=-9, latsz<0))
    ind = {'lat':range(50,110)}
elif(posidx==10):
    idx = np.where(np.logical_and(latsz>=0, latsz<9))
    ind = {'lat':range(60,120)}
elif(posidx==11):
    idx = np.where(np.logical_and(latsz>=9, latsz<18))
    ind = {'lat':range(70,130)}
elif(posidx==12):
    idx = np.where(np.logical_and(latsz>=18, latsz<27))
    ind = {'lat':range(80,140)}
elif(posidx==13):
    idx = np.where(np.logical_and(latsz>=27, latsz<36))
    ind = {'lat':range(90,150)}
elif(posidx==14):
    idx = np.where(np.logical_and(latsz>=36, latsz<45))
    ind = {'lat':range(100,160)}
elif(posidx==15):
    idx = np.where(np.logical_and(latsz>=45, latsz<54))
    ind = {'lat':range(110,170)}
elif(posidx==16):
    idx = np.where(np.logical_and(latsz>=54, latsz<63))
    ind = {'lat':range(120,180)}
elif(posidx==17):
    idx = np.where(np.logical_and(latsz>=63, latsz<72))
    ind = {'lat':range(130,180)}
elif(posidx==18):
    idx = np.where(np.logical_and(latsz>=72, latsz<81))
    ind = {'lat':range(140,180)}
elif(posidx==19):
    idx = np.where(latsz>=81)
    ind = {'lat':range(150,180)}

lonsz = lonsz[idx]
latsz = latsz[idx]

if(not lonsz.size):
    sys.exit("Only land in the run with this idx")

dep = dd * np.ones(latsz.shape)

times = np.array([datetime(2005, 11, 25) - delta(days=x) for x in range(0,int(365*5+1+30*10),9)])
time = np.empty(shape=(0));lons = np.empty(shape=(0));lats = np.empty(shape=(0));
for i in range(len(times)):
    lons = np.append(lons,lonsz)
    lats = np.append(lats, latsz)
    time = np.append(time, np.full(len(lonsz),times[i])) 
#%%
def set_fieldset(ufiles, wfiles, bfile):#dirread_top+'domain/coordinates.nc'):#
    filenames = { 'U': {'lon': [ufiles[0]],
                        'lat': [ufiles[0]],
                        'depth': [ufiles[0]],
                        'data':ufiles},
                'V' : {'lon': [ufiles[0]],
                        'lat': [ufiles[0]],
                        'depth': [ufiles[0]],
                        'data':ufiles},
                'W' : {'lon': [wfiles[0]],
                        'lat': [wfiles[0]],
                        'depth': [wfiles[0]],
                        'data':wfiles},  
                'S' : {'lon': [ufiles[0]],
                        'lat': [ufiles[0]],
                        'depth': [ufiles[0]],
                        'data':ufiles},   
                'T' : {'lon': [ufiles[0]],
                        'lat': [ufiles[0]],
                        'depth': [ufiles[0]],
                        'data':ufiles}   ,     
                'B' : {'lon': [bfile],
                        'lat': [bfile],
                        'depth': [ufiles[0]],
                        'data':bfile}#,      
                }

    variables = {'U': 'UVEL',
                 'V': 'VVEL',
                 'W': 'WVEL',
                 'T': 'TEMP',
                 'S': 'SALT',
                 'B':'depth'}

    dimensions = {'U':{'lon': 'u_lon', 'lat': 'u_lat', 'depth': 'w_dep', 'time': 'Time'},#
                  'V': {'lon': 'u_lon', 'lat': 'u_lat', 'depth': 'w_dep', 'time': 'Time'},#
                    'W': {'lon': 'u_lon', 'lat': 'u_lat', 'depth': 'w_dep', 'time': 'Time'},#
                    'T': {'lon': 'u_lon', 'lat': 'u_lat', 'depth': 'w_dep', 'time': 'Time'},
                    'S': {'lon': 'u_lon', 'lat': 'u_lat', 'depth': 'w_dep', 'time': 'Time'},
                    'B': {'lon': 'u_lon', 'lat': 'u_lat'} }
    fieldset = FieldSet.from_pop(filenames, variables, dimensions, indices = ind, allow_time_extrapolation=False)
    fieldset.U.vmax = 10    # set max of flow to 10 m/s
    fieldset.V.vmax = 10
    fieldset.W.vmax = 10; fieldset.W.vmin = -10;
    fieldset.T.vmin = -5
    return fieldset
        

def periodicBC(particle, fieldSet, time):
    if particle.lon > 360:
        particle.lon -= 360        
    if particle.lon < 0:
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

def SampleSurf(particle, fieldset, time):
    particle.temp = fieldset.T[time+particle.dt, fieldset.surface, particle.lat, particle.lon]
    particle.salin = fieldset.S[time+particle.dt, fieldset.surface, particle.lat, particle.lon]               

def Age(particle, fieldset, time):
    particle.age = particle.age + math.fabs(particle.dt)  

def DeleteParticle(particle, fieldset, time):
    particle.delete()

def initials(particle, fieldset, time):
    if particle.age==0.:
        particle.depth = fieldset.B[time+particle.dt, fieldset.surface, particle.lat, particle.lon]
        if(particle.depth  > 5370.):
            particle.age = (particle.depth - 5370.)*fieldset.sinkspeed
            particle.depth = 5370.        
        particle.lon0 = particle.lon
        particle.lat0 = particle.lat
        particle.depth0 = particle.depth
        
def FirstParticle(particle, fieldset, time):
    particle.lon = particle.lon0
    particle.lat = particle.lat0
    particle.depth = fieldset.dwellingdepth   

def run_corefootprintparticles(dirwrite,outfile,lonss,latss,dep):
    files = sorted(glob(dirread_pop+'control_PD_1egree/t.x1_SAMOC_flux.160???.interp.nc'))
    wfiles = sorted(glob(dirread_pop+'control_PD_1egree/WVEL_x1_SAMOC_flux_160???.nc'))   
    bfile = dirread_pop+'bathymetry_POP_1res.nc'

    fieldset = set_fieldset(files, wfiles, bfile)    
    fieldset.B.allow_time_extrapolation = True
    fieldset.add_periodic_halo(zonal=True)
       
    fieldset.add_constant('dwellingdepth', np.float(dd))
    fieldset.add_constant('sinkspeed', sp/86400.)
    fieldset.add_constant('maxage', 300000.*86400)
    fieldset.add_constant('surface', 5.00622)

    smagorinski.prepare(fieldset, Cs=Cs)
    fieldset.add_constant('Cs', Cs)

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

    kernels = pset.Kernel(initials) + Sink  + pset.Kernel(AdvectionRK4_3D) + smagorinski.kernels.default + Age + periodicBC  

    pset.execute(kernels, runtime=delta(days=2170), dt=delta(minutes=-5), output_file=pfile, verbose_progress=False, recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})

    print('Execution finished')

outfile = "grid_smagorinski_Cs"+str(Cs)+"_id"+str(posidx)+'_dd'+str(int(dd)) +'_sp'+str(int(sp))+"_res"+str(res)
run_corefootprintparticles(dirwrite,outfile,lons,lats,dep)

