# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 15:31:22 2017

@author: nooteboom
"""

from parcels import (ParticleSet, JITParticle, AdvectionRK4_3D,
                     ErrorCode, ParticleFile, Variable, Field)
from datetime import timedelta as delta
from datetime import  datetime
import numpy as np
from glob import glob
import sys
import setfieldsets as setf
import popkernels as popk

dirread_top = '/projects/0/topios/hydrodynamic_data/NEMO-MEDUSA/ORCA0083-N006/'
dirread_pop = '/projects/0/palaeo-parcels/POP/POPdata/'

sp = 6. #The sinkspeed m/day
dd = 10. #The dwelling depth
res = 1 #resolution in degrees
Cs = float(sys.argv[1]) #Diffusion paramete
gmbol = False # Put to true for GM parameterisation/bolus speed
posidx = int(sys.argv[2]) #ID of the file to define latitude and longitude ranges

dirwrite = '/projects/0/palaeo-parcels/POP/POPres/particlefiles/sp%d_dd%d/'%(int(sp),int(dd))

# determine the release grid and the indices to load
if(posidx==1):
    lons, lats = np.meshgrid(np.arange(0,360)+0.5,np.arange(-79, -62)+0.5)
    ind = {'lat':range(0, 83)}
elif(posidx==2):
    lons, lats = np.meshgrid(np.arange(0,360)+0.5,np.arange(-62, -45)+0.5)
    ind = {'lat':range(0, 120)}
elif(posidx==3):
    lons, lats = np.meshgrid(np.arange(0,360)+0.5,np.arange(-45, -27)+0.5)
    ind = {'lat':range(33, 140)}
elif(posidx==4):
    lons, lats = np.meshgrid(np.arange(0,360)+0.5,np.arange(-27, -9)+0.5)
    ind = {'lat':range(70, 200)}
elif(posidx==5):
    lons, lats = np.meshgrid(np.arange(0,360)+0.5,np.arange(-9, 9)+0.5)
    ind = {'lat':range(90, 300)}
elif(posidx==6):
    lons, lats = np.meshgrid(np.arange(0,360)+0.5,np.arange(9, 27)+0.5)
    ind = {'lat':range(170, 320)}
elif(posidx==7):
    lons, lats = np.meshgrid(np.arange(0,360)+0.5,np.arange(27, 45)+0.5)
    ind = {'lat':range(220, 350)}
elif(posidx==8):
    lons, lats = np.meshgrid(np.arange(0,360)+0.5,np.arange(45, 62)+0.5)
    ind = {'lat':range(270, 384)}
elif(posidx==9):
    lons, lats = np.meshgrid(np.arange(0,360)+0.5,np.arange(62, 71)+0.5)
    ind = {'lat':range(300, 384)}

# reshape longitudes and latitudes
lons = lons.flatten(); lats = lats.flatten();
lonsz = []; latsz = [];
lons[lons>320.01] -= 360

# Delete the particles on land
bathy = Field.from_netcdf([dirread_pop+'bathymetry_POP_lowres_320t384.nc'], 'Bathymetry', 
                         {'lon':'ULONG','lat':'ULAT'}, interp_method='bgrid_tracer')

grid = bathy.grid
grid.add_periodic_halo(zonal=True, meridional=False, halosize=20)
bathy.grid = grid

bathy.add_periodic_halo(zonal=True, meridional=False, halosize=20)
for i in range(lons.shape[0]):
    if(bathy[0,0,lats[i], lons[i]]>0):
        lonsz.append(lons[i])
        latsz.append(lats[i])

# The boundary of the grid is at 320 degrees
lonsz = np.array(lonsz); latsz = np.array(latsz);
lonsz[lonsz>320.01] -= 360

if(not lonsz.size):
    sys.exit("Only land in the run with this idx %d"%(posidx))

dep = dd * np.ones(latsz.shape)
times = np.array([datetime(2004, 12, 30) - delta(days=x) for x in range(0,int(365*4),3)])[:140]
time = np.empty(shape=(0));lons = np.empty(shape=(0));lats = np.empty(shape=(0));
for i in range(len(times)):
    lons = np.append(lons,lonsz)
    lats = np.append(lats, latsz)
    time = np.append(time, np.full(len(lonsz),times[i])) 


def run_corefootprintparticles(dirwrite,outfile,lonss,latss,dep):
    files = sorted(glob(dirread_pop+'control_PD_1egree_extra_BOLUS/tavg/t.*'))
    dfile = [dirread_pop+'control_PD_1egree/t.x1_SAMOC_flux.160001.interp.nc']
    bfile = [dirread_pop+'bathymetry_POP_lowres_320t384.nc']
    dimfile = [dirread_pop+'coordinates_curvilinear_pop_grid_320x384.nc']
    afile = [dirread_pop+'spinup_B_2000_cam5_f09_g16.pop.h.1000-01.nc']

    if(gmbol):
        fieldset = setf.set_pop_fieldset_bolus(files, dimfile, dfile, bfile,afile, indices = ind)
    else:
        fieldset = setf.set_pop_fieldset(files, dimfile, dfile, bfile, afile, indices = ind)

    fieldset.add_periodic_halo(zonal=True, halosize=20)
       
    fieldset.add_constant('dwellingdepth', np.float(dd))
    fieldset.add_constant('sinkspeed', sp/86400.)
    fieldset.add_constant('maxage', 300000.*86400)
    fieldset.add_constant('surface', 5.00622)
    fieldset.add_constant('gmbol', gmbol)

    fieldset.add_constant('Cs', Cs)

    class DinoParticle(JITParticle):
        temp = Variable('temp', dtype=np.float32, initial=np.nan)
        age = Variable('age', dtype=np.float32, initial=0.)
        salin = Variable('salin', dtype=np.float32, initial=np.nan)
        lon0 = Variable('lon0', dtype=np.float32, initial=0.)
        lat0 = Variable('lat0', dtype=np.float32, initial=0.)
        depth0 = Variable('depth0',dtype=np.float32, initial=0.) 
        beached = Variable('beached',dtype=np.float32, initial=0.) 

    pset = ParticleSet.from_list(fieldset=fieldset, pclass=DinoParticle, lon=lonss.tolist(), lat=latss.tolist(), 
                       time = time)

    pfile = ParticleFile(dirwrite + outfile, pset, write_ondelete=True)

    if(gmbol):
        advectionkernel = pset.Kernel(popk.AdvectionRK4_3D_addbolus)
        if(Cs>0):
            advectionkernel += popk.smagorinsky_bolus
    else:
        advectionkernel = pset.Kernel(AdvectionRK4_3D)
        if(Cs>0):
            advectionkernel += popk.smagorinsky

    kernels = pset.Kernel(popk.initials) + popk.Sink  + advectionkernel + popk.Age + popk.periodicBC  

    pset.execute(kernels, runtime=delta(days=5*365), dt=delta(minutes=-15), output_file=pfile, verbose_progress=False, recovery={ErrorCode.ErrorOutOfBounds: popk.DeleteParticle})

    print('Execution finished')

if(gmbol):
    outfile = "grid_smagorinskiwn_gm_Cs"+str(Cs)+"_id"+str(posidx)+'_dd'+str(int(dd)) +'_sp'+str(int(sp))+"_res"+str(res)
else:
    outfile = "grid_smagorinskiwn_Cs"+str(Cs)+"_id"+str(posidx)+'_dd'+str(int(dd)) +'_sp'+str(int(sp))+"_res"+str(res)
run_corefootprintparticles(dirwrite,outfile,lons,lats,dep)

