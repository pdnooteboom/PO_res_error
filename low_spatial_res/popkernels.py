import math
from random import random

def periodicBC(particle, fieldSet, time):
    if particle.lon > 320:
        particle.lon -= 360
    if particle.lon < -40:
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
        particle.depth = fieldset.B[time, particle.depth, particle.lat, particle.lon] - 0.01
        if(particle.depth  > 5470.):
            particle.age = (particle.depth - 5470.)*fieldset.sinkspeed
            particle.depth = 5470.
        particle.lon0 = particle.lon
        particle.lat0 = particle.lat
        particle.depth0 = particle.depth

def FirstParticle(particle, fieldset, time):
    particle.lon = particle.lon0
    particle.lat = particle.lat0
    particle.depth = fieldset.dwellingdepth

def AdvectionRK4_3D_addbolus(particle, fieldset, time):
    """Advection of particles using fourth-order Runge-Kutta integration including vertical velocity.
    Function needs to be converted to Kernel object before execution"""
    m_to_degzon = 1 / (1852*60*math.cos(particle.lat*math.pi/180))
    m_to_degmer = 1 / (1852*60)
    (u1, v1, w1) = fieldset.UVW[time, particle.depth, particle.lat, particle.lon]
    u1 += fieldset.Ubolus[time, particle.depth, particle.lat, particle.lon] * m_to_degzon
    v1 += fieldset.Vbolus[time, particle.depth, particle.lat, particle.lon] * m_to_degmer
    lon1 = particle.lon + u1*.5*particle.dt
    lat1 = particle.lat + v1*.5*particle.dt
    dep1 = particle.depth + w1*.5*particle.dt
    (u2, v2, w2) = fieldset.UVW[time + .5 * particle.dt, dep1, lat1, lon1]
    u2 += fieldset.Ubolus[time + .5 * particle.dt, dep1, lat1, lon1] * m_to_degzon
    v2 += fieldset.Vbolus[time + .5 * particle.dt, dep1, lat1, lon1] * m_to_degmer
    lon2 = particle.lon + u2*.5*particle.dt
    lat2 = particle.lat + v2*.5*particle.dt
    dep2 = particle.depth + w2*.5*particle.dt
    (u3, v3, w3) = fieldset.UVW[time + .5 * particle.dt, dep2, lat2, lon2]
    u3 += fieldset.Ubolus[time + .5 * particle.dt, dep2, lat2, lon2] * m_to_degzon
    v3 += fieldset.Vbolus[time + .5 * particle.dt, dep2, lat2, lon2] * m_to_degmer
    lon3 = particle.lon + u3*particle.dt
    lat3 = particle.lat + v3*particle.dt
    dep3 = particle.depth + w3*particle.dt
    (u4, v4, w4) = fieldset.UVW[time + particle.dt, dep3, lat3, lon3]
    u4 += fieldset.Ubolus[time + particle.dt, dep3, lat3, lon3] * m_to_degzon
    v4 += fieldset.Vbolus[time + particle.dt, dep3, lat3, lon3] * m_to_degmer
    particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
    particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt
    particle.depth += (w1 + 2*w2 + 2*w3 + w4) / 6. * particle.dt



def smagorinsky_bolus(particle, fieldset, time):
    m_to_degzon = 1 / (1852*60*math.cos(particle.lat*math.pi/180))
    m_to_degmer = 1 / (1852*60)
    dx = 0.1;
    dudx = (fieldset.U[time, particle.depth, particle.lat, particle.lon+dx]+fieldset.Ubolus[time, particle.depth, particle.lat, particle.lon+dx]*m_to_degzon-
            (fieldset.U[time, particle.depth, particle.lat, particle.lon-dx]+fieldset.Ubolus[time, particle.depth, particle.lat, particle.lon-dx]*m_to_degzon)) / (2*dx)

    dudy = (fieldset.U[time, particle.depth, particle.lat+dx, particle.lon]+fieldset.Ubolus[time, particle.depth, particle.lat+dx, particle.lon]*m_to_degzon-
            (fieldset.U[time, particle.depth, particle.lat-dx, particle.lon]+fieldset.Ubolus[time, particle.depth, particle.lat-dx, particle.lon]*m_to_degzon)) / (2*dx)

    dvdx = (fieldset.V[time, particle.depth, particle.lat, particle.lon+dx]+fieldset.Vbolus[time, particle.depth, particle.lat, particle.lon+dx]*m_to_degmer-
            (fieldset.V[time, particle.depth, particle.lat, particle.lon-dx]+fieldset.Vbolus[time, particle.depth, particle.lat, particle.lon-dx]*m_to_degmer)) / (2*dx)

    dvdy = (fieldset.V[time, particle.depth, particle.lat+dx, particle.lon]+ fieldset.Vbolus[time, particle.depth, particle.lat+dx, particle.lon]*m_to_degmer-
            (fieldset.V[time, particle.depth, particle.lat-dx, particle.lon]+fieldset.Vbolus[time, particle.depth, particle.lat-dx, particle.lon]*m_to_degmer)) / (2*dx)

    A = fieldset.cell_areas[time, 0, particle.lat, particle.lon]
    deg_to_m = (1852*60)**2*math.cos(particle.lat*math.pi/180)
    A = A / deg_to_m
    Vh = fieldset.Cs * A * math.sqrt(dudx**2 + 0.5*(dudy + dvdx)**2 + dvdy**2)

    xres = 1 # [degrees] 
    yres = 1

    hitBoundary = True
    tries = 0

    while hitBoundary and tries <5:
        dlat = yres * random.normalvariate(0., 1.) * math.sqrt(2*math.fabs(particle.dt)* Vh) #/ dy_cell
        dlon = xres * random.normalvariate(0., 1.) * math.sqrt(2*math.fabs(particle.dt)* Vh) #/ dx_cell
        if(fieldset.U[time, particle.depth, particle.lat+dlat, particle.lon+dlon] > 0):
            if(fieldset.V[time, particle.depth, particle.lat+dlat, particle.lon+dlon] > 0):
                hitBoundary = False
            elif(fieldset.V[time, particle.depth, particle.lat+dlat, particle.lon+dlon] < 0):
                hitBoundary = False
        elif(fieldset.U[time, particle.depth, particle.lat+dlat, particle.lon+dlon] < 0):
            if(fieldset.V[time, particle.depth, particle.lat+dlat, particle.lon+dlon] > 0):
                hitBoundary = False
            elif(fieldset.V[time, particle.depth, particle.lat+dlat, particle.lon+dlon] < 0):
                hitBoundary = False
        tries += 1
        if tries == 5:
            particle.beached = 1
            dlat = 0
            dlon = 0

    particle.lat += dlat
    particle.lon += dlon

def smagorinsky(particle, fieldset, time):
    dx = 0.01;
    dudx = (fieldset.U[time, particle.depth, particle.lat, particle.lon+dx]-fieldset.U[time, particle.depth, particle.lat, particle.lon-dx]) / (2*dx)
    dudy = (fieldset.U[time, particle.depth, particle.lat+dx, particle.lon]-fieldset.U[time, particle.depth, particle.lat-dx, particle.lon]) / (2*dx)
    dvdx = (fieldset.V[time, particle.depth, particle.lat, particle.lon+dx]-fieldset.V[time, particle.depth, particle.lat, particle.lon-dx]) / (2*dx)
    dvdy = (fieldset.V[time, particle.depth, particle.lat+dx, particle.lon]-fieldset.V[time, particle.depth, particle.lat-dx, particle.lon]) / (2*dx)

    A = fieldset.cell_areas[time, 0, particle.lat, particle.lon]
    deg_to_m = (1852*60)**2*math.cos(particle.lat*math.pi/180)
    A = A / deg_to_m
    Vh = fieldset.Cs * A * math.sqrt(dudx**2 + 0.5*(dudy + dvdx)**2 + dvdy**2)

    xres = 1# [degrees] 
    yres = 1

    hitBoundary = True
    tries = 0

    while hitBoundary and tries <5:
        dlat = yres * random.normalvariate(0., 1.) * math.sqrt(2*math.fabs(particle.dt)* Vh) #/ dy_cell
        dlon = xres * random.normalvariate(0., 1.) * math.sqrt(2*math.fabs(particle.dt)* Vh) #/ dx_cell
        if(fieldset.U[time, particle.depth, particle.lat+dlat, particle.lon+dlon] > 0):
            if(fieldset.V[time, particle.depth, particle.lat+dlat, particle.lon+dlon] > 0):
                hitBoundary = False
            elif(fieldset.V[time, particle.depth, particle.lat+dlat, particle.lon+dlon] < 0):
                hitBoundary = False
        elif(fieldset.U[time, particle.depth, particle.lat+dlat, particle.lon+dlon] < 0):
            if(fieldset.V[time, particle.depth, particle.lat+dlat, particle.lon+dlon] > 0):
                hitBoundary = False
            elif(fieldset.V[time, particle.depth, particle.lat+dlat, particle.lon+dlon] < 0):
                hitBoundary = False
        tries += 1
        if tries == 5:
            particle.beached = 1
            dlat = 0
            dlon = 0

    particle.lat += dlat
    particle.lon += dlon
