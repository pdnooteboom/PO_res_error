#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 15:38:49 2019

@author: nooteboom
"""

from parcels import Field
import math
from random import random
import numpy as np

def prepare(fieldset, Cs=0.1):
    """
    Add cell_areas field and gradients of U and V that are both necessary for Smagorinsky parametrization.
        
    :param fieldset: mod:`parcels.fieldset.FieldSet` object to add necessary fields to
    """
    fieldset.add_constant('Cs', Cs)
    fieldset.add_constant('resolutionx',fieldset.U.grid.lon[1]-fieldset.U.grid.lon[0])
    fieldset.add_constant('resolutiony',fieldset.U.grid.lat[1]-fieldset.U.grid.lat[0])
    
#    (dU_dx, dU_dy) = fieldset.U.gradient()
#    fieldset.add_field(dU_dx)
#    fieldset.add_field(dU_dy)

#    (dV_dx, dV_dy) = fieldset.V.gradient()
#    fieldset.add_field(dV_dx)
#    fieldset.add_field(dV_dy)
    
    x = fieldset.U.grid.lon
    y = fieldset.U.grid.lat
    
    cell_areas = Field(name='cell_areas', data=fieldset.U.cell_areas(), lon=x, lat=y)
    fieldset.add_field(cell_areas)
    fieldset.U.calc_cell_edge_sizes()
   
    cell_edge_sizes_x = Field(name='cell_edge_sizes_x', data=fieldset.U.cell_edge_sizes['x'], lon=x, lat=y)
    cell_edge_sizes_y = Field(name='cell_edge_sizes_y', data=fieldset.U.cell_edge_sizes['y'], lon=x, lat=y)
    fieldset.add_field(cell_edge_sizes_x)
    fieldset.add_field(cell_edge_sizes_y)

class kernels:
    def __init__(self):
        self.members = ['uniform', 'spatially_varying']
 
#    def default(particle, fieldset, time): 
#        """Smagorinski parametrization kernel. 
#        """
#        dudx = fieldset.dU_dx[time, 0., particle.lat, particle.lon]
#        dudy = fieldset.dU_dy[time, 0., particle.lat, particle.lon]
#        dvdx = fieldset.dV_dx[time, 0., particle.lat, particle.lon]
#        dvdy = fieldset.dV_dy[time, 0., particle.lat, particle.lon]
#
#        A = fieldset.cell_areas[time, 0., particle.lat, particle.lon]
#        Vh = fieldset.Cs * A * math.sqrt(dudx**2 + 0.5*(dudy + dvdx)**2 + dvdy**2)
#
#        r = 1/3.
#        xres = fieldset.resolutionx # [degrees] 
#        yres = fieldset.resolutiony
#        dx_cell = fieldset.cell_edge_sizes_x[0, 0, particle.lat, particle.lon] # [meters]
#        dy_cell = fieldset.cell_edge_sizes_y[0, 0, particle.lat, particle.lon]
#
#        hitBoundary = True
#        tries = 0
#        while hitBoundary and tries <10:
#            dlat = yres * random.uniform(-1., 1.) * math.sqrt(2*math.fabs(particle.dt)* Vh/r) / dy_cell
#            dlon = xres * random.uniform(-1., 1.) * math.sqrt(2*math.fabs(particle.dt)* Vh/r) / dx_cell
#            if fieldset.UV[time, 0., particle.lat+dlat, particle.lon+dlon] != (0, 0):
#                hitBoundary = False
#            tries += 1
#            if tries == 10:
#                dlat = 0
#                dlon = 0
#
#        particle.lat += dlat
#        particle.lon += dlon    
        
 
    def default(particle, fieldset, time): 
        """Smagorinski parametrization kernel. 
        """
        dx = 0.01;
        dudx = (fieldset.U[time, particle.depth, particle.lat, particle.lon+dx]-fieldset.U[time, particle.depth, particle.lat, particle.lon-dx]) / (2*dx)
        dudy = (fieldset.U[time, particle.depth, particle.lat+dx, particle.lon]-fieldset.U[time, particle.depth, particle.lat-dx, particle.lon]) / (2*dx)
        dvdx = (fieldset.V[time, particle.depth, particle.lat, particle.lon+dx]-fieldset.V[time, particle.depth, particle.lat, particle.lon-dx]) / (2*dx)
        dvdy = (fieldset.V[time, particle.depth, particle.lat+dx, particle.lon]-fieldset.V[time, particle.depth, particle.lat-dx, particle.lon]) / (2*dx)

#        dudx = fieldset.dU_dx[time, particle.depth, particle.lat, particle.lon]
#        dudy = fieldset.dU_dy[time, particle.depth, particle.lat, particle.lon]
#        dvdx = fieldset.dV_dx[time, particle.depth, particle.lat, particle.lon]
#        dvdy = fieldset.dV_dy[time, particle.depth, particle.lat, particle.lon]

        A = fieldset.cell_areas[time, 0, particle.lat, particle.lon]
        Vh = fieldset.Cs * A * math.sqrt(dudx**2 + 0.5*(dudy + dvdx)**2 + dvdy**2)

        r = 1/3.
        xres = fieldset.resolutionx # [degrees] 
        yres = fieldset.resolutiony
        dx_cell = fieldset.cell_edge_sizes_x[0, 0, particle.lat, particle.lon] # [meters]
        dy_cell = fieldset.cell_edge_sizes_y[0, 0, particle.lat, particle.lon]

        hitBoundary = True
        tries = 0
        while hitBoundary and tries <10:
            dlat = yres * random.uniform(-1., 1.) * math.sqrt(2*math.fabs(particle.dt)* Vh/r) / dy_cell
            dlon = xres * random.uniform(-1., 1.) * math.sqrt(2*math.fabs(particle.dt)* Vh/r) / dx_cell
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
            if tries == 10:
                dlat = 0
                dlon = 0

        particle.lat += dlat
        particle.lon += dlon            
        
        
