#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 13:29:05 2019

@author: tommyesse
"""

import numpy as np 
from decimal import Decimal
import math as mt 

kb = 1.380649e-23 #Boltzamnn Constant
mass_elec = 9.1909e-31 # kg 
elem_charge = 1.6e-19 # elementary charge
mKr = 83.798 #u of Krypton
conv = 1.66054e-27# 1u in kg
eps0 = 8.851e-12


# from ev to K 
def ev_to_K(T):
    T_K = T * 1.160451812e4
    return T_K


## Eletrons simulated Per second. 


def mdotelec_to_curr(mdot_elec):
    mdot_elec *= 1e-6 # conversion in kg from mg 
    num_elec = mdot_elec / mass_elec
    curr = num_elec * elem_charge
    return curr 


# From mdot of neutral gas Kr, particles per dt. IT WORKS 
def mdot_to_kr_dt(mdot,spwt,dt):
    mKr_kg = mKr * conv 
    mdot *= 1e-6 # conversion in kg from mg 
    part_s = mdot / mKr_kg
    part_s_sim = part_s / spwt
    part_dt = part_s_sim * dt
    return part_dt


# current --> number of electrons inside per dt. IT WORKS 
def curr_to_el_dt(curr,spwt_el,dt):
    num_el_s = curr / elem_charge
    num_el_sim_s = num_el_s / spwt_el
    num_el_sim_dt = str(num_el_sim_s * dt)
    return num_el_sim_dt


def num_nodes(domain_size, spacing):
    nodesx= domain_size[0] / spacing[0]
    nodesy = domain_size[1] / spacing[1]
    return (nodesx, nodesy)


# Domain 
domain_size = np.array([140e-3,15e-3]) # domain size in mm 
spacing = np.array([5e-4, 5e-4]) # size of one cell
min_nodesx, min_nodesy = num_nodes(domain_size, spacing)


dt = 5e-10
mdot_kr = 1 # mg/s
curr_el = 2 # A 
spwt_kr = 1e10 # Particles of Kr in one macroparticle 
spwt_el = 5e10 # Particles of electrons in one macroparticle 
num_it = 3000000 # number of iterations
time_sim = dt * num_it # Time simulation
 

# Number of electrons and neutrlas injected per dt 
def num_neut_ele():
    print(f"Number of neutrals per dt= {dt} : " '%.2E' % Decimal(mdot_to_kr_dt(mdot_kr,spwt_kr)))
    print(f"Number of electrons per dt= {dt} : " '%.2E' % Decimal(curr_to_el_dt(curr_el,spwt_el)))
    



 
def drift_ele(ne,I,r_or): # ne: density electrons in the orifice
    A = mt.pi * r_or**2
    drift_v = I / (ne* elem_charge * A) # drift velocity of electrons 
    print(f" Drift velocity of electrons in the orifice: {drift_v}")


# Thermal velocity 
def temp_to_vel(temperature,mass_particle):
    v_th = np.sqrt(2*kb*temperature/mass_particle)
    return v_th

# Density from pressure. # TBC
def pressure_to_particles(pressure,kb,T_neutr,spwt): 
    density = pressure/(kb*T_neutr) # #/m^3
    part = density/spwt 
    return part, density 



# Debye lenght 
def debye(Te,n0):
    debye = np.sqrt((eps0*kb*ev_to_K(Te))/(n0*elem_charge**2))
    return debye













    


    
    