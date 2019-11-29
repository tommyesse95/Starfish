#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 13:29:05 2019

@author: tommyesse
"""

import numpy as np 
from decimal import Decimal
import math as mt 

# from ev to K 
def ev_to_K(T):
    kb = 8.611733e-5 #Boltzamnn Constant
    T_K = T / kb
    return T_K


## Eletrons simulated Per second. 


def mdotelec_to_curr(mdot_elec):
    mass_elec = 9.1909e-31 # kg 
    mdot_elec *= 1e-6 # conversion in kg from mg 
    num_elec = mdot_elec / mass_elec
    elem_charge = 1.6e-19
    curr = num_elec * elem_charge
    return curr 


# From mdot of neutral gas Kr, particles per dt 
def mdot_to_kr_dt(mdot,spwt):
    mKr = 83.798 #u of Krypton
    conv = 1.66054e-27# 1u in kg
    mKr_kg = mKr * conv 
    mdot *= 1e-6 # conversion in kg from mg 
    part_s = mdot / mKr_kg
    part_s_sim = part_s / spwt
    part_dt = part_s_sim * dt
    return part_dt

# Fot the neutral gas I can alos fix the pressure. See discharge code for 
# computation 
    

# current --> number of electrons inside per dt. 
def curr_to_el_dt(curr,spwt):
    elem_charge = 1.61e-19
    num_el_s = curr / elem_charge
    num_el_sim_s = num_el_s / spwt
    num_el_sim_dt = str(num_el_sim_s * dt)
    return num_el_sim_dt


def num_nodes(domain_size, spacing):
    nodesx= domain_size[0] / spacing[0]
    nodesy = domain_size[1] / spacing[1]
    return (nodesx, nodesy)

domain_size = np.array([140e-3,15e-3]) # domain size in mm 
spacing = np.array([5e-4, 5e-4]) # size of one cell
min_nodesx, min_nodesy = num_nodes(domain_size, spacing)


mdot_kr = 1 # mg/s
curr_el = 2 # A 
spwt_kr = 1e10
spwt_el = 5e10
dt = 5e-10
num_it = 3000000
time_sim = dt * num_it
 

print(f"Number of neutrals per dt= {dt} : " '%.2E' % Decimal(mdot_to_kr_dt(mdot_kr,spwt_kr)))
print(f"Number of electrons per dt= {dt} : " '%.2E' % Decimal(curr_to_el_dt(curr_el,spwt_el)))




# Drift of electrons in the orifice

ne = 1e20 # density electrons in the orifice in m -3 
I = 10 # discharge current 
r_orifice = 2e-3 # mm 
A = mt.pi * r_orifice**2
elem_charge = 1.61e-19
drift_v = I / (ne* elem_charge * A) # drift velocity of electrons 

# THe thermal velocity is much bigger. 

print(f" Drift velocity of electrons in the orifice: {drift_v}")






    


    
    