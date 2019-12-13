#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 13:29:05 2019

@author: tommyesse
"""

import numpy as np 
from decimal import Decimal
import math as mt 


# NOTES 
# Max number of particles should be 10^5. 
# ni = 10 ^18, ne = 10^18, nn related to pressure 

kb = 1.380649e-23 #Boltzamnn Constant
mass_elec = 9.1909e-31 # kg 
elem_charge = 1.6e-19 # elementary charge
mKr = 83.798 #u of Krypton
conv = 1.66054e-27 # 1u in kg
mKr_kg = mKr * conv
eps0 = 8.851e-12
T_neutr = 273
r_cyl = 15e-3
h_cyl = 38e-3
d_anode = 34e-3
vol_cyl = mt.pi * r_cyl**2 * h_cyl

dist_inl_anode = mt.sqrt(r_cyl**2 + d_anode**2)



# from ev to K 
def ev_to_K(T):
    T_K = T * 1.160451812e4
    return T_K

# from pressure (pa) to density (# / m^3)
def pressure_to_density(pressure,kb,T_neutr):
    density = pressure/(kb*T_neutr)
    return density

def density_to_pressure(density,kb,T_neutr):
    pressure = density * kb * T_neutr
    return pressure 


def therm_velocity(T, mass):
    vel = mt.sqrt((8 * kb * T)/(mass*mt.pi))
    return vel 


def neutral_density(mdot, r_out):
    v_therm = therm_velocity(T_neutr, mKr_kg)
    A = mt.pi * r_out**2
    n = (4 * mdot) / (A * v_therm * mKr_kg)
    part = n * vol_cyl
    pres = density_to_pressure(n,kb,T_neutr)
    pd = pres * dist_inl_anode
    pd_Torr = pd /(133.322 * 1e-2)
    print(f"Density: {n}")
    print("Particles:") 
    print('%.2E' % Decimal(str(part))) 
    print(f"Pressure: {pres} Pa")
    print(f"pd (Torr cm): {pd_Torr}")
    return n, part, pres, pd_Torr

    
# Debye length
def debye(Te,n0):
    debye = np.sqrt((eps0*kb*ev_to_K(Te))/(n0*elem_charge**2))
    return debye


def num_nodes(domain_size, spacing):
    nodesx= domain_size[0] / spacing[0] + 1
    nodesy = domain_size[1] / spacing[1] + 1 
    return (nodesx, nodesy)


def drift_ele(ne,I,r_or): # ne: density electrons in the orifice
    A = mt.pi * r_or**2
    drift_v = I / (ne* elem_charge * A) # drift velocity of electrons 
    print(f" Drift velocity of electrons in the orifice: {drift_v}")


# time simulation 
def num_it(dt,time_sim):
    num_it = time_sim/dt # Time simulation
    return num_it

def time_sim(dt, n_it):
    time_sim = dt * n_it
    return time_sim


# Try maybe to invert the process to extrapolare the mdot necessary
def neutral_dt(mdot,dt,num_it):
    print("Part per dt", "Tot part", "Average pressure", "Average density")
    mKr_kg = mKr * conv 
    mdot *= 1e-6 # conversion in kg from mg 
    part_s = mdot / mKr_kg
    part_dt = part_s * dt
    return part_dt


# current --> number of electrons inside per dt. IT WORKS 
def el(curr,dt):
    num_el_s = curr / elem_charge
    num_el_sim_dt = num_el_s * dt  
    return num_el_sim_dt
    
def dt_el(ne):
    omega = mt.sqrt((ne*elem_charge**2)/(mass_elec * eps0))
    dt = 1/omega
    return dt











    


    
    