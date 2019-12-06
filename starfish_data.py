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
eps0 = 8.851e-12
T_neutr = 293


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


## Eletrons simulated Per second. (from mdot)
def mdotelec_to_curr(mdot_elec):
    mdot_elec *= 1e-6 # conversion in kg from mg 
    num_elec = mdot_elec / mass_elec
    curr = num_elec * elem_charge
    return curr 





# Definition of domain 
    
# Debye lenght 
def debye(Te,n0):
    debye = np.sqrt((eps0*kb*ev_to_K(Te))/(n0*elem_charge**2))
    return debye

def domain():
    def num_nodes(domain_size, spacing):
        nodesx= domain_size[0] / spacing[0]
        nodesy = domain_size[1] / spacing[1]
        return (nodesx, nodesy)
    
    domain_size = np.array([140e-3,15e-3]) # domain size in mm 
    spacing = np.array([5e-4, 5e-4]) # size of one cell
    min_nodesx, min_nodesy = num_nodes(domain_size, spacing)






 
def drift_ele(ne,I,r_or): # ne: density electrons in the orifice
    A = mt.pi * r_or**2
    drift_v = I / (ne* elem_charge * A) # drift velocity of electrons 
    print(f" Drift velocity of electrons in the orifice: {drift_v}")


# Thermal velocity 
def temp_to_vel(temperature,mass_particle):
    v_th = np.sqrt(2*kb*temperature/mass_particle)
    return v_th



# NUmber of particles 
# Density from pressure. # TBC
    
def neutrals():
    print("Insert pressure in Pa:")
    pressure = int(input())
       
    print("\nInsert r cylinder (mm):")
    r_cyl = int(input())
    
    print("\nInsert h cylinder (mm):")
    h = int(input())
    
    print("\nInsert r inlet (mm):")
    r_inl = int(input())
    
    print("\nInsert vel injection (m/s):")
    u = int(input())
    
    vol = mt.pi * r_cyl**2 * h * 1e-9 # m^3 
    n = pressure/(kb*T_neutr) # #/m^3
    part = n * vol
    #A = mt.pi * r_inl**2 * 1e-6 
    A = 2 * mt.pi * r_cyl * r_inl * 1e-6
    mdot = mKr * n * u * A * conv
    
    print(f"Number of particles: {part}, spwt?")
    
    return part, n, mdot 


# time simulation 
def num_it(dt,time_sim):
    
    num_it = time_sim/dt # Time simulation
    return num_it



r_cyl = 15 
h = 45
vol = mt.pi * r_cyl**2 * h * 1e-9 # m^3 

# From mdot of neutral gas Kr, particles per dt. IT WORKS 

# Try maybe to invert the process to extrapolare the mdot necessary
def neutral(mdot,dt,num_it):
    print("Part per dt", "Tot part", "Average pressure", "Average density")
    mKr_kg = mKr * conv 
    mdot *= 1e-6 # conversion in kg from mg 
    part_s = mdot / mKr_kg
    part_dt = part_s * dt
    part_tot = part_dt * num_it 
    n = part_tot / vol
    p = density_to_pressure(n,kb,T_neutr)
    return part_dt, part_tot, p, n


# current --> number of electrons inside per dt. IT WORKS 
def el(curr,dt,num_it):
    print("tot_num_el", "density")
    num_el_s = curr / elem_charge
    num_el_sim_dt = num_el_s * dt
    tot_num_el = num_it * num_el_sim_dt
    n_el = tot_num_el / vol 
    return tot_num_el, n_el

# Density of electrons change as well as the one for the neutrals 

# For definition dt 
    
# def plasma_freq_el(ne,)










    


    
    