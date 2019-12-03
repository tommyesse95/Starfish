# -*- coding: utf-8 -*-
"""
Created on Fri Nov  15 12:47:27 2019

@author: Tommaso
"""

# Test code. Constant Electric Field between two plates

import numpy as np
from numpy import linalg as LA
from scipy import interpolate
import random
import cross_sections as cr

# Definition of constants

epsilon_0 = 8.8e-12 
mass_Kr = 84*1.6e-27 # kg 
elementary_charge = 1.6e-19
electron_max_cx_tot = 3.7e-19 #Maximum total cross section for the interactions involved.
E_exc = 9.96*elementary_charge
E_ion = 13.9996*elementary_charge
T_neutr = 300
kb = 1.380649e-23
electron_mass = 9.1e-31
electron_max_cx_tot = 3.7e-19
T_rod = 2000 # tempreature wire
A_kr = 16
B_kr = 240
sec_coeff = 0.1 # number of secondary electrons produced per incident positive ion




def temp_to_vel(temperature,mass_particle):
    v_th = np.sqrt(2*k_B*temperature/mass_particle)
    return v_th



class particle:
        position = np.array([0.0,0.0,0.0])
        velocity = np.array([0.0,0.0,0.0])

        def __init__(self,position):
            self.position = position

        def new_velocity(self,new_velocity):
            self.velocity = new_velocity

        def new_position(self,new_position):
            self.position = new_position

        def velSquared(self):
            return self.velocity.dot(self.velocity)

        def energy(self):
            return 0.5*self.velocity.dot(self.velocity)*electron_mass/elementary_charge


class swarm:

    particle_Array = []

    def __init__(self,name):
        self.name = name

    def add_particle(self,particle):
        self.particle_Array.append(particle)

    def size(self):
        return len(self.particle_Array)


def Maxwel_sample(v_th):
    sum_array = [random.random(),random.random(),random.random()]
    f_M = 2*(np.sum(sum_array)-1.5)
    return v_th*f_M


def position_checker(particle, d_plate):
    if np.abs(particle.position[1]) > d_plate:
        return True
    elif (particle.position[2] < 0.):
        return True
    elif particle.position[1] < 0:
        return True
    else:
        return False

def pusher(particle,delta_t,mass_electron,charge_electron,E_field):
    vel_old = particle.velocity
    pos_old = particle.position
    vel_new = np.zeros(3)
    pos_new = np.zeros(3)
    vel_new[0] = -elementary_charge*E_field[0]/mass_electron*delta_t+vel_old[0]
    vel_new[1] = -elementary_charge*E_field[1]/mass_electron*delta_t+vel_old[1]
    vel_new[2] = -elementary_charge*E_field[2]/mass_electron*delta_t+vel_old[2]
    pos_new[0] = pos_old[0]+vel_new[0]*delta_t
    pos_new[1] = pos_old[1]+vel_new[1]*delta_t
    pos_new[2] = pos_old[2]+vel_new[2]*delta_t
    particle.new_position(pos_new)
    particle.new_velocity(vel_new)
    

def does_it_collide(neutral_density,delta_x,electron_energy):
    probability = 1 - np.exp(-neutral_density*electron_max_cx_tot*delta_x)
    collides = (random.random() < probability)
    return collides


def pressure_to_density(pressure,kb,T_neutr):
    density = pressure/(kb*T_neutr)
    return density


def elastic_vel_update(particle1,rel_vel_mod,v_center_of_mass):
    new_vel = np.zeros(3)
    new_vel2 = np.zeros(3)
    cosXi = 2*random.random()-1
    sinXi = np.sqrt(1-cosXi**2)
    eta = (2*np.pi*random.random())
    reduced_mass = mass_Kr/(electron_mass+mass_Kr)
    new_vel[0] = rel_vel_mod*cosXi*reduced_mass + v_center_of_mass[0]
    new_vel[1] = rel_vel_mod*sinXi*np.cos(eta)*reduced_mass + v_center_of_mass[1]
    new_vel[2] = rel_vel_mod*sinXi*np.sin(eta)*reduced_mass + v_center_of_mass[2]
    particle1.new_velocity(new_vel)
    new_vel2[0] = v_center_of_mass[0] - rel_vel_mod*cosXi*reduced_mass
    new_vel2[1] = v_center_of_mass[1] - rel_vel_mod*sinXi*np.cos(eta)*reduced_mass
    new_vel2[2] = v_center_of_mass[2] - rel_vel_mod*sinXi*np.sin(eta)*reduced_mass
    return new_vel2


def e_excitation(particle1):
    v_the = temp_to_vel(300,mass_Kr)
    v2 = np.array([Maxwel_sample(v_the),np.abs(Maxwel_sample(v_the)),Maxwel_sample(v_the)])
    red_mass = electron_mass*mass_Kr/(electron_mass+mass_Kr)
    rel_vel_mod = np.sqrt(np.abs(LA.norm((particle1.velocity-v2))**2-2*E_exc/red_mass))
    v_center_of_mass = (particle1.velocity*electron_mass+v2*mass_Kr)/(electron_mass+mass_Kr)
    elastic_vel_update(particle1,rel_vel_mod,v_center_of_mass)


def e_elastic(particle1):
    v_the = temp_to_vel(300,mass_Kr)
    v2 = np.array([Maxwel_sample(v_the),np.abs(Maxwel_sample(v_the)),Maxwel_sample(v_the)])
    rel_vel = particle1.velocity-v2
    v_center_of_mass = (particle1.velocity*electron_mass+v2*mass_Kr)/(electron_mass+mass_Kr)
    rel_vel_mod = LA.norm(rel_vel)
    elastic_vel_update(particle1,rel_vel_mod,v_center_of_mass)



def e_ionization(particle1,swarm1):
    v_the = temp_to_vel(300,mass_Kr)
    v2 = np.array([Maxwel_sample(v_the),np.abs(Maxwel_sample(v_the)),Maxwel_sample(v_the)])
    rel_vel = particle1.velocity-v2
    red_mass = electron_mass*mass_Kr/(electron_mass+mass_Kr)
    col_energy = 0.5*red_mass*LA.norm(rel_vel)**2-E_ion
    X = 0.5+np.sqrt(0.25-0.25*random.random())
    col_energy_new = col_energy*X
    col_energy_rest = col_energy - col_energy_new
    rel_vel_new = np.sqrt(2*col_energy_new/red_mass)
    rel_vel_rest = np.sqrt(2*col_energy_rest/red_mass)
    v_center_of_mass = (particle1.velocity*electron_mass+v2*mass_Kr)/(electron_mass+mass_Kr)
    particle2 = particle(particle1.position)
    particle2.new_velocity(elastic_vel_update(particle1,rel_vel_new,v_center_of_mass))
    v_center_of_mass2 = particle2.velocity
    elastic_vel_update(particle2,rel_vel_rest,v_center_of_mass2)
    swarm1.add_particle(particle2)



def energy_tester(test_particle):
    v_the = temp_to_vel(2000,electron_mass)
    vel = np.array([Maxwel_sample(v_the),np.abs(Maxwel_sample(v_the)),Maxwel_sample(v_the)])
    test_particle.new_velocity(vel)
    return 0.5*test_particle.velSquared()*electron_mass/elementary_charge

def pressure_pd(pd,distance):
    pressure = pd_si / distance
    return pressure

def distance_pd(pd,pressure):
    distance = pd_si / pressure
    return distance



number_of_Electrons = 100  
num_of_iter = 10 #tba  
v_thermal = temp_to_vel(T_rod,electron_mass)
delta_t = 1e-12 # tba I set this because with those data we were in that part of curve. 
d_plate = 2e-3 # tba I set this, but maybe with lower voltage it takes a while
pd = np.array([1.5, 2, 2.5, 3.0, 3.5, 4.0, 5]) # tba 
# Voltage = np.linspace(250,550,10)
Voltage = np.arange(250,600,10) # tba 

for i in range(number_of_Electrons): # tba same swarm of elec for pd?
    
    # Initialization swarm electrons 
    random_angle = (random.random()*np.pi) # Angle of departure for electorns
    starting_position = np.array([0.0, 0.0, 0.0])
    electron_swarm = swarm("electrons")
    electron = particle(starting_position)
    vel = np.array([Maxwel_sample(v_thermal),np.abs(Maxwel_sample(v_thermal)),Maxwel_sample(v_thermal)])
    electron.new_velocity(vel)
    electron_swarm.add_particle(electron)

for j in range(len(pd)):
    print(f"\n\n Simulation {j+1} with pd = {pd[j]} Torr cm \n",file=open("output.txt", "a"))
    pd_si = pd[j] * 133.322 * 1e-2 # conversion in Pa*m
    pressure = pressure_pd(pd_si, d_plate) # Pressure in Pascal considering pd
    n0 = pressure_to_density(pressure,kb,T_neutr) #density of the neutral gas
    sec_emission = 0

#  V_br = (B_kr*pd_si) / (np.log(A_kr*pd_si) - np.log(np.log(1+ 1/sec_coeff)))

# Check this before simulation 

#    lambda_free = 1 / (electron_max_cx_tot * n0) # mean free path
#    est_in_vel = np.sqrt((2*Voltage[j]*elementary_charge)/(electron_mass)) #Velocity for conservation energy
#    deltat_free = lambda_free/est_in_vel # delta t needed for distance lamda free
#    factor_distance = d_plate[j] / lambda_free
#    print(f"dt for mean free path: {deltat_free}")
                
            
    for h in range(len(Voltage)):
        print(f"Voltage: {Voltage[h]}\n",file=open("output.txt", "a"))
        if sec_emission > number_of_Electrons :
            break
        Ey = -Voltage[j] / d_plate[j]
        E_field = np.array([0, Ey, 0])
        # Initialization 
        num_excitations = 0 
        num_ionizations = 0 
        num_scatterings = 0
        no_collisions = 0
        nothing = 0
        time_in_domain = 0
        out_of_bounds = False
    
        
    
    
        for k in range(num_of_iter):
            sec_emission = round(sec_coeff * num_ionizations)
            if sec_emission > number_of_Electrons[j]:
                time_break = time_in_domain * delta_t
                V_break = Voltage[h]
                print("Breakdown! after {time_break} s\n")
                break 
            for p in range(electron_swarm.size()):
                out_of_bounds = position_checker(electron_swarm.particle_Array[p], d_plate)
                while not out_of_bounds:
                    pusher(electron_swarm.particle_Array[p],delta_t,electron_mass,-1.0*elementary_charge,E_field)
                    Energy_part = 0.5*electron_swarm.particle_Array[p].velSquared()*electron_mass
                    Energy_electron = Energy_part / elementary_charge
                    out_of_bounds = position_checker(electron_swarm.particle_Array[p], d_plate)
                    delta_x = delta_t * np.sqrt(electron_swarm.particle_Array[p].velSquared())
                    if(does_it_collide(n0,delta_x,Energy_part)):
                        rand = random.random()
                        if(rand < cr.electron_exc_rel_cx(Energy_electron)):
                            e_excitation(electron_swarm.particle_Array[p])
                            num_excitations += 1
                            print(f"Electron {p+1}: excitation")
                        elif(rand <(cr.electron_exc_rel_cx(Energy_electron)+cr.electron_ion_rel_cx(Energy_electron))):
                            e_ionization(electron_swarm.particle_Array[p],electron_swarm)
                            num_ionizations += 1
                            print(f"Electron {p+1}: IONIZATION}")
                        elif(rand <(cr.electron_exc_rel_cx(Energy_electron)+cr.electron_ion_rel_cx(Energy_electron)+cr.electron_elas_scat_cx(Energy_electron))):
                            e_elastic(electron_swarm.particle_Array[p])
                            num_scatterings += 1
        #                   print(f"Electron {p+1}: scattering with energy: {Energy_part}")
                        else:
                            nothing += 1
                    else:
                        print(f"Electron {p}: No collisions")
                        no_collisions += 1
    
                    time_in_domain +=1
    
        sec_emission = round(sec_coeff * num_ionizations)
    
        if sec_emission > number_of_Electrons:
            print(f"Breakdown! at {time_break} s with {V_break} V", file=open("output.txt", "a"))
            print(f"Pressure neutral gas: {pressure} Pa", file=open("output.txt", "a"))
            print(f"Density neutral gas: {n0} per m\u00b3", file=open("output.txt", "a"))
            print(f"pd = {pd_si} Pa m", file=open("output.txt", "a"))
            
        else:
            print("No Breakdown", file=open("output.txt", "a"))
            
        
                
       
        
        
       
        



