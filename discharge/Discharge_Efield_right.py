#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 17:33:33 2019

@author: tommyesse
"""

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
from general_data import * 
import cross_sections as cr


def csv_reader(filename):
    data = np.genfromtxt(filename,delimiter=',')
    return data

def E_field_interpolator(E_field,coord):
    f = interpolate.interp2d(E_field[:,1],E_field[:,2],E_field[:,(3+coord)],kind='linear')
    return f

def E_field_interpolator3(E_field,coord,position):
    position_y = position[0]
    position_z = position[1]
    if(np.abs(position_y) > 0.0144):
        if (position_y > 0):
            position_y = 0.0144
        else:
            position_y = -0.0144
    f = interpolate.griddata(E_field[:,1:3],E_field[:,(3+coord)],(np.abs(position_y),position_z),method='linear')
    if(position_y < 0.0 and coord == 1):
        f = -f       
    if(f > 1e5):
        f = 1e5
    elif(f < -1e5):
        f = -1e5
    return f


def temp_to_vel(temperature,mass_particle):
    v_th = np.sqrt(0.5*k_B*temperature/mass_particle)
    return v_th


def E_field_interpolator2(Voltage_diff,coord):
        anode_position = np.array([0.0,1.5e-3,0.0])
        if (coord == 0):
            f = 0
        elif(coord == 1):
            f = Voltage_diff/(anode_position[1])
        elif(coord == 2):
            f = 0
        return f

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

#NOte this needs to be updated in case of actual electric field
def position_checker(particle):
    if np.abs(particle.position[1]) > 0.015:
        return True
    elif (particle.position[2] < 0.00 or particle.position[2] > 0.1):
        return True
    elif (np.abs(particle.position[1]) > 0.04):
        return True
    else:
        return False

def pusher(particle,delta_t,mass_electron,charge_electron,E_field):
    vel_old = particle.velocity
    pos_old = particle.position
    vel_new = np.zeros(3)
    pos_new = np.zeros(3)
    vel_new[0] = -elementary_charge*E_field_interpolator3(E_field,0,(particle.position[1],particle.position[2]))/mass_electron*delta_t+vel_old[0]
    vel_new[1] = -elementary_charge*E_field_interpolator3(E_field,1,(particle.position[1],particle.position[2]))/mass_electron*delta_t+vel_old[1] 
    vel_new[2] = -elementary_charge*E_field_interpolator3(E_field,2,(particle.position[1],particle.position[2]))/mass_electron*delta_t+vel_old[2]

    pos_new[0] = pos_old[0]+vel_new[0]*delta_t
    pos_new[1] = pos_old[1]+vel_new[1]*delta_t
    pos_new[2] = pos_old[2]+vel_new[2]*delta_t
    particle.new_position(pos_new)
    particle.new_velocity(vel_new)
    
def does_it_collide(neutral_density,delta_t,electron_energy):
    # probability = 1 - np.exp(-neutral_density*electron_tot_cx(electron_energy)*delta_t)
    probability = 1 - np.exp(-neutral_density*electron_max_cx_tot*delta_t)
    collides = (random.random() > probability)
    return collides

def pressure_to_density(pressure):
    T_neutr = 300
    density = pressure/T_neutr/8.314*avogadro
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
    rel_vel_mod = LA.norm((particle1.velocity-v2))**2-2*E_exc/red_mass
    if rel_vel_mod < 0.0:
        print("PROBLEMS")
    rel_vel_mod = np.sqrt(rel_vel_mod)
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
    
def main():
    print("Simulation Started \n")
    number_Of_Electrons = 10 
    num_of_iter = 10
    rod_position = np.array([0.0,0.0,0.066]) #  xyz
    rod_radius = 0.0005 # radius wire 
    T_rod = 2000 # tempreature wire 
    v_thermal = temp_to_vel(T_rod,electron_mass)
    delta_t = 1e-11 # time step simulation 
    #Voltage = 200
    #fx = E_field_interpolator2(Voltage,0)
    #fy = E_field_interpolator2(Voltage,1)
    #fz = E_field_interpolator2(Voltage,2)
    
#    number_Of_Iter = 10 # Number of iterations of the swarm 
    
    # Definition of the electric field 
    
    
    
    E_field = csv_reader(test_filename)
    E_field = E_field.clip(-1e5,1e5)  # Trick to limit the electric fiel 

    num_excitations = 0 # Inizialitation number of excitations 
    num_ionizations = 0 # Initialization number of ionizations
    no_collisions = 0
    nothing = 0
    n0 = pressure_to_density(1)
    out_of_bounds = False
    
    for i in range(number_Of_Electrons):
        
        print("Random Initialization of the swarm \n")
        random_angle = (random.random()*np.pi) # Angle of departure for electorns
        shift_of_particle = np.array([0.0,rod_radius*(1.01)*np.sin(random_angle),rod_radius*(1.1)*np.cos(random_angle)]) # tba! why those numbers?
        starting_position = rod_position+shift_of_particle
        electron_swarm = swarm("electrons") 
        electron = particle(starting_position) 
        v_thermal = temp_to_vel(T_rod,electron_mass) # Thermal velocity electrons
        vel = np.array([Maxwel_sample(v_thermal),np.abs(Maxwel_sample(v_thermal)),Maxwel_sample(v_thermal)])
        electron.new_velocity(vel)
        electron_swarm.add_particle(electron)
#        mainfor k in range(number_Of_Iter):
#            print(f"Number of iter: {k} out of bonds? {out_of_bounds}")
        
            
    print(f"Electron swarm size: {electron_swarm.size()}\n")
    for k in range(num_of_iter):
        for p in range(electron_swarm.size()):
            
            out_of_bounds = position_checker(electron_swarm.particle_Array[p])
            while not out_of_bounds: 
                pusher(electron_swarm.particle_Array[p],delta_t,electron_mass,-1.0*elementary_charge,E_field)
                Energy_part = 0.5*electron_swarm.particle_Array[p].velSquared()*electron_mass
                out_of_bounds = position_checker(electron_swarm.particle_Array[p])#                    if(does_it_collide(n0,delta_t,Energy_part)):
                if(does_it_collide(n0,delta_t,Energy_part)):
                    rand = random.random()
                    if(rand < cr.electron_exc_rel_cx(Energy_part)):
                        e_excitation(electron_swarm.particle_Array[p])
                        num_excitations += 1
                        print(f"Electron {p+1}: excitation")                            
                    elif(rand <(cr.electron_exc_rel_cx(Energy_part)+cr.electron_ion_rel_cx(Energy_part))):    
                        e_ionization(electron_swarm.particle_Array[p],electron_swarm)
                        # Here you wanna a particle back but actually nothing will change
                        
                        num_ionizations += 1
                        print(f"Electron {p+1}: ionization, Swarm increased? {electron_swarm.size()}")    
                            #elif(rand <(electron_exc_rel_cx(Energy_part)
                            #+electron_ion_rel_cx(Energy_part)
                            #+electron_elas_scat_cx(Energy_part))):
                                #e_elastic(electron_swarm.particle_Array[p])
                                #print("scatterings")
                                #print(time_in_domain)
                                #print(electron_swarm.particle_Array[p].energy())
                    else: 
                        nothing += 1
                else:
                    print(f"Electron {p}: No collisions")
                    no_collisions += 1
            
            
            
            
    print(f"\nTotal number of excitations: {num_excitations}")
    print(f"\nTotal number of ionizations: {num_ionizations}") 
    print(f"\nTotal number of collissions wihout any phenomena: {nothing}")
    print(f"\nTotal number of electrons without collisions: {no_collisions}")       
    

main()

#def I_th(Temperature,surface):
#    ans =surface*A_0*Lambda_R*Temperature**2.0*np.exp(-W_emitter*elementary_charge/(k_B*Temperature))
#    return ans
#
#def surface(length, radius):
#    ans = 2*np.pi*radius*length
#    return ans
#
#def wire_resistance(length,radius):
#    ans = Omega*length/(np.pi*radius**2.0)
#    return ans
#
#def heat_deposit(length,radius,wire_current):
#    ans = wire_resistance(length,radius)*wire_current**2.0
#    return ans
#
#def surface_temperature(length,radius,wire_current):
#    ans = np.sqrt(np.sqrt(heat_deposit(length,radius,wire_current)/(surface(length, radius)*SB*epsilon_emitter)+T_inf**4.0))
#    return ans