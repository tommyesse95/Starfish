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


#def csv_reader(filename):
#    data = np.genfromtxt(filename,delimiter=',')
#    return data


#def E_field_interpolator(E_field,coord):
#    f = interpolate.interp2d(E_field[:,1],E_field[:,2],E_field[:,(3+coord)],kind='linear')
#    return f


#def E_field_interpolator3(E_field,coord,position):
#    position_y = position[0]
#    position_z = position[1]
#    if(np.abs(position_y) > 0.0144):
#        if (position_y > 0):
#            position_y = 0.0144
#        else:
#            position_y = -0.0144
#    f = interpolate.griddata(E_field[:,1:3],E_field[:,(3+coord)],(np.abs(position_y),position_z),method='linear')
#    if(position_y < 0.0 and coord == 1):
#        f = -f
#    if(f > 1e5):
#        f = 1e5
#    elif(f < -1e5):
#        f = -1e5
#    return f



def temp_to_vel(temperature,mass_particle):
    v_th = np.sqrt(0.5*k_B*temperature/mass_particle)
    return v_th


#def E_field_interpolator2(Voltage_diff,coord):
#        anode_position = np.array([0.0,1.5e-3,0.0])
#        if (coord == 0):
#            f = 0
#        elif(coord == 1):
#            f = Voltage_diff/(anode_position[1])
#        elif(coord == 2):
#            f = 0
#        return f


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
#def position_checker(particle):
#    if np.abs(particle.position[1]) > 0.015:
#        return True
#    elif (particle.position[2] < 0.00 or particle.position[2] > 0.1):
#        return True
#    elif (np.abs(particle.position[1]) > 0.04):
#        return True
#    else:
#        return False


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
#    vel_new[0] = -elementary_charge*E_field_interpolator3(E_field,0,(particle.position[1],particle.position[2]))/mass_electron*delta_t+vel_old[0]
#    vel_new[1] = -elementary_charge*E_field_interpolator3(E_field,1,(particle.position[1],particle.position[2]))/mass_electron*delta_t+vel_old[1]
#    vel_new[2] = -elementary_charge*E_field_interpolator3(E_field,2,(particle.position[1],particle.position[2]))/mass_electron*delta_t+vel_old[2]
    vel_new[0] = -elementary_charge*E_field[0]/mass_electron*delta_t+vel_old[0]
    vel_new[1] = -elementary_charge*E_field[1]/mass_electron*delta_t+vel_old[1]
    vel_new[2] = -elementary_charge*E_field[2]/mass_electron*delta_t+vel_old[2]

    pos_new[0] = pos_old[0]+vel_new[0]*delta_t
    pos_new[1] = pos_old[1]+vel_new[1]*delta_t
    pos_new[2] = pos_old[2]+vel_new[2]*delta_t
    particle.new_position(pos_new)
    particle.new_velocity(vel_new)

def does_it_collide(neutral_density,delta_x,electron_energy):
    # probability = 1 - np.exp(-neutral_density*electron_tot_cx(electron_energy)*delta_t)
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
    rel_vel_mod = np.abs(LA.norm((particle1.velocity-v2))**2-2*E_exc/red_mass)
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

def pressure_pd(pd,distance):
    pressure = pd_si / distance
    return pressure

def distance_pd(pd,pressure):
    distance = pd_si / pressure
    return distance


num_sim = 3

for j in range(num_sim):
    print(f"Simulation {j+1} Started \n")
    number_of_Electrons =  np.array([500,700,1000])
    num_of_iter = 10
    T_neutr = 300
    kb = 1.380649e-23
    electron_mass = 9.1e-31
    elementary_charge = 1.6e-19
    electron_max_cx_tot = 3.7e-19
    T_rod = 2000 # tempreature wire
    v_thermal = temp_to_vel(T_rod,electron_mass)
    delta_t = np.array([1e-12,1e-12,1e-12]) # time step simulation
    Voltage = ([250,400,500])
    pd = np.array([5, 5, 5]) # Torr*cm
    pd_si = pd[j] * 133.322 * 1e-2 # conversion in Pa*m
    d_plate = np.array([5e-4, 2e-3, 2e-2])  # distance between parallel plates
    pressure = pressure_pd(pd_si, d_plate[j]) # Pressure in Pascal considering pd
    n0 = pressure_to_density(pressure,kb,T_neutr) #density of the neutral gas

    lambda_free = 1 / (electron_max_cx_tot * n0) # mean free path
    est_in_vel = np.sqrt((2*Voltage[j]*elementary_charge)/(electron_mass)) #Velocity for conservation energy
    deltat_free = lambda_free/est_in_vel # delta t needed for distance lamda free
    factor_distance = d_plate[j] / lambda_free
    print(f"dt for mean free path: {deltat_free}")

    Ey = -Voltage[j] / d_plate[j]
    E_field = np.array([0, Ey, 0])

    A_kr = 16
    B_kr = 240
    sec_coeff = 0.1 # number of secondary electrons produced per incident positive ion
    M = 1 + 1/sec_coeff


    # V_br= (pd_si * B_kr**(2)) / (np.log(M/(A_kr*pd_si)))**(2) source

    # wikipedia source
    V_br = (B_kr*pd_si) / (np.log(A_kr*pd_si) - np.log(np.log(1+ 1/sec_coeff)))


    num_excitations = 0 # Inizialitation number of excitations
    num_ionizations = 0 # Initialization number of ionizations
    num_scatterings = 0
    no_collisions = 0
    nothing = 0
    time_in_domain = 0
    out_of_bounds = False

    for i in range(number_of_Electrons[j]):

        print("Random Initialization of the swarm \n")
        random_angle = (random.random()*np.pi) # Angle of departure for electorns
    #        shift_of_particle = np.array([0.0,rod_radius*(1.01)*np.sin(random_angle),rod_radius*(1.1)*np.cos(random_angle)]) # tba! why those numbers?
    #        starting_position = rod_position
        starting_position = np.array([0.0, 0.0, 0.0])
        electron_swarm = swarm("electrons")
        electron = particle(starting_position)
        v_thermal = temp_to_vel(T_rod,electron_mass) # Thermal velocity electrons
        vel = np.array([Maxwel_sample(v_thermal),np.abs(Maxwel_sample(v_thermal)),Maxwel_sample(v_thermal)])
        electron.new_velocity(vel)
        electron_swarm.add_particle(electron)


    for k in range(num_of_iter):
        sec_emission = round(sec_coeff * num_ionizations)
        if sec_emission > number_of_Electrons[j]:
            time = time_in_domain
            print("Breakdown! after {time_in_domain * delta_t[j]} s\n")
            break
        for p in range(electron_swarm.size()):
            out_of_bounds = position_checker(electron_swarm.particle_Array[p], d_plate[j])
            while not out_of_bounds:
                pusher(electron_swarm.particle_Array[p],delta_t[j],electron_mass,-1.0*elementary_charge,E_field)
                Energy_part = 0.5*electron_swarm.particle_Array[p].velSquared()*electron_mass
                Energy_electron = Energy_part / elementary_charge
                out_of_bounds = position_checker(electron_swarm.particle_Array[p], d_plate[j])
                delta_x = delta_t[j] * np.sqrt(electron_swarm.particle_Array[p].velSquared())
                if(does_it_collide(n0,delta_x,Energy_part)):
                    rand = random.random()
                    if(rand < cr.electron_exc_rel_cx(Energy_electron)):
                        e_excitation(electron_swarm.particle_Array[p])
                        num_excitations += 1
                        print(f"Electron {p+1}: excitation, with energy: {Energy_part}")
                    elif(rand <(cr.electron_exc_rel_cx(Energy_electron)+cr.electron_ion_rel_cx(Energy_electron))):
                        e_ionization(electron_swarm.particle_Array[p],electron_swarm)
                        num_ionizations += 1
                        print(f"Electron {p+1}: IONIZATION, with energy: {Energy_part}")
                    elif(rand <(cr.electron_exc_rel_cx(Energy_electron)+cr.electron_ion_rel_cx(Energy_electron)+cr.electron_elas_scat_cx(Energy_electron))):
                        e_elastic(electron_swarm.particle_Array[p])
                        num_scatterings += 1
    #                       print(f"Electron {p+1}: scattering with energy: {Energy_part}")
                    else:
                        nothing += 1
                else:
                    print(f"Electron {p}: No collisions")
                    no_collisions += 1

                time_in_domain +=1

    sec_emission = round(sec_coeff * num_ionizations)

    time_break = time_in_domain * delta_t[j]
    print(f"\n\nSimulation {j+1} with {number_of_Electrons[j]} electrons", file=open("output.txt", "a"))

    if sec_emission > number_of_Electrons[j]:
            print("Breakdown! at {time_break} s", file=open("output.txt", "a"))

    print(f"Pressure neutral gas: {pressure} Pa", file=open("output.txt", "a"))
    print(f"Density neutral gas: {n0} per m\u00b3", file=open("output.txt", "a"))
    print(f"Distance between plates: {d_plate[j]} m", file=open("output.txt", "a"))

    print(f"pd = {pd_si} Pa m", file=open("output.txt", "a"))

    
    print(f"\nTotal number of excitations: {num_excitations}", file=open("output.txt", "a"))
    print(f"Total number of IONIZATIONS: {num_ionizations}", file=open("output.txt", "a"))
    print(f"Total number of scatterings: {num_scatterings}", file=open("output.txt", "a"))
    print(f"Total number of collissions wihout any phenomena: {nothing}", file=open("output.txt", "a"))
    print(f"Total number of electrons without collisions: {no_collisions}", file=open("output.txt", "a"))
    print(f"Total number of secondary electrons: {sec_emission}", file=open("output.txt", "a"))
    print(f"Analytical breakdown voltage: {V_br} V", file=open("output.txt", "a"))
    print(f"Voltage used: {Voltage[j]} V", file=open("output.txt", "a"))








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
