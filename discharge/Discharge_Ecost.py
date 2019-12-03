# -*- coding: utf-8 -*-
"""
Created on Fri Nov  15 12:47:27 2019

@author: Tommaso
"""

# Make sure to have also the file "cross_sections" 

# Test code. Constant Electric Field between two plates

import numpy as np
import random
from matplotlib import pyplot as plt 
import csv

# Definition of constants

epsilon_0 = 8.8e-12 
mass_Kr = 84*1.6e-27 # kg 
elementary_charge = 1.6e-19
electron_max_cx_tot = 3.7e-19 #Maximum total cross section for the interactions involved.
E_exc = 9.96*elementary_charge
E_ion = 13.9996*elementary_charge
T_neutr = 300
k_b = 1.380649e-23
electron_mass = 9.1e-31
electron_max_cx_tot = 3.7e-19
T_rod = 2000 # tempreature wire
A_kr = 16
B_kr = 240
sec_coeff = 0.1 # number of secondary electrons produced per incident positive ion
factor = 1.1 # For definition breakdown voltage 



def electron_elas_scat_cx(electron_energy):
    ans = 0.0
    cx = np.array([[0.01336406983303936, 3.181396515358613e-19],
                   [0.018270087215254787, 2.923940874194074e-19],
                   [0.028013768681598505, 2.4698476044385235e-19],
                   [0.037510023752653165, 2.1893446594321645e-19],
                   [0.05291504402947375, 1.849335501313526e-19],
                   [0.06868208597082541, 1.4885892182773155e-19],
                   [0.0910210354887281, 1.227454028542252e-19],
                   [0.1231903547660483, 8.864290474836367e-20],
                   [0.17938272211706383, 5.539307791537611e-20],
                   [0.21875937572384052, 4.0003154620811525e-20],
                   [0.31528330216067424, 2.2975054195115398e-20],
                   [0.39265920041240704, 1.4707513883285014e-20],
                   [0.47889612979455004, 1.0000000000000082e-20],
                   [0.5719757530891328, 7.221688363647997e-21],
                   [0.661895166661915, 6.5577359661024676e-21],
                   [0.7422313175233057, 6.79924566405827e-21],
                   [0.8409129724948268, 7.763407021146366e-21],
                   [1.1364808869109344, 1.4531273098790106e-20],
                   [0.9525240818926235, 1.0121283787935762e-20],
                   [1.3142129751297982, 2.1115789787161924e-20],
                   [1.7579208969728017, 3.676588514755182e-20],
                   [2.692014953486649, 7.221688363647998e-20],
                   [4.038025043534839, 1.2880942571685e-19],
                   [5.995678266734137, 1.964236230044633e-19],
                   [7.538332264704461, 2.325370435519632e-19],
                   [9.091848033840687, 2.5918661240103032e-19],
                   [11.673950992523647, 2.7199127560420315e-19],
                   [14.08332863822641, 2.560807678468362e-19],
                   [19.867221622875878, 2.163109646270179e-19],
                   [28.03109289406893, 1.6393047069070092e-19],
                   [40.803911015503004, 1.272658967139995e-19],
                   [58.169855573717584, 1.0121283787935762e-19],
                   [72.41280520634076, 8.758069292951211e-20],
                   [89.20883576358081, 7.578463300432568e-20],
                   [124.54849896355462, 6.174132084957152e-20],
                   [192.98252877276298, 5.030031220223962e-20],
                   [286.7967661319877, 4.2488583367628814e-20],
                   [408.81146099696934, 3.632531793187723e-20]])

    ans = np.interp(electron_energy,cx[:,0],cx[:,1],left=0,right=0)/electron_max_cx_tot 
    return ans

def electron_exc_rel_cx(electron_energy):
    ans = 0.0
    cx = np.array([[0,0],
                  [9.968803734054186, 1.2113071338026534e-22],
                  [9.975255644084465, 4.326748936771704e-22],
                  [9.996879949184952, 6.309371745042097e-22],
                  [10.184971794196398, 2.018561554192476e-21],
                  [10.369773189076172, 4.982876900465995e-22],
                  [10.685826239717441, 1.0105513611436902e-21],
                  [11.017087007198649, 1.5203454937762005e-21],
                  [11.454334188743996, 1.6185669425990547e-21],
                  [11.629901767702616, 2.048682469343214e-21],
                  [12.182648727313248, 2.5925436666318627e-21],
                  [12.760986772830162, 3.385552430688689e-21],
                  [13.263205383106026, 4.352413306118517e-21],
                  [14.55546632986207, 6.5452903527589604e-21],
                  [16.222361517354674, 1.08152741235494e-20],
                  [17.80625921355134, 1.4569820331078074e-20],
                  [19.399558272779824, 1.6001318583751074e-20],
                  [21.471389368401315, 1.599023934804686e-20],
                  [26.51141501787568, 1.499429807995171e-20],
                  [34.84372896501458, 1.4054396562213607e-20],
                  [47.247172608140104, 1.3170604566523676e-20],
                  [73.73400361583263, 1.1948975758556937e-20],
                  [105.60740702280589, 1.051134144672416e-20],
                  [144.34219134826114, 9.105387318438538e-21],
                  [197.28941207367293, 7.764485569375252e-21]])
    
    ans = np.interp(electron_energy,cx[:,0],cx[:,1],left=0,right=0)/electron_max_cx_tot 
    return ans

def electron_ion_rel_cx(electron_energy):
    ans = 0.0
    cx = np.array([[15.849, 0.045e-20],
                  [19.623, 1.264e-20],
                  [22.642, 1.614e-20],
                  [27.17, 2.246e-20],
                  [33.208, 2.867e-20],
                  [39.245, 3.183e-20],
                  [49.811, 3.431e-20],
                  [80, 3.442e-20],
                  [114.717, 3.194e-20],
                  [141.887, 3.002e-20],
                  [176.604, 2.788e-20],
                  [209.811, 2.596e-20],
                  [243.019, 2.415e-20],
                  [273.208, 2.269e-20],
                  [309.434, 2.111e-20],
                  [348.679, 1.986e-20],
                  [398.491, 1.828e-20],
                  [496.604, 1.591e-20],
                  [596.226, 1.4e-20],
                  [694.34, 1.264e-20],
                  [793.962, 1.14e-20],
                  [893.585, 1.061e-20],
                  [992.453, 0.971e-20]])
   
    ans = np.interp(electron_energy,cx[:,0],cx[:,1],left=0,right=0)/electron_max_cx_tot 
    return ans

def electron_tot_cx(electron_energy):
    ans = 0.0
    cx = np.array([[0,0],
                   [2.7878e-20,0],
                   [2.8118e-20,5.0174e-20],
                   [3.6049e-20,3.77e-20],
                   [4.422e-20,2.6556e-20],
                   [4.8065e-20,2.31e-20],
                   [5.5996e-20,1.7556e-20],
                   [6.4167e-20,1.2974e-20],
                   [7.2098e-20,1.04e-20],
                   [8.796e-20,6.5293e-21],
                   [1.0406e-19,4.9577e-21],
                   [1.1872e-19,4.4156e-21],
                   [1.442e-19,5.3e-21],
                   [1.6006e-19,6.7059e-21],
                   [1.603e-19,9.5845e-21],
                   [2.2422e-19,1.9274e-20],
                   [2.4033e-19,2.2375e-20],
                   [2.8046e-19,2.6408e-20],
                   [3.6049e-19,4.2349e-20],
                   [4.8065e-19,6.7374e-20],
                   [6.4071e-19,1.0435e-19],
                   [7.2098e-19,1.2361e-19],
                   [8.0101e-19,1.4384e-19],
                   [9.6131e-19,1.8308e-19],
                   [1.1216e-18,2.1927e-19],
                   [1.2817e-18,2.4885e-19],
                   [1.442e-18,2.642e-19],
                   [1.602e-18,2.7246e-19],
                   [1.7623e-18,2.7464e-19],
                   [1.9226e-18,2.7474e-19],
                   [2.0829e-18,2.6744e-19],
                   [2.2432e-18,2.6305e-19],
                   [2.4033e-18,2.5559e-19],
                   [2.5633e-18,2.4743e-19],
                   [2.6181e-18,2.4526e-19],
                   [2.8839e-18,2.348e-19],
                   [3.2043e-18,2.242e-19],
                   [3.2045e-18,1.3004e-21],
                   [9.4206e-18,5.8322e-22],
                   [2.243e-17,2.9466e-22],
                   [2.2431e-17,3.2375e-20],
                   [3.6049e-17,2.6771e-20],
                   [4.8066e-17,2.307e-20],
                   [6.4087e-17,1.964e-20],
                   [8.812e-17,1.608e-20],
                   [1.1215e-16,1.367e-20],
                   [1.442e-16,1.146e-20],
                   [1.7624e-16,9.87e-21],
                   [2.4033e-16,7.85e-21]])
    
    ans = np.interp(electron_energy,cx[:,0],cx[:,1],left=0,right=0) 
    return ans





def temp_to_vel(temperature,mass_particle):
    v_th = np.sqrt(2*k_b*temperature/mass_particle)
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
    
    def reinitialize(self):
        self.particle_Array = []


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
    rel_vel_mod = np.sqrt(np.abs(np.linalg.norm((particle1.velocity-v2))**2-2*E_exc/red_mass))
    v_center_of_mass = (particle1.velocity*electron_mass+v2*mass_Kr)/(electron_mass+mass_Kr)
    elastic_vel_update(particle1,rel_vel_mod,v_center_of_mass)


def e_elastic(particle1):
    v_the = temp_to_vel(300,mass_Kr)
    v2 = np.array([Maxwel_sample(v_the),np.abs(Maxwel_sample(v_the)),Maxwel_sample(v_the)])
    rel_vel = particle1.velocity-v2
    v_center_of_mass = (particle1.velocity*electron_mass+v2*mass_Kr)/(electron_mass+mass_Kr)
    rel_vel_mod = np.linalg.norm(rel_vel)
    elastic_vel_update(particle1,rel_vel_mod,v_center_of_mass)



def e_ionization(particle1,swarm1):
    v_the = temp_to_vel(300,mass_Kr)
    v2 = np.array([Maxwel_sample(v_the),np.abs(Maxwel_sample(v_the)),Maxwel_sample(v_the)])
    rel_vel = particle1.velocity-v2
    red_mass = electron_mass*mass_Kr/(electron_mass+mass_Kr)
    col_energy = 0.5*red_mass*np.linalg.norm(rel_vel)**2-E_ion
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



# SIMULATION 

# Notes! With 0.1 we don't reach the breakdown voltage. with 0.15 we reach 240
# V. The minimum value should be somewhere in between. 

number_of_Electrons = 100  
num_of_iter = 10 
v_thermal = temp_to_vel(T_rod,electron_mass) 
d_plate = 2e-3 # tba
pd = np.array([0.12, 0.15])# Torr cm  TBC
# Voltage = np.linspace(250,550,10)
Voltage = np.arange(200,800,10) # tba  
starting_position = np.array([0.0, 0.0, 0.0])
V_break = []




for j in range(len(pd)):
    print(f"\n\nSimulation {j+1} with pd = {pd[j]} Torr cm \n",file=open("output.txt", "a"))
    pd_si = pd[j] * 133.322 * 1e-2 # conversion in Pa*m
    pressure = pressure_pd(pd_si, d_plate) # Pressure in Pascal considering pd
    n0 = pressure_to_density(pressure,k_b,T_neutr) #density of the neutral gas
    breakdown = False

#  V_br = (B_kr*pd_si) / (np.log(A_kr*pd_si) - np.log(np.log(1+ 1/sec_coeff)))

# Check this before simulation 

   
                          
    for h in range(len(Voltage)):
        if breakdown: 
            break
        lambda_free = 1 / (electron_max_cx_tot * n0) # mean free path
        max_vel = np.sqrt((2*Voltage[h]*elementary_charge)/(electron_mass)) # Velocity for conservation energy
        deltat_free = lambda_free/max_vel # delta t needed for distance lamda free
        delta_t = 1e-1 * deltat_free
        print(f"Voltage: {Voltage[h]} with dt = {delta_t} s",file=open("output.txt", "a"))
        Ey = -Voltage[h] / d_plate
        E_field = np.array([0, Ey, 0])
        # Initialization variables for siimulation 
        num_excitations = 0 
        num_ionizations = 0 
        num_scatterings = 0
        no_collisions = 0
        nothing = 0
        time_in_domain = 0
        out_of_bounds = False
        electron_swarm = swarm("electrons")
        electron_swarm.reinitialize()
        
        for i in range(number_of_Electrons): 
            # Initialization swarm electrons 
            random_angle = (random.random()*np.pi) 
            electron = particle(starting_position)
            vel = np.array([Maxwel_sample(v_thermal),np.abs(Maxwel_sample(v_thermal)),Maxwel_sample(v_thermal)])
            electron.new_velocity(vel)
            electron_swarm.add_particle(electron)
    
        
        for k in range(num_of_iter):
            sec_emission = round(sec_coeff * num_ionizations)
            if sec_emission > factor * number_of_Electrons:
                time_break = time_in_domain * delta_t
                V_break.append(Voltage[h])
                breakdown = True 
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
                        if(rand < electron_exc_rel_cx(Energy_electron)):
                            e_excitation(electron_swarm.particle_Array[p])
                            num_excitations += 1
                            print(f"Electron {p+1}: excitation")
                        elif(rand <(electron_exc_rel_cx(Energy_electron)+electron_ion_rel_cx(Energy_electron))):
                            e_ionization(electron_swarm.particle_Array[p],electron_swarm)
                            num_ionizations += 1
                            print(f"Electron {p+1}: IONIZATION")
                        elif(rand <(electron_exc_rel_cx(Energy_electron)+electron_ion_rel_cx(Energy_electron)+electron_elas_scat_cx(Energy_electron))):
                            e_elastic(electron_swarm.particle_Array[p])
                            num_scatterings += 1
                        else:
                            nothing += 1
                    else:
                        no_collisions += 1
    
                    time_in_domain +=1
    
               
        if breakdown:
            print(f"Breakdown! at {time_break} s with {V_break[j]} V", file=open("output.txt", "a"))
            print(f"Pressure neutral gas: {pressure} Pa", file=open("output.txt", "a"))
            print(f"Density neutral gas: {n0} per m\u00b3", file=open("output.txt", "a"))
            print(f"pd = {pd[j]} Torr cm \n", file=open("output.txt", "a"))
            
        else:
            V_break.append('None')
            print(f"No Breakdown. num ionizations: {num_ionizations}\n", file=open("output.txt", "a"))
                
            


with open('data_sparc.csv', 'w', newline='') as f: # TBC CHange number file 
    thewriter = csv.writer(f)
    
    for i in range(len(pd)):
        thewriter.writerow([pd[i], V_break[i]])
        
## Plotting 
#            
#plt.figure(figsize=(8,5), dpi=100)
#
#plt.plot (pd, V_break, 'b^--')
#
#plt.xlabel('Pd (Torr cm)')
#
#plt.ylabel('Voltage (V)')
#
#plt.title('Breakdown voltage for Kr', fontdict={'fontname': 'Comic Sans MS', 'fontsize': 20})
#
#plt.grid(True)
#
## X, Y axis Tickmarks (scale of your graph)
#plt.xticks([0,1.5,2.5,3.5,4.5,5.5,10])
#
#plt.yticks([200,250,300,450,500,700])
#
## Save figure (dpi 300 is good when saving so graph has high resolution)
#plt.savefig('breakdown_voltage.png', dpi=300)
#
#plt.show()
       
        
        
       
        



