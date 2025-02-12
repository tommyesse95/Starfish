# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 12:47:27 2019

@author: Sander
"""

import numpy as np
from numpy import linalg as LA
from numpy import genfromtxt
from scipy import interpolate
import random

epsilon_0 = 8.8e-12
A_0 = 1.2e6
Lambda_R = 0.5
W_emitter = 4.5 #Tungsten
electron_mass = 9.1e-31
elementary_charge = 1.6e-19
k_B = 1.38e-23
k_B_eV = 8.617e-5
A_material = 15000 #to be checked
phi_s = 7
avogadro = 6.022e23
U_plus = 14.1 # to be checked 
neutral_density = 10e22
work_function = 1.56
Omega = 5.6e-8 #Tungsten\
T_inf = 300
SB = 5.67e-8
epsilon_emitter = 0.35
V_drop = 60
mfr = 1e-6 #kg/s
mass_Kr = 84*1.6e-27
test_filename = "/home/sander/Desktop/PythoN_Data/0_ygrid_electric_field_20.csv"
A_Kr = 16#to be checked
B_Kr = 240 #to be checked
A = A_Kr
B = B_Kr
mass_Atom = mass_Kr
electron_max_cx_tot = 3.7e-19 #Maximum total cross section for the interactions involved.
E_exc = 9.96*elementary_charge
E_ion = 13.9996*elementary_charge

def csv_reader(filename):
    data = genfromtxt(filename,delimiter=',')
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

def electron_elas_scat_cx(electron_energy):
    ans = 0.0
    electron_energy = electron_energy/elementary_charge
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
    if (electron_energy <= 0.01337):
        ans = 0.0
    elif (electron_energy >= 408.8):
        ans = 0.0
    else:
       ans = np.interp(electron_energy,cx[:,0],cx[:,1])/electron_max_cx_tot 
    return ans

def electron_exc_rel_cx(electron_energy):
    ans = 0.0
    electron_energy = electron_energy/elementary_charge
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
    if (electron_energy <= 9.969):
        ans = 0.0
    elif (electron_energy >= 197.28):
        ans = 0.0
    else:
       ans = np.interp(electron_energy,cx[:,0],cx[:,1])/electron_max_cx_tot 
    return ans

def electron_ion_rel_cx(electron_energy):
    ans = 0.0
    electron_energy = electron_energy/elementary_charge
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
    if (electron_energy <= 15.8491):
        ans = 0.0
    elif (electron_energy >= 992.450):
        ans = 0.0
    else:
       ans = np.interp(electron_energy,cx[:,0],cx[:,1])/electron_max_cx_tot 
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
    if (electron_energy <= 2.7879e-20):
        ans = 0.0
    elif (electron_energy >= 2.4033e-16):
        ans = 0.0
    else:
        ans = np.interp(electron_energy,cx[:,0],cx[:,1]) 
    return ans

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
    density = pressure/300/8.314*avogadro
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
    

def e_ionization(particle1,swarm1): #To be Checked
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
    return particle2

def energy_tester(test_particle):
    v_the = temp_to_vel(2000,electron_mass)
    vel = np.array([Maxwel_sample(v_the),np.abs(Maxwel_sample(v_the)),Maxwel_sample(v_the)])
    test_particle.new_velocity(vel)
    return 0.5*test_particle.velSquared()*electron_mass/elementary_charge
    
def main():
    number_Of_Sims = 10
    rod_position = np.array([0.0,0.0,0.066])#xyz
    rod_radius = 0.0005
    T_rod = 2000
    v_thermal = temp_to_vel(T_rod,electron_mass)
    delta_t = 1e-11
    #Voltage = 200
    #E_field = csv_reader(test_filename)
    #fx = E_field_interpolator2(Voltage,0)
    #fy = E_field_interpolator2(Voltage,1)
    #fz = E_field_interpolator2(Voltage,2)
    number_Of_Iter = 10
    
    E_field = csv_reader(test_filename)
    E_field = E_field.clip(-1e5,1e5)
    
    for i in range(number_Of_Sims):
        #Determining the position of the particle and initialisation
        n0 = pressure_to_density(1)
        new_part = 0
        random_angle = (random.random()*np.pi)
        shift_of_particle = np.array([0.0,rod_radius*(1.01)*np.sin(random_angle),rod_radius*(1.1)*np.cos(random_angle)])
        starting_position = rod_position+shift_of_particle
        electron_swarm = swarm("electrons")
        electron = particle(starting_position)
        v_thermal = temp_to_vel(2000,electron_mass)    
        #Determine Particle Velocity 
        vel = np.array([Maxwel_sample(v_thermal),np.abs(Maxwel_sample(v_thermal)),Maxwel_sample(v_thermal)])
        electron.new_velocity(vel)
        electron_swarm.add_particle(electron)
        out_of_bounds = False  
        time_in_domain = 0.0
        for k in range(number_Of_Iter):#Sloppy
            for p in range(electron_swarm.size()):
                time_in_domain = 0.0
                out_of_bounds = False
                out_of_bounds = position_checker(electron_swarm.particle_Array[p])
                while not out_of_bounds: 
                    pusher(electron_swarm.particle_Array[p],delta_t,electron_mass,-1.0*elementary_charge,E_field)
                    Energy_part = 0.5*electron_swarm.particle_Array[p].velSquared()*electron_mass
                    out_of_bounds = position_checker(electron_swarm.particle_Array[p])#                    if(does_it_collide(n0,delta_t,Energy_part)):
                    if(does_it_collide(n0,delta_t,Energy_part)):
                        rand = random.random()
                        if(rand < electron_exc_rel_cx(Energy_part)):
                            #e_excitation(electron_swarm.particle_Array[p])
                            print("excitations")                            
                        elif(rand <(electron_exc_rel_cx(Energy_part)+electron_ion_rel_cx(Energy_part))):    
                            #e_ionization(electron_swarm.particle_Array[p],electron_swarm)
                            new_part += 1
                            print("ionizations")    
                        #elif(rand <(electron_exc_rel_cx(Energy_part)
                        #+electron_ion_rel_cx(Energy_part)
                        #+electron_elas_scat_cx(Energy_part))):
                            #e_elastic(electron_swarm.particle_Array[p])
                            #print("scatterings")
                            #print(time_in_domain)
                            #print(electron_swarm.particle_Array[p].energy())
                    time_in_domain += delta_t  
            out_of_bounds = False
            print(new_part)
    return electron_swarm

def I_th(Temperature,surface):
    ans =surface*A_0*Lambda_R*Temperature**2.0*np.exp(-W_emitter*elementary_charge/(k_B*Temperature))
    return ans

def surface(length, radius):
    ans = 2*np.pi*radius*length
    return ans

def wire_resistance(length,radius):
    ans = Omega*length/(np.pi*radius**2.0)
    return ans

def heat_deposit(length,radius,wire_current):
    ans = wire_resistance(length,radius)*wire_current**2.0
    return ans

def surface_temperature(length,radius,wire_current):
    ans = np.sqrt(np.sqrt(heat_deposit(length,radius,wire_current)/(surface(length, radius)*SB*epsilon_emitter)+T_inf**4.0))
    return ans