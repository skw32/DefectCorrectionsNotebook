# SUZY EDITS, PLAN:
## Take functions from the top as separate script
## Incorporate 'main' section below into notebook using notebook parameters as function args
## (Replace command line args for bottom of script with notebook parameters also)
## Generate plot in notebook and also save .png into defects_output_dir
#################################################

#!/usr/bin/pythonw
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 13:29:30 2018

@author: tong
"""

import matplotlib.pyplot as plt
#import matplotlib as mpl
import random
import time
from scipy import optimize
#import matplotlib.colors
import numpy as np
#from scipy.optimize import curve_fit
#from atomic_V import atomic_3d   



import os
import sys
import math
import csv




def sys_basic(inputfile): 
    info = open(inputfile,'r')
    res = info.readline().split()
    defect_info = res
    res = info.readline()
    charge = float(res)
    res = info.readline()
    shift_H = float(res)
    res = info.readline()
    shift_D = float(res)
    return defect_info,charge, shift_H, shift_D   

def py_read(filename):
    data = np.load(filename)
    grid = np.shape(data)
    return grid,data   

# Suzy edits: use this function as is, but use notebook parameters as function args?
def geo_basic(hostF,defectF): 
    host = open(hostF, 'r') 
    defect = open(defectF,'r')
    host_len = len(host.readlines())
    defect_len = len(defect.readlines())
    atom_no = host_len - 3 
    defect_no = defect_len -3 
    return atom_no, defect_no 
    
# Suzy edits: use this function as is, but use notebook parameters as function args?
def geo_compare(atom_no,defect_no, hostF,defectF):
    host = open(hostF, 'r')
    defect = open(defectF,'r')
    lattice_vectors = np.zeros([3,3])
    host_corr = np.zeros([atom_no,3])
    defect_corr = np.zeros([defect_no,3])
    host_name = ['']*atom_no
    defect_name = ['']*defect_no
    for i in range(3): 
        res = host.readline().split()
        for j in range(3): 
            k=j+1
            lattice_vectors[i][j] = float(res[k]) 
    for i in range(atom_no): 
        res = host.readline().split()
        for j in range(3): 
            k = j+1
            host_corr[i][j] = float(res[k])
        host_name[i] = res[4]
    for i in range(3):
        defect.readline()
    for i in range(defect_no): 
        res = defect.readline().split()
        for j in range(3):
            k = j+1
            defect_corr[i][j] = float(res[k])
        defect_name[i] = res[4]
    if (defect_no == atom_no):
        for i in range(atom_no): 
            if (defect_name[i] != host_name[i]): 
                defect_line = i
                break
    else: 
        for i in range(atom_no): 
            dist = np.linalg.norm(host_corr[i,:]-defect_corr[i,:])
            if (dist > 0.2):
                defect_line = i
                break
    return defect_line, lattice_vectors, host_corr  


# Suzy edits: use this function as is, but use notebook parameters as function args?
#obtain the atomic potential from the calculations  
def atomic_pot_fhiaims(atom_no,defect_no,defect_line,latvec,host_corr,host,defect,shift_H,shift_D): 
    host_pot = open(host,'r')
    host_pot.readline()
    defect_pot = open(defect,'r')
    defect_pot.readline()
    h_pot = np.zeros([atom_no-1])
    D_pot = np.zeros([atom_no-1])
    Final = np.zeros([atom_no-1,2]) 
#    host_pot.readline()
    if (defect_line == 0):
        host_pot.readline()
        for i in range(1,atom_no): 
            res = host_pot.readline().split()
            h_pot[i-1] = float(res[1]) 
    elif (defect_line == atom_no -1): 
        for i in range(0,atom_no-1): 
            res = host_pot.readline().split()
            h_pot[i] = float(res[1])
    else:
        for i in range(0,defect_line): 
            res = host_pot.readline().split()
            h_pot[i] = float(res[1])
        host_pot.readline()    
        for i in range(defect_line+1,atom_no): 
            res = host_pot.readline().split()
            h_pot[i-1] = float(res[1])
#    defect_pot.readline() 
#    print defect_line 
    if (defect_no == atom_no):
        if (defect_line == 0):
            defect_pot.readline()
            for i in range(1,atom_no): 
                res2 = defect_pot.readline().split()
                D_pot[i-1] = float(res2[1]) 
        elif (defect_line == atom_no -1): 
            for i in range(0,atom_no-1): 
                res2 = defect_pot.readline().split()
                D_pot[i] = float(res2[1])
        else:
            for i in range(0,defect_line): 
                res2 = defect_pot.readline().split()
                D_pot[i] = float(res2[1])
            defect_pot.readline()
            for i in range(defect_line+1,atom_no): 
                res2 = defect_pot.readline().split()
                D_pot[i-1] = float(res2[1])
    else: 
        for i in range(defect_no): 
            res2 = defect_pot.readline().split()
            D_pot[i] = float(res2[1])
    print "H_pot:",h_pot 
    print "D_pot:",D_pot
    result = D_pot -shift_D -(h_pot - shift_H)  
    rel_defect_corr = np.zeros([atom_no-1,3])
    distance = np.zeros([atom_no-1])
    for i in range(0,defect_line): 
        rel_defect_corr[i,:] = host_corr[i,:] - host_corr[defect_line,:]
    for i in range(defect_line+1,atom_no):
        rel_defect_corr[i-1,:] = host_corr[i,:] - host_corr[defect_line,:]
    for i in range(atom_no-1): 
        a = rel_defect_corr[i,0]*latvec[0,:] + rel_defect_corr[i,1]*latvec[1,:] +rel_defect_corr[i,2]*latvec[2,:]
        distance[i] = np.linalg.norm(a)
    
    for i in range(0,defect_line): 
        Final[i,0] =  distance[i]
        Final[i,1] = result[i]
    for i in range(defect_line+1,atom_no):
        Final[i-1,0] = distance[i-1]
        Final[i-1,1] = result[i-1]
    return Final
        
        
'''   

#compute the atomic potential from the model   
def compute_atomic_pot(atom_no,defect_line,grid,latvec,host_corr,V_G,beta,G1,G2,G3): 
    latvec = latvec*bohr
    volume = np.dot(latvec[0,:],np.cross(latvec[1,:],latvec[2,:])) 
    rlatvec= np.zeros([3,3])
    step = np.zeros([3,3])
#    rlatvec[0,:] = 2*np.pi*np.cross(latvec[1,:],latvec[2,:])/volume
#    rlatvec[1,:] = 2*np.pi*np.cross(latvec[2,:],latvec[0,:])/volume
#    rlatvec[2,:] = 2*np.pi*np.cross(latvec[0,:],latvec[1,:])/volume
#    step[0,:] = rlatvec[0,:]/float(grid[0])
#    step[1,:] = rlatvec[1,:]/float(grid[1])
#    step[2,:] = rlatvec[2,:]/float(grid[2])
#    print rlatvec 
    #calculate the atom distance to defect atoms in realspace 
    rel_defect_corr = np.zeros([atom_no-1,3])
    host_corr_rel = np.zeros([atom_no,3])
    distance = np.zeros([atom_no-1])
    result = np.zeros([atom_no-1,2])
    for i in range(0,defect_line): 
        rel_defect_corr[i,:] = host_corr[i,:] - host_corr[defect_line,:]
    for i in range(defect_line+1,atom_no):
        rel_defect_corr[i-1,:] = host_corr[i,:] - host_corr[defect_line,:]
    for i in range(atom_no-1): 
        a = rel_defect_corr[i,0]*latvec[0,:] + rel_defect_corr[i,1]*latvec[1,:] +rel_defect_corr[i,2]*latvec[2,:]
        distance[i] = np.linalg.norm(a)
    #calculate the atomic site potential  
    for i in range(atom_no): 
        host_corr_rel[i,:] = host_corr[i,:] - host_corr[defect_line,:]
    V_atomic = atomic_3d(defect_line,atom_no,beta,grid,latvec,host_corr,V_G,G1,G2,G3)
    
    
#    for i in range(0,defect_line): 
#        R = host_corr[i,0]*latvec[0,:] + host_corr[i,1]*latvec[1,:] + host_corr[i,2]*latvec[2,:]
#        for j in range(1,grid[0]):
#            for k in range(1,grid[1]):
#                for l in range(1,grid[2]): 
#                    G = j*step[0,:]+k*step[1,:]+l*step[2,:] 
#                    #second = np.exp(iGR)
#                    second = np.exp(1j*np.dot(G,R))
#                    second = np.real(np.dot(second,V_G[j][k][l]))
#                    print second
#                    #third = np.exp(-0.5 G**2*beta**2)
#                    third = np.exp(-0.5*np.dot(G,G)*beta*beta)
#                    print third
                    #first = V(G) 
#                    V_G[j][k][l] = second*third 
#                    V_atomic[i] += second*third 
#    for i in range(defect_line+1,atom_no): 
#        R = host_corr[i,0]*latvec[0,:] + host_corr[i,1]*latvec[1,:] + host_corr[i,2]*latvec[2,:]
#        for j in range(1,grid[0]):
#            for k in range(1,grid[1]):
#                for l in range(1,grid[2]): 
#                    G = j*step[0,:]+k*step[1,:]+l*step[2,:] 
#                    #second = np.exp(iGR)
#                    second = np.exp(1j*G*R)
#                    #third = np.exp(-0.5 G**2*beta**2)
#                    third = np.exp(-0.5*G*G*beta*beta)
#                    #first = V(G) 
#                    V_atomic[i-1] = V_atomic[i-1] + (1/np.sqrt(volume))*V_G[j][k][l]*second*third 
#    V_atomic = V_atomic*grid[0]*grid[1]*grid[2]*hartree
    result[:,0] = distance/bohr
    result[:,1] = V_atomic 
    return result
    
'''


#main 
#charge = float(sys.argv[1])
#shift_H = float(sys.argv[2]) 
#shift_D = float(sys.argv[3])
#shift_H =  -15.02300189   
#shift_D =  -15.24318524 

defect_info, charge, shift_H, shift_D = sys_basic('system.in')
print "Charge:", charge
print "Shift_Host:", shift_H 
print "Shift_Defect:", shift_D 
atom_no,defect_no = geo_basic('geometry-host.in','geometry-defect.in')  
#print "atom_No, Defect_No, defect_line:", atom_no, defect_no 
defect_line,lattice_vec, host_corr = geo_compare(atom_no,defect_no,'geometry-host.in','geometry-defect.in')
print "atom_No, Defect_No, defect_line:", atom_no, defect_no, defect_line 
result = atomic_pot_fhiaims(atom_no,defect_no,defect_line,lattice_vec,host_corr,'pot-host.out','pot-defect.out',shift_H,shift_D)
np.savetxt('atom_potential_FHI-aims',result)

#lattice_constant_a = 0.5*np.sqrt(lattice_vec[0][0]**2 + lattice_vec[0][1]**2 + lattice_vec[0][2]**2) 
#lattice_constant_b = 0.5*np.sqrt(lattice_vec[1][0]**2 + lattice_vec[1][1]**2 + lattice_vec[1][2]**2)
#lattice_constant_c = 0.5*np.sqrt(lattice_vec[2][0]**2 + lattice_vec[2][1]**2 + lattice_vec[2][2]**2)

#lattice_ab = 0.5*np.sqrt ( (lattice_vec[0][0] -lattice_vec[1][0])**2 + (lattice_vec[0][1] -lattice_vec[1][1])**2+ (lattice_vec[0][2] -lattice_vec[1][2])**2 )

#lattice_constant_z = lattice_vec[2][2]
volume = np.dot(lattice_vec[0,:],np.cross(lattice_vec[1,:],lattice_vec[2,:]))  
Rws = (volume*3.0/(4.0*np.pi))**(1.0/3.0)

#print lattice_constant_a, lattice_constant_b, lattice_constant_c, lattice_ab
#Rws = min(lattice_constant_a,lattice_constant_b,lattice_constant_c,lattice_ab)
print "Rws:",Rws 

temp_a = (lattice_vec[0][0]+lattice_vec[1][0]+lattice_vec[2][0])**2 
temp_b = (lattice_vec[0][1]+lattice_vec[1][1]+lattice_vec[2][1])**2 
temp_c = (lattice_vec[0][2]+lattice_vec[1][2]+lattice_vec[2][2])**2 

Rwl  = 0.5*np.sqrt(temp_a + temp_b + temp_c)
print "Rwl:",Rwl 
#Rws= 0.5*lattice_constant_z 
#Rwl=0.5*np.sqrt(3)*lattice_constant_z

Model = np.loadtxt('atom_potential_model')
#Model[:,0] = -1.0*Model[:,0]

Diff = result[:,1] - Model[:,1] 

Dim = len(Diff)

number = 0
value = 0
for i in range(Dim): 
    if ((Model[i,0] > Rws) and (Model[i,0] < Rwl)): 
        value += Diff[i]
        number += 1.0 

value = charge*value/number 
print "Correction:",value 

value = -1.0*value 
value = '%.2f' % value 


print defect_info
title = defect_info[0] + ': '+r'$-\Delta q$'+'='+str(value) 


y_min = min(Model[:,1])
y_max = max(Model[:,1])

Ymin = y_min - 0.5*(y_max-y_min)
Ymax = y_max + 0.5*(y_max-y_min)      

plt.clf()
plt.xlim((0,Rwl))
plt.ylim((Ymin,Ymax))
plt.axvline(x=Rws,linestyle='--',color='orange') 
plt.title(title,fontsize=15,fontname = "Times New Roman")
#plt.plot(X1,Y1,'or',label='V(Defect,q)-V(Defect,0)')
plt.plot(result[:,0],result[:,1],'o',ms=5,markerfacecolor='none', markeredgecolor='red',label='V(Defect,q)-V(Host)')
plt.plot(Model[:,0],Model[:,1],'o',ms=5,markerfacecolor='none', markeredgecolor='blue',label='V(Model)')
plt.plot(Model[:,0],Diff,'ok',ms=5,label='(V(Defect,q)-V(Host))-V(Model)')
plt.legend()
plt.xlabel(r'Distance from Defect ($\AA$)',fontsize=15,fontname = "Times New Roman")
plt.ylabel('Potential (eV)',fontsize=15,fontname = "Times New Roman")
plt.savefig("PA",dpi=600)
#result = np.fft.ifftn(model*grid[0]*grid[1]*grid[2])
#result = compute_atomic_pot(atom_no,defect_line,grid,lattice_vec,host_corr,model,1.5,G1,G2,G3)
#result = -1.0*result 
#print result






   
