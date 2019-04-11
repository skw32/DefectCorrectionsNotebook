# SUZY EDITS, PLAN:
## Take functions from the top as separate script
## Incorporate 'main' section below into notebook using notebook parameters as function args
## (Replace command line args for bottom of script with notebook parameters also)
## Store outputs from 'main' section as files in defects_output_dir but also as notebook params to use for plot
#################################################

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 13:29:30 2018

@author: tong
"""


#import matplotlib.pyplot as plt
#import matplotlib as mpl
import random
import time
from scipy import optimize
#import matplotlib.colors
import numpy as np
#from scipy.optimize import curve_fit
from APA_script.atomic_V import atomic_3d   # Suzy edits: added APA_script path for calling from top dir with notebook

import os
import sys
import math
import csv

# Suzy edits
def suzy_testing_import(word):
    print("Look, it's working, here's a word: "+str(word))

hartree = 27.2116  
bohr = 1.8897259886     

# Suzy edits: leave as it is?
def py_read(filename):
    data = np.load(filename)
    grid = np.shape(data)
    return grid,data   

# Suzy edits: replace this with my notebook parameters/ functions for counting species, etc.?
def geo_basic(hostF,defectF): 
    host = open(hostF, 'r') 
    defect = open(defectF,'r')
    host_len = len(host.readlines())
    defect_len = len(defect.readlines())
    atom_no = host_len - 3 
    defect_no = defect_len -3 
    return atom_no, defect_no 
    
# Suzy edits: replace this with my notebook parameters/ functions for counting species, etc.?
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
def atomic_pot_fhiaims(atom_no,defect_no,defect_line,host,defect): 
    host_pot = open(host,'r')
    host_pot.readline()
    defect_pot = open(defect,'r')
    defect_pot.readline()
    h_pot = np.zeros([atom_no-1])
    D_pot = np.zeros([atom_no-1])
    for i in range(0,defect_line): 
        res = host_pot.readline().split()
        h_pot[i] = float(res[1])
    for i in range(defect_line+1,atom_no): 
        res = host_pot.readline().split()
        h_pot[i-1] = float(res[1])
    if (defect_no == atom_no):
        for i in range(0,defect_line): 
            res2 = defect_pot.readline().split()
            D_pot[i] = float(res2[1])
        for i in range(defect_line+1,atom_no): 
            res2 = defect_pot.readline().split()
            D_pot[i-1] = float(res2[1])
    else: 
        for i in range(defect_no): 
            res2 = defect_pot.readline().split()
            D_pot[i] = float(res2[1])
    result = D_pot -h_pot  
    return result 
        
        

    
# Suzy edits: use this function as is, but use notebook parameters as function args?
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
    



#main 
beta = float(sys.argv[1]) # Suzy edits: replace with sigma parameter from notebook?
#dir_temp = sys.argv[2]


atom_no,defect_no = geo_basic('geometry-host.in','geometry-defect.in')  
defect_line,lattice_vec, host_corr = geo_compare(atom_no,defect_no,'geometry-host.in','geometry-defect.in')
grid_t, model = py_read('V_G-model.npy')
grid_g1, G1 = py_read('G1.npy')
grid_g2, G2 = py_read('G2.npy')
grid_g3, G3 = py_read('G3.npy')
#print grid 
grid = np.array(grid_t)
grid.astype(int)
#print grid

model = model*hartree

#result = atomic_pot_fhiaims(atom_no,defect_no,defect_line,'host_pot.out','m1_pot.out')

#result = np.fft.ifftn(model*grid[0]*grid[1]*grid[2])
result = compute_atomic_pot(atom_no,defect_line,grid,lattice_vec,host_corr,model,beta,G1,G2,G3)
result[:,1] = -1.0*result[:,1] 
#print result
np.savetxt('atom_potential_model',result)





   
