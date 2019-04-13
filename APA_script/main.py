# SUZY EDITS, PLAN:
## Take functions from the top as separate script
## Incorporate 'main' section below into notebook using notebook parameters as function args
## (Replace command line args for bottom of script with notebook parameters also)
## Store outputs from 'main' section as files in defects_output_dir but also as notebook params to use for plot
#################################################

"""
Created on Thu Nov 29 13:29:30 2018

@author: tong

Modified by Suzanne K. Wallace (with permission) for incorporation into DefectCorrectionsNotebook
"""

import random
import time
from scipy import optimize
import numpy as np
from PoissonSolver import atomic_3d   # Suzy edits: after updating coffee_poisson_solver_ko conda pkg
import os
import sys
import math
import csv

# Suzy edits
def suzy_testing_import(word):
    print("Look, it's working, here's a word: "+str(word))
      

'''
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
'''

        
# Suzy edits: use this function as is, but use notebook parameters as function args?
#compute the atomic potential from the model   
def compute_atomic_pot(host_atom_num,defect_line,grid,latvec,host_corr,V_G,sigma,G1,G2,G3): 
    latvec = latvec*bohr
    volume = np.dot(latvec[0,:],np.cross(latvec[1,:],latvec[2,:])) 
    rlatvec= np.zeros([3,3])
    step = np.zeros([3,3])
    #calculate the atom distance to defect atoms in realspace 
    rel_defect_corr = np.zeros([host_atom_num-1,3])
    host_corr_rel = np.zeros([host_atom_num,3])
    distance = np.zeros([host_atom_num-1])
    result = np.zeros([host_atom_num-1,2])
    for i in range(0,defect_line): 
        rel_defect_corr[i,:] = host_corr[i,:] - host_corr[defect_line,:]
    for i in range(defect_line+1,host_atom_num):
        rel_defect_corr[i-1,:] = host_corr[i,:] - host_corr[defect_line,:]
    for i in range(host_atom_num-1): 
        a = rel_defect_corr[i,0]*latvec[0,:] + rel_defect_corr[i,1]*latvec[1,:] +rel_defect_corr[i,2]*latvec[2,:]
        distance[i] = np.linalg.norm(a)
    #calculate the atomic site potential  
    for i in range(host_atom_num): 
        host_corr_rel[i,:] = host_corr[i,:] - host_corr[defect_line,:]
    V_atomic = atomic_3d(defect_line,host_atom_num,sigma,grid,latvec,host_corr,V_G,G1,G2,G3)
    result[:,0] = distance/bohr
    result[:,1] = V_atomic 
    return result


#obtain the atomic potential from the calculations  
def atomic_pot_fhiaims_plot(host_atom_num,defect_atom_num,defect_line,latvec,host_corr,host_geom,defect_geom,shift_H,shift_D): 
    host_pot = open(host_geom,'r')
    host_pot.readline()
    defect_pot = open(defect_geom,'r')
    defect_pot.readline()
    h_pot = np.zeros([host_atom_num-1])
    D_pot = np.zeros([host_atom_num-1])
    Final = np.zeros([host_atom_num-1,2]) 
#    host_pot.readline()
    if (defect_line == 0):
        host_pot.readline()
        for i in range(1,host_atom_num): 
            res = host_pot.readline().split()
            h_pot[i-1] = float(res[1]) 
    elif (defect_line == host_atom_num -1): 
        for i in range(0,host_atom_num-1): 
            res = host_pot.readline().split()
            h_pot[i] = float(res[1])
    else:
        for i in range(0,defect_line): 
            res = host_pot.readline().split()
            h_pot[i] = float(res[1])
        host_pot.readline()    
        for i in range(defect_line+1,host_atom_num): 
            res = host_pot.readline().split()
            h_pot[i-1] = float(res[1])
#    defect_pot.readline() 
#    print defect_line 
    if (defect_atom_num == host_atom_num):
        if (defect_line == 0):
            defect_pot.readline()
            for i in range(1,host_atom_num): 
                res2 = defect_pot.readline().split()
                D_pot[i-1] = float(res2[1]) 
        elif (defect_line == host_atom_num -1): 
            for i in range(0,host_atom_num-1): 
                res2 = defect_pot.readline().split()
                D_pot[i] = float(res2[1])
        else:
            for i in range(0,defect_line): 
                res2 = defect_pot.readline().split()
                D_pot[i] = float(res2[1])
            defect_pot.readline()
            for i in range(defect_line+1,host_atom_num): 
                res2 = defect_pot.readline().split()
                D_pot[i-1] = float(res2[1])
    else: 
        for i in range(defect_atom_num): 
            res2 = defect_pot.readline().split()
            D_pot[i] = float(res2[1])
    print "H_pot:",h_pot 
    print "D_pot:",D_pot
    result = D_pot -shift_D -(h_pot - shift_H)  
    rel_defect_corr = np.zeros([host_atom_num-1,3])
    distance = np.zeros([host_atom_num-1])
    for i in range(0,defect_line): 
        rel_defect_corr[i,:] = host_corr[i,:] - host_corr[defect_line,:]
    for i in range(defect_line+1,host_atom_num):
        rel_defect_corr[i-1,:] = host_corr[i,:] - host_corr[defect_line,:]
    for i in range(host_atom_num-1): 
        a = rel_defect_corr[i,0]*latvec[0,:] + rel_defect_corr[i,1]*latvec[1,:] +rel_defect_corr[i,2]*latvec[2,:]
        distance[i] = np.linalg.norm(a)
    
    for i in range(0,defect_line): 
        Final[i,0] =  distance[i]
        Final[i,1] = result[i]
    for i in range(defect_line+1,host_atom_num):
        Final[i-1,0] = distance[i-1]
        Final[i-1,1] = result[i-1]
    return Final
    






   
