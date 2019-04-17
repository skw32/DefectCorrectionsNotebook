"""
Created on Thu Nov 29 13:29:30 2018

@author: tong

Modified by Suzanne K. Wallace (with permission) for incorporation into DefectCorrectionsNotebook
"""

import numpy as np
from PoissonSolver import atomic_3d   # Suzy edits: after updating coffee_poisson_solver_ko conda pkg
import logging
import re
logger = logging.getLogger()


def read_free_atom_pot(planar_pot_file):
    '''
    Reads free atom potential from top line of 'plane_average_realspace_ESP.out' FHI-aims output file
    '''
    try:
        with open(planar_pot_file, 'r') as f:
            for line in f:
                if re.search('#Numerical', line):
                    words = line.split()
                    free_atom_pot = float(words[7])
                    if line == None:
                        logger.info('Warning! - No average free-atom electrostatic potential found in '+str(planar_pot_file))
    except IOError:
        logger.info("Could not open "+str(planar_pot_file))
    return free_atom_pot 
    
       
# Compute the atomic potentials for the CoFFEE charge model   
def model_atomic_pot(host_atom_num,defect_line,grid,latvec,host_coords,V_G,sigma,G1,G2,G3):
    bohr = 1.8897259886 
    latvec = latvec*bohr
    #volume = np.dot(latvec[0,:],np.cross(latvec[1,:],latvec[2,:])) #Redundant lines?
    #rlatvec= np.zeros([3,3])
    #step = np.zeros([3,3])
    # Calculate the atom distance to defect atoms in realspace 
    rel_defect_corr = np.zeros([host_atom_num-1,3])
    host_corr_rel = np.zeros([host_atom_num,3])
    distance = np.zeros([host_atom_num-1])
    model_atom_pots = np.zeros([host_atom_num-1,2])
    for i in range(0,defect_line): 
        rel_defect_corr[i,:] = host_coords[i,:] - host_coords[defect_line,:]
    for i in range(defect_line+1,host_atom_num):
        rel_defect_corr[i-1,:] = host_coords[i,:] - host_coords[defect_line,:]
    for i in range(host_atom_num-1): 
        a = rel_defect_corr[i,0]*latvec[0,:] + rel_defect_corr[i,1]*latvec[1,:] +rel_defect_corr[i,2]*latvec[2,:]
        distance[i] = np.linalg.norm(a)
    # Calculate the atomic site potential  
    for i in range(host_atom_num): 
        host_corr_rel[i,:] = host_coords[i,:] - host_coords[defect_line,:]
    V_atomic = atomic_3d(defect_line,host_atom_num,sigma,grid,latvec,host_coords,V_G,G1,G2,G3)
    model_atom_pots[:,0] = distance/bohr
    model_atom_pots[:,1] = V_atomic 
    return model_atom_pots


# Obtain the atomic potential from outputs of FHI-aims calculations  
def fhiaims_atomic_pot(host_atom_num,defect_atom_num,defect_line,latvec,host_coords,host_atom_pot,defect_atom_pot,shift_H,shift_D): 
    host_pot = open(host_atom_pot,'r')
    host_pot.readline()
    defect_pot = open(defect_atom_pot,'r')
    defect_pot.readline()
    h_pot = np.zeros([host_atom_num-1])
    D_pot = np.zeros([host_atom_num-1])
    Final = np.zeros([host_atom_num-1,2]) 
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
    result = D_pot -shift_D -(h_pot - shift_H)  
    rel_defect_corr = np.zeros([host_atom_num-1,3])
    distance = np.zeros([host_atom_num-1])
    for i in range(0,defect_line): 
        rel_defect_corr[i,:] = host_coords[i,:] - host_coords[defect_line,:]
    for i in range(defect_line+1,host_atom_num):
        rel_defect_corr[i-1,:] = host_coords[i,:] - host_coords[defect_line,:]
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