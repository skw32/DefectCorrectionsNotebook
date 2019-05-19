"""
Created on Thu Nov 29 13:29:30 2018

@author: tong

Modified by Suzanne K. Wallace (with permission) for incorporation into DefectCorrectionsNotebook and generalisation to include interstitial defects
"""

import numpy as np
from PoissonSolver import atomic_3d   # SKW edits: after updating coffee_poisson_solver_ko conda pkg
import logging
import re
logger = logging.getLogger()


def read_free_atom_pot(planar_pot_file):
    '''
    Arguments: planar average of electrostatic potential file from FHI-aims calculation
    Returns: Free atom potential from top line of 'plane_average_realspace_ESP.out' FHI-aims output file
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
    

def host_coords_skip_defect(defect_type, defect_line, host_atom_num, lattice_vec_array, host_coords_array, defect_coords_array):
    '''
    Arguments: 
    - defect_type and defect_line are parameters from notebook workflow with default names
    - host_coords_array and defect_coords_array are numpy arrays for fractional coordinates obtained e.g. with 'host_coords_array = dsa.coords_to_array(host_coords_frac)'
    Returns: coordinates of atoms in host supercell relative to defect position, skipping coordinates of defect in the case of antisite or vacancy
    '''
    if (defect_type == 'interstitial'):
        distance = np.zeros([host_atom_num]) # All host atoms plotted for interstitials
        rel_defect_corr = np.zeros([host_atom_num,3])
    else:
        distance = np.zeros([host_atom_num-1]) # For antisite or vacancy, atom corresponding to defect is not plotted
        rel_defect_corr = np.zeros([host_atom_num-1,3])

    if (defect_type == 'interstitial'): # All host coordinates are plotted, relative to defect location in defect supercell
        for i in range(host_atom_num):
            rel_defect_corr[i,:] = host_coords_array[i,:] - defect_coords_array[defect_line,:]
            a = rel_defect_corr[i,0]*lattice_vec_array[0,:] + rel_defect_corr[i,1]*lattice_vec_array[1,:] +rel_defect_corr[i,2]*lattice_vec_array[2,:]
            distance[i] = np.linalg.norm(a)

    if (defect_type == 'antisite' or defect_type == 'vacancy'): # Defect is skipped in host supercell and all other atoms coordinates are plotted relative to defect position
        if (defect_line == 0):
            for i in range(1,host_atom_num): # Omit first set of coordinates for defect 
                rel_defect_corr[i,:] = host_coords_array[i,:] - host_coords_array[defect_line,:]
        elif (defect_line > 0 and defect_line < host_atom_num-1): # Omit defect_line
            for i in range(0,defect_line): 
                rel_defect_corr[i,:] = host_coords_array[i,:] - host_coords_array[defect_line,:]
            for i in range(defect_line+1,host_atom_num):
                rel_defect_corr[i-1,:] = host_coords_array[i,:] - host_coords_array[defect_line,:]
        else: # Defect line is final line of host supercell
            for i in range(0,host_atom_num-1): # Omit final set of coordinates for defect 
                rel_defect_corr[i,:] = host_coords_array[i,:] - host_coords_array[defect_line,:]
        for i in range(host_atom_num-1): 
            a = rel_defect_corr[i,0]*lattice_vec_array[0,:] + rel_defect_corr[i,1]*lattice_vec_array[1,:] +rel_defect_corr[i,2]*lattice_vec_array[2,:]
            distance[i] = np.linalg.norm(a)

    return distance


# Compute the atomic potentials for the CoFFEE charge model   
def model_atomic_pot(defect_type,host_atom_num,defect_line,grid,lattice_vec_array,host_coords_array,defect_coords_array,V_G,sigma,G1,G2,G3):
    '''
    Arguments:
    Returns:
    '''
    bohr = 1.8897259886 
    lattice_vec_array = lattice_vec_array*bohr
    # Calculate coordinates of atoms in host lattice relative to defect coordinates, but omitting coordinates of defect
    distance = host_coords_skip_defect(defect_type, defect_line, host_atom_num, lattice_vec_array, host_coords_array, defect_coords_array)
    # SKW: **CHECK WITH TONG:** Looks like V_atomic always skips defect_line, so gives broadcasting error for interstitials when filling model_atom_pots array
    V_atomic = atomic_3d(defect_type,defect_line,host_atom_num,sigma,grid,lattice_vec_array,host_coords_array,V_G,G1,G2,G3)
    if (defect_type == 'interstitial'): # All atoms included when calculating 'distance' above
        model_atom_pots = np.zeros([host_atom_num,2])
    else: # Defect species skipped when calculating 'distance' above for antisites or vacancies
        model_atom_pots = np.zeros([host_atom_num-1,2])
    model_atom_pots[:,0] = distance/bohr
    model_atom_pots[:,1] = V_atomic 

    return model_atom_pots


# Obtain the atomic potential from outputs of FHI-aims calculations  
def fhiaims_atomic_pot(defect_type, host_atom_num,defect_atom_num,defect_line,lattice_vec_array,host_coords_array, defect_coords_array, host_atom_pot,defect_atom_pot,shift_H,shift_D): 
    '''
    Arguments:
    Returns:
    '''
    # Data for atom potentials from FHI-aims outputs
    host_pot = open(host_atom_pot,'r')
    host_pot.readline() # Used to skip header of file
    defect_pot = open(defect_atom_pot,'r')
    defect_pot.readline() # Used to skip header of file
    if (defect_type == 'antisite' or defect_type == 'vacancy'):
        h_pot = np.zeros([host_atom_num-1]) # subtract 1 because of skipping defect line
        D_pot = np.zeros([host_atom_num-1])
        Final = np.zeros([host_atom_num-1,2])
    if (defect_type == 'interstitial'): # no need to subtract 1 as line to skip is additional in defect supercell
        h_pot = np.zeros([host_atom_num])
        D_pot = np.zeros([host_atom_num])
        Final = np.zeros([host_atom_num,2])
    # Read in potentials from FHI-aims output for perfect host supercell and skip reading in defect_line
    if (defect_type == 'interstitial'): # Read all lines for atom potentials in host
        for i in range(0,host_atom_num):
            res = host_pot.readline().split()
            h_pot[i] = float(res[1])
    if (defect_type == 'vacancy' or defect_type == 'antisite'): # Skip the line in the host that corresponds to the defect
        if (defect_line == 0):
            host_pot.readline() # Skip first line
            for i in range(0,host_atom_num-1): 
                res = host_pot.readline().split()
                h_pot[i] = float(res[1]) 
        elif (defect_line == host_atom_num -1): 
            for i in range(0,host_atom_num-1): # Skip last line 
                res = host_pot.readline().split()
                h_pot[i] = float(res[1])
        else:
            for i in range(0,defect_line): 
                res = host_pot.readline().split()
                h_pot[i] = float(res[1])
            host_pot.readline() # Skip defect_line   
            for i in range(defect_line,host_atom_num-1): 
                res = host_pot.readline().split()
                h_pot[i] = float(res[1])          
    # Read in potentials from FHI-aims output for defect supercell and skip reading in defect_line
    if (defect_type == 'vacancy'): # Read all atom potentials from defect supercell
        for i in range(0,host_atom_num-1): 
            res2 = defect_pot.readline().split()
            D_pot[i] = float(res2[1])
    if (defect_type == 'antisite'): # Read all atom potentials, but skip line for the defect
        if (defect_line == 0):
            defect_pot.readline() # Skip first line
            for i in range(0,host_atom_num-1): 
                res2 = defect_pot.readline().split()
                D_pot[i] = float(res2[1]) 
        elif (defect_line == host_atom_num -1): 
            for i in range(0,host_atom_num-1): # Skip last line
                res2 = defect_pot.readline().split()
                D_pot[i] = float(res2[1])
        else:
            for i in range(0,defect_line): 
                res2 = defect_pot.readline().split()
                D_pot[i] = float(res2[1])
            defect_pot.readline() # Skip defect_line
            for i in range(defect_line,host_atom_num-1): 
                res2 = defect_pot.readline().split()
                D_pot[i] = float(res2[1])
    if (defect_type == 'interstitial'): # Read all atom potentials, but skip line for the defect
        if (defect_line == 0):
            defect_pot.readline() # Skip first line
            for i in range(0,host_atom_num): 
                res2 = defect_pot.readline().split()
                D_pot[i] = float(res2[1]) 
        elif (defect_line > 0 and defect_line < defect_atom_num-1):
            for i in range(0,defect_line): 
                res2 = defect_pot.readline().split()
                D_pot[i] = float(res2[1])
            defect_pot.readline() # Skip defect_line
            for i in range(defect_line,host_atom_num): 
                res2 = defect_pot.readline().split()
                D_pot[i] = float(res2[1])
        else: # Interstitial is final line of defect supercell file, so only last line needs to be skipped
            for i in range(0,host_atom_num-1): # Skip last line
                res2 = defect_pot.readline().split()
                D_pot[i] = float(res2[1])
                        
    # Preparing pa plot with potentials read in above and free atom potential shifts read in from planar average FHI-aims output file
    result = D_pot -shift_D -(h_pot - shift_H)  
    # Calculate coordinates of atoms in host lattice relative to defect coordinates, but omitting coordinates of defect
    distance = host_coords_skip_defect(defect_type, defect_line, host_atom_num, lattice_vec_array, host_coords_array, defect_coords_array)

    if (defect_type == 'interstitial'): # All host coordinates are plotted, relative to defect location in defect supercell
        Final = np.zeros([host_atom_num,2])
        for i in range(0,host_atom_num): 
            Final[i,0] =  distance[i]
            Final[i,1] = result[i]    
            
    if (defect_type == 'antisite' or defect_type == 'vacancy'): # Defect was skipped in host supercell and all other atoms coordinates are plotted relative to defect position
        Final = np.zeros([host_atom_num-1,2])
        for i in range(0,host_atom_num-1): 
            Final[i,0] =  distance[i]
            Final[i,1] = result[i]
    return Final