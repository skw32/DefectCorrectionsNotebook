import numpy as np
import re
from math import sqrt

def read_lattice_vectors(geom_file):
    '''
    Input crystal geometry file in format for FHI-aims (geometry.in)
    Function searches for lattice vectors using string 'lattice_vector'
    Returns x, y and z components of a1, a2 and a3 lattice vectors
    E.g. x_vecs[1], y_vecs[1], z_vecs[1] would be x, y, z components of a2
    '''
    x_vecs = []
    y_vecs = []
    z_vecs = []
    try:
        with open(geom_file, 'r') as f:
            for line in f:
                if re.search('lattice_vector', line):
                    words = line.split()
                    x_vecs.append(float(words[1]))
                    y_vecs.append(float(words[2]))
                    z_vecs.append(float(words[3]))
                    if line == None:
                        logger.info('Warning! - No lattice vectors found in '+str(geom_file))
    except IOError:
        logger.info("Could not open "+str(geom_file))
    return x_vecs, y_vecs, z_vecs


def get_supercell_dimensions(geom_file):
    '''
    Take maximum of each direction to be supercell dimension for orthogonal unit cells
    (allowing for some numerical noise in off-diagonals)
    Return list 'supercell_dims' where x = supercell_dims[0], y = supercell_dims[0], z = supercell_dims[2]
    '''
    x_vecs, y_vecs, z_vecs = read_lattice_vectors(geom_file)
    supercell_dims = []
    supercell_dims.append(max(x_vecs))
    supercell_dims.append(max(y_vecs))
    supercell_dims.append(max(z_vecs))
    return supercell_dims


def read_atom_coords(geom_file):
    '''
    Input crystal geometry file in format for FHI-aims (geometry.in)
    Function searches for atom using string 'atom' to allow for either 'atom' or 'atom_frac' in the file format
    Returns list of lists for all atom coordinates where atom_coords[row][col]
    Columns are: x, y, z, species
    '''
    atom_coords = []
    try:
        with open(geom_file, 'r') as f:
            for line in f:
                if re.search('atom', line):
                    words = line.split()
                    atom_coords.append((float(words[1]), float(words[2]), float(words[3]), str(words[4])))
                    if line == None:
                        logger.info('Warning! - No atom coordinates found in '+str(geom_file))
    except IOError:
        logger.info("Could not open "+str(geom_file))
    return atom_coords


def find_defect_type(host_coords, defect_coords):
    '''
    Compares number of atoms in defect and host supercells to determine type of defect
    host_atom_num == defect_atom_num+1 --> vacancy
    host_atom_num == defect_atom_num-1 --> interstitial
    host_atom_num == defect_atom_num --> antisite
    '''
    host_atom_num = len(host_coords)
    defect_atom_num = len(defect_coords)
    if (host_atom_num == defect_atom_num+1):
        defect_type = 'vacancy'
    elif (host_atom_num == defect_atom_num-1):
        defect_type = 'interstitial'
    elif (host_atom_num == defect_atom_num):
        defect_type = 'antisite'
    else:
        logger.info('Error finding defect type')
    return defect_type


def count_species(host_coords, defect_coords):
    '''
    Read through species in atom_coords[row][3] for host and defect supercells
    First function output is a list of all different species present in the host supercell
    Next two outputs are the number of each of these species in the host and defect supercell, in the same order
    Assumption is made that only intrinsic defects are present, hence same species are present in host and defect supercells
    '''
    # Obtain list of all species contained in host supercell
    species = []
    current_species = host_coords[0][3]
    species.append(host_coords[0][3])
    for i in range(0, len(host_coords)):
        if (host_coords[i][3] != current_species):
            species.append(host_coords[i][3])
            current_species = host_coords[i][3]         
    # Count number of each species in host supercell
    host_species_nums = []
    for j in range(0, len(species)):
        species_count = 0
        for i in range(0, len(host_coords)):
            if (host_coords[i][3] == species[j]):
                species_count += 1
        host_species_nums.append(int(species_count))
    # Count number of each species in defect supercell
    defect_species_nums = []
    for j in range(0, len(species)):
        species_count = 0
        for i in range(0, len(defect_coords)):
            if (defect_coords[i][3] == species[j]):
                species_count += 1
        defect_species_nums.append(int(species_count))       
    return species, host_species_nums, defect_species_nums
    
       
def find_vacancy(host_coords, defect_coords):
    '''
    Find species where count is one less in defect supercell than in host supercell
    '''
    species, host_species_nums, defect_species_nums = count_species(host_coords, defect_coords)
    species_vac = 'no vacancy'
    for i in range (0, len(species)):
        if (host_species_nums[i] == defect_species_nums[i]+1):
            species_vac = species[i]
    if (species_vac == 'no vacancy'):
        logger.info('Error finding vacancy')
    return species_vac


def find_interstitial(host_coords, defect_coords):
    '''
    Find species where count is one more in defect supercell than in host supercell
    '''
    species, host_species_nums, defect_species_nums = count_species(host_coords, defect_coords)
    species_int = 'no interstitial'
    for i in range (0, len(species)):
        if (host_species_nums[i] == defect_species_nums[i]-1):
            species_int = species[i]
    if (species_int == 'no interstitial'):
        logger.info('Error finding interstitial')
    return species_int


def find_antisite(host_coords, defect_coords):
    '''
    Find species where count is one less in defect supercell than in host (species_out)
    Find species where count is one more in defect supercell than in host (species_in) 
    '''
    species, host_species_nums, defect_species_nums = count_species(host_coords, defect_coords)
    species_in = 'no species in'
    species_out = 'no species out'
    for i in range (0, len(species)):
        if (host_species_nums[i] == defect_species_nums[i]-1):
            species_in = species[i]
        if (host_species_nums[i] == defect_species_nums[i]+1):
            species_out = species[i]
    if (species_in == 'no species in' or species_out == 'no species out'):
        logger.info('Error finding antisite')
    return species_in, species_out


def vacancy_coords(host_coords, defect_coords):
    '''
    Define vacancy coordinates in perfect host supercell
    '''
    species_vac = find_vacancy(host_coords, defect_coords)  
    # Read in coordinates of vacancy species in perfect host supercell
    host_vac_coords = []
    for i in range (0, len(host_coords)):
        if (host_coords[i][3] == species_vac):
            host_vac_coords.append(host_coords[i][:3])
    # Read in coordinates of vacancy species in defect supercell
    defect_vac_coords = []
    for i in range (0, len(defect_coords)):
        if (defect_coords[i][3] == species_vac):
            defect_vac_coords.append(defect_coords[i][:3]) 
    # Find closest vacancy species in defect supercell for each one in host supercell
    all_closest_species = []
    for x_host, y_host, z_host in host_vac_coords:
        # USE np.argmin??
        closest_species = None
        min_distance = None
        for x_defect, y_defect, z_defect in defect_vac_coords:
            distance_to_defect = sqrt( (abs(x_host-x_defect)*abs(x_host-x_defect)) + (abs(y_host-y_defect)*abs(y_host-y_defect)) + (abs(z_host-z_defect)*abs(z_host-z_defect)))
            if min_distance is None or distance_to_defect < min_distance:
                min_distance = distance_to_defect
                closest_species = [x_defect, y_defect, z_defect]
        all_closest_species.append(closest_species + [min_distance])
    # Find which species in host where the 'closest distance' to a species in the defect supercell is largest
    # This is identified as the vacancy in te host supercell
    x_vac, y_vac, z_vac = host_coords[np.argmax([i[3] for i in all_closest_species])][:3]
    # Above in one-liner version of code commented out below
    #max_dist = 0
    #for i in range(0, len(all_closest_species)):
    #    if (all_closest_species[i][3] > max_dist):
    #        x_vac, y_vac, z_vac = host_coords[i][:3]
    #        max_dist = all_closest_species[i][3] 
    return species_vac, x_vac, y_vac, z_vac


def interstitial_coords(host_coords, defect_coords):
    '''
    Define interstitial coordinates in defect supercell
    '''
    species_int = find_interstitial(host_coords, defect_coords)  
    # Read in coordinates of interstitial species in perfect host supercell
    host_int_coords = []
    for i in range (0, len(host_coords)):
        if (host_coords[i][3] == species_int):
            host_int_coords.append(host_coords[i][:3])
    # Read in coordinates of interstitial species in defect supercell
    defect_int_coords = []
    for i in range (0, len(defect_coords)):
        if (defect_coords[i][3] == species_int):
            defect_int_coords.append(defect_coords[i][:3]) 
    # Find closest interstitial species in host supercell for each one in defect supercell
    all_closest_species = []
    for x_defect, y_defect, z_defect in defect_int_coords:
        # USE np.argmin??
        closest_species = None
        min_distance = None
        for x_host, y_host, z_host in host_int_coords:
            distance_to_defect = sqrt( (abs(x_host-x_defect)*abs(x_host-x_defect)) + (abs(y_host-y_defect)*abs(y_host-y_defect)) + (abs(z_host-z_defect)*abs(z_host-z_defect)))
            if min_distance is None or distance_to_defect < min_distance:
                min_distance = distance_to_defect
                closest_species = [x_host, y_host, z_host]
        all_closest_species.append(closest_species + [min_distance])
    # Find which species in defect where the 'closest distance' to a species in the host supercell is largest
    # This is identified as the interstitial in the defect supercell
    x_int, y_int, z_int = defect_coords[np.argmax([i[3] for i in all_closest_species])][:3]
    return species_int, x_int, y_int, z_int


def antisite_coords(host_coords, defect_coords):
    '''
    Define antisite coordinates in defect supercell
    '''
    species_in, species_out = find_antisite(host_coords, defect_coords)
    # Find species_in in defect supercell mostly using function for finding interstitial
    # Read in coordinates of antisite_in species in perfect host supercell
    host_in_coords = []
    for i in range (0, len(host_coords)):
        if (host_coords[i][3] == species_in):
            host_in_coords.append(host_coords[i][:3])
    # Read in coordinates of antisite_in species in defect supercell
    defect_in_coords = []
    for i in range (0, len(defect_coords)):
        if (defect_coords[i][3] == species_in):
            defect_in_coords.append(defect_coords[i][:3]) 
    # Find closest interstitial species in host supercell for each one in defect supercell
    all_closest_species = []
    for x_defect, y_defect, z_defect in defect_in_coords:
        # USE np.argmin??
        closest_species = None
        min_distance = None
        for x_host, y_host, z_host in host_in_coords:
            distance_to_defect = sqrt( (abs(x_host-x_defect)*abs(x_host-x_defect)) + (abs(y_host-y_defect)*abs(y_host-y_defect)) + (abs(z_host-z_defect)*abs(z_host-z_defect)))
            if min_distance is None or distance_to_defect < min_distance:
                min_distance = distance_to_defect
                closest_species = [x_host, y_host, z_host]
        all_closest_species.append(closest_species + [min_distance])
    # Find which species in defect where the 'closest distance' to a species in the perfect supercell is largest
    # This is identified as the species added into in the defect supercell
    x_in, y_in, z_in = defect_coords[np.argmax([i[3] for i in all_closest_species])][:3]
    return species_in, species_out, x_in, y_in, z_in


def defect_to_boundary(x_defect, y_defect, z_defect, supercell_x, supercell_y, supercell_z):
    # Finding minimum x, y, z distance of defect from supercell boundaries
    x_min = x_defect if (x_defect <= supercell_x/2.0) else supercell_x - x_defect
    y_min = y_defect if (y_defect <= supercell_y/2.0) else supercell_y - y_defect
    z_min = z_defect if (z_defect <= supercell_z/2.0) else supercell_z - z_defect
    return x_min, y_min, z_min