import numpy as np
import re
from math import sqrt
import logging
logger = logging.getLogger()


def count_atoms(geom_file: str) -> int:
    ''' Counts number of lines in file starting with 'atom' to allow for use of 'atom' or 'atom_frac'.
    NOTE: It is important lines containing atom coordinates have been deleted (not commented out) to create the defect.

    Args: 
        geom_file: input crystal geometry file in format for FHI-aims (geometry.in)
    
    Returns: 
        Number of atoms identified in file (int)
    '''
    atom_num = 0
    try:
        with open(geom_file, 'r') as f:
            for line in f:
                if re.search('atom', line):
                    atom_num += 1
    except IOError:
        logger.info("Could not open "+str(geom_file))
    return atom_num


def read_lattice_vectors(geom_file: str) -> list:
    ''' Function searches for lattice vectors using string 'lattice_vector'

    Args: 
        geom_file: input crystal geometry file in format for FHI-aims (geometry.in)

    Returns: 
        lists for x, y and z components of a1, a2 and a3 lattice vectors
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


def lattice_vectors_array(geom_file: str) -> tuple:
    ''' Function searches for lattice vectors using string 'lattice_vector' and returns them as numpy array.

    Args: 
        geom_file: input crystal geometry file in format for FHI-aims (geometry.in)

    Returns: 
        Each component of the lattice vectors as elements of a 3x3 numpy array
    '''
    latt_vec_array = np.zeros([3,3])
    try:
        with open(geom_file, 'r') as f:
            i = 0
            for line in f:
                if re.search('lattice_vector', line):
                    words = line.split()
                    latt_vec_array[i][0] = float(words[1])
                    latt_vec_array[i][1] = float(words[2])
                    latt_vec_array[i][2] = float(words[3])
                    i += 1
                    if line == None:
                        logger.info('Warning! - No lattice vectors found in '+str(geom_file))
    except IOError:
        logger.info("Could not open "+str(geom_file))      
    return latt_vec_array


def get_supercell_dimensions(geom_file: str) -> list:
    ''' Take maximum of each direction to be supercell dimension for orthogonal unit cells
    (allowing for some numerical noise in off-diagonals)

    Args: 
        geom_file: input crystal geometry file in format for FHI-aims (geometry.in)

    Returns: 
        List 'supercell_dims' where x = supercell_dims[0], y = supercell_dims[0], z = supercell_dims[2]
    '''
    x_vecs, y_vecs, z_vecs = read_lattice_vectors(geom_file)
    supercell_dims = []
    supercell_dims.append(max(x_vecs))
    supercell_dims.append(max(y_vecs))
    supercell_dims.append(max(z_vecs))
    return supercell_dims


def read_atom_coords(geom_file: str) -> list:
    ''' Function searches for atom using string 'atom' to allow for either 'atom' or 'atom_frac' in the file format
    If coordinates are fractional, these are converted to Cartesian coordinates

    Args: 
        geom_file: input crystal geometry file in format for FHI-aims (geometry.in)

    Returns: 
        List of lists for all atom coordinates where atom_coords[row][col]
        Columns are: x, y, z, species and each row is a different atom
    '''
    atom_coords = []
    # For converting from fractional to Cartesian coordinates
    latvec = lattice_vectors_array(geom_file)
    try:
        with open(geom_file, 'r') as f:
            for line in f:
                if re.search('atom', line):
                    words = line.split()
                    if (words[0] == 'atom_frac'):
                        # Convert from fractional coordinates to Cartesian coordinates using lattice vectors
                        cart_coords = float(words[1])*latvec[0,:] + float(words[2])*latvec[1,:] + float(words[3])*latvec[2,:]
                        atom_coords.append((cart_coords[0], cart_coords[1], cart_coords[2], str(words[4])))
                    else:
                        atom_coords.append((float(words[1]), float(words[2]), float(words[3]), str(words[4])))
                    if line == None:
                        logger.info('Warning! - No atom coordinates found in '+str(geom_file))
    except IOError:
        logger.info("Could not open "+str(geom_file))
    return atom_coords


def coords_to_array(coord_list: list) -> tuple:
    '''
    Args: 
        coords_list: list of atomic coordinates outputted from 'read_atom_coords' function

    Returns: 
        Only the coordinates (not also species type as in read_atom_coords) as a numpy array
    '''
    coords_array = np.zeros([len(coord_list),3])
    for i in range(len(coord_list)):
        coords_array[i][0], coords_array[i][1], coords_array[i][2] = coord_list[i][0], coord_list[i][1], coord_list[i][2]
    return coords_array


# To-DO: Add test for this function
def frac_coords_convert(geom_file: str) -> tuple:
    '''
    Args: 
        geom_file: input crystal geometry file in format for FHI-aims (geometry.in)
    
    Returns: 
        Inverts lattice vectors to use for conversion from Cartesian to fractional coordinates.
    '''
    latvec = lattice_vectors_array(geom_file)
    #lattice_inv_mat = np.linalg.inv(latvec)
    super_mat = latvec.transpose()
    super_invmat = np.linalg.inv(super_mat)
    def wrap_vec(v):
        relvec = np.dot(super_invmat,v)
        wrapvec = (relvec+1e-5) % 1.0 - 1e-5
        return np.dot(super_mat,wrapvec)
    return super_invmat, wrap_vec


# TO-DO: Add test for this function
def read_atom_coords_frac(geom_file: str) -> tuple:
    ''' Function searches for atom using string 'atom' to allow for either 'atom' or 'atom_frac' in the file format
    If coordinates are not already fractional, they are converted to fractional for compatibility with CoFFEE code routines

    Args: 
        geom_file: input crystal geometry file in format for FHI-aims (geometry.in)
    
    Returns: 
        Numpy array of fractional coordinates
    '''
    # For converting from Cartesian to fractional coordinates
    super_invmat, wrap_vec = frac_coords_convert(geom_file)
    # Initialise np array for coordinates
    coord_num = count_atoms(geom_file)
    coords = np.zeros([coord_num,3])
    frac_coords = np.zeros([coord_num,3])
    try:
        with open(geom_file, 'r') as f:
            i = 0
            for line in f:
                if re.search('atom', line):
                    words = line.split()
                    if (words[0] != 'atom_frac'):
                        # Convert from Cartesian coordinates to fractional coordinates using lattice vectors
                        coords[i,:] = float(words[1]), float(words[2]), float(words[3])
                        frac_coords[i,:] =  wrap_vec(coords[i,:])
                        frac_coords[i,:] = np.dot(super_invmat,coords[i,:])
                    else:
                        frac_coords[i,:] = float(words[1]), float(words[2]), float(words[3])
                    if line == None:
                        logger.info('Warning! - No atom coordinates found in '+str(geom_file))
                    i +=1 
    except IOError:
        logger.info("Could not open "+str(geom_file))   
    return frac_coords


def find_defect_type(host_coords: list, defect_coords: list) -> str:
    ''' Compares number of atoms in defect and host supercells to determine type of defect
    host_atom_num == defect_atom_num+1 --> vacancy
    host_atom_num == defect_atom_num-1 --> interstitial
    host_atom_num == defect_atom_num --> antisite

    Args: 
        host_coords: lists of coordinates of host supercell obtained with 'read_atom_coords' function 
        defect_coords: lists of coordinates of defect supercell obtained with 'read_atom_coords' function 
    
    Returns: 
        defect_type (vacancy, interstitial, antisite) as a string
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


def count_species(host_coords: list, defect_coords: list) -> list:
    ''' Reads through species in atom_coords[row][3] for host and defect supercells
    Assumption is made that only intrinsic defects are present, hence same atom types are present in host and defect supercells
    TZ: This assumption is not always work, and extrinsic defects is normal in the ``antisite" case, and ``interstitials". 
        As a result, I change the species count based on the defect supercell: 
        For intrinsic defects, it will be the same as before: 
        However, for extrinsic defects, it will count the extrinsic species, and the number in the host would be zero. 

    Args: 
        host_coords: lists of coordinates of host supercell obtained with 'read_atom_coords' function 
        defect_coords: lists of coordinates of defect supercell obtained with 'read_atom_coords' function

    Returns: 
        First function output is a list of all different species present in the host supercell
        Next two outputs are the number of each of these species in the host and defect supercell, in the same order
    '''
    species = [] 
    current_species = defect_coords[0][3]
    species.append(defect_coords[0][3])
    for i in range(0,len(defect_coords)):
        if (defect_coords[i][3] != current_species): 
            species.append(defect_coords[i][3])
            current_species = defect_coords[i][3]
    #eliminate the duplicate species 
    species = list(set(species))       
    #Suzy original code, keep it here. 
    #Obtain list of all species contained in host supercell
    #species = []
    #current_species = host_coords[0][3]
    #species.append(host_coords[0][3])
    #for i in range(0, len(host_coords)):
    #    if (host_coords[i][3] != current_species):
    #        species.append(host_coords[i][3])
    #        current_species = host_coords[i][3]         
    
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
    
       
def find_vacancy(host_coords: list, defect_coords: list) -> str:
    ''' Find species where count is one less in defect supercell than in host supercell.

    Args: 
        host_coords: lists of coordinates of host supercell obtained with 'read_atom_coords' function 
        defect_coords: lists of coordinates of defect supercell obtained with 'read_atom_coords' function

    Returns: 
        vacancy species as a string
    '''
    species, host_species_nums, defect_species_nums = count_species(host_coords, defect_coords)
    species_vac = 'no vacancy'
    for i in range (0, len(species)):
        if (host_species_nums[i] == defect_species_nums[i]+1):
            species_vac = species[i]
    if (species_vac == 'no vacancy'):
        logger.info('Error finding vacancy')
    return species_vac


def find_interstitial(host_coords: list, defect_coords: list) -> str:
    ''' Find species where count is one more in defect supercell than in host supercell.

    Args: 
        host_coords: lists of coordinates of host supercell obtained with 'read_atom_coords' function 
        defect_coords: lists of coordinates of defect supercell obtained with 'read_atom_coords' function 

    Returns: 
        interstitial species as a string
    '''
    species, host_species_nums, defect_species_nums = count_species(host_coords, defect_coords)
    species_int = 'no interstitial'
    for i in range (0, len(species)):
        if (host_species_nums[i] == defect_species_nums[i]-1):
            species_int = species[i]
    if (species_int == 'no interstitial'):
        logger.info('Error finding interstitial')
    return species_int


def find_antisite(host_coords: list, defect_coords: list) -> str:
    ''' Find species where count is one less in defect supercell than in host (species_out).
    Find species where count is one more in defect supercell than in host (species_in).

    Args: 
        host_coords: lists of coordinates of host supercell obtained with 'read_atom_coords' function 
        defect_coords: lists of coordinates of defect supercell obtained with 'read_atom_coords' function

    Returns: 
        Two strings, the first is the species added into the defect supercell to make the antisite defect 
        and the second is the species removed from the host
    
    '''
    species, host_species_nums, defect_species_nums = count_species(host_coords, defect_coords)
    species_in = 'no species in'
    species_out = 'no species out'     
    #If this is only used to find antisite, it's more easy to do the selections. 
    #the following code is not always work for the reason that people may change the antisite defects at any positions they want, as a result, the same species may not be in the same positions.   for example: 

    for i in range (0, len(species)):
        if (host_species_nums[i] == defect_species_nums[i]-1):
            species_in = species[i]
        if (host_species_nums[i] == defect_species_nums[i]+1):
            species_out = species[i]
    if (species_in == 'no species in' or species_out == 'no species out'):
        logger.info('Error finding antisite')
    return species_in, species_out


def vacancy_coords(host_coords: list, defect_coords: list) -> tuple:
    '''
    Args: 
        host_coords: lists of coordinates of host supercell obtained with 'read_atom_coords' function 
        defect_coords: lists of coordinates of defect supercell obtained with 'read_atom_coords' function

    Returns: 
        Vacancy species as string, vacancy coordinates in the perfect host supercell and the line in the geometry file for the defect
        defect_line for a vacancy is defined as the line number in the perfect host supercell of the atom missing in the vacancy supercell
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
        closest_species = None
        min_distance = None
        for x_defect, y_defect, z_defect in defect_vac_coords:
            distance_to_defect = sqrt( (abs(x_host-x_defect)*abs(x_host-x_defect)) + (abs(y_host-y_defect)*abs(y_host-y_defect)) + (abs(z_host-z_defect)*abs(z_host-z_defect)))
            if min_distance is None or distance_to_defect < min_distance:
                min_distance = distance_to_defect
                closest_species = [x_defect, y_defect, z_defect]
        all_closest_species.append(closest_species + [min_distance])
    # Find which species in host where the 'closest distance' to a species in the defect supercell is largest
    # This is identified as the vacancy in the host supercell
    x_vac, y_vac, z_vac = host_vac_coords[np.argmax([i[3] for i in all_closest_species])][:3]
    # Above in one-liner version of code commented out below
    '''
    max_dist = 0
    for i in range(0, len(all_closest_species)):
        if (all_closest_species[i][3] > max_dist):
            x_vac, y_vac, z_vac = host_coords[i][:3]
            max_dist = all_closest_species[i][3]
    '''
    # Find line number in coordinates list that corresponds to the defect (in host supercell for vacancy)
    tmp = host_vac_coords[np.argmax([i[3] for i in all_closest_species])][:3]
    for i in range (0, len(host_coords)):
        if (host_coords[i][0:3] == tmp):
            defect_line = i
    return species_vac, x_vac, y_vac, z_vac, defect_line


def interstitial_coords(host_coords: list, defect_coords: list) -> tuple:
    '''
     Args: 
        host_coords: lists of coordinates of host supercell obtained with 'read_atom_coords' function 
        defect_coords: lists of coordinates of defect supercell obtained with 'read_atom_coords' function

    Returns: 
        Vacancy species as string, vacancy coordinates in the perfect host supercell and the line in the geometry file for the defect
        defect_line for an interstitial is defined as the line number in the defect supercell of the atom not present in the host supercell
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
    x_int, y_int, z_int = defect_int_coords[np.argmax([i[3] for i in all_closest_species])][:3]
    # Find line number in coordinates list that corresponds to the defect (in defect supercell for interstitial)
    tmp = defect_int_coords[np.argmax([i[3] for i in all_closest_species])][:3]
    for i in range (0, len(defect_coords)):
        if (defect_coords[i][0:3] == tmp):
            defect_line = i
    return species_int, x_int, y_int, z_int, defect_line


def antisite_coords(host_coords: list, defect_coords: list) -> tuple:
    ''' NOTE: TZ:  The way to find the defect is worng for some specific structures, we should find a better way to find the defect atoms. 
    I try to fixed that within the antisite_coords first, probably similar job should be done in the other routines, i.e. vacancies, and interstitial.
    
    Args: 
        host_coords: lists of coordinates of host supercell obtained with 'read_atom_coords' function 
        defect_coords: lists of coordinates of defect supercell obtained with 'read_atom_coords' function

    Returns: 
        Vacancy species as string, vacancy coordinates in the perfect host supercell and the line in the geometry file for the defect
        defect_line for an antisite is defined as the line number in the defect supercell of the atom not present in the host supercell   
    '''
    species_in, species_out = find_antisite(host_coords, defect_coords)

    # Find species_in in defect supercell mostly using function for finding interstitial
    # Read in coordinates of antisite_out species in perfect host supercell
    # TZ: To make it work for all cases, this should be species_out in the host, and species_in in the defects. 
    # The defect atoms can be find to be the smallest distance between this coordinates.  the other difference existed in the crystal 
    # should always be larger than this small changes.   

    host_in_coords = []
    for i in range (0, len(host_coords)):
        if (host_coords[i][3] == species_in):
            host_in_coords.append(host_coords[i][:3]) 
    # Read in coordinates of antisite_out species in defect supercell
    defect_in_coords = []
    for i in range (0, len(defect_coords)):
        if (defect_coords[i][3] == species_in):
            defect_in_coords.append(defect_coords[i][:3])  
    # Find closest antisite_in species in host supercell for each one in defect supercell
# Original code from suzy 
#    all_closest_species = []
#    for x_defect, y_defect, z_defect in defect_in_coords:
#        closest_species = None
#        min_distance = None
#        for x_host, y_host, z_host in host_in_coords:
#            distance_to_defect = sqrt( (abs(x_host-x_defect)*abs(x_host-x_defect)) + (abs(y_host-y_defect)*abs(y_host-y_defect)) + (abs(z_host-z_defect)*abs(z_host-z_defect)))
#            if min_distance is None or distance_to_defect < min_distance:
#                min_distance = distance_to_defect
#                closest_species = [x_defect, y_defect, z_defect]
#        all_closest_species.append(closest_species + [min_distance])
    if not host_in_coords: 
        #if the species_in coords in the host cell are an empty list, which means the species in is extrinsic defects. 
        x_in, y_in, z_in = defect_in_coords[0][:3] 
    else: 
        #if the species_in coords are not belong to a extrinsic defect type, we compare the min_distance value of the corrd of species_in in the host and defect 
        #atoms has the largest min distance value is the defect atoms
        all_closest_species = []
        for x_defect, y_defect, z_defect in defect_in_coords:
            closest_species = None
            min_distance = None
            for x_host, y_host, z_host in host_in_coords:
                distance_to_defect = sqrt( (abs(x_host-x_defect)*abs(x_host-x_defect)) + (abs(y_host-y_defect)*abs(y_host-y_defect)) + (abs(z_host-z_defect)*abs(z_host-z_defect)))
                if min_distance is None or distance_to_defect < min_distance:
                    min_distance = distance_to_defect
                    closest_species = [x_defect, y_defect, z_defect]
            all_closest_species.append(closest_species + [min_distance])
    # Find which species in defect where the 'closest distance' to a species in the perfect supercell is largest
    # This is identified as the species added into in the defect supercell
        x_in, y_in, z_in = defect_in_coords[np.argmax([i[3] for i in all_closest_species])][:3]
    # Find line number in coordinates list that corresponds to the defect (in defect supercell for antisite)
    tmp = (x_in, y_in, z_in)
    #tmp = defect_in_coords[np.argmin([i[3] for i in all_closest_species])][:3]
    for i in range (0, len(defect_coords)):
        if (defect_coords[i][0:3] == tmp):
            defect_line = i
    return species_in, species_out, x_in, y_in, z_in, defect_line


def defect_to_boundary(x_defect: float, y_defect: float, z_defect: float, supercell_x: float, supercell_y: float, supercell_z: float) -> float:
    '''
    Args: 
        x_defect: x-coordinate of defect within the supercell
        y_defect: y-coordinate of defect within the supercell
        z_defect: z-coordinate of defect within the supercell
        supercell_x: dimension of supercell along the x-direction
        supercell_y: dimension of supercell along the y-direction
        supercell_z: dimension of supercell along the z-direction
        (All are extracted as part of the notebook workflow and above are their default names in the notebook)
    
    Returns: 
        Closest distances in the x, y and z-directions of the defect to any of the supercell boundaries
    '''
    # Finding minimum x, y, z distance of defect from supercell boundaries
    x_min = x_defect if (x_defect <= supercell_x/2.0) else supercell_x - x_defect
    y_min = y_defect if (y_defect <= supercell_y/2.0) else supercell_y - y_defect
    z_min = z_defect if (z_defect <= supercell_z/2.0) else supercell_z - z_defect
    return x_min, y_min, z_min
