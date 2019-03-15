# Tests plans:

# Unit tests
## Check that values read for FHI-aims dummy structure by functions are a given value (within set tolerance)
### For: 'read_lattice_vectors', 'read_atom_coords', 'find_defect_type', 'count_species', 

# Integration tests
## Check that functions using other functions give a certain value within set tolerance
### For: 'get_supercell_dimensions', 'find_vacancy', 'find_interstitial', 'find_antisite', 
### 'vacancy_coords', 'interstitial_coords', 'antisite_coords', 'defect_to_boundary'

# -----------------------------------------------------------------------------------------------------------

# Decide how to import custom functions from main dir???