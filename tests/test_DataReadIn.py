# Tests plans:

# Unit tests
## Check that values read for FHI-aims dummy structure by functions are a given value (within set tolerance)
### For: 'read_lattice_vectors', 'read_atom_coords', 'find_defect_type', 'count_species', 

# Integration tests
## Check that functions using other functions give a certain value within set tolerance
### For: 'get_supercell_dimensions', 'find_vacancy', 'find_interstitial', 'find_antisite', 
### 'vacancy_coords', 'interstitial_coords', 'antisite_coords', 'defect_to_boundary'

# -----------------------------------------------------------------------------------------------------------

import DefectSupercellAnalyses as dsa
import pytest

#def test_test():
#    dummy_value = dsa.get_supercell_dimensions("WorkflowTests/TestData/vacancy/geometry.in")
#    assert dummy_value == 10

# Unit tests for functions

# Testing functions that use other functions
def test_get_supercell_dimensions():
    dims_test = dsa.get_supercell_dimensions("tests/TestData/vacancy/geometry.in")
    dims_verified = [14.8588928, 12.9147523, 12.33466882]
    assert dims_test == pytest.approx(dims_verified)