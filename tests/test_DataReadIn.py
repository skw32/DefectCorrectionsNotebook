import DefectSupercellAnalyses as dsa
import pytest

# Unit tests for functions
def test_read_lattice_vectors():
    vecs_test_x, vecs_test_y, vecs_test_z = dsa.read_lattice_vectors("tests/TestData/interstitial/geometry.in")
    vecs_verified_x, vecs_verified_y, vecs_verified_z = [14.8588928, 0.00046026, -0.00010556], [-0.00093786, 12.9147523, -0.00021786], [-0.00064498, -0.0010864, 12.33466882]
    assert vecs_test_x == pytest.approx(vecs_verified_x)
    assert vecs_test_y == pytest.approx(vecs_verified_y)
    assert vecs_test_z == pytest.approx(vecs_verified_z)
'''
def test_read_atom_coords():

    assert == pytest.approx()
'''

# Testing functions that use other functions
def test_get_supercell_dimensions():
    dims_test = dsa.get_supercell_dimensions("tests/TestData/vacancy/geometry.in")
    dims_verified = [14.8588928, 12.9147523, 12.33466882]
    assert dims_test == pytest.approx(dims_verified)
'''
def test_count_species():
    
    assert == pytest.approx()
'''
def test_find_defect_type():
    host_coords = dsa.read_atom_coords("tests/TestData/perfect/geometry.in")
    vacancy_coords = dsa.read_atom_coords("tests/TestData/vacancy/geometry.in")
    antisite_coords = dsa.read_atom_coords("tests/TestData/antisite/geometry.in")
    interstitial_coords = dsa.read_atom_coords("tests/TestData/interstitial/geometry.in")
    vacancy_test = dsa.find_defect_type(host_coords, vacancy_coords)
    antisite_test = dsa.find_defect_type(host_coords, antisite_coords)
    interstitial_test = dsa.find_defect_type(host_coords, interstitial_coords)
    assert vacancy_test == 'vacancy'
    assert antisite_test == 'antisite'
    assert interstitial_test == 'interstitial'
'''
def test_find_vacancy():

    assert == pytest.approx()
def test_find_interstitial():

    assert == pytest.approx()
def test_find_antisite():

    assert == pytest.approx()
def test_vacancy_coords():

    assert == pytest.approx()
def test_interstitial_coords():

    assert == pytest.approx()
def test_antisite_coords():

    assert == pytest.approx()
def test_defect_to_boundary():

    assert == pytest.approx()
'''
    