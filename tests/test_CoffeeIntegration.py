import CoffeeConvenienceFunctions as ccf
import re
import os
import pytest

'''
def test_coffee_solver:
    # Known energy output for test system
    expectedEnergy = 0.4008
    coffee_solver_output = 'coffee_solver_test.out'
    # Parameters for known test system
    sigma = 2.108966292
    # Run coffee Poission solver with test system (for quick convergence)
    run_CoFFEE_solver(path_to_coffee_dir, defect_outputs_dir, 2, 2, 2, "tests/TestData/antisite/geometry.in", sigma, cutoff, 1, defect_x, defect_y, defect_z, dielectric_xx, dielectric_yy, dielectric_zz)

    # CHECK FOR NO ERROR MESSAGES!


    try:
        with open(coffee_solver_output), 'r') as f:
            for line in f:
                if re.search('!', line):
                    words = line.split()
                    E_q_per_m = float(words[4])
                    if line == None:
                        print('Warning! - Error finding E_q^{per,m} from '+str(coffee_solver_output))
    except IOError:
        print("Could not open "+str(coffee_solver_output)) 

    assert  E_q_per_m == pytest.approx(expectedEnergy)
    '''