# For 'run_CoFFEE_solver'
## Check that strings written to put into sub-process are what are expected
### def check_subprocess_strings, assert command_for_subprocess == "correct string"

# For 'write_CoFFEE_in_file' and 'run_CoFFEE_solver'
## Write dummy 'in' coffee input files that is expected to run very quickly and compare to expected result 
### def check_basic_solver, assert 


import re
from ??/CoffeeConvenienceFunctions import *



def check_subprocess_strings:

    assert command_for_subprocess == "correct string"




def test_coffee_solver:
    coffee_solver_output = 'coffee_solver_test.out'
    # Known energy output for test system
    expectedOutput = 6.5
    # Run coffee Poission solver with test system (for quick convergence)
    run_CoFFEE_solver(path_to_coffee_dir, defect_outputs_dir, super_x, super_y, super_z, defect_geom, sigma, cutoff, defect_charge, defect_x, defect_y, defect_z, dielectric_xx, dielectric_yy, dielectric_zz)

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

    assert  E_q_per_m < expectedOutput+0.1 and E_q_per_m > expectedOutput+0.1

