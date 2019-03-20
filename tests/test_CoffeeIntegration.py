import CoffeeConvenienceFunctions as ccf
import re
import os
import pytest

def test_coffee_solver():
    # Known energy output for test system (2x2x2 supercell of antisite test defect supercell)
    expectedEnergy = 0.4008
    coffee_solver_output = 'cm_2x2x2.out'
    defect_outputs_dir = 'ProcessedDefects'
    os.makedirs(defect_outputs_dir, exist_ok=True)
    coffee_exe = "coffee.py"
    # Parameters for known test system
    super_x, super_y, super_z = 2, 2, 2
    defect_geom = 'tests/TestData/antisite/geometry.in'
    sigma = 2.108966292
    cutoff = 6.773799249378561
    defect_charge = 1
    defect_x, defect_y, defect_z = 5.58437198, 8.56614992, 6.21005598
    dielectric_xx, dielectric_yy, dielectric_zz = 7.49, 6.92, 7.19
    # Run coffee Poission solver with test system (for quick convergence)
    ccf.run_CoFFEE_solver(coffee_exe, defect_outputs_dir, super_x, super_y, super_z, defect_geom, sigma, cutoff, defect_charge, defect_x, defect_y, defect_z, dielectric_xx, dielectric_yy, dielectric_zz)
    try:
        with open(os.path.join(defect_outputs_dir, coffee_solver_output), 'r') as f:
            for line in f:
                if re.search('!', line):
                    words = line.split()
                    if  words[1:4] == ['Total', 'Energy', '(eV):']:
                        testEnergy = float(words[4])
                    if line == None:
                        print('Warning! - Error finding E_q^{per,m} from '+str(os.path.join(defect_outputs_dir, coffee_solver_output)))
                    else:
                        print('Check '+str(os.path.join(defect_outputs_dir, coffee_solver_output))+' for error messages')
    except IOError:
        print('Could not open '+str(os.path.join(defect_outputs_dir, coffee_solver_output))) 

    assert  testEnergy == pytest.approx(expectedEnergy)