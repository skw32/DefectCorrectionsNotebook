'''
User instructions:

This script has been written to allow DefectCorrectionsNotebook.ipynb to be run for a set of defect structures.
The user needs to add inputs from the 'User inputs' cell that were in the notebook file.
Below that, the user must then add in a list for the location of data for the different defects to be processed.

After making these modification, run this script from command line using 'python DefectCorrectionsDataset.py'.

'''

from NotebookScripter import run_notebook_in_jupyter
import os
import shutil
from ast import literal_eval
import PyladaDefectsImageCharge as pdic
import numpy as np

############## USER INPUTS START BELOW HERE #################

# Option to speed up analysis when calculating LZ image-charge corrections for dataset.
# The term can be pre-computed if you have several defects in your dataset with the same charge state.
# If you would like to do this, set 'precompute_LZ' to True, uncomment the next set of lines adding in your dielectric constants, lattice vectors and list of defect charge states.
precompute_LZ = False
'''
dielectric_xx = 7.49
dielectric_yy = 6.92
dielectric_zz = 7.19
latt_vecs = np.array([[1.48588928e+01, -9.37860000e-04, -6.44980000e-04], [4.60260000e-04, 1.29147523e+01, -1.08640000e-03], [-1.05560000e-04, -2.17860000e-04, 1.23346688e+01]])
dielectric_av = (dielectric_xx+dielectric_yy+dielectric_zz)/3.0
charge_list = [-2, -1, 1, 2]
precomputed_IC = []
for precompute in charge_list:
    # Using modified function from pylada-defects to compute LZ image charge correction
    Madelung_energy, ThirdOrder, csh, scaling_f, E_corr_LZ = pdic.get_imagecharge(np.transpose(latt_vecs), charge=precompute, epsilon=dielectric_av,cutoff=50.)
    precomputed_IC.append(E_corr_LZ)
assert (len(charge_list) == len(precomputed_IC)), 'Error pre-computing LZ IC term.'
'''

# Define directory containing all defects data
base_dir = "./sample_data/relaxed_defects/Cu3AsS4"
# Info from 'User inputs' cell of notebook that apply to all defects in dataset
global_configuration = {
    "dielectric_xx": 7.49,
    "dielectric_yy": 6.92,
    "dielectric_zz": 7.19,
    "path_to_all_defects": base_dir,
    "path_to_host": './sample_data/relaxed_defects/Cu3AsS4/perfect',
    "LZ": False,
    "FNV": True,
    "atom_centered_pa": True,
    "planar_av_pa": False
}

# Dictionary storing location of data for each charged defect to be analysed 
# base_dir path (defined above) will be joined to the defect paths below to allow the dictionary to be more compact
# Dictionary key defines the name of the directory output data will be stored in for each defect and then dictionary entries are:
# location of neutral defect data, location of charged defect data, defect chargse state
defect_dataset = {
    "As-on-Cu_q=+1_sg=1": [
        # neutral_dir, charge_dir, charge_state
        "neutral_antisite", "charged_antisite", 1
    ],
    "V-Cu_q=-1_sg=1": [
        # neutral_dir, charge_dir, charge_state
        "neutral_vacancy", "charged_vacancy", -1
    ],
    "Cu_i_q=+1_pore1": [
        # neutral_dir, charge_dir, charge_state
        "neutral_interstitial", "charged_interstitial", 1
    ]
}

################# END OF USER INPUTS #####################


# Notebook parameters for specific charged defect to be processed
configurations = []
for name, (neutral_dir, charge_dir, charge_state) in defect_dataset.items():
    # make sure required inputs exists
    path_to_defect = os.path.join(base_dir, charge_dir)
    path_to_neutral = os.path.join(base_dir, neutral_dir)  
    assert os.path.isdir(path_to_defect), 'required input directory is missing {0}'.format(path_to_defect)
    assert os.path.isdir(path_to_neutral), 'required input directory is missing {0}'.format(path_to_neutral)
    if (precompute_LZ == True):
        for i in charge_list:
            if (charge_state == charge_list[i]):
                precomputed_val = precomputed_IC[i]
        config = {
            "defect_outputs_dir": name,
            "path_to_defect": path_to_defect,
            "path_to_neutral": path_to_neutral,
            "defect_charge": charge_state,
            "IC": precomputed_val
        }
    else:
        config = {
            "defect_outputs_dir": name,
            "path_to_defect": path_to_defect,
            "path_to_neutral": path_to_neutral,
            "defect_charge": charge_state,
            "IC": None
        }
    configurations.append(config)
#print(configurations)

# Use NotebookScripter to run notebook for each charged defect in turn
for config in configurations: 
    try:
        notebook = run_notebook_in_jupyter("./DefectCorrectionsNotebook.ipynb", **global_configuration, **config,timeout=1200)("defect_outputs_dir", save_output_notebook="output_DefectCorrectionsNotebook.ipynb")
        shutil.move("./output_DefectCorrectionsNotebook.ipynb", os.path.join(notebook.defect_outputs_dir, "output_DefectCorrectionsNotebook.ipynb"))
    except Exception as err:
        print("Caught error when executing notebook: {0}".format(err))

print("")
print("FINISHED")
print("Outputs for all defects processed can be found in directories ProcessedDefects/: "+ ", ".join(defect_dataset.keys()) )
print("See log.info files in each subdirectory for an overview of the analysis and corrections_summary.dat for final correction terms.")
print("Have a nice day!")