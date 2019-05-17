'''
User instructions:

This script has been written to allow DefectCorrectionsNotebook.ipynb to be run for a set of defect structures.
The user needs to add inputs from the 'User inputs' cell that were in the notebook file.
Below that, the user must then add in a list of the location of data for the different defects to be processed.

After making these modification, run this script from command line using 'python DefectCorrectionsDataset.py'.

'''

from NotebookScripter import run_notebook
import os

############## USER INPUTS START BELOW HERE #################

# Define directory containing all defects data
base_dir = "./sample_data"
# Info from 'User inputs' cell of notebook that apply to all defects in dataset
global_configuration = {
    "dielectric_xx": 7.49,
    "dielectric_yy": 6.92,
    "dielectric_zz": 7.19,
    "path_to_all_defects": base_dir,
    "path_to_host": './sample_data/perfect',
    "manual_cutoff": None,
    "LZ": True,
    "FNV": True,
    "atom_centered_pa": True,
    "planar_av_pa": True
}
# Dictionary storing location of data for each charged defect to be analysed 
# base_dir (defined above) path will be joined to the defect paths below to allow the dictionary to be more compact
# Dictionary key defines the name of the directory output data will be stored in for each defect and then dictionary entries are:
# location of neutral defect data, location of charged defect data, defect charge state
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


# Notebook parameters for specific charged defect being processed
configurations = []
for name, (neutral_dir, charge_dir, charge_state) in defect_dataset.items():
    # make sure required inputs exists
    path_to_defect = os.path.join(base_dir, charge_dir)
    path_to_neutral = os.path.join(base_dir, neutral_dir)  
    assert os.path.isdir(path_to_defect), 'required input directory is missing {0}'.format(path_to_defect)
    assert os.path.isdir(path_to_neutral), 'required input directory is missing {0}'.format(path_to_neutral)
    config = {
        "defect_outputs_dir": name,
        "path_to_defect": path_to_defect,
        "path_to_neutral": path_to_neutral,
        "defect_charge": charge_state
    }
    configurations.append(config)
#print(configurations)

# Use NotebookScripter to run notebook for each charged defect in turn
for config in configurations:
    try:
        notebook = run_notebook("./DefectCorrectionsNotebook.ipynb", **global_configuration, **config)
    except Exception as err:
        print("Caught error when executing notebook: {0}".format(err))

print("")
print("FINISHED")
print("Outputs for all defects processed can be found in directories ProcessedDefects/: "+ ", ".join(defect_dataset.keys()) )
print("See log.info files in each subdirectory for an overview of the analysis and corrections_summary.dat for final correction terms")
print("Have a nice day!")