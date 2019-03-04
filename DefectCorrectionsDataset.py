# Script must be run from command line as 'ipython DefectCorrectionsDataset.py'

from NotebookScripter import run_notebook
import os

############## USER INPUTS ##############

# Define directory containing all defects data
base_dir = "/Users/suzannewallace/PhD/PhD_year3/PhDFinalProjects/DefectAnalysis/EnargiteDefects/fromLandau/Results/data/final_one_shots"
# Info from 'User inputs' cell of notebook that apply to all defects in dataset
global_configuration = {
    "dielectric_xx": 7.49,
    "dielectric_yy": 6.92,
    "dielectric_zz": 7.19,
    "path_to_coffee_dir": '/Users/suzannewallace/PhD/PhD_year3/PhDFinalProjects/DefectAnalysis/CoFFEE_1.1',
    "path_to_all_defects": base_dir,
    "path_to_host": '/Users/suzannewallace/PhD/PhD_year3/PhDFinalProjects/DefectAnalysis/EnargiteDefects/fromLandau/Results/data/final_one_shots/PerfectReference',
    "charge_model_file": 'charge_model.dat',
    "pa_plot_file": 'pa_plot.png',
    "manual_cutoff": None
}
# Dictionary storing location of data for each charged defect to be analysed
# Dictionary key defines name of directory output data will be stored in for each defect
# Location of neutral defect data, charged defect data and defect charge state must be inputted for each defect that is to be processed
defect_dataset = {
    "V-S_q=+1_sg=1": [
        # neutral_dir, charge_dir, charge_state
        "VacancySupercells/V-S/neutral/DefectSpacegroup1", "VacancySupercells/V-S/charged/+1/DefectSpacegroup1", 1
    ],
    "V-S_q=+2_sg=1": [
        # neutral_dir, charge_dir, charge_state
        "VacancySupercells/V-S/neutral/DefectSpacegroup1", "VacancySupercells/V-S/charged/+2/DefectSpacegroup1", 2
    ],
    "V-S_q=+1_sg=6_s1": [
        # neutral_dir, charge_dir, charge_state
        "VacancySupercells/V-S/neutral/DefectSpacegroup6/structure1", "VacancySupercells/V-S/charged/+1/DefectSpacegroup6/structure1", 1
    ],
    "V-S_q=+1_sg=6_s2": [
        # neutral_dir, charge_dir, charge_state
        "VacancySupercells/V-S/neutral/DefectSpacegroup6/structure2", "VacancySupercells/V-S/charged/+1/DefectSpacegroup6/structure2", 1
    ],
    "V-S_q=+2_sg=6_s1": [
        # neutral_dir, charge_dir, charge_state
        "VacancySupercells/V-S/neutral/DefectSpacegroup6/structure1", "VacancySupercells/V-S/charged/+2/DefectSpacegroup6/structure1", 2
    ],
    "V-S_q=+2_sg=6_s2": [
        # neutral_dir, charge_dir, charge_state
        "VacancySupercells/V-S/neutral/DefectSpacegroup6/structure2", "VacancySupercells/V-S/charged/+2/DefectSpacegroup6/structure2", 2
    ],
    "V-Cu_q=-1_sg=1": [
        # neutral_dir, charge_dir, charge_state
        "VacancySupercells/V-Cu/neutral/DefectSpacegroup1", "VacancySupercells/V-Cu/charged/-1/DefectSpacegroup1", -1
    ],
    "V-Cu_q=-1_sg=6": [
        # neutral_dir, charge_dir, charge_state
        "VacancySupercells/V-Cu/neutral/DefectSpacegroup6", "VacancySupercells/V-Cu/charged/-1/DefectSpacegroup6", -1
    ],
    "V-As_q=-1_sg=6": [
        # neutral_dir, charge_dir, charge_state
        "VacancySupercells/V-As/neutral/DefectSpacegroup6", "VacancySupercells/V-As/charged/-1/DefectSpacegroup6", -1
    ],
    "V-As_q=-2_sg=6": [
        # neutral_dir, charge_dir, charge_state
        "VacancySupercells/V-As/neutral/DefectSpacegroup6", "VacancySupercells/V-As/charged/-2/DefectSpacegroup6", -2
    ],
    "V-As_q=-3_sg=6": [
        # neutral_dir, charge_dir, charge_state
        "VacancySupercells/V-As/neutral/DefectSpacegroup6", "VacancySupercells/V-As/charged/-3/DefectSpacegroup6", -3
    ],
    "V-As_q=-4_sg=6": [
        # neutral_dir, charge_dir, charge_state
        "VacancySupercells/V-As/neutral/DefectSpacegroup6", "VacancySupercells/V-As/charged/-4/DefectSpacegroup6", -4
    ],
    "V-As_q=-5_sg=6": [
        # neutral_dir, charge_dir, charge_state
        "VacancySupercells/V-As/neutral/DefectSpacegroup6", "VacancySupercells/V-As/charged/-5/DefectSpacegroup6", -5
    ],
    "As_Cu_q=+1_sg=1": [
        # neutral_dir, charge_dir, charge_state
        "AntisiteSupercells/As-Cu/neutral/DefectSpacegroup1", "AntisiteSupercells/As-Cu/charged/+1/DefectSpacegroup1", 1
    ],
    "As_Cu_q=+2_sg=1": [
        # neutral_dir, charge_dir, charge_state
        "AntisiteSupercells/As-Cu/neutral/DefectSpacegroup1", "AntisiteSupercells/As-Cu/charged/+2/DefectSpacegroup1", 2
    ],
    "As_Cu_q=+3_sg=1": [
        # neutral_dir, charge_dir, charge_state
        "AntisiteSupercells/As-Cu/neutral/DefectSpacegroup1", "AntisiteSupercells/As-Cu/charged/+3/DefectSpacegroup1", 3
    ],
    "As_Cu_q=+4_sg=1": [
        # neutral_dir, charge_dir, charge_state
        "AntisiteSupercells/As-Cu/neutral/DefectSpacegroup1", "AntisiteSupercells/As-Cu/charged/+4/DefectSpacegroup1", 4
    ],
    "As_Cu_q=+1_sg=6": [
        # neutral_dir, charge_dir, charge_state
        "AntisiteSupercells/As-Cu/neutral/DefectSpacegroup6", "AntisiteSupercells/As-Cu/charged/+1/DefectSpacegroup6", 1
    ],
    "As_Cu_q=+2_sg=6": [
        # neutral_dir, charge_dir, charge_state
        "AntisiteSupercells/As-Cu/neutral/DefectSpacegroup6", "AntisiteSupercells/As-Cu/charged/+2/DefectSpacegroup6", 2
    ],
    "As_Cu_q=+3_sg=6": [
        # neutral_dir, charge_dir, charge_state
        "AntisiteSupercells/As-Cu/neutral/DefectSpacegroup6", "AntisiteSupercells/As-Cu/charged/+3/DefectSpacegroup6", 3
    ],
    "As_Cu_q=+4_sg=6": [
        # neutral_dir, charge_dir, charge_state
        "AntisiteSupercells/As-Cu/neutral/DefectSpacegroup6", "AntisiteSupercells/As-Cu/charged/+4/DefectSpacegroup6", 4
    ],

    "Cu_As_q=-1_sg=6": [
        # neutral_dir, charge_dir, charge_state
        "AntisiteSupercells/Cu-As/neutral/DefectSpacegroup6", "AntisiteSupercells/Cu-As/charged/-1/DefectSpacegroup6", -1
    ],
    "Cu_As_q=-2_sg=6": [
        # neutral_dir, charge_dir, charge_state
        "AntisiteSupercells/Cu-As/neutral/DefectSpacegroup6", "AntisiteSupercells/Cu-As/charged/-2/DefectSpacegroup6", -2
    ],
    "Cu_As_q=-3_sg=6": [
        # neutral_dir, charge_dir, charge_state
        "AntisiteSupercells/Cu-As/neutral/DefectSpacegroup6", "AntisiteSupercells/Cu-As/charged/-3/DefectSpacegroup6", -3
    ],
    "Cu_As_q=-4_sg=6": [
        # neutral_dir, charge_dir, charge_state
        "AntisiteSupercells/Cu-As/neutral/DefectSpacegroup6", "AntisiteSupercells/Cu-As/charged/-4/DefectSpacegroup6", -4
    ]
}

############## END OF USER INPUTS ##############


# Notebook parameters for specific charged defect being processed
configurations = []
for name, (neutral_dir, charge_dir, charge_state) in defect_dataset.items():
    # make sure required inputs exists
    path_to_defect = os.path.join(base_dir, charge_dir)
    path_to_neutral = os.path.join(base_dir, neutral_dir)  
    assert os.path.exists(path_to_defect), 'required input directory is missing {0}'.format(path_to_defect)
    assert os.path.exists(path_to_neutral), 'required input directory is missing {0}'.format(path_to_neutral)
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

print("FINISHED")
print("Outputs for all defects processed can be found in directories ProcessedDefects/: "+ ", ".join(defect_dataset.keys()) )