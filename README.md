** UNDER CONSTRUCTION **

# DefectCorrectionsNotebook: Overview

This project contains a workflow which aims to allow for convenient and explanatory post-processing of charged defect supercells from electronic structure calculations with the all-electron electronic structure software package [FHI-aims](https://aimsclub.fhi-berlin.mpg.de/). There are two major components to the workflow:

1. **DefectCorrectionsNotebook.ipynb** contains a tutorial notebook for performing a finite-size correction scheme for charged defect supercells for one defect at a time. More information on the correction scheme and processing steps is contained in the notebook. 

2. **DefectCorrectionsDataset.py** allows for running the notebook from the command line for a set of defect supercells data.

Other components of this project include:
- DefectSupercellAnalysis.py: This contains functions used in the notebook workflow that were  written for reading in structural information of defect supercells in the geometry file format of FHI-aims ('geometry.in').
- CoffeeConvenienceFunctions.py: This contains functions used in the notebook workflow for easy integration with the [CoFFEE](https://www.sciencedirect.com/science/article/pii/S0010465518300158) software package for performing processing steps to charged defect supercells.
- LogFileSetup.py: This file contains the default format of the log file used to store intermediate processing results from the notebook.
- WorkflowTests: This directory contains tests for functions written for this workflow and sample data to use with the tests.

## Installation instructions 
This information is contained within the notebook (DefectCorrectionsNotebook.ipynb) but it repeated here for completeness.

This workflow uses python3. The most convenient way to setup the python environment for this workflow is to use [Anaconda](https://www.anaconda.com/distribution/). All dependencies present when testing this workflow are listed in DefectCorrections_conda_env.yml. This environment can be re-created using conda with 

`conda env create -n chooseYourEnvName --file DefectCorrections_conda_env.yml` 

To use this workflow you must then activate this conda environment with 

`conda activate YourChosenEnvName`

## License
**To add** + refer to use of CoFFEE which has a BSD 3-clause license

## Code contributions
Contributions to extend the functionality of this workflow are very much welcomed! For this, we welcome contributors to fork their own copy of this repository for local developments before submitting a pull request to merge with this repository. See [here](https://guides.github.com/activities/forking/) for more details on this procedure.

## Authors
- Suzanne K. Wallace
- Tong Zhu

## Contact
For any queries or bugs to report, please contact suzywallace501@gmail.com

## Acknowledgements
- Volker Blum (Duke University)
- Nathaniel Cohen for support in the development of this workflow and integration with the [NotebookScripter](https://github.com/breathe/NotebookScripter) library
- Extensive use has been made of the [CoFFEE](https://www.sciencedirect.com/science/article/pii/S0010465518300158) software package for implementing finite-size corrections to defect supercells
- Many useful discussions on defect correction techniques were provided by [Stephan Lany, Anuj Goyal and Prashun Gorai](https://github.com/pylada/pylada-defects) (NREL)

