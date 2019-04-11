[![Build Status](https://travis-ci.org/skw32/DefectCorrectionsNotebook.svg?branch=master)](https://travis-ci.org/skw32/DefectCorrectionsNotebook)
[![Coverage Status](https://coveralls.io/repos/github/skw32/DefectCorrectionsNotebook/badge.png?branch=master)](https://coveralls.io/github/skw32/DefectCorrectionsNotebook?branch=master)

** UNDER CONSTRUCTION **

To-do:

- Add new FNV alignment step, option for LZ Ecor and final KO pa step
- Finish writing paper.md (with sample outputs)
- Add tests for newer potential alignment steps and LZ Ecor
- Compare workflow results to Tong's defects processed by hand

# DefectCorrectionsNotebook: Overview

This project contains a workflow which aims to allow for convenient and explanatory post-processing of charged defect supercells from electronic structure calculations with the all-electron electronic structure software package [FHI-aims](https://aimsclub.fhi-berlin.mpg.de/). There are two major components to the workflow:

1. **DefectCorrectionsNotebook.ipynb** contains a tutorial notebook for performing a finite-size correction scheme for charged defect supercells for one defect at a time. More information on the correction scheme and processing steps is contained in the notebook.

2. **DefectCorrectionsDataset.py** allows for running the notebook from the command line for a set of defect supercells data.

Other components of this project include:

- DefectSupercellAnalysis.py: This contains functions used in the notebook workflow that were written for reading in structural information of defect supercells in the geometry file format of FHI-aims ('geometry.in').
- CoffeeConvenienceFunctions.py: This contains functions used in the notebook workflow to automatically generate input files for the [CoFFEE](https://www.sciencedirect.com/science/article/pii/S0010465518300158) software package for performing processing steps for applying corrections to charged defect supercells.
- LogFileSetup.py: This file contains the default format of the log file used to store intermediate processing results from the notebook.
- coffee.py: This is the main executable for the [CoFFEE](https://www.sciencedirect.com/science/article/pii/S0010465518300158) package (from version 1.1) that is used in this workflow.
- tests: This directory contains tests for functions written for this workflow and sample data to use with the tests.
- conftest.py is a configuration file for pytest to define the root directory for the testsuite.
- .travis.yml and .coveragerc are files required for Travis-CI and test coverage respectively.

## Installation instructions

First download a copy of the repository. All analysis should be performed from the directory containing the notebook file (DefectCorrectionsNotebook.ipynb). The following information is also contained within the notebook, but is repeated here for completeness.

This workflow uses python3. The most convenient way to setup the python environment for this workflow is to use [Anaconda](https://www.anaconda.com/distribution/). 


### Option 1: Create conda env directly from .yml file
All dependencies present when testing this workflow are listed in DefectCorrections_conda_env.yml. This environment can be re-created using conda with 

`conda env create -n chooseYourEnvName --file DefectCorrections_conda_env.yml` 

To use this workflow you must then activate this conda environment with 

`conda activate YourChosenEnvName`

### Option 2: Create conda env manually

Alternatively, a conda environment can be created for this workflow and the necessary packages installed one at a time, as outlined below

`conda create -n chooseYourEnvName python=3`

`conda activate YourChosenEnvName`

`conda install -c suzannekwallace coffee_poisson_solver_ko`

`conda install jupyter pytest pytest-cov`

`pip install NotebookScripter coveralls`

## License

This software has been made available under a 3-Clause BSD License.

## Code contributions

Contributions to extend the functionality of this workflow are very much welcomed! For this, we welcome contributors to fork their own copy of this repository for local developments before submitting a pull request to merge with this repository. See [here](https://guides.github.com/activities/forking/) for more details on this procedure.

## Relevant citations

Methods used in this workflow are taken from a number of publications and so these should be cited when making use of this workflow. These include doi: 10.1088/0965-0393/17/8/084002, doi: 10.1103/PhysRevLett.102.016402, doi: 10.1016/j.cpc.2018.01.011 and doi: 10.1103/PhysRevB.89.195205.

## Authors

- Tong Zhu
- Suzanne K. Wallace

## Contact

For any queries or bugs to report, please contact suzywallace501@gmail.com

## Acknowledgements

- Volker Blum (Duke University)
- Nathaniel Cohen for support in the development of this workflow and for the [NotebookScripter](https://github.com/breathe/NotebookScripter) library
- Extensive use has been made of the [CoFFEE](https://www.sciencedirect.com/science/article/pii/S0010465518300158) software package for implementing finite-size corrections to defect supercells
- Many useful discussions on defect correction techniques were provided by [Stephan Lany, Anuj Goyal and Prashun Gorai](https://github.com/pylada/pylada-defects) (NREL)
