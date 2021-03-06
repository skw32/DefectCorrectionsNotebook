Tests in this directory can be run with `pytest` from the top level directory of this repository, or, with `pytest --cov-config .coveragerc --cov=. tests/` for test coverage information.

The different test scripts include:
- test_DataReadIn.py: this script contains tests for functions in 'DefectSupercellAnalyses.py' in the previous directory which check that data from FHI-aims structure files is being read in as expected using sample data contained in the 'TestData' directory and also tests that the function used to configure the log file written to during the notebook workflow is being called without error.
- test_CoffeeIntegration.py: this script contains tests for functions in 'CoffeeConvenienceFunctions.py' that automatically generate inputs for coffee.py and that instructions to run coffee.py are implemented correctly.
- test_Notebook.py: this involves integration tests for multiple operations performed within the 'DefectCorrectionsNotebook.ipynb' for energy correction terms for sample defects. Sample data is inputted into the notebook and parameters (such as energy correction terms) are extracted using the NotebookScripter library.


The 'TestData' directory contains sample FHI-aims structure files (geometry.in) for a perfect supercell of Cu3AsS4 and defect supercells: S vacancy, Cu interstitial, As-on-Cu antisite, all with a +1 charge state. The dataset is fixed to ensure it is valid to compare new test results to the pre-computed values.