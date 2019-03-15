**To-do: hook up with Travis-CI or circle-CI?? and link tests to functions in main dir**

Tests in this directory can be run with `pytest ChoosenScriptName.py`

The different test scripts include:
- testDataReadIn.py: this script contains tests for functions in 'DefectSupercellAnalyses.py' in the previous directory which check that data from FHI-aims structure files is being read in as expected using sample data contained in the 'TestData' directory.
- testCoffeeIntegration.py: this script contains tests for functions in 'CoffeeConvenienceFunctions.py' that automatically generate inputs for coffee.py and that instructions to run coffee.py are implemented correctly.
- testNotebook.ipynb: this involves integration tests for multiple operations performed within the 'DefectCorrectionsNotebook.ipynb' and an end-to-end test for final correction energy outputted for the charge model for a sample defect. **check compatibility of pytest with notebooks** - or extract functions from notebook with NotebookScripter and run tests this way??


The 'TestData' directory contains sample FHI-aims structure files (geometry.in) for a perfect supercell of Cu3AsS4 and defect supercells: S vacancy, Cu interstitial, As-on-Cu antisite, all with a +1 charge state.