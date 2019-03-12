# Add functions to check sample FHI-aims structure data is being read correctly by DefectSupercellAnalyses.py functions??

# Ideas
## Add dummy data: structure files
### Sample perfect host, sample: vacancy, antisite and interstitial

## Unit tests
### Check that values read for aims dummy structure by functions are a given value (within set tolerance)

## Integration tests
### Check that functions using other functions gives a certain value within set tolerance??

## End-to-end tests
### Check Ecor outputted from workflow for dummy data is a given value (within set tolerance)
### Check that plots are generated --> pytest with raise exception errors for this when trying to run notebook for above 
### (check that this is not suppressed in scripter, if so change implementation for test)
#### for E_iso
#### for each potential alignment step