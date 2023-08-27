# System requirements:
* MathWorks MATLAB 2020a or newer

# Installation guide: 
* Please copy folder +SupFun in the same directory as the matlab functions POPSIM_PreProcessing.m and POPSIM_SIMreconstruction.m  Do not change the name of that folder.

# Instructions for use:
* Download Processing_ExampleData.tar.gz from https://zenodo.org/record/8277661
* For registering the POPSIM data, run POPSIM_PreProcessing.m on the raw data file 1_CH00_000001.tif to register the different SIM direction. 
* The program will output registered datasets reg2_SIM00-1.tif to reg2_SIM22-1.tif
* For a SIM reconstruction, run POPSIM_SIMreconstruction.m  on the pre-registered data. It will create a new folder "SIM" which contains a widefield and a SIM reconstruction.

