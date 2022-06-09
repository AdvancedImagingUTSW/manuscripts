System requirements:

MathWors Matlab 2020a or newer

Installation guide:

please copy folder +SupFun in the same directory as the matlab functions OPMSIM_PreProcessing.m and SimReconRL.m  Do not change the name of that folder.

Demo:

Download

Instructions for use:

download OPSIM_ExampleData.tar.gz from https://zenodo.org/record/6481084#.YmVM-7lOmHs 

For registering the OPSIM data, run OPMSIM_PreProcessing.m on the raw data file 1CH00_00000.tif to register the different SIM views. At one point, the program will prompt you to select a region of interest 
within that dataset. The program will output registered datasets reg2_SIM00-0.tif to reg2_SIM22-0tif

For a SIM reconstruction, run SimReconRL.m on the pre-registered data. It will create a new folder "SIM2" which contains maximum intensity projections and 3D stacks of the reconstruction.

Expected run time for demo on a "normal" desktop computer:
