# A Single-Objective Light-Sheet Microscope with 200 nm-Scale Resolution.

## Abstract
We present a single-objective light-sheet microscope, also known as an oblique-plane microscope, that uses a bespoke glass-tipped tertiary objective and improves the resolution, field of view, usability, and stability over previous variants. Owing to its high numerical aperture optics, this microscope achieves the highest lateral resolution in light-sheet fluorescence microscopy, and its axial resolution is similar to that of Lattice Light-Sheet Microscopy. Given this performance, we demonstrate high-resolution imaging of clathrin-mediated endocytosis, vimentin, the endoplasmic reticulum, membrane dynamics, and natural killer cell-mediated cell death. Furthermore, we image biological phenomena that would be otherwise challenging or impossible to perform in a traditional light-sheet microscope geometry, including cell migration through a confined space within a microfluidic device, photoactivation of PI3K, and diffusion of cytoplasmic rheological tracers at a volumetric rate of 14 Hz.

## BioRxiv Preprint

## Authors
Etai Sapoznik (1,2), Bo-Jui Chang (1), Robert J. Ju (3), Erik S. Welf (1,2), David Broadbent (4), Alexandre F. Carisey (5), Samantha J. Stehbens (3), Kyung-min Lee (6), Arnaldo Marín (6), Ariella B. Hanker (6), Jens C. Schmidt (4,7), Carlos L. Arteaga (6), Bin Yang (8), Rory Kruithoff (9), Doug P. Shepherd (9), Alfred Millett-Sikking (10), Andrew G. York (10), Kevin M. Dean (1*), Reto Fiolka (1,2*)

## Affiliations
* 1 – Department of Cell Biology, University of Texas Southwestern Medical Center, Dallas, TX, USA.
* 2 – Department of Bioinformatics, University of Texas Southwestern Medical Center, Dallas, TX, USA.
* 3 – Institute for Molecular Bioscience, University of Queensland, St Lucia, Queensland, Australia.
* 4 – Institute for Quantitative Health Sciences and Engineering, Michigan State University, East Lansing, MI, USA.
* 5 – William T. Shearer Center for Human Immunobiology, Baylor College of Medicine and Texas Children’s Hospital, Houston, TX, USA.
* 6 – Harold C. Simmons Comprehensive Cancer Center and the Department of Internal Medicine, University of Texas Southwestern Medical Center, Dallas, TX, USA.
* 7 - Department of Obstetrics, Gynecology, and Reproductive Biology, Michigan State University, East Lansing, MI, USA.
* 8 – Chan Zuckerberg Biohub, San Francisco, CA, USA.
* 9 – Department of Physics and the Center for Biological Physics, Arizona State University, Tempe, AZ, USA.
* 10 – Calico Life Sciences LLC, South San Francisco, CA, USA.

## Correspondence
* Kevin.Dean@utsouthwestern.edu
* Reto.Fiolka@utsouthwestern.edu


## Software
Shearing and Deconvolution Routines for OPM (Snouty)

Python and MATLAB routines for shearing and deconvolving data, respectively.  

The shearing program was developed in Python 3.6 operating in an Anaconda environment on a Linux operating system.  Changes would need to be made in order to use it on a Windows Machine.  Use the code at your own risk, and expect bugs.  Changes will be made in the future to improve readability and reduce redundancy, but time is of the essence during manuscript submission.

Our data acquisition software saves data according to the following structure:  .../cellType/label/date/cell#/1_CH0#_######.tif.  Thus, all of the information is provided in the file path itself.  A hypothetical example is MCF7/AktPH-GFP/200218/Cell5/1_CH00_000003.tif.

The deskewing software is designed to be operated at the date level of the path.  It then goes through that directory, identifies the number of Cell# subdirectories, and then within those directories, the number of channels, and timepoints.  It then deskews all of these files in a parallel process.  You execute the software as follows: python deskewDirectory.py MCF7/AktPH-GFP/200218/

The lateral pixel size, z step size, oblique illumination angle, and number of parallel threads, is all specified within the function itself.

The deconvolution software is used on the raw, non-sheared data.  It uses an experimentally measured PSF as a prior, and updates this PSF with a double blind deconvolution process.  This code is written in MATLAB.  The experimentally measured PSF needs to have the same voxel dimensions as the data that was acquired, or be appropriately scaled.  After deconvolution, the image is then sheared using the aforementioned Python code.

- Kevin Dean
