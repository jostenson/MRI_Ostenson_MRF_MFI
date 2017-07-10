# SETUP #
MRI Ostenson et al., MFI for Spiral Deblurring
=============================

Requirements:
-------------
* MATLAB  - for reconstruction, processing, and figure generation
* Python - for file download and SHA1 checksum verification of data files and publicly available contributed code
* C compiler - for generating MEX files from contributing code

Tested Configuration:
---------------------
* Mac OS X 10511.5 (El Capitan)
* MATLAB R2015a
* Python 2.7.13

Installation Options:
---------------------
* Click the `Download ZIP` button on the lower right hand side of the repository to download the code to your local machine
* OR clone the git repository to your local machine
* OR fork to your own repository and then clone the git repository to your local machine

Usage:
------
* Go to `./code/` and run `batch_proc.m`
* Code is not optimized for speed, the total processing time is ~48 hours

Folder Structure:
--------

`./code/` - (with downloaded contributors) contains all code necessary to reconstruct, process, and generate figures
`./data_in/` - the data input directory
`./data_out/` - the reconstruction and processing output directory
`./figures/` - the figure output directory

