What is SSH-aerosol
===================

SSH-aerosol model represents the physico chemical transformation undergone by
aerosols in the troposphere.

Installation
============

Building
--------

> compile

The arguement for the chemical mechanism can be used. 
Default mechanism is cb05 and the available mechanisms are cb05, racm and racm2.
    
-c or --chemisty
    
> compile -c=cb05

or 

> compile --chemistry=cb05

scons is needed to compile the program. Makefile for make is also provided. It
is a preliminary version and therefore is not fully tested. 

Run
---

> ./ssh-aerosol INIT/namelist_xxx.ssh


Input files
-----------

- namelist_xxx.ssh

This file is the main configuration file.

- species-list-cb05.dat

This file contains the list of gas-phase species.

- species-list-aer.dat

This file contains the list of aerosol species.

- init_gas.dat

This file contains initial concentrations for gas-phase species.

- init_aero.dat

This file contains initial mass concentrations for aerosol species.

- init_num.dat

This file contains initial number concentrations for aerosol species.

- coef_s1_f1_b50.nc

This file contains an example coagulation coefficient data.


Download input data
===================

Photolysis
----------

Please follow README in the "photolysis" sub-directory to download input data.
photolysis-cb05.dat is provided for the indices of photolysis reactions.
It should not be modified.


Implementation examples
=======================

Using shared library: see ./example_using_so

Using Fortran90 module: see ./example_fortran
