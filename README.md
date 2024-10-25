What is SSH-aerosol
===================

SSH-aerosol model represents the physico chemical transformation undergone by
aerosols in the troposphere.

Installation
============

Building
--------

> ./compile

scons is needed to compile the program. Makefile for make is also provided. It
is a preliminary version and therefore is not fully tested. 

Run
---

> ./ssh-aerosol namelist.ssh


Input files
-----------

- namelist.ssh

  This file is the main configuration file.

- species-list/species-list.dat

  This file contains the list of gas-phase species.

- species-list/species-list-aer.dat

  This file contains the list of aerosol species.

- inputs/inputs-XXX/init_gas_XXX.dat

  This file contains initial concentrations for gas-phase species.

- inputs/inputs-XXX/init_aero_XXX.dat

  This file contains initial mass concentrations for aerosol species.

- inputs/inputs-XXX/init_num_XXX.dat

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
