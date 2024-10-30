

# Instructions to prepare files for a new chemical kinetic mechanism


## Prepare the input files

- species list

it is recommended to use the postfix ".species"

- reactions list

it is recommended to use the postfix ".reactions"

## Location of the input files :

- create a folder

> mkdir src/include/CHEMISTRY/new_mechanism

- put the file for species list and the file for reactions list into the new folder
 
- modify namelist file using the path to the new input files

&gas_phase_species
species_list_file = "./src/include/CHEMISTRY/new_folder/new.species",
reaction_list_file =    "./src/include/CHEMISTRY/new_folder/new.reactions",



