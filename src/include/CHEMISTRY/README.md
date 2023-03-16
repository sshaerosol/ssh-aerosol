

# Instructions to prepare files for a new chemical kinetic mechanism


## Prepare the input files

- species list

it should have the postfix ".species"

- reactions list

it should have the postfix ".reactions"

## Location of the input files :

- create a folder

> mkdir src/include/CHEMISTRY/new_mechanism

- put the file for species list and the file for reactions list into the new folder
 

## Compile with the new chemical kinetic mechanism

> compile -c=new_mechanism

