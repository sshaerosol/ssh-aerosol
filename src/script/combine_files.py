#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (C) 2024, CEREA
#     Author(s): Youngseob Kim
#
from __future__ import with_statement

import os
import string
import sys
from replacement_species import *
from user_defined_scheme import *
import configparser, argparse

config = configparser.ConfigParser()

parser = argparse.ArgumentParser()
parser.add_argument('--usrcfg', type = str)
parser.add_argument('--species', type = str)
parser.add_argument('--reactions', type = str)
parser.add_argument('--scheme', type = str)
parser.add_argument('--matching', type = str)
args = parser.parse_args()

Ns = None
Ns_gas = 0
Ns_aerosol = 0
species_name=[]
species_weight=[]

""" Input files

species_gas = sys.argv[1]
reactions_gas = sys.argv[2]

species_h2o = "../h2o/h2o.species"
reactions_h2o = "../h2o/h2o.reactions"

"""

species_gas = args.species # sys.argv[1]
reactions_gas = args.reactions # sys.argv[2]
scheme = args.scheme # sys.argv[3]
matching_file = args.matching

"""
Selection of schemes
"""
if scheme == "h2o":
    user_defined = False
elif scheme == "user":
    user_defined = True
    
if user_defined == True:

    config_inputfile = args.usrcfg
    config_outputfile = 'user_defined_scheme_tmp.cfg'
    with open(config_inputfile, 'r') as file:
        lines = file.readlines()

    with open(config_outputfile, 'w') as file:
        for line in lines:
            # Remove everything after the first '#' character
            new_line = line.split('#', 1)[0]
            # Remove any trailing whitespace and
            # write the line back to the file if it is not empty
            new_line = new_line.rstrip()
            
            if new_line:
                file.write(new_line + '\n')
    
    config.read(config_outputfile)

    list_schemes = config['schemes']

    options_list = {}
    for key in list_schemes:
        options_list[key] = list_schemes.getboolean(key)
    
    reactions_precursor, species_precursor, species_aerosol = \
        user_defined_scheme(options_list, config_outputfile)


else:
    species_precursor = "../h2o/h2o.species"
    reactions_precursor = "../h2o/h2o.reactions"
    species_aerosol = "../h2o/species-list-aer-h2o.dat"


species_inorganic_aerosol = "../inorganic/species-list-aer-inorg.dat"
species_water_aerosol = "../inorganic/species-list-aer-water.dat"    
print("Ozone species:", species_gas)
print("SOA species:", species_precursor)
print("SOA reactions:", reactions_precursor)

species_out = "combined_species.dat"
reactions_out = "combined_reactions.dat"

with open(species_gas,'r') as f:
    # Jumps the title line.
    # f.readline()
    for line in f:
        line = line.strip()
        if not line or line[0] in "#%!-":
            continue # A comment,
        line=line.split()
        Ns_gas = Ns_gas + 1
        # if Ns_gas is None:
        #     Ns_gas = int(line[0])
        #     continue; # The first value is the number of species.
        species_name.append(line[0])
        species_weight.append(line[1])

# Checks Ns even if it is not needed in this script.
if (Ns_gas != len(species_name)):
    print("Error " + Ns_gas + " " + len(species_name))
    sys.exit(1)
        

with open(species_precursor,'r') as f:
    # Jumps the title line.
    # f.readline()
    for line in f:
        line = line.strip()
        if not line or line[0] in "#%!-":
            continue # A comment,
        line=line.split()
        Ns_aerosol = Ns_aerosol + 1
        # if Ns_aerosol is None:
        #     Ns_aerosol = int(line[0])
        #     continue; # The first value is the number of species.
        species_name.append(line[0])
        species_weight.append(line[1])        

# Checks Ns even if it is not needed in this script.
Ns = Ns_gas + Ns_aerosol
if (Ns != len(species_name)):
    print("Error " + str(Ns) + " " + str(len(species_name)))
    sys.exit(1)


# """ Write a combined list of the species
# """
# header = "File for chemical species CB05 (name and molar mass)\n" + \
#          "# gaseous species # aqueous species\n" + \
#          str(Ns) + "  0\n" + \
#          "---Gas-phase----"

# with open(species_out, "w") as f:
#     f.write(header)
#     for k, v in zip(species_name, species_weight):
#         k += " "
#         f.write("\n")
#         f.write(k.ljust(8))
#         f.write(v)

""" Read reactions list of the mechanism.
"""

with open(reactions_gas,'r', encoding="UTF-8") as f1, open(reactions_out, 'w') as f2:
    # Jumps the title line.
    f1.readline()
    for line in f1:
        line = line.strip()
        if not line or line[0] in "#%!-":
            continue # A comment,

        if (line != "END"):
            f2.write(line)
            f2.write("\n")

with open(reactions_precursor,'r', encoding="UTF-8") as f1, open(reactions_out, 'a') as f2:

    species_matching = read_matching_file(matching_file)
    
    # Jumps the title line.
    f1.readline()
    for line in f1:
        line = line.strip()
        if not line or line[0] in "#%!-":
            continue # A comment,


        # Replace the species names
        # from those of the h2o mechanism
        # to those of the given mechanism.
        line_output = replacement_species(line,
                                          species_matching)
        
        f2.write(line_output)
        f2.write("\n")

# End of file        
with open(reactions_out, 'a') as f:
    f.write("END")


""" Write a combined list of the species as the input format of ssh-aerosol
"""
with open("species-list-user.dat","w") as f:
    f.write("# This file is automatically generated by combine_files.py.")
    for k, v in zip(species_name, species_weight):
        k += " "
        f.write("\n")
        f.write(k.ljust(8))
        f.write(v)

""" Write header
"""
output_filename = "species-list-aer-user.dat"
header = "# Name Type group ID	MW       Precursor   coll_fac        " + \
    "mole_diam       surf_tens accomod  mass_dens non-volatile  " + \
    "partitioning  smiles psat_torr  dHvap	Henry	Tref\n"
with open(output_filename,"w") as f:
    f.write(header)
       

""" Write inorganic aerosol species list. This should be called after 
the organic aerosol species.
"""
with open(output_filename,"a") as f, \
     open(species_inorganic_aerosol, "r") as f2:
    f2.readline()
    for line in f2:
        line = line.strip()
        if not line or line[0] in "#%!-":
            continue # comment
        
        f.write(line)
        f.write("\n")

        

""" Write organic aerosol species list.
"""
with open(output_filename,"a") as f, \
     open(species_aerosol, "r") as f2:
    f2.readline()
    for line in f2:
        line = line.strip()
        if not line or line[0] in "#%!-":
            continue # comment
        
        f.write(line)
        f.write("\n")

""" Write water aerosol species list. This should be called last.
"""
with open(output_filename,"a") as f, \
     open(species_water_aerosol, "r") as f2:
    f2.readline()
    for line in f2:
        line = line.strip()
        if not line or line[0] in "#%!-":
            continue # comment
        
        f.write(line)
        f.write("\n")        
        

        
