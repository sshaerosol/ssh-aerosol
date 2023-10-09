#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (C) 2023, CEREA
#     Author(s): Youngseob Kim
#
from __future__ import with_statement

import os
import string
import sys
from replacement_species import *

Ns = None
Ns_gas = None
Ns_h2o = None
species_name=[]
species_weight=[]

""" Input files

species_gas = sys.argv[1]
reactions_gas = sys.argv[2]

species_h2o = "../h2o/h2o.species"
reactions_h2o = "../h2o/h2o.reactions"

"""

species_gas = sys.argv[1]
reactions_gas = sys.argv[2]

species_h2o = "../h2o/h2o.species"
reactions_h2o = "../h2o/h2o.reactions"

species_out = "combined_species.dat"
reactions_out = "combined_reactions.dat"

with open(species_gas,'r') as f:
    # Jumps the title line.
    f.readline()
    for line in f:
        line = line.strip()
        if not line or line[0] in "#%!-":
            continue # A comment,
        line=line.split()
        if Ns_gas is None:
            Ns_gas = int(line[0])
            continue; # The first value is the number of species.
        species_name.append(line[0])
        species_weight.append(line[1])

# Checks Ns even if it is not needed in this script.
if (Ns_gas != len(species_name)):
    print("Error " + Ns_gas + " " + len(species_name))
    sys.exit(1)
        

with open(species_h2o,'r') as f:
    # Jumps the title line.
    f.readline()
    for line in f:
        line = line.strip()
        if not line or line[0] in "#%!-":
            continue # A comment,
        line=line.split()
        if Ns_h2o is None:
            Ns_h2o = int(line[0])
            continue; # The first value is the number of species.
        species_name.append(line[0])
        species_weight.append(line[1])        

# Checks Ns even if it is not needed in this script.
Ns = Ns_gas + Ns_h2o
if (Ns != len(species_name)):
    print("Error " + str(Ns) + " " + str(len(species_name)))
    sys.exit(1)

# ####################
# # TODO For debugging purpose, an output of the photolysis reactions is done,
# # but a mapping should be done with the real photolysis reactions names
# # to generate automatically the section [photolysis_reaction_index].
# def format_reaction(reaction_str):
#     lhs, rhs = reaction_str.split('->')
#     lhs = sorted([ s.strip() for s in lhs.split('+')])
#     rhs = sorted([ s.strip() for s in rhs.split('+')])
#     return ' + '.join(lhs) + ' -> ' + ' + '.join(rhs)

# photolysis_reactions = []
# current_reaction = ""
# with open(sys.argv[2],'r', encoding="UTF-8") as f:
#     # Jumps the title line.
#     f.readline()
#     for line in f:
#         line = line.strip()
#         if not line or line[0] in "#%!-":
#             continue # A comment.
#         if '->' in line:
#             current_reaction = line
#             continue # A reaction.
#         line=line.split()
#         if line[0:2] != ["KINETIC", "PHOTOLYSIS"]:
#             continue # current_reaction is not a photolysis reaction.

#         photolysis_reactions.append(format_reaction(current_reaction))
# ####################


""" Write a combined list of the species
"""
header = "File for chemical species CB05 (name and molar mass)\n" + \
         "# gaseous species # aqueous species\n" + \
         str(Ns) + "  0\n" + \
         "---Gas-phase----"

with open(species_out, "w") as f:
    f.write(header)
    for k, v in zip(species_name, species_weight):
        k += " "
        f.write("\n")
        f.write(k.ljust(8))
        f.write(v)



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

with open(reactions_h2o,'r', encoding="UTF-8") as f1, open(reactions_out, 'a') as f2:

    species_matching = read_matching_file()
    
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
        


#with open("species.spack.dat","w") as f:
#    f.write("[species]\n")
#    f.write("\n".join(species_name))
#    f.write("\n\n[molecular_weight]\n")
#    for k, v in zip(species_name, species_weight):
#        k += " "
#        v += "\n"
#        f.write(k.ljust(8))
#        f.write(v)
#    f.write("""

# Here is the list of photolysis reactions in good order.
# This section is generated to help debugging
# the section '[photolysis_reaction_index]'
#""")
#    f.write("\n#".join(photolysis_reactions))
#    f.write("\n")
    
#write species-list.dat automaticaly
with open("species-list.dat","w") as f:
    f.write("# The order of species in this list should not be changed.")
    for k, v in zip(species_name, species_weight):
        k += " "
        f.write("\n")
        f.write(k.ljust(8))
        f.write(v)
