#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (C) 2015, ENPC
#     Author(s): Sylvain Doré
#
# This file is part of the air quality modeling system Polyphemus.
#
# Polyphemus is developed in the INRIA project-team CLIME and in
# the ENPC - EDF R&D joint laboratory CEREA.
#
# Polyphemus is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Polyphemus is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.
#
# For more information, visit the Polyphemus web site:
#      http://cerea.enpc.fr/polyphemus/
from __future__ import with_statement

import os
import string
import sys


Ns=None
species_name=[]
species_weight=[]

with open(sys.argv[1],'r') as f:
    # Jumps the title line.
    f.readline()
    for line in f:
        line = line.strip()
        if not line or line[0] in "#%!-":
            continue # A comment,
        line=line.split()
        if Ns is None:
            Ns=int(line[0])
            continue; # The first value is the number of species.
        species_name.append(line[0])
        species_weight.append(line[1])

# Checks Ns even if it is not needed in this script.
assert Ns == len(species_name)


####################
# TODO For debugging purpose, an output of the photolysis reactions is done,
# but a mapping should be done with the real photolysis reactions names
# to generate automatically the section [photolysis_reaction_index].
def format_reaction(reaction_str):
    lhs, rhs = reaction_str.split('->')
    lhs = sorted([ s.strip() for s in lhs.split('+')])
    rhs = sorted([ s.strip() for s in rhs.split('+')])
    return ' + '.join(lhs) + ' -> ' + ' + '.join(rhs)

photolysis_reactions = []
current_reaction = ""
with open(sys.argv[2],'r', encoding="UTF-8") as f:
    # Jumps the title line.
    f.readline()
    for line in f:
        line = line.strip()
        if not line or line[0] in "#%!-":
            continue # A comment.
        if '->' in line:
            current_reaction = line
            continue # A reaction.
        line=line.split()
        if line[0:2] != ["KINETIC", "PHOTOLYSIS"]:
            continue # current_reaction is not a photolysis reaction.

        photolysis_reactions.append(format_reaction(current_reaction))
####################


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
