#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (C) 2023, CEREA
#     Author(s): Youngseob Kim
#
from __future__ import with_statement

import os
import string
import sys

def read_matching_file(matching_file):

    old = []
    new = []
    with open(matching_file, 'r', encoding="UTF-8") as f_replace:
        for line in f_replace:  
            line = line.strip()
            if not line or line[0] in "#%!-":
                continue # A comment
            old.append(line.split()[0])
            new.append(line.split()[1])
        
    species_matching = dict(zip(old, new))
    print(species_matching)
    
    return species_matching
    
def replacement_species(input_line, species_matching):


    list_key = list(species_matching.keys())

    output = input_line
    for key in list_key:
        value = species_matching.get(key)
        output = output.replace(key, value) 

    return output

# line = "alpha_pinene + hydroxyl = alpha_pinene + hydroxyl"
# print(line)
# line_out = replacement_species(line)
# print(line_out)
                                    
