#!/usr/bin/env python3
import sys, os
import configparser

config = configparser.ConfigParser()



def append_aerosol_species(new_species, file_out, opt_write):

    with open(file_out, opt_write) as f2:
        for line in new_species:
            f2.write('\t'.join(line))
            f2.write("\n")

def copy_reactions(file_in, file_out, opt_write):

    with open(file_in, 'r', encoding = "UTF-8") as f1, \
        open(file_out, opt_write) as f2:
        f1.readline()
        for line in f1:
            line = line.strip()
            if not line or line[0] in "#%!-":
                continue # comment
            if (line != "END"):
                f2.write(line)
                f2.write("\n")
                

def read_species(file_name):

    species_name = []
    species_weight = []
    ns = None

    with open(file_name, 'r') as f:

        # Jump the title line
        f.readline()
        for line in f:
            line = line.strip()
            if not line or line[0] in "#%!-":
                continue # comment
            line = line.split()
            if ns is None:
                ns = int(line[0])
                continue # First value is the number of species.
            species_name.append(line[0])
            species_weight.append(line[1])

    return species_name, species_weight, ns


def add_scheme(section,
               species_name, species_weight):

    section = section.lower()
    species = config[section]['species']
    reactions = config[section]['reactions']

    species_name_scheme, species_wt, ns = \
        read_species(species)

    for spe, wt in zip(species_name_scheme, \
                       species_wt):
        if spe in species_name:
            print(spe + " exists")
        else:
            species_name.append(spe)
            species_weight.append(wt)

    return species_name, species_weight, reactions


def read_aerosol_species(file_name):

    aerosol_species_name = []
    aerosol_properties = []
    ns_aer = 0
    
    with open(file_name, 'r') as f:

        # Jump the title line
        f.readline()
        for line in f:
            ns_aer = ns_aer + 1
            line = line.strip()
            if not line or line[0] in "#%!-":
                continue # comment
            line = line.split()
            aerosol_species_name.append(line[0])
            aerosol_properties.append(line)
            
    return aerosol_species_name, aerosol_properties


def add_aerosol_scheme(section,
                       aerosol_species_name):

    section = section.lower()

    aerosol_properties = []
    
    # Read aerosol species
    aerosol_species = config[section]['aerosol_species']

    aerosol_species_name_scheme, aerosol_properties_scheme = \
        read_aerosol_species(aerosol_species)

    for spe, properties in zip(aerosol_species_name_scheme, \
                               aerosol_properties_scheme):
        if spe in aerosol_species_name:
            print(spe + " exists")
        else:
            aerosol_species_name.append(spe)
            aerosol_properties.append(properties)
    

    return aerosol_species_name, aerosol_properties


""" Read the list of RO2 reactions from the given schemes,
and make a new file listing all the RO2 reactions.
"""
def add_ro2_species(section,
                    ro2_species_name,
                    ro2_group_id):

    section = section.lower()

    try:
        # Read RO2 reactions
        ro2_reactions = config[section]['ro2_reactions']

        ro2_species_name_scheme, ro2_group_id_scheme = \
            read_ro2_species(ro2_reactions)

        for spe, gid in zip(ro2_species_name_scheme, \
                           ro2_group_id_scheme):
            if spe in ro2_species_name:
                print(spe + " exists in RO2 reactions.")
            else:
                ro2_species_name.append(spe)
                ro2_group_id.append(gid)
                
        
    except Exception:
        print("RO2 reaction not used for ", section)
        pass

    
def read_ro2_species(file_name):

    ro2_species_name = []
    ro2_group_id = []
    
    with open(file_name, 'r') as f:
        
        # Jump the title line
        f.readline()
        for line in f:
            line = line.strip()
            if not line or line[0] in "#%!-":
                continue # comment
            line = line.split()
            ro2_species_name.append(line[0])
            ro2_group_id.append(line[1])
    return ro2_species_name, ro2_group_id


def user_defined_scheme(options_list, configfile):
    
    config.read(configfile)

    species_name = []
    species_weight = []
    aerosol_species_name = []

    print(options_list)

    # Remove old files
    species_precursor = "../user_defined/user_defined.species"
    reactions_precursor = "../user_defined/user_defined.reactions"
    species_aerosol = "../user_defined/species-list-aer.dat"
    ro2_reactions = "../user_defined/ro2_reactions.dat"
    if os.path.isfile(species_precursor):
        os.remove(species_precursor)
    if os.path.isfile(reactions_precursor):
        os.remove(reactions_precursor)
    if os.path.isfile(species_aerosol):
        os.remove(species_aerosol)
    if os.path.isfile(ro2_reactions):
        os.remove(ro2_reactions)                

    # Write header for reaction file
    header = "#==== User-definded reactions ====\n"
    with open(reactions_precursor, 'w') as f:
        f.write(header)

    # Write header for aerosol species file
    header = "# Name Type group ID	MW       Precursor   coll_fac        mole_diam       surf_tens accomod  mass_dens non-volatile  partitioning  smiles psat_torr  dHvap	Henry	Tref\n"
    with open(species_aerosol, 'w') as f:
        f.write(header)        

    # RO2 variables
    ro2_species_name = []
    ro2_group_id = []

    # Loop of the schemes.
    for opt in options_list:

        section = opt
        if (options_list.get(opt)):
            print ("===== " + opt + " ====")

            
            
            species_name, species_weight, reactions = \
                add_scheme(section,
                           species_name,
                           species_weight)

            aerosol_species_name, properties = \
                add_aerosol_scheme(section,
                                   aerosol_species_name)

            add_ro2_species(section,
                            ro2_species_name,
                            ro2_group_id)
        
            copy_reactions(reactions, reactions_precursor, "a")
            append_aerosol_species(properties, species_aerosol, "a")
            

    # Write species file
    species_header = "File for chemical species user-defined\n" + \
        "# gaseous species # aqueous species\n" + \
        str(len(species_name)) + " 0\n" + \
        "---Gas-phase---"

    with open(species_precursor, 'w') as f:
        f.write(species_header)
        for k, v in zip(species_name, species_weight):
            k += " "
            f.write("\n")
            f.write(k.ljust(8))
            f.write(v)


    # Write RO2 reactions file
    header = "#==== User-definded RO2 reactions ===="
    with open(ro2_reactions, 'w') as f:
        f.write(header)
        for k, v in zip(ro2_species_name, ro2_group_id):
            k += " "
            f.write("\n")
            f.write(k.ljust(8))
            f.write(v)


    return reactions_precursor, species_precursor, species_aerosol
