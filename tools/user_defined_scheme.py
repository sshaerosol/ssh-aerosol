#!/usr/bin/env python3
import sys, os
import configparser

config = configparser.ConfigParser()


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

    
def user_defined_scheme(options_list):

    config.read('user_defined_scheme.cfg')


    species_name = []
    species_weight = []

    # Toluene
    section = "Toluene"
    opt_tol = options_list.get(section)
    if opt_tol:
        species_name, species_weight, toluene_reactions = \
            add_scheme(section,
                       species_name,
                       species_weight)
        
    # Xylene
    section = "Xylene"
    opt_xyl = options_list.get(section)
    if opt_xyl:
        species_name, species_weight, xylene_reactions = \
            add_scheme(section,
                       species_name,
                       species_weight)
    
    # Alpha-pinene
    section = "Alpha-pinene"
    opt_api = options_list.get(section)
    if opt_api:
        species_name, species_weight, api_reactions = \
            add_scheme(section,
                       species_name,
                       species_weight)

    # Beta-pinene
    section = "Beta-pinene"
    opt_bpi = options_list.get(section)
    if opt_bpi:
        species_name, species_weight, bpi_reactions = \
            add_scheme(section,
                       species_name,
                       species_weight)

    # Limonene
    section = "Limonene"
    opt_lim = options_list.get(section)
    if opt_lim:
        species_name, species_weight, lim_reactions = \
            add_scheme(section,
                       species_name,
                       species_weight)

    # Isoprene
    section = "Isoprene"
    opt_iso = options_list.get(section)
    if opt_iso:
        species_name, species_weight, iso_reactions = \
            add_scheme(section,
                       species_name,
                       species_weight)

    # Humulene
    section = "Humulene"
    opt_hum = options_list.get(section)
    if opt_hum:
        species_name, species_weight, hum_reactions = \
            add_scheme(section,
                       species_name,
                       species_weight)

    # POA
    section = "POA"
    opt_poa = options_list.get(section)
    if opt_poa:
        species_name, species_weight, poa_reactions = \
            add_scheme(section,
                       species_name,
                       species_weight)

    # Benzene
    section = "Benzene"
    opt_ben = options_list.get(section)
    if opt_ben:
        species_name, species_weight, ben_reactions = \
            add_scheme(section,
                       species_name,
                       species_weight)

    # Cresol
    section = "Cresol" 
    opt_cre = options_list.get(section)
    if opt_cre:
        species_name, species_weight, cre_reactions = \
            add_scheme(section,
                       species_name,
                       species_weight)

    # Syringol
    section = "Syringol"
    opt_syr = options_list.get(section)
    if opt_syr:
        species_name, species_weight, syr_reactions = \
            add_scheme(section,
                       species_name,
                       species_weight)

    # Guaiacol
    section = "Guaiacol"
    opt_gua = options_list.get(section)
    if opt_gua:
        species_name, species_weight, gua_reactions = \
            add_scheme(section,
                       species_name,
                       species_weight)

    # Naphthalene
    section = "Naphthalene"
    opt_nap = options_list.get(section)
    if opt_nap:
        species_name, species_weight, nap_reactions = \
            add_scheme(section,
                       species_name,
                       species_weight)
        
        # section = section.lower()
        # nap_species = config[section]['species']
        # nap_reactions = config[section]['reactions']

        # species_name_nap, species_wt_nap, ns_nap = \
        #     read_species(nap_species)                

    # Methyl-naphthalene
    section = "Methyl-naphthalene"
    opt_mna = options_list.get(section)
    if opt_mna:
        species_name, species_weight, mna_reactions = \
            add_scheme(section,
                       species_name,
                       species_weight)
            
    # Remove old files
    species_aerosol = "../user_defined/user_defined.species"
    reactions_aerosol = "../user_defined/user_defined.reactions"
    if os.path.isfile(species_aerosol):
        os.remove(species_aerosol)
    if os.path.isfile(reactions_aerosol):
        os.remove(reactions_aerosol)

    # Write species file
    header = "File for chemical species user-defined\n" + \
              "# gaseous species # aqueous species\n" + \
              str(len(species_name)) + " 0\n" + \
              "---Gas-phase---"

    with open(species_aerosol, 'w') as f:
        f.write(header)
        for k, v in zip(species_name, species_weight):
            k += " "
            f.write("\n")
            f.write(k.ljust(8))
            f.write(v)


    # Write reaction file
    header = "#==== User-definded reactions ====\n"
    with open(reactions_aerosol, 'w') as f:
        f.write(header)
        
    if opt_tol:        
        copy_reactions(toluene_reactions, reactions_aerosol, "a")
    if opt_xyl:                
        copy_reactions(xylene_reactions, reactions_aerosol, "a")
    if opt_api:                
        copy_reactions(api_reactions, reactions_aerosol, "a")    
    if opt_bpi:
        copy_reactions(bpi_reactions, reactions_aerosol, "a")    
    if opt_lim:
        copy_reactions(lim_reactions, reactions_aerosol, "a")
    if opt_iso:
        copy_reactions(iso_reactions, reactions_aerosol, "a")
    if opt_hum:
        copy_reactions(hum_reactions, reactions_aerosol, "a")
    if opt_poa:
        copy_reactions(poa_reactions, reactions_aerosol, "a")            
    if opt_ben:
        copy_reactions(ben_reactions, reactions_aerosol, "a")
    if opt_cre:
        copy_reactions(cre_reactions, reactions_aerosol, "a")
    if opt_syr:
        copy_reactions(syr_reactions, reactions_aerosol, "a")
    if opt_gua:
        copy_reactions(gua_reactions, reactions_aerosol, "a")
    if opt_nap:
        copy_reactions(nap_reactions, reactions_aerosol, "a")
    if opt_mna:
        copy_reactions(mna_reactions, reactions_aerosol, "a")    
    
    print(options_list)

    return reactions_aerosol, species_aerosol

