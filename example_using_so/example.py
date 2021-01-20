#! /usr/bin/env python3

#
# Import modules
#
import ctypes as ct
import numpy as np
import matplotlib.pyplot as plt

#
# Load ssh-aerosol shared library generated with the command ./compile --sharedlib=yes
#
try :
    print("Loading shared library ssh-aerosol")
    libssh = ct.CDLL("../src/libssh-aerosol.so")
except OSError:
    print("Shared library ssh-aerosol could not be found")
    exit()

#
# Python assumes every function returns a int
#
libssh.api_sshaerosol_set_standalone_.restype=None
libssh.api_sshaerosol_set_logger_.restype=None
libssh.api_sshaerosol_initialize_.restype=None
libssh.api_sshaerosol_init_distributions_.restype=None
libssh.api_sshaerosol_get_initial_t_.restype=ct.c_double
libssh.api_sshaerosol_get_dt_.restype=ct.c_double
libssh.api_sshaerosol_get_current_t_.restype=ct.c_double
libssh.api_sshaerosol_get_gas_.restype=None
libssh.api_sshaerosol_get_aero_.restype=None
libssh.api_sshaerosol_get_aero_num_.restype=None
libssh.api_sshaerosol_get_aero_name_.restype=None
libssh.api_sshaerosol_initoutput_.restype=None
libssh.api_sshaerosol_report_.restype=None
libssh.api_sshaerosol_output_.restype=None
libssh.api_sshaerosol_initphoto_.restype=None
libssh.api_sshaerosol_set_current_t_.restype=None
libssh.api_sshaerosol_updatephoto_.restype=None
libssh.api_sshaerosol_emission_.restype=None
libssh.api_sshaerosol_gaschemistry_.restype=None
libssh.api_sshaerosol_aerodyn_.restype=None
libssh.api_sshaerosol_output_.restype=None
libssh.api_sshaerosol_finalize_.restype=None

#
# Set the standalone flag to false
# Set the logger flag to true
#
libssh.api_sshaerosol_set_standalone_(ct.byref(ct.c_bool(False)))
libssh.api_sshaerosol_set_logger_(ct.byref(ct.c_bool(True)))

#
# Initialize with a namelist
#
print("Initialize ssh-aerosol.")
libssh.api_sshaerosol_initialize_(ct.create_string_buffer(b"namelist.ssh"))
libssh.api_sshaerosol_init_distributions_()
print("Initialization completed.")

#
# Get number of gas species, aerosols and aerosols layers
#
ngas = libssh.api_sshaerosol_get_ngas_()
n_aero_layer = libssh.api_sshaerosol_get_n_aerosol_layers_()
n_size = libssh.api_sshaerosol_get_nsize_()
print("Number of gas species, aerosols and aerosols layers : " + str(ngas) + " " + str(n_size) + " " + str(n_aero_layer))

#
# Read the initial time and the time step
#
libssh.api_sshaerosol_get_initial_t_.restype=ct.c_double
libssh.api_sshaerosol_get_dt_.restype=ct.c_double
libssh.api_sshaerosol_get_current_t_.restype=ct.c_double
#
t_init = libssh.api_sshaerosol_get_initial_t_()
dt = libssh.api_sshaerosol_get_dt_()

#
# Read the gas concentration
#
cgas = (ct.c_double * ngas)()
libssh.api_sshaerosol_get_gas_(cgas)
#print("cgas : " + str(cgas[:]))

#
# Read the aerosol concentration
#
caero = (ct.c_double * (n_aero_layer * n_size))()
libssh.api_sshaerosol_get_aero_(caero)
#print("caero : " + str(caero[:]))

#
# Read the aerosol numbers
#
caeronum = (ct.c_double * n_size)()
libssh.api_sshaerosol_get_aero_num_(caeronum)
#print("caeronum : " + str(caeronum[:]))

#
# Read the name of aerosol number 3
#
cname = ct.create_string_buffer(81)
libssh.api_sshaerosol_get_aero_name_(ct.byref(ct.c_int(3)),cname)
#print("Name of aerosol number 3 : " + str(cname.value.decode()))


#
# Call ssh-aerosol built-in reporting
#
libssh.api_sshaerosol_initoutput_()
libssh.api_sshaerosol_report_()
libssh.api_sshaerosol_output_()

#
# Initialize photo-chemistry
#
libssh.api_sshaerosol_initphoto_()

#
# Perform 60 time steps with chemistry + aerosol dynamic
#
time = t_init
for it in range(60):
   # Update time in ssh-aerosol
   if it != 0:
      time += dt
      libssh.api_sshaerosol_set_current_t_(ct.byref(ct.c_double(time)))
   # Update photolysis in ssh-aerosol
   libssh.api_sshaerosol_updatephoto_()
   # Emission
   libssh.api_sshaerosol_emission_()
   # Gas chemistry
   libssh.api_sshaerosol_gaschemistry_()
   # Aerosol dynamic
   libssh.api_sshaerosol_aerodyn_()
   # ssh-aerosol built-in reporting
   libssh.api_sshaerosol_output_()

#
# Finalize
#
libssh.api_sshaerosol_finalize_()

