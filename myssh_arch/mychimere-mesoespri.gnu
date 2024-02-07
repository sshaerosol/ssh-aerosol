#!/bin/bash

#---------------------------------------------------------------------------------
#	Architecture file for compiling and running CHIMERE	
#	Specify path to libraries, compilers and utilities 
#---------------------------------------------------------------------------------

## Modules
module purge
module load gnu/9.3.0
module load openmpi/4.0.5

#---------------------------------------------------------------------------------
# 	Compilers
#---------------------------------------------------------------------------------
export my_compilerF90=gfortran		# Path to Fortran 90 compiler
export my_compilerC=gcc			# Path to C compiler
export my_compilerCpp=g++			# Path to C++ compiler


#---------------------------------------------------------------------------------
# 	MPI - parallel execution of chimere
#---------------------------------------------------------------------------------
export  my_mpiframe=openmpi  		                            	# implementaion of MPI norm [ ompi / ccrt ] TO REMOVE
export  my_mpibin=/net/nfs/tools/meso-sl6/openmpi/4.0.5-gcc-9.3.0/bin    			# Path to MPI binary directory
export  my_mpirun=${my_mpibin}/mpirun    		# Path to mpirun to execute parallel job in MPI
export  my_mpif90=${my_mpibin}/mpif90    		# Wrapper to my_compilerF90 to link with MPI library
export  my_mpicc=${my_mpibin}/mpicc     		# Wrapper to my_compilerC to link with MPI library
export  my_mpilib=${my_mpibin}/../lib    			# Path to MPI libraries directory
export  my_mpiinc=${my_mpibin}/../include    		# Path to MPI include files directory


#---------------------------------------------------------------------------------
# 	HDF5  - parallel version	
#---------------------------------------------------------------------------------
dirhdf=/net/nfs/tools/PrgEnv/linux-scientific6-x86_64/gcc-9.3.0/hdf5-1.10.7-nnvab7d5rtajgbzj23cdrxr25ejyqp6o
export my_hdflib=${dirhdf}/lib		# Path to HDF5 parallel library directory
export my_hdfinc=${dirhdf}/include	# Path to HDF5 parallel include files directory


#---------------------------------------------------------------------------------
# 	NETCDF-C  - link with HDF5 parallel 
#---------------------------------------------------------------------------------
dirnetcdfC=/net/nfs/tools/PrgEnv/linux-scientific6-x86_64/gcc-9.3.0/netcdf-c-4.7.4-vnbhiojcfwagpd3udmndopc6o4vw22zj
export my_netcdfCbin=${dirnetcdfC}/bin 		# Path to NETCDF-C (linked with HDF5 parallel) binaries directory 
export my_netcdfClib=${dirnetcdfC}/lib		# Path to NETCDF-C (linked with HDF5 parallel) library directory
export my_netcdfCinc=${dirnetcdfC}/include		# Path to NETCDF-C (linked with HDF5 parallel) library directory


#---------------------------------------------------------------------------------
# 	NETCDF-Fortran  - link with HDF5 parallel and NETCDF-C
#---------------------------------------------------------------------------------
dirnetcdfF90=/net/nfs/tools/PrgEnv/linux-scientific6-x86_64/gcc-9.3.0/netcdf-fortran-4.5.3-4he3oodq53qw4yhwz6sngsvzt7dekcnv
export my_netcdfF90bin=${dirnetcdfF90}/bin      # PATH to NETCDF-Fortran (linked with HDF5 parallel and NETCDF-C) binaries  directory
export my_netcdfF90lib=${dirnetcdfF90}/lib		# Path to NETCDF-Fortran (linked with HDF5 parallel and NETCDF-C) library  directory
export my_netcdfF90inc=${dirnetcdfF90}/include	# Path to NETCDF-Fortran (linked with HDF5 parallel and NETCDF-C) include files  directory


#---------------------------------------------------------------------------------
# 	GRIB  - link with jasper 
#---------------------------------------------------------------------------------
dirgrib=/net/nfs/tools/PrgEnv/linux-scientific6-x86_64/gcc-9.3.0/eccodes-2.20.0-wn3jopdqng6a2trbjp5imdauzefjuniv/
export my_griblib=${dirgrib}/lib64     	# Path to GRIB library directory
export my_gribinc=${dirgrib}/include 	# Path to GRIB include files directory
dirjasper=/net/nfs/tools/PrgEnv/linux-scientific6-x86_64/gcc-9.3.0/jasper-2.0.16-gb2xfrgrl4w6fu2kgkjmo4ukrhr4yupu
export my_jasperlib=${dirjasper}/lib64 		                	                # Path to JASPER library directory
export my_jasperinc=${dirjasper}/include/jasper			                        # Path to JASPER include files directory


#---------------------------------------------------------------------------------
# 	BLITZ
#---------------------------------------------------------------------------------
#export my_blitzinc=/home/rpennel/local/blitz/include/blitz-0.10		 # Path to BLITZ include files 
export my_blitzinc=/net/nfs/tools/PrgEnv/linux-scientific6-x86_64/gcc-9.3.0/blitz-1.0.1-rpk4svg7s4d7jwbshhheevkyyn2djd2a/include


#---------------------------------------------------------------------------------
# 	Utilities	
#---------------------------------------------------------------------------------
export my_make=make 	                                        	# Path to make 
export my_awk=awk		                                        	# Path to awk
export my_python3=/opt/anaconda3/bin/python3
export my_ncdump=${dirnetcdfC}/bin/ncdump		# Path to ncdump


#---------------------------------------------------------------------------------
# 	Makefile header needed to compile CHIMERE and WRF 
#	     - with this architecture configuration - 	
#---------------------------------------------------------------------------------
export my_hdr=Makefile.hdr.gfortran-64-ompi   		            	# Makefile header to compile CHIMERE in makefiles.hdr directory
export configure_wrf_file_name=configure.wrf.gfortran             	# Makefile header to compile WRF in config_wrf directory
export configure_wps_file_name=configure_gfortran.wps          	# Makefile header to compile WPS in config_wps directory
export configure_xios_file_name=arch-gfortran.fcm       # Makefile header (*.fcm) to compile XIOS in config_xios directory


#---------------------------------------------------------------------------------
#	Export of Shared Library to be available at run time 	
#---------------------------------------------------------------------------------
export LD_LIBRARY_PATH=${my_hdflib}:${my_netcdfF90lib}:${my_netcdfClib}:${my_griblib}:${my_mpilib}:${my_mpilib}/openmpi:${my_jasperlib}:$LD_LIBRARY_PATH

#---------------------------------------------------------------------------------
# 	SSH-aerosol	
#---------------------------------------------------------------------------------
export my_ssh=${chimere_root}/ssh-aerosol/src/

