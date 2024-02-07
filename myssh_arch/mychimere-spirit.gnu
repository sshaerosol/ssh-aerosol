#!/bin/bash

#---------------------------------------------------------------------------------
#	Architecture file for compiling and running CHIMERE	
#	Specify path to libraries, compilers and utilities 
#---------------------------------------------------------------------------------

## Modules
module purge
module load gcc/9.4.0
module load openmpi/4.0.7
module load hdf5/1.10.7-mpi
module load netcdf-c/4.7.4-mpi
module load netcdf-fortran/4.5.3-mpi
module load eccodes/2.21.0
module load jasper/2.0.32
module load blitz/1.0.2

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
export  my_mpibin=${OPENMPI_ROOT}/bin    			# Path to MPI binary directory
export  my_mpirun=${my_mpibin}/mpirun    		# Path to mpirun to execute parallel job in MPI
export  my_mpif90=$MPIF90    		# Wrapper to my_compilerF90 to link with MPI library
export  my_mpicc=$MPICC     		# Wrapper to my_compilerC to link with MPI library
export  my_mpilib=${OPENMPI_ROOT}/lib    			# Path to MPI libraries directory
export  my_mpiinc=${OPENMPI_ROOT}/include    		# Path to MPI include files directory


#---------------------------------------------------------------------------------
# 	HDF5  - parallel version	
#---------------------------------------------------------------------------------
dirhdf=${HDF5_ROOT}
export my_hdflib=${dirhdf}/lib		# Path to HDF5 parallel library directory
export my_hdfinc=${dirhdf}/include	# Path to HDF5 parallel include files directory


#---------------------------------------------------------------------------------
# 	NETCDF-C  - link with HDF5 parallel 
#---------------------------------------------------------------------------------
dirnetcdfC=${NETCDF_C_ROOT}
export my_netcdfCbin=${dirnetcdfC}/bin 		# Path to NETCDF-C (linked with HDF5 parallel) binaries directory 
export my_netcdfClib=${dirnetcdfC}/lib		# Path to NETCDF-C (linked with HDF5 parallel) library directory
export my_netcdfCinc=${dirnetcdfC}/include		# Path to NETCDF-C (linked with HDF5 parallel) library directory


#---------------------------------------------------------------------------------
# 	NETCDF-Fortran  - link with HDF5 parallel and NETCDF-C
#---------------------------------------------------------------------------------
dirnetcdfF90=${NETCDF_FORTRAN_ROOT}
export my_netcdfF90bin=${dirnetcdfF90}/bin      # PATH to NETCDF-Fortran (linked with HDF5 parallel and NETCDF-C) binaries  directory
export my_netcdfF90lib=${dirnetcdfF90}/lib		# Path to NETCDF-Fortran (linked with HDF5 parallel and NETCDF-C) library  directory
export my_netcdfF90inc=${dirnetcdfF90}/include	# Path to NETCDF-Fortran (linked with HDF5 parallel and NETCDF-C) include files  directory


#---------------------------------------------------------------------------------
# 	GRIB  - link with jasper 
#---------------------------------------------------------------------------------
dirgrib=${ECCODES_ROOT}
export my_griblib=${dirgrib}/lib     	# Path to GRIB library directory
export my_gribinc=${dirgrib}/include 	# Path to GRIB include files directory
dirjasper=${JASPER_ROOT}
export my_jasperlib=${dirjasper}/lib 	 # Path to JASPER library directory
export my_jasperinc=${dirjasper}/include # Path to JASPER include files directory


#---------------------------------------------------------------------------------
# 	BLITZ
#---------------------------------------------------------------------------------
export my_blitzinc=${BLITZ_ROOT}/include


#---------------------------------------------------------------------------------
# 	Utilities	
#---------------------------------------------------------------------------------
export my_make=make 	                                        	# Path to make 
export my_awk=awk		                                        	# Path to awk
export my_python3=/usr/bin/python3
export my_ncdump=${dirnetcdfC}/bin/ncdump		# Path to ncdump


#---------------------------------------------------------------------------------
# 	Makefile header needed to compile CHIMERE and WRF 
#	     - with this architecture configuration - 	
#---------------------------------------------------------------------------------
export my_hdr=Makefile.hdr.gfortran-64-ompi   		            # Makefile header to compile CHIMERE in makefiles.hdr directory
export configure_wrf_file_name=configure.wrf.gfortran.spirit  	# Makefile header to compile WRF in config_wrf directory
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
