#!/bin/bash

module purge
#module unload blitz
module load flavor/hdf5/parallel
module load hdf5
module load netcdf-fortran
module load grib
module load jasper
module load blitz/1.0.2

source ${FORTRAN_INTEL_ROOT}/bin/compilervars.sh intel64
source ${C_INTEL_ROOT}/bin/compilervars.sh intel64

#---------------------------------------------------------------------------------
#	Architecture file for compiling and running CHIMERE	
#	Specify path to libraries, compilers and utilities 
#---------------------------------------------------------------------------------


#---------------------------------------------------------------------------------
# 	Compilers
#---------------------------------------------------------------------------------
export my_compilerF90=${FORTRAN_INTEL_ROOT}/bin/intel64/ifort # Path to Fortran 90 compiler
export my_compilerC=${C_INTEL_EXEDIR}/icc                     # Path to C compiler
export my_compilerCpp=${C_INTEL_EXEDIR}/icpc                  # Path to C++ compiler


#---------------------------------------------------------------------------------
# 	MPI - parallel execution of chimere
#---------------------------------------------------------------------------------
export  my_mpiframe=ccrt              # implementaion of MPI norm [ ompi / ccrt ] TO REMOVE
export  my_mpibin=${MPI_EXEDIR}       # Path to MPI binary directory
export  my_mpirun=ccc_mprun           # Path to mpirun to execute parallel job in MPI
export  my_mpif90=${my_mpibin}/mpif90 # Wrapper to my_compilerF90 to link with MPI library
export  my_mpicc=${my_mpibin}/mpicc   # Wrapper to my_compilerC to link with MPI library
export  my_mpilib=${MPI_ROOT}/lib     # Path to MPI libraries directory
export  my_mpiinc=${MPI_ROOT}/include # Path to MPI include files directory


#---------------------------------------------------------------------------------
# 	HDF5  - parallel version	
#---------------------------------------------------------------------------------
export my_hdflib=${HDF5_LIBDIR}		# Path to HDF5 parallel library directory
export my_hdfinc=${HDF5_INCDIR}		# Path to HDF5 parallel include files directory


#---------------------------------------------------------------------------------
# 	NETCDF-C  - link with HDF5 parallel 
#---------------------------------------------------------------------------------
export my_netcdfCbin=${NETCDFC_EXEDIR} 		# Path to NETCDF-C (linked with HDF5 parallel) binaries directory 
export my_netcdfClib=${NETCDFC_LIBDIR}		# Path to NETCDF-C (linked with HDF5 parallel) library directory


#-------------------------------------------------------------------------------
# 	NETCDF-Fortran  - link with HDF5 parallel and NETCDF-C
#---------------------------------------------------------------------------------
export my_netcdfF90bin=${NETCDFFORTRAN_EXEDIR}  # PATH to NETCDF-Fortran (linked with HDF5 parallel and NETCDF-C) binaries  directory
export my_netcdfF90lib=${NETCDFFORTRAN_LIBDIR}	# Path to NETCDF-Fortran (linked with HDF5 parallel and NETCDF-C) library  directory
export my_netcdfF90inc=${NETCDFFORTRAN_INCDIR}	# Path to NETCDF-Fortran (linked with HDF5 parallel and NETCDF-C) include files  directory


#---------------------------------------------------------------------------------
# 	GRIB  - link with jasper 
#---------------------------------------------------------------------------------
export my_griblib=${GRIB_LIBDIR}     			# Path to GRIB library directory
export my_gribinc=${GRIB_INCDIR} 			# Path to GRIB include files directory
export my_jasperlib=${JASPER_LIBDIR} 			# Path to JASPER library directory
export my_jasperinc=${JASPER_INCDIR} 			# Path to JASPER include files directory


#---------------------------------------------------------------------------------
# 	BLITZ
#---------------------------------------------------------------------------------
export my_blitzinc=${BLITZ_INCDIR}		 # Path to BLITZ include files 


#---------------------------------------------------------------------------------
# 	Utilities	
#---------------------------------------------------------------------------------
export my_make=make                      # Path to make
export my_awk=awk                        # Path to awk
export my_ncdump=${NETCDF_EXEDIR}/ncdump # Path to ncdump


#---------------------------------------------------------------------------------
# 	Makefile header needed to compile CHIMERE and WRF 
#	     - with this architecture configuration - 	
#---------------------------------------------------------------------------------
export my_hdr=Makefile.hdr.ifort-64-ompi   			# Makefile header to compile CHIMERE in makefiles.hdr directory
export configure_wrf_file_name="configure.wrf.ifort_irene"  	# Makefile header to compile WRF in config_wrf directory
export configure_wps_file_name=configure_ifort.wps  		# Makefile header to compile WPS in config_wps directory


#---------------------------------------------------------------------------------
#	Export of Shared Library to be available at run time 	
#---------------------------------------------------------------------------------
export LD_LIBRARY_PATH=${my_hdflib}:${my_netcdfF90lib}:${my_netcdfClib}:${my_griblib}:${my_mpilib}:${my_mpilib}/openmpi:$LD_LIBRARY_PATH

#export my_ssh=${CCCWORKDIR}/SSH/ssh-aerosol-master-b417da263350a38a98a57e4847af23a0742f7823/src/
export my_ssh=${CCCWORKDIR}/SSH/ssh-aerosol/src/

#---------------------------------------------------------------------------------
# 	Monitoring	
#---------------------------------------------------------------------------------

export my_prof=vtune
export my_prof=none

if [ $my_prof == "vtune" ] ; then
   module load vtune
   export my_ompiparams="$VTUNE_EXEDIR/amplxe-cl -collect hotspots -r ${PWD}/vtune_results"
   export my_mode=PROD
fi
if [ $my_prof == "intel" ] ; then
   module load advisor
   export my_ompiparams="$ADVISOR_EXEDIR/advixe-cl -collect survey --search-dir src:r=${chimere_root}/src -project-dir ${chimere_root}/advise_results --"
   export my_mode=PROD
fi
if [ $my_prof == "ipm" ] ; then
   module load ipm
   export LD_PRELOAD=$IPM_LIBDIR/libipm.so
   export my_mode=PROD
   #imp_parse -html $tmpdir/*.xml
fi

