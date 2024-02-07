module purge
module load slurm
module load openmpi/gnu/64/4.0.5
module load hdf5/parallel/gnu/1.12.0
module load netcdf/netcdf-c/gnu/4.7.4
module load netcdf/netcdf-f/gnu/4.5.3
module load blitz/gnu/blitz
module load eccodes/gnu/2.22.1
module load git/2.18.4
module load cdo/gnu/1.9.10rc1
module load nco/gnu/4.9.5
#---------------------------------------------------------------------------------
#	Architecture file for compiling and running CHIMERE	
#	Specify path to libraries, compilers and utilities 
#---------------------------------------------------------------------------------


#---------------------------------------------------------------------------------
# 	Compilers
#---------------------------------------------------------------------------------
export my_compilerF90=gfortran # Path to Fortran 90 compiler
export my_compilerC=gcc                     # Path to C compiler
export my_compilerCpp=g++                  # Path to C++ compiler


#---------------------------------------------------------------------------------
# 	MPI - parallel execution of chimere
#---------------------------------------------------------------------------------
export  my_mpiframe=openmpi              # implementaion of MPI norm [ ompi / ccrt ] TO REMOVE
export  my_mpibin=${MPI_EXEDIR}       # Path to MPI binary directory
export  my_mpirun=${MPI_EXEDIR}/mpirun           # Path to mpirun to execute parallel job in MPI
export  my_ompiparams="--use-hwthread-cpus"
export  my_mpif90=${my_mpibin}/mpif90 # Wrapper to my_compilerF90 to link with MPI library
export  my_mpicc=${my_mpibin}/mpicc   # Wrapper to my_compilerC to link with MPI library
export  my_mpilib=${MPI_ROOT}/lib     # Path to MPI libraries directory
export  my_mpiinc=${MPI_ROOT}/include # Path to MPI include files directory

# Specific environment variable to avoid flood by warnings on irene AMD

export UCX_LOG_LEVEL=error

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
export my_netcdfCinc=${NETCDFC_INCDIR}		# Path to NETCDF-C (linked with HDF5 parallel) library directory


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
export my_blitzlib=${BLITZ_LIBDIR}		 # Path to BLITZ include files 


#---------------------------------------------------------------------------------
# 	Utilities	
#---------------------------------------------------------------------------------
export my_make=make                      # Path to make
export my_awk=awk                        # Path to awk
export my_python3=/usr/bin/python3
export my_ncdump=${NETCDFC_EXEDIR}/ncdump # Path to ncdump


#---------------------------------------------------------------------------------
# 	Makefile header needed to compile CHIMERE and WRF 
#	     - with this architecture configuration - 	
#---------------------------------------------------------------------------------
export my_hdr=Makefile.hdr.gfortran-64-ompi   			# Makefile header to compile CHIMERE in makefiles.hdr directory
export configure_wrf_file_name="configure.wrf.gfortran.chimatmos2"  	# Makefile header to compile WRF in config_wrf directory
export configure_wps_file_name=configure_gfortran.wps  		# Makefile header to compile WPS in config_wps directory
export configure_xios_file_name=arch-gfortran.fcm           # Makefile header (*.fcm) to compile XIOS in config_xios directory

#---------------------------------------------------------------------------------
#	Export of Shared Library to be available at run time 	
#---------------------------------------------------------------------------------
export LD_LIBRARY_PATH=${my_hdflib}:${my_netcdfF90lib}:${my_netcdfClib}:${my_griblib}:${my_mpilib}:${my_mpilib}/openmpi:${my_jasperlib}:$LD_LIBRARY_PATH


#---------------------------------------------------------------------------------
# 	Monitoring	
#---------------------------------------------------------------------------------

export my_prof=vtune
export my_prof=none

if [ $my_prof == "vtune" ] ; then
   module load vtune
   export my_ompiparams="$my_ompiparams $VTUNE_EXEDIR/amplxe-cl -collect hotspots -r ${PWD}/vtune_results"
   export my_mode=PROF
fi
if [ $my_prof == "intel" ] ; then
   module load advisor
   export my_ompiparams="$my_ompiparams $ADVISOR_EXEDIR/advixe-cl -collect survey --search-dir src:r=${chimere_root}/src -project-dir ${chimere_root}/advise_results --"
   export my_mode=PROF
fi
if [ $my_prof == "ipm" ] ; then
   module load ipm
   export LD_PRELOAD=$IPM_LIBDIR/libipm.so
   export my_mode=PROF
   #imp_parse -html $tmpdir/*.xml
fi

#---------------------------------------------------------------------------------
# 	SSH-aerosol	
#---------------------------------------------------------------------------------
export my_ssh=${chimere_root}/ssh-aerosol/src/