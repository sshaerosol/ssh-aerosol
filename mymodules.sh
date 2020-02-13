#!/bin/bash

module purge
module load intel/17.0.4.196
module load flavor/hdf5/parallel
module load hdf5
module load netcdf-fortran
module load pnetcdf
module load grib
module load jasper
module load blitz
#module switch mpi/openmpi/1.8.8

## Codes
export my_mpiframe=openmpi
export my_ompiparams="-mca btl self,sm,tcp -mca btl_tcp_eager_limit 4095 -x LD_LIBRARY_PATH"
export my_ompiparams=""
export my_lamparams="-ssi rpi sysv -ssi rpi_tcp_short 4095 -ssi rpi_sysv_short 4095"
export my_hostfile=
export my_make="make"
export my_awk=awk
export my_netcdfdir=${NETCDFFORTRAN_ROOT}
export my_netcdfdirc=${NETCDFC_ROOT}
export my_ncdump=${NETCDF_EXEDIR}/ncdump
export my_netcdflib=${NETCDFFORTRAN_LIBDIR}
export my_pnetcdflib=${PNETCDF_LIBDIR}
export my_pnetcdfinc=${PNETCDF_INCDIR}
export my_netcdfinc=${NETCDFFORTRAN_INCDIR}
export my_netcdflibc=${NETCDFC_LIBDIR}
export my_hdfdir=${HDF5_ROOT}
export my_hdflib=${HDF5_LIBDIR}
export my_hdfinc=${HDF5_INCDIR}
export my_mpidir=${MPI_ROOT}
export my_mpilib=${my_mpidir}/lib
export my_mpiinc=${my_mpidir}/include
export my_mpibin=${MPI_EXEDIR}
export my_mpirun=${my_mpibin}/mpirun
export my_mpirun=ccc_mprun


export my_hdr=Makefile.hdr.ifort-64-ompi
export my_compile_ecm=no
export my_clean_source=yes
export my_super_clean=no

export my_inteldir=${FORTRAN_INTEL_ROOT}/bin/intel64
source ${FORTRAN_INTEL_ROOT}/bin/compilervars.sh intel64
source ${C_INTEL_ROOT}/bin/compilervars.sh intel64

# only for ifort
export my_ifort=${my_inteldir}/ifort
export my_gfortran=ifort


# C++ Compiler
export my_icc=${C_INTEL_EXEDIR}/icc

# only if you use OpenMPI
export my_mpif90=${my_mpibin}/mpif90

# only if you use LAM
export my_mpif77=${my_mpibin}/mpif77

# Needed for mct compilation
export my_mpicc=mpicc

# only if you use the GRIB_API library
export my_gribapi=${GRIB_ROOT}
export my_griblib=${GRIB_LIBDIR}
export my_gribinc=${GRIB_INCDIR}
export my_jasperlib=${JASPER_LIBDIR}
#export my_blitzdir=/usr/lib64
export my_blitzinc=$BLITZ_ROOT/include/

# dynamic libraries for MPI
#export LD_LIBRARY_PATH=${my_netcdflib}:${my_netcdflibc}:${my_mpilib}:${my_mpilib}/openmpi:${my_griblib}:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${my_netcdflib}:${my_netcdflibc}:${my_mpilib}:${my_mpilib}/openmpi:${my_griblib}:$LD_LIBRARY_PATH

# Name of the configure.wrf file. The file must be located within the "${chimere_root}/config_wrf" directory.
export configure_wrf_file_name="configure.wrf.ifort_curie"

## Mode
export my_mode=PROD
export my_bigarray=Yes
