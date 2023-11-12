#!/bin/bash

export SSH_dir=`pwd`

# Deal with prompt 
while (($# > 0))
do
    case $1 in
        "-h"|"--h"|"--help"|"-help") 
            echo "build-ssh.sh - compiles SSH-aerosol to use with CHIMERE"
            echo "build-ssh.sh [options]"
            echo "      [ -h | --h | -help | --help ] : this help message"
            echo "      [--arch arch] : to choose target architecture (file ${chimere-folder}/mychimere/mychimere-arch must exist - default : last architecture used)"
            exit ;;
	"--arch") arch=$2 ; arch_defined=true ; shift ; shift ;;
        *) code=$1 ; shift ;;
    esac
done
compil=${arch##*.}

if test -f ./makefiles.hdr/Makefile.hdr.${compil}-64-ompi; then
    ln -vfns ./makefiles.hdr/Makefile.hdr.${compil}-64-ompi Makefile.hdr
else
    echo "Error in architecture definition. The file ./makefiles.hdr/Makefile.hdr.${compil}-64-ompi does not exist. Please check --arch [architecture choice]"
    exit 1
fi

if test -f ../mychimere/mychimere-${arch}; then
    ln -vfns ../mychimere/mychimere-${arch} mymodules.sh
else
    echo "Error in architecture definition. The file ../mychimere/mychimere-${arch} does not exist. Please run ./make --arch [architecture choice]"
    exit 1
fi

source mymodules.sh

# Read the option for Chemistry kinetic mechanism.
# Accepted values: racm, racm2, cb05
# Default: chemistry=cb05
for i in "$@"
do
case $i in
    -c=*|--chemistry=*)
    chemistry="${i#*=}"
    shift
    ;;
esac
done
if [ -z $chemistry ]; then
    chemistry="cb05"
fi
cd src/include/spack/src/
make clean ; make
cp ../../CHEMISTRY/${chemistry}/spack_config ./
cp ../../CHEMISTRY/${chemistry}/*.species ./species
cp ../../CHEMISTRY/${chemistry}/*.reactions ./reactions
./generate_chem_files.sh species reactions
cp *.f90 ../../CHEMISTRY/common/

cd ${SSH_dir}
make


