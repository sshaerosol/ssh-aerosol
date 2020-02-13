#!/bin/bash

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
make

cp ../../CHEMISTRY/${chemistry}/spack_config ./
cp ../../CHEMISTRY/${chemistry}/*.species ./species
cp ../../CHEMISTRY/${chemistry}/*.reactions ./reactions
./generate_chem_files.sh species reactions
cp *.f90 ../../CHEMISTRY/common/

cd ../../../../
make


