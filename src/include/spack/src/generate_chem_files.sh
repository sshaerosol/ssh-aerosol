#!/usr/bin/env bash

echo""
echo "**********************************************************"
echo "* This script generates chemical mechanism fortran code"
echo "* from a list of species and reactions."
echo "**********************************************************"
echo""

script_dir=$(dirname $(readlink -f "$0"))
script_name=$(basename $0)


#
# Command line parsing
#

function usage {
  echo
  echo "Usage:"
  echo "$script_name [species file] [reactions file]"
  echo "Arguments:"
  echo "  [species file]: file containing names of chemical species to be used (\"ciMechanism\" located in src/mechanism/)."
  echo "  [reactions file]: file containing related chemical mechanism (\"Mechanism\" located in src/mechanism/)."
  echo
  exit 1
}

if [[ $# != 2 ]]; then
  echo "Error on the number of arguments."
  usage
fi

SPECIES=$1
REACT=$2

if [[ ! -f $SPECIES ]]; then
  echo "Error: the species file \"$SPECIES\" does not exist."
  has_error=1
fi
if [[ ! -f $REACT ]]; then
  echo "Error: the reactions file \"$REACT\" does not exist."
  has_error=1
fi

[[ $has_error ]] && usage


#
# Chemical files generation
#

tmp_files="inSPACK perm.dat ${SPECIES}_preprocessed non_zero.dat"
output_file="dratedc.f fexchem.f fexloss.f fexprod.f jacdchemdc.f kinetic.f rates.f dimensions.f LU_decompose.f LU_solve.f non_zero.dat species.spack.dat"

# Automatically delete temporary files on exit.
function on_exit {
  [[ $? != 0 ]] && echo "[ERROR] ***** Failed to generate Spack files. *****"
  rm -f $tmp_files
}
trap on_exit EXIT

# Removes tabulations in the input files
sed -i 's/\t/    /g' $1 $2

echo "======= 1 - Spack..."
set -e
rm -f $output_file
$script_dir/spack_generator spack_config $SPECIES $REACT

$script_dir/find_perm $SPECIES non_zero.dat perm.dat
$script_dir/spack_input $SPECIES ${SPECIES}_preprocessed perm.dat

rm -f $output_file
$script_dir/spack_generator spack_config ${SPECIES}_preprocessed $REACT

# Those unused files don't even compile:
rm -f fexloss.f fexprod.f

# Replaces strings "E+" by "D+" and "E-" by "D-" in the output to get double precision.
sed -i 's/E\([+-]\)/D\1/g' *.f

echo "======= 2 - LU decomposition and solver(s)..."

# Reads the Spack config (to get the function suffix).
. spack_config

headers_dir=$script_dir/headers
$script_dir/LU_generator ${SPECIES}_preprocessed non_zero.dat $headers_dir/header_license.txt $headers_dir/header_decompose.txt $headers_dir/header_solve.txt $headers_dir/header_solve_tr.txt $function_suffix

echo "======= 3 - Generating Polyphemus configuration files..."
$script_dir/make_config.py ${SPECIES}_preprocessed $REACT 
set +e

echo "Done."
