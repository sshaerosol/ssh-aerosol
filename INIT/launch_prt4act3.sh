#!/bin/bash

# This script modifies namelist_prt4act3.sh and runs generated namelists with ./ssh-aerosol

namelist="INIT/namelist_prt4act3.ssh"
nml="INIT/prt4act3.nml"
resf="results/prt4act3/"
activity_models="1 2 3"
aerosol_list_suffixes="HPHI HPHO BOTH ''"

mkdir -p ${resf}

echo "Cleaning up previous results in ${resf} ..."
rm -rf ${resf}/*

for i in $activity_models; do
    for j in $aerosol_list_suffixes; do
        cp $namelist $nml
        
        sed -i "s/activity_model = [1-3]/activity_model = $i/g" $nml
        
        if [ "$j" = "''" ]; then
            output_dir="default${i}"
        else
            sed -i "s/species-list-aer.dat/species-list-aer-${j}.dat/g" $nml
            output_dir="${j}${i}"
        fi
        
        # Replace the line contains "output_directory = " with "output_directory = ${resf}${output_dir}"
        sed -i "/output_directory = /c\output_directory = \"${resf}${output_dir}\"," $nml
        
        echo "Running ${namelist} with activity_model = ${i} and aerosol_list_suffix = ${j} ..."
        
        ./ssh-aerosol $nml
        if [ $? -eq 0 ]; then
            mv $nml ${resf}${output_dir}/
        else
            exit 1
        fi
    done
done
