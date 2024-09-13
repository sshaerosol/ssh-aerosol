#!/bin/sh

# namelist
namelist="INIT/namelist_mcm-wSMILES.ssh"

# change tag_RO2 and result forder
for i in 0 1 2 3
do
    # change tag_RO2
    sed "s/tag_RO2 = 3,/tag_RO2 = $i,/1" $namelist > $namelist.$i
    # change result folder
    sed -i "s/output_directory = \"results\/mcm-wSMILES\",/output_directory = \"results\/mcm-wSMILES-tag_RO2-$i\",/1" $namelist.$i
    # run simulation
    ./ssh-aerosol $namelist.$i >> toto$i
    # move
    resultdir="results/mcm-wSMILES-tag_RO2-$i/"
    mkdir -p $resultdir
    mv $namelist.$i results/mcm-wSMILES-tag_RO2-$i/
    mv toto$i results/mcm-wSMILES-tag_RO2-$i/
    # print
    echo "Finish mcm simulation with tag_RO2 = $i. Results see: results/mcm-wSMILES-tag_RO2-$i/"
done

