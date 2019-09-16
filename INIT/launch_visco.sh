for ivisco in `seq 1 8`
do
    ./ssh-aerosol INIT/namelist_visco${ivisco}_poamp.ssh
    ./ssh-aerosol INIT/namelist_visco${ivisco}_poalp.ssh
    ./ssh-aerosol INIT/namelist_visco${ivisco}_soalp.ssh
done
