for ivisco in "0","1" "2" "3" "4" "5" 
do
    ./ssh-aerosol INIT/namelist_visco${ivisco}_poamp.ssh
    #../git-ssh/ssh-aerosol/ssh-aerosol namelist_visco${ivisco}_poalp.ssh
    #../git-ssh/ssh-aerosol/ssh-aerosol namelist_visco${ivisco}_soalp.ssh
done
