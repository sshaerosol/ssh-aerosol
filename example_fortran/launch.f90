
PROGRAM main
  USE SSHaerosol
  IMPLICIT NONE

  INTEGER :: nbins
  INTEGER :: nspecies
  INTEGER :: naero
  INTEGER :: itest
  INTEGER :: nlayers
  INTEGER :: iemis
  DOUBLE PRECISION :: duration,lat,lon,temp,pres,relh
  double precision,allocatable,dimension(:) :: Vlayers
  double precision,allocatable,dimension(:) :: number
  double precision,allocatable,dimension(:) :: gas_conc
  double precision,allocatable,dimension(:,:) :: aero_conc
  character (len=40) :: namelist_ssh  

  itest=1

  if (itest==0) then
     ! Nucleation example
     print*,"Nucleation test case"    

     namelist_ssh="./namelist_modssh.ssh"     

     duration=3600.d0
     pres=1.01325e05
     temp=288.15d0
     relh=0.6d0

     iemis=1

     call initialisation_modssh(namelist_ssh,iemis,naero,nspecies,nbins,nlayers,Vlayers)
     print*,naero,nspecies,nbins,nlayers
     print*,Vlayers

     allocate(number(nbins))
     allocate(gas_conc(nspecies))
     allocate(aero_conc(nbins,naero))

     number=0.d0
     gas_conc=0.d0
     aero_conc=0.d0
     gas_conc(5)=8.0d0
     aero_conc(:,4)=(/7.222506E-19, 1.191796E-17, 1.783646E-16, 2.421107E-15, 2.980732E-14, 3.328443E-13, 3.371117E-12, 3.096894E-11, 2.580500E-10, 1.950345E-09, 1.337071E-08, 8.314528E-08, 4.689924E-07, 2.399630E-06, 1.113844E-05, 4.714408E-05, 1.987428E-04, 1.190799E-03, 8.254623E-03, 3.246224E-02, 5.914383E-02, 6.289222E-02, 7.917661E-02, 1.392658E-01, 2.426590E-01, 3.854985E-01, 5.556413E-01, 7.265657E-01, 8.619621E-01, 9.278976E-01, 9.067632E-01, 8.053857E-01, 6.525903E-01, 4.879282E-01, 3.484007E-01, 2.601101E-01, 2.374289E-01, 2.873215E-01, 4.142062E-01, 6.217591E-01, 9.106726E-01, 1.273746E+00, 1.690945E+00, 2.127184E+00, 2.534756E+00, 2.860762E+00, 3.057965E+00, 3.095896E+00, 2.968541E+00, 2.695891E+00/)
     aero_conc(:,5)=(/3.095360E-19, 5.107696E-18, 7.644198E-17, 1.037617E-15, 1.277457E-14, 1.426476E-13, 1.444764E-12, 1.327240E-11, 1.105928E-10, 8.358622E-10, 5.730304E-09, 3.563369E-08, 2.009968E-07, 1.028413E-06, 4.773617E-06, 2.020461E-05, 8.517548E-05, 5.103422E-04, 3.537695E-03, 1.391239E-02, 2.534736E-02, 2.695381E-02, 3.393283E-02, 5.968534E-02, 1.039967E-01, 1.652136E-01, 2.381320E-01, 3.113853E-01, 3.694123E-01, 3.976704E-01, 3.886128E-01, 3.451653E-01, 2.796816E-01, 2.091121E-01, 1.493146E-01, 1.114758E-01, 1.017553E-01, 1.231378E-01, 1.775169E-01, 2.664682E-01, 3.902882E-01, 5.458911E-01, 7.246907E-01, 9.116503E-01, 1.086324E+00, 1.226041E+00, 1.310556E+00, 1.326812E+00, 1.272232E+00, 1.155382E+00/)

     !To create the coefficient repartition file
     ! call initialisation_repart_coeff(namelist_ssh)

     !Launch SSH
     print*,"Launch SSH-aerosol for nucleation test case"
     call launch_ssh_aerosolonly(namelist_ssh,iemis,naero,nbins,nspecies,duration,temp,pres,relh,gas_conc,aero_conc,number)
     print*,"SSH-aerosol suceeded"
     !call launch_ssh_aerosolonly(namelist_ssh,iemis,naero,nbins,nspecies,duration,temp,pres,relh,gas_conc,aero_conc,number)
     !print*,"C"

     print*,"Number concentrations: ",number
     print*,"Ammoniac concentration: ",gas_conc(5)
     print*,"Ammonium: ",aero_conc(:,5)

  else ! Visco 5 POAmP test case     
     
     namelist_ssh="./namelist_modssh_visco.ssh"
     
     duration=2419200.0d0 !3600.d0
     pres=1.01325e05
     temp=298.0
     relh=0.1d0

     iemis=0
     
     call initialisation_modssh(namelist_ssh,iemis,naero,nspecies,nbins,nlayers,Vlayers)
     print*,naero,nspecies,nbins,nlayers
     print*,Vlayers
         
     allocate(number(nbins))
     allocate(gas_conc(nspecies))
     allocate(aero_conc(nbins,naero))

     number=0.d0
     gas_conc=0.d0
     aero_conc=0.d0
     
     gas_conc(29)=5.0d0
     gas_conc(25)=0.01d0
     aero_conc(1,93:97)=5.0d0*Vlayers(:)
     print*,aero_conc(1,93:97)

     print*,"Launch SSH-aerosol for nucleation test case"
     call launch_ssh_aerosolonly(namelist_ssh,iemis,naero,nbins,nspecies,duration,temp,pres,relh,gas_conc,aero_conc,number)
     print*,"SSH-aerosol suceeded"

     print*,"POAmP gas concentration: ",gas_conc(29)
     print*,"POAmP aerosol concentration: ",sum(aero_conc(1,113:117))
     print*,"SOAlP aerosol concentration: ",sum(aero_conc(1,93:97))
     print*,"Number: ",number
  endif
  
END PROGRAM main
