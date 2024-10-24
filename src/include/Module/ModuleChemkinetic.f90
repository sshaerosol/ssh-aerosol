module mod_sshchemkinetic
  
  use aInitialization, only : attenuation, n_gas, &
       relative_humidity, tag_genoa, &
                          humidity, temperature, pressure, &
                          photo_rcn, TB_rcn, fall_rcn, extra_rcn, &
                          index_RCT, index_PDT, index_PDT_ratio, &
                          fall_coeff, extra_coeff, &
                          photo_ratio, Arrhenius, &
                          photo_ratio_read, &
                          ratio_PDT, kinetic_rate, chem_prod, chem_loss, &
                          drv_knt_rate, rcn_rate, gas_yield, species_name, &
                          concentration_gas_all, SumMc, YlH2O, &
                          szas, nsza, nRO2_group, RO2index, nRO2_chem, &
                          hetero_rcn, hetero_ind, &
                          concentration_mass, concentration_number, &
                          diam_bound, fixed_density, &
                          lwc_cloud_threshold, mass_density, molecular_weight, n_aerosol, &
                          n_fracmax, n_size, wet_diameter, &
                          with_fixed_density, with_heterogeneous, &
                          n_sizebin, &
                          wall_rcn, wall_coeff, & ! wall loss
                          irdi_rcn, irdi_ind ! irreversible dicarbonyl

  use mod_meteo
  
  implicit none
  
contains

subroutine ssh_compute_ro2(RO2s)

! compute ro2 concentrations

  implicit none

  integer i,j,s
  double precision,INTENT(OUT) :: RO2s(nRO2_group)
  
  RO2s = 0d0 ! init
  do i = 1, nRO2_chem
    s = RO2index(i,1)! isps
    j = RO2index(i,2)! igroup
    RO2s(j) = RO2s(j) + gas_yield(s)
  enddo

end subroutine ssh_compute_ro2

SUBROUTINE ssh_gck_compute_gas_phase_water(temp0, rh0, water)
  
  ! compute water in the gas phase - from GECKO BOXMODEL
    IMPLICIT NONE

    double precision, INTENT(IN) :: temp0       ! input temperature
    double precision, INTENT(IN) :: rh0        ! input relative humidity
    double precision, INTENT(OUT) :: water     ! output water concentration (molec/cm3)

    ! Constants
    double precision, PARAMETER :: avogadro=6.02214d23    ! avogadro number
    double precision, PARAMETER :: Rgas=8.3144621d0 ! gas constant (J.K-1.mol-1)
    double precision, PARAMETER :: c1 = 610.94d0, c2 = 17.625d0, c3 = 243.04d0  ! Magnus formula parameters

    double precision :: TC, psat_H2O, p_H2O

    ! Convert temperature to Celsius
    TC = temp0 - 273.16d0

    ! Calculate saturation water pressure
    psat_H2O = c1 * dEXP(c2 * TC / (c3 + TC))

    ! Calculate water pressure
    p_H2O = psat_H2O * rh0 / 1d2

    ! Calculate water concentration (molec/cm3)
    water = p_H2O / (Rgas * temp0 * 1.d6 / avogadro)

END SUBROUTINE ssh_gck_compute_gas_phase_water


subroutine ssh_dratedc()
 
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This routine computes the derivative of reaction  rates.    
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     kinetic_rate: kinetic rates.
!     gas_yield: chemical concentrations.
!
!     -- INPUT/OUTPUT VARIABLES
!
!     -- OUTPUT VARIABLES
!
!     drv_knt_rate: derivative of reaction rates wrt Y.
!
!------------------------------------------------------------------------
 
  implicit none

  integer i,j,k,s,js,ntot
  integer, parameter :: nrct = 1 ! no. possible other reactants
  
  drv_knt_rate = 0.d0
  ntot = size(drv_knt_rate)
  
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,s,js)
  do i=1, ntot
    j = index_RCT(i,1) ! ircn
    k = index_RCT(i,2) ! isps
    
    drv_knt_rate(i) = kinetic_rate(j)

    ! check for another reactants - now check only one
    do s = max(1,i-nrct), min(i+nrct,ntot)
       if (s.ne.i .and. index_RCT(s,1).eq.j) then ! find
         js = index_RCT(s,2) ! isps of s
         drv_knt_rate(i) = drv_knt_rate(i) * gas_yield(js)
       endif
    enddo
  enddo
!$OMP END PARALLEL DO

END subroutine ssh_dratedc

! =================================================================

subroutine ssh_fexloss()
 
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This routine computes the chemical loss  term L in a P-Lc formulation.     
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     drv_knt_rate: derivative of reaction rates wrt Y.
!
!     -- INPUT/OUTPUT VARIABLES
!
!     -- OUTPUT VARIABLES
!
!     chem_loss: array of chemical loss terms.
!
!------------------------------------------------------------------------
 
  implicit none
 
  integer i,k

! Chemical loss terms.
  chem_loss = 0.d0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,k)
  do i=1, size(index_RCT,1)
    k = index_RCT(i,2) ! isps
    chem_loss(k) = chem_loss(k) + drv_knt_rate(i)
  enddo
!$OMP END PARALLEL DO

END subroutine ssh_fexloss

subroutine ssh_fexprod()
 
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This routine computes the production  term P in a P-Lc formulation.      
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     rcn_rate: reaction rates.
!
!     -- INPUT/OUTPUT VARIABLES
!
!     -- OUTPUT VARIABLES
!
!     chem_prod: array of chemical production terms.
!
!------------------------------------------------------------------------
 
  implicit none

  integer i,j,k,s,s0,s1,s2 ! s for ratios
  double precision :: ratio

! Chemical production terms.
  chem_prod = 0.d0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,s,s0,s1,s2)
  do i=1, size(index_PDT,1)
    j = index_PDT(i,1) ! index in reaction list
    k = index_PDT(i,2) ! index in species list
    s = index_PDT(i,3) ! index in ratio list

    if (k == 0) cycle   ! empty product
    if (s == 0) then    ! no read ratio
        ratio = 1.d0
    elseif (s > 0) then ! Read ratio
        ratio = ratio_PDT(s)
    else ! Read ratio as a function
        s = abs(s)
        s0 = index_PDT_ratio(s,1) ! ratio index
        s1 = index_PDT_ratio(s,2) ! k1
        s2 = index_PDT_ratio(s,3) ! k2
        ! k1 + k2
        ratio = kinetic_rate(s1) + kinetic_rate(s2)
        if (ratio <= 0.d0) then
            ratio = 5d-1 * ratio_PDT(s0)
        else ! ratio = k1/(k1+k2+...)
            ratio = (kinetic_rate(s1)/ratio) * ratio_PDT(s0)
        endif
        !print*, "ratio: ",k,ratio_PDT(s0),s1,s2,j,ratio
    endif
    
    ! Get production
    chem_prod(k) = chem_prod(k) + rcn_rate(j) * ratio

  enddo
!$OMP END PARALLEL DO

END subroutine ssh_fexprod

! =================================================================

subroutine ssh_rates()
 
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This routine computes the reaction rates.    
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     kinetic_rate: kinetic rates.
!     gas_yield: chemical concentrations.
!
!     -- INPUT/OUTPUT VARIABLES
!
!     -- OUTPUT VARIABLES
!
!     rcn_rate: reaction rates.
!
!------------------------------------------------------------------------
 
  implicit none

  integer i,j,k 

  rcn_rate = kinetic_rate ! init
  
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k)
  do i=1, size(index_RCT,1) ! number of reactants
    j = index_RCT(i,1) ! ircn
    k = index_RCT(i,2) ! index in species list
    rcn_rate(j) = rcn_rate(j) * gas_yield(k)
  enddo
!$OMP END PARALLEL DO

END subroutine ssh_rates

! =================================================================

subroutine ssh_basic_kinetic(ro2_basic_rate, nro2_s)
     
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This routine computes the kinetic rates for the gas-phase.
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     temperature: temperature ([K]).
!     humidity: water massic fraction.
!     pressure: pressure ([Pa]).
!
!     -- INPUT/OUTPUT VARIABLES
!
!     -- OUTPUT VARIABLES
!
!     kinetic_rate: kinetic rates.
!
!------------------------------------------------------------------------
 
  implicit none

  integer, intent(in) :: nro2_s
  double precision, INTENT(OUT) :: ro2_basic_rate(nro2_s)

  integer :: i,j,k,s


  ! Arrhenius ! k = C1 * T**C2 * dEXP(-C3/T)
  ! Define labels for extra reactions
  integer, PARAMETER ::geckoLabel(6) = (/100, 200, 500, 501, 502, 550/)
  integer, PARAMETER ::mcmLabel(4) = (/91, 92, 93, 94/)

  ! Heterogeneous reactions
  DOUBLE PRECISION dsf_aero(n_size)

  ! Wall loss
  double precision :: Psat, Masmol, cstar, aw, cbar, awc, denom

  ! Irreversible dicarbonyl
  double precision :: facteur, xlw

  if (humidity == -999.d0) then
     call ssh_compute_sh(relative_humidity, temperature, &
          pressure, humidity)
  endif
  
  xlw = humidity
  
  do i=1, size(kinetic_rate) ! reactions
    kinetic_rate(i) = Arrhenius(i,1) * (temperature**Arrhenius(i,2)) &
                      * dEXP(-Arrhenius(i,3)/temperature)
  enddo

  ! falloff reactions
  do i=1, size(fall_rcn)
    j = fall_rcn(i) ! reaction index
    call ssh_gck_forate(j,i)
  enddo
  
  ! extra reactions
  do i=1, size(extra_rcn)
  
    j = extra_rcn(i) ! reaction index
    k = INT(extra_coeff(i,1)) ! get kinetic label

    ! GECKO
    if (any(k == geckoLabel)) then ! 100,200,500,501,502,550
        call ssh_gck_extrarate(j,i,k) ! GECKO EXTRA reactions & ISOM

    ! MCM
    else if (any(k == mcmLabel)) then ! 91,92,93,94

      if (k.eq.94) then ! KRO2: qfor = 1.26D-12 * RO2
        ! in the format: KINETIC RO2 [n] EXTRA 94 1. 0. 0.
        kinetic_rate(j) = 1.26D-12
      else if (k.ne.91) then ! No photolysis
        call ssh_MCM_rate(j,i,k,9d1)
      endif
    
    ! Additional kinetic
    else if (k .eq. 99) then
        call ssh_genoa_spec(j,i)

    ! From SPACK - need to be completed !!!
    else if ((k .eq. 10) .or. (k .eq. 20)) then
        call ssh_spack_spec(j,i,k)

    else ! print label, ircn, iex
        print*, "EXTRA type unknown: ",k,j,i
        stop
    endif
  enddo
  
  ! need to be computed LAST !!!
  ! TBs(5) = ["O2 ", "H2O", "M  ", "N2 ", "H2 "]

  s = 0 ! count no.ro2 reaction
  do i=1,size(TB_rcn,1) ! no.TB reactions
  
    j = TB_rcn(i,1) ! reaction index
    k = TB_rcn(i,2) ! TB (+) or RO2s (-) index

    if (k.gt.0) then ! TB
      select case (k)
        case(1) ! O2
            kinetic_rate(j) = kinetic_rate(j) * SumMc * 0.2d0
        case(2) ! H2O
            kinetic_rate(j) = kinetic_rate(j) * YlH2O
        case(3) ! M
            kinetic_rate(j) = kinetic_rate(j) * SumMc
        case(4) ! N2
            kinetic_rate(j) = kinetic_rate(j) * SumMc * 0.8d0
        case(5) ! H2
            kinetic_rate(j) = kinetic_rate(j) * SumMc * 5.8d-7
      end select
    else ! RO2-RO2 reactions
      s = s + 1
      ro2_basic_rate(s) = kinetic_rate(j)
      kinetic_rate(j) = 0 ! Reset
    endif
  enddo

  !==============================!
  !=== heterogenous reactions ===!
  !==============================!
  do i=1, size(hetero_rcn)

     j = hetero_rcn(i) ! reaction index
     
     if (with_heterogeneous == 1) then

        ! Compute particle diameter.
        call ssh_compute_granulo(fixed_density, not(with_fixed_density), &
             n_size, n_aerosol, concentration_mass, mass_density, &
             concentration_number, n_fracmax, diam_bound, dsf_aero)
        

        ! liquid_water_content in kg/kg (QCLOUD from WRF)
        call ssh_hetero_rate(i, j, hetero_ind(i), &
             n_gas, n_sizebin, temperature, pressure, &
             wet_diameter, concentration_number, &
             dsf_aero, molecular_weight, lwc_cloud_threshold, &
             kinetic_rate(j))
        
     else
        kinetic_rate(j) = 0.d0
     endif
  
  enddo

  !===========================!
  !=== Wall loss reactions ===!
  !===========================!
  do i = 1, size(wall_rcn)

     j = wall_rcn(i) ! reaction index

     ! Check coefficients
     if (wall_coeff(i, 5) <= 0.d0) then
        write(*,*) "Error in WALL coefficient", wall_coeff(i, :)
        stop
     endif
     
     Psat =  wall_coeff(i, 2) / 760. * dexp(wall_coeff(i, 3) * 1000. &
          / 8.314*(1./298.-1./temperature))
     Masmol =  wall_coeff(i, 5)
     cstar = Masmol * 1d+6 * Psat / (8.205d-5 * temperature)
     aw = 10**(-0.27440D+01)*cstar**(-0.14066D+01)
     cbar = (8000*8.314*temperature/(3.14*Masmol))**0.5
     awc = (aw*cbar/4.)
     denom = awc / (( wall_coeff(i, 4) * 0.50000D-05)**0.5)
     kinetic_rate(j) =  0.33330D+02 * awc / (1. + 1.5708 * denom)
  
  enddo

  !===============================!
  !=== Irreversible dicarbonyl ===!
  !===============================!
  do i=1, size(irdi_rcn)

     j = irdi_rcn(i) ! reaction index
     
     !if (with_heterogeneous == 1) then

     ! Compute particle diameter.
     call ssh_compute_granulo(fixed_density, not(with_fixed_density), &
              n_size, n_aerosol, concentration_mass, mass_density, &
              concentration_number, n_fracmax, diam_bound, dsf_aero)
        
     ! liquid_water_content in kg/kg (QCLOUD from WRF)
     call ssh_dicarb_rate(i, j, irdi_ind(i), &
!              n_sizebin, temperature, pressure, option_cloud, &
              n_sizebin, temperature, pressure, &
              wet_diameter, concentration_number, &
!              dsf_aero, lwc_cloud_threshold, &
!              cloud_water, concentration_gas_all(idOH),&
              concentration_gas_all(idOH),&
              kinetic_rate(j))
        
     !else
     !   kinetic_rate(j) = 0.d0
     !endif
     
  
  enddo
  

  
  
END subroutine ssh_basic_kinetic

subroutine ssh_update_kinetic_pho(azi)
     
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This routine update photolysis kinetic rate.
!
!     Need to have basic rate first.
!
!------------------------------------------------------------------------
 
  implicit none

  integer :: i,j,k,ik,s,tag
  double precision,INTENT(IN) :: azi
  double precision :: photo, zp,aziloc

  ! Photolysis

  aziloc=min(azi,90.d0)
  do i=1, size(photo_rcn,1)
     j = photo_rcn(i,1) ! reaction index
     k = photo_rcn(i,2) ! pholysis index
     if (option_photolysis == 1) then
        ! Get zone s in the range of [1, nsza-1]
        if (aziloc .ge. szas(nsza)) then
           s = nsza - 1
        else
           do s=1, nsza-1
              if (aziloc.ge.szas(s) .and. aziloc.lt.szas(s+1)) then
                 exit
              endif
           enddo
        endif
        ! Get photolysis rate in zone s
        if (k.lt.0) then ! spack
           ik = abs(k)
           photo = photo_ratio(ik,s,4)
           photo = photo_ratio(ik,s,3) + (aziloc-szas(s))*photo
           photo = photo_ratio(ik,s,2) + (aziloc-szas(s))*photo
           photo = photo_ratio(ik,s,1) + (aziloc-szas(s))*photo
        else ! gck
           zp = (aziloc-szas(s))/(szas(s+1)-szas(s)) !zp=(xin-xp(i-1))/xp(i)-xp(i-1)
           photo =((zp*photo_ratio(k,s,1) + photo_ratio(k,s,2)) &
                *zp + photo_ratio(k,s,3))*zp + photo_ratio(k,s,4)
           photo = dabs(photo)
        endif
        ! Read from binary files     
     else if (option_photolysis == 2) then
        photo = photolysis_rate(abs(k))
     endif

     ! Add cloud attenuation and a ratio (Arrhenius(j, 1)
     kinetic_rate(j) = Arrhenius(j,1) * max(photo*attenuation, 0.d0)
  enddo
  !else ! no photolysis
  !  do i=1, size(photo_rcn,1) ! no.TB reactions
  !      j = photo_rcn(i,1) ! reaction index
   !     kinetic_rate(j) = 0.d0 ! set to zero
  !  enddo
  !endif
  
  ! extra reactions - MCM photolysis
  do i=1, size(extra_rcn)
    j = extra_rcn(i) ! reaction index
    k = INT(extra_coeff(i,1)) ! get kinetic label

    ! MCM photolysis
    if (k == 91) then
        call ssh_MCM_rate(j,i,k,aziloc)
    endif
  enddo

END subroutine ssh_update_kinetic_pho 

subroutine ssh_kinetic(azi,RO2s)
     
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This routine computes the kinetic rates for the gas-phase.
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     temperature: temperature ([K]).
!     humidity: water massic fraction.
!     pressure: pressure ([Pa]).
!     AZI: zenithal angle ([degree]).
!     attenuation: attenuation variable.
!
!     -- INPUT/OUTPUT VARIABLES
!
!     -- OUTPUT VARIABLES
!
!     kinetic_rate: kinetic rates.
!
!------------------------------------------------------------------------
 
  implicit none

  integer :: i,j,k,ik,s,tag
  double precision,INTENT(IN) :: RO2s(nRO2_group),azi

  !integer,INTENT(IN) :: photo_rcn(:,:), TB_rcn(:,:)
  !integer,INTENT(IN) :: fall_rcn(:), extra_rcn(:)

  ! for photolysis
  ! nsza, szas(11)
  double precision :: photo, zp

  ! Arrhenius ! k = C1 * T**C2 * dEXP(-C3/T)
  ! Define labels for extra reactions
  integer, PARAMETER ::geckoLabel(6) = (/100, 200, 500, 501, 502, 550/)
  integer, PARAMETER ::mcmLabel(4) = (/91, 92, 93, 94/)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
  do i=1, size(kinetic_rate) ! reactions
    kinetic_rate(i) = Arrhenius(i,1) * (temperature**Arrhenius(i,2)) &
                      * dEXP(-Arrhenius(i,3)/temperature)
  enddo
!$OMP END PARALLEL DO

  ! Photolysis
  if (azi .lt. 9d1) then ! with photolysis

    do i=1, size(photo_rcn,1)
        j = photo_rcn(i,1) ! reaction index
        k = photo_rcn(i,2) ! pholysis index
        
        tag = 0 ! tag for computing photolysis rate

        do s=1, nsza-1 ! get photolysis rate
            if (azi.ge.szas(s) .and. azi.lt.szas(s+1)) then ! find zone
                if (k.lt.0) then ! spack
                    ik = abs(k)
                    photo = photo_ratio(ik,s,4)
                    photo = photo_ratio(ik,s,3) + (azi-szas(s))*photo
                    photo = photo_ratio(ik,s,2) + (azi-szas(s))*photo
                    photo = photo_ratio(ik,s,1) + (azi-szas(s))*photo
                else ! gck
                    zp = (azi-szas(s))/(szas(s+1)-szas(s)) !zp=(xin-xp(i-1))/xp(i)-xp(i-1)
                    photo =((zp*photo_ratio(k,s,1) + photo_ratio(k,s,2)) &
                            *zp + photo_ratio(k,s,3))*zp + photo_ratio(k,s,4)
                    photo = dabs(photo)
                endif
                tag = 1
                exit
            endif
        enddo
        
        if (tag.eq.0) then
            print*, "Not find photolysis info: ",j,k,azi
            stop 
        endif

       ! Cloud attenuation.
       kinetic_rate(j) = kinetic_rate(j) * max(photo*attenuation, 0.d0) ! add a ratio
    enddo

  else ! no photolysis

    do i=1, size(photo_rcn,1) ! no.TB reactions
        j = photo_rcn(i,1) ! reaction index
        kinetic_rate(j) = 0.d0 ! set to zero
    enddo

  endif

  ! falloff reactions
  do i=1, size(fall_rcn)
  
    j = fall_rcn(i) ! reaction index
    call ssh_gck_forate(j,i)
  enddo
  
  ! extra reactions
  do i=1, size(extra_rcn)
  
    j = extra_rcn(i) ! reaction index
    k = INT(extra_coeff(i,1)) ! get kinetic label
    
    ! GECKO
    if (any(k == geckoLabel)) then ! 100,200,500,501,502,550

        call ssh_gck_extrarate(j,i,k) ! GECKO EXTRA reactions & ISOM

    ! MCM
    else if (any(k == mcmLabel)) then ! 91,92,93,94

      if (k.eq.94) then ! KRO2: qfor = 1.26D-12 * RO2
        ! in the format: KINETIC RO2 [n] EXTRA 94 1. 0. 0.
        kinetic_rate(j) = 1.26D-12
      else 
        call ssh_MCM_rate(j,i,k,azi)
      endif
    
    ! Additional kinetic
    else if (k .eq. 99) then
        call ssh_genoa_spec(j,i)

    ! From SPACK - need to be completed !!!
    else if (k .eq. 10) then
    
        call ssh_spack_spec(j,i,k)

    else ! print label, ircn, iex
        print*, "EXTRA type unknown: ",k,j,i
        stop
    endif

  enddo
  
  ! need to be computed LAST !!!
  ! TBs(5) = ["O2 ", "H2O", "M  ", "N2 ", "H2 "]
  do i=1,size(TB_rcn,1) ! no.TB reactions
  
    j = TB_rcn(i,1) ! reaction index
    k = TB_rcn(i,2) ! TB (+) or RO2s (-) index

    if (k.gt.0) then ! TB
      select case (k)
        case(1) ! O2
            kinetic_rate(j) = kinetic_rate(j) * SumMc * 0.2d0
        case(2) ! H2O
            kinetic_rate(j) = kinetic_rate(j) * YlH2O
        case(3) ! M
            kinetic_rate(j) = kinetic_rate(j) * SumMc
        case(4) ! N2
            kinetic_rate(j) = kinetic_rate(j) * SumMc * 0.8d0
        case(5) ! H2
            kinetic_rate(j) = kinetic_rate(j) * SumMc * 5.8d-7
      end select
    else ! RO2-RO2 reaction & !KRO2
        kinetic_rate(j) = kinetic_rate(j) * RO2s(abs(k))
    endif
  enddo

END subroutine ssh_kinetic 

! =================================================================

SUBROUTINE ssh_gck_forate(ire,ifo)

!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This routine computes kinetic rates for FALL OFF reactions.
!     Extract from GECKO BOX MODEL    
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     ire: index in reaction list
!     ifo: index in falloff reaction list
!     fall_coeff: 
!     Arrhenius: 
!     SumMc
!
!
!     -- OUTPUT VARIABLES
!
!     qfor: kinetic_rate.
!
!------------------------------------------------------------------------

    IMPLICIT NONE 
                                     
    ! INPUT                                                                 
    INTEGER, INTENT(in) :: ire,ifo 
                                                                    
    ! LOCAL                                                                 
    double precision :: xkom,xki,kratio,factor
    ! patch MCM - see http://mcm.york.ac.uk/parameters/complex.htt for      
    ! details                                                               
    double precision :: fcent,factor2,qfor
                                                                        
! expression fall off type TROE                                         
                                                                        
! First IF uses the usual Troe equation to compute rate constant        
! xkom: ko*M  (effective low pressure rate coef.) /// xki: k_infiny     
! JMLT: NOTE: Arrhenius(1) is NO LONGER used in log form by this point

     xki = Arrhenius(ire,1)*((temperature/3d2)**Arrhenius(ire,2)) &
    &      *dEXP(-Arrhenius(ire,3)/temperature)                                   
     xkom = fall_coeff(ifo,1)*((temperature/3d2)**fall_coeff(ifo,2)) &
    &     * dEXP(-fall_coeff(ifo,3)/temperature)*SumMc
      
     if (xki .ne. 0.d0) then
       kratio = xkom/xki
     else
       print*,ire,ifo,"MCk: xki is zero",Arrhenius(ire,:), fall_coeff(ifo,:)
       stop
     endif 
      
     factor = 1./(1.+(dLOG10(kratio))**2.) 
                                                                        
     IF ((fall_coeff(ifo,4)/=0.d0) .AND. (fall_coeff(ifo,5)==0.d0)) THEN 
       qfor = (xkom/(1.+kratio)) * (fall_coeff(ifo,4)**factor) 
                                                                        
! The following IF use a different equation for MCM rates constant      
     ELSE IF ((fall_coeff(ifo,4)/=0.d0) .AND. (fall_coeff(ifo,5)==1.d0)) THEN 
       factor2 = 1./(1.+(dLOG10(kratio)/(0.75-1.27*dLOG10(fall_coeff(ifo,4))))**2.)                                              
       qfor = (xkom/(1.+kratio)) * (fall_coeff(ifo,4)**factor2) 
                                                                        
     ELSE IF ((fall_coeff(ifo,4)==0.d0) .AND. (fall_coeff(ifo,5)==204.) .AND. &
    &          (fall_coeff(ifo,6)==0.17) .AND. (fall_coeff(ifo,7)==51.)) THEN      
       fcent = dEXP(-temperature/fall_coeff(ifo,5)) + fall_coeff(ifo,6) &
    &        + dEXP(-fall_coeff(ifo,7)/temperature)                                  
       factor2 = 1./(1.+(dLOG10(kratio)/(0.75-1.27*dLOG10(fcent)))**2.) 
       qfor = (xkom/(1.+kratio)) * (fcent**factor2) 
     ELSE
       print*, fall_coeff(ifo,:)
       print*, "MCk: unexpected setup of the auxiliary fall_coeff in falloff rxn",ire,ifo
       STOP 
     ENDIF 

    kinetic_rate(ire) = qfor * fall_coeff(ifo,8) ! ADD A RATIO
    
    !qln = log(qfor) 
                                                                        
END SUBROUTINE ssh_gck_forate                                 

subroutine ssh_gck_extrarate(ire, iex, label)

!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This routine computes kinetic rates for FALL OFF reactions.
!     Extract from GECKO BOX MODEL    
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     ire: index in reaction list
!     iex: index in extra reaction list
!     fall_coeff: 
!     extra_coeff: 
!     SumMc
!     YlH2O
!
!
!     -- OUTPUT VARIABLES
!
!     qfor: kinetic_rate.
!
!------------------------------------------------------------------------

!=======================================================================
! Purpose: Compute the rate coefficient for the reactions using the
! keyword EXTRA in the mechanism, i.e.
!      A + B + EXTRA => X +Y
!              EXTRA /label data1 data2 .../
! where data_x are used to compute these "special" reaction rate coef.
!=======================================================================

    IMPLICIT NONE 

    double precision :: kd, xk1, xk2, xk3, qfor
    double precision :: water_dimer

    ! INPUT                                                                 
    INTEGER, INTENT(in) :: ire, iex, label

! CASE SELECTOR FOR LABEL
! -----------------------------
!  label=INT(extra_coeff(iex,1))
  qfor = kinetic_rate(ire) ! init
    
  SELECT CASE(label)

!**********************************************************
!* identificateur=100 --> 0 + O2 -> O3                    *
!* qln en entree est la bonne valeur de k pour la         *
!* recation d'ordre 3. A multiplier par M et O2 (0.2*M)   *
!**********************************************************
! label 100: O+O2+M=>O3 // 3rd order reaction 
    CASE (100)
      qfor = qfor * SumMc * SumMc * 2d-1
 
!*********************************************************              
! identificateur=500 --> reaction avec H2O               *              
! identificateur=501 --> 2 reaction :                    *              
!                          1 avec H2O (qfor)             *              
!                          2 avec H2O et M (q2)          *              
!            les parametres de 2 sont dans extra_coeff   *              
!*********************************************************              
! label 500: reaction with water vapor  
    CASE (500)
      qfor = qfor * YlH2O 

! label 501: reaction with H20 (xk1) and H2O+M (xk2)  // specific to HO2+HO2 => H2O2
    CASE (501) 
      xk1 = qfor * YlH2O 
      xk2 = extra_coeff(iex,2)*dEXP(-extra_coeff(iex,4)/temperature)*  &
            YlH2O * SumMc 
      qfor = xk1 + xk2

!*********************************************************              
! identificateur=502 --> reaction avec H2O dimer          *    
!*********************************************************                
! label 502: reaction with water dimer (Kd according to Scribano et al., 2006)
    CASE (502)
      ! calculate water monomers vs dimers concentration                      
      ! dimerization constant [atm-1] according to Scribano et al., 2006      
      ! atm^(-1)          
      kd = 4.7856E-4*dEXP(1851/temperature-5.10485E-3*temperature) 
      !convert atm-1 to molec-1 cm3                                           
      kd = kd*8.314*temperature*1E6/(1.01325E5*6.02E23)
      water_dimer = kd * (YlH2O**2.) 
      qfor = qfor * water_dimer

                                                                    
!*********************************************************              
! identificateur=550 --> reaction OH + HNO3              *              
! k=k0 + k3M/(1+K3M/K2)                                  *              
!    k0 : donnees par les parametres de Arrhenius (qfor)  *              
!    k2 : parametres 2,3,4 de extra_coeff                    *              
!    k3 : parametres 5,6,7 de extra_coeff                    *              
!*********************************************************              
! label 550: OH+HNO3 reaction: k=k0 + k3*M/(1+K3*M/K2)
    CASE (550)

      xk2 = extra_coeff(iex,2)*(temperature**extra_coeff(iex,3))* &
 &       dEXP(-extra_coeff(iex,4)/temperature)
      xk3 = extra_coeff(iex,5)*(temperature**extra_coeff(iex,6))* &
 &       dEXP(-extra_coeff(iex,7)/temperature)
      qfor = qfor + xk3*SumMc / (1.+ ((xk3*SumMc)/xk2) )

! Add: ISOM for isomerisation reaction
    CASE (200) 
      qfor = qfor*(extra_coeff(iex,2)*(temperature**4) +  &
 &      extra_coeff(iex,3) * (temperature**3) + &
 &      extra_coeff(iex,4) * (temperature**2) + &
 &      extra_coeff(iex,5) * temperature + extra_coeff(iex,6))

!********************************************************               
! identificateur inconnu                                *               
!********************************************************               
! unidentified label
    CASE DEFAULT
      print*, '--error-- in akkextra. Type unknown: ', &
                label,ire,iex
      STOP 
  END SELECT
  
  ! add a ratio
  kinetic_rate(ire) = qfor * extra_coeff(iex,8)
  
end subroutine ssh_gck_extrarate

subroutine ssh_mcm_rate(ire,iex,label,azi)

!C------------------------------------------------------------------------
!C
!C     -- DESCRIPTION: all types of MCM rate
!C
!C     with Arrihenis(ire,1) the ratio
!C
!C  Extracted from mcm_3-3-1_unix/mcm_3-3-1_fortran_complete.txt
!C
!C     WZZ 11/06/2023
!C
!C------------------------------------------------------------------------

    IMPLICIT NONE 

    ! INPUT                                                                 
    INTEGER, INTENT(in) :: ire, iex, label
    double precision, intent(in) :: azi
    
    integer :: ind
    double precision :: ka, kb, kd, qfor

  SELECT CASE(label)

    CASE (91) ! Photolysis Rates
!C------------------------------------------------------------------------
!C
!C  IN THE FORMAT J = l* cosX ** m *dEXP(-n*secX)
!C
!C------------------------------------------------------------------------
        ka = 0.d0 !init cosX
        kb = 0.d0 !init secsX

        if (azi .lt. 9d1) then ! with photolysis
          ka = max(0d0,dcos(azi/1.8D2* 3.14159265358979323846D0)) ! cosX
          if (ka.gt.0d0) kb=1.0d+0/(ka+1.0D-30) !secX
        endif
      
        if (ka .lt. 1.d-10) then
            qfor = 1.0d-30
        else
            qfor = extra_coeff(iex,2)* ka**(extra_coeff(iex,3)) &
                           * dEXP(-extra_coeff(iex,4)*kb)
        endif
        
    CASE (92) ! Generic Rate Coefficients
!C------------------------------------------------------------------------
!C
!C  https://mcm.york.ac.uk/MCM/rates/generic
!C
!C------------------------------------------------------------------------
        ind = INT(extra_coeff(iex,2))
        SELECT CASE(ind)
          CASE (1) ! KRO2NO
            qfor = 2.7D-12*dEXP(360/temperature)
          CASE (2) ! KRO2HO2
            qfor = 2.91D-13*dEXP(1300/temperature)
          CASE (3) ! KAPHO2
            qfor = 5.2D-13*dEXP(980/temperature)
          CASE (4) ! KAPNO
            qfor = 7.5D-12*dEXP(290/temperature)
          CASE (5) ! KRO2NO3
            qfor = 2.3D-12
          CASE (6) ! KNO3AL
            qfor = 1.44D-12*dEXP(-1862/temperature)
          CASE (7) ! KDEC
            qfor = 1.00D+06
          CASE (8) ! KROPRIM
            qfor = 2.50D-14*dEXP(-300/temperature)
          CASE (9) ! KROSEC
            qfor = 2.50D-14*dEXP(-300/temperature)
          CASE (10) ! KCH3O2
            qfor = 1.03D-13*dEXP(365/temperature)
          CASE (11) ! K298CH3O2
            qfor = 3.5D-13
          CASE (12) ! K14ISOM1
            qfor = 3.00D7*dEXP(-5300/temperature)
        CASE DEFAULT
          print*, 'Error: MCM Generic Rate Coefficients id non known: ', &
                    label,ind,ire,iex
          STOP
        END SELECT
        
    CASE (93) ! Complex Rate Coefficients
!C------------------------------------------------------------------------
!C
!C  https://mcm.york.ac.uk/MCM/rates/complex
!C
!C------------------------------------------------------------------------

        ind = INT(extra_coeff(iex,2))
        SELECT CASE(ind)

        CASE(1)  ! KMT01
            ka = 1.0D-31*SumMc*(temperature/3d2)**(-1.6d0)
            kb = 5.0D-11*(temperature/3d2)**(-0.3)
            kd = 1d1**(dLOG10(0.85d0)/(1.0D0+(dLOG10(ka/kb) &
                    /(0.75D0-1.27D0*dLOG10(0.85d0)))**2D0))
            qfor = (ka*kb)*kd/(ka+kb)
        CASE(2)  ! KMT02
            ka = 1.3D-31*SumMc*(temperature/3d2)**(-1.5d0)
            kb = 2.3D-11*(temperature/3d2)**0.24
            kd = 1d1**(dLOG10(0.6d0)/(1.0D0+(dLOG10(ka/kb) &
                    /(0.75D0-1.27D0*dLOG10(0.6d0)))**2D0))
            qfor = (ka*kb)*kd/(ka+kb)
        CASE(3)  ! KMT03
            ka = 3.6D-30*SumMc*(temperature/3d2)**(-4.1d0)
            kb = 1.9D-12*(temperature/3d2)**0.2
            kd = 1d1**(dLOG10(0.35d0)/(1.0D0+(dLOG10(ka/kb) &
                    /(0.75D0-1.27D0*dLOG10(0.35d0)))**2D0))
            qfor = (ka*kb)*kd/(ka+kb)
        CASE(4)  ! KMT04
            ka = 1.3D-3*SumMc*(temperature/3d2)**(-3.5d0) &
                    *dEXP(-11d3/temperature)
            kb = 9.7D+14*(temperature/3d2)**0.1 &
                    *dEXP(-11080/temperature)
            kd = 1d1**(dLOG10(0.35d0)/(1.0D0+(dLOG10(ka/kb) &
                    /(0.75D0-1.27D0*dLOG10(0.35d0)))**2D0))
            qfor = (ka*kb)*kd/(ka+kb)
        CASE(5)  ! KMT05
            qfor = 1.44D-13 * (1+(SumMc/4.2D+19))
        CASE(6)  ! KMT06
            qfor = 1 + (1.40D-21*dEXP(2200/temperature)*YlH2O)
        CASE(7)  ! KMT07
            ka = 7.4D-31*SumMc*(temperature/3d2)**(-2.4d0)
            kb = 3.3D-11*(temperature/3d2)**(-0.3d0)
            kd = 1d1**(dLOG10(0.81d0)/(1.0D0+(dLOG10(ka/kb) &
                    /(0.75D0-1.27D0*dLOG10(0.81d0)))**2D0))
            qfor = (ka*kb)*kd/(ka+kb)
        CASE(8)  ! KMT08
            ka = 3.2D-30*SumMc*(temperature/3d2)**(-4.5d0)
            kb = 3.0D-11
            kd = 1d1**(dLOG10(0.41d0)/(1+(dLOG10(ka/kb) &
                    /(0.75-1.27*(dLOG10(0.41d0))))**2))
            qfor = (ka*kb)*kd/(ka+kb)
        CASE(9)  ! KMT09
            ka = 1.4D-31*SumMc*(temperature/3d2)**(-3.1d0)
            kb = 4.0D-12
            kd = 1d1**(dLOG10(0.4d0)/(1+(dLOG10(ka/kb) &
                    /(0.75-1.27*(dLOG10(0.4d0))))**2))
            qfor = (ka*kb)*kd/(ka+kb)
        CASE(10)  ! KMT10
            ka = 4.10D-05*SumMc*dEXP(-10650/temperature)
            kb = 6.0D+15*dEXP(-11170/temperature)
            kd = 1d1**(dLOG10(0.4d0)/(1+(dLOG10(ka/kb) &
                    /(0.75-1.27*(dLOG10(0.4d0))))**2))
            qfor = (ka*kb)*kd/(ka+kb)
        CASE(11)  ! KMT11
            qfor = 2.40D-14*dEXP(460/temperature) + &
                    (6.50D-34*dEXP(1335/temperature)*SumMc) &
                    /(1+(6.50D-34*dEXP(1335/temperature)*SumMc) &
                    /(2.70D-17*dEXP(2199/temperature)))
        CASE(12)  ! KMT12
            ka = 2.5D-31*SumMc*(temperature/3d2)**(-2.6d0)
            kb = 2.0D-12
            kd = 1d1**(dLOG10(0.53d0)/(1+(dLOG10(ka/kb) &
                    /(0.75-1.27*(dLOG10(0.53d0))))**2))
            qfor = (ka*kb)*kd/(ka+kb)
        CASE(13)  ! KMT13
            ka = 2.5D-30*SumMc*(temperature/3d2)**(-5.5d0)
            kb = 1.8D-11
            kd = 1d1**(dLOG10(0.36d0)/(1+(dLOG10(ka/kb) &
                    /(0.75-1.27*(dLOG10(0.36d0))))**2))
            qfor = (ka*kb)*kd/(ka+kb)
        CASE(14)  ! KMT14
            ka = 9.0D-5*dEXP(-9690/temperature)*SumMc
            kb = 1.1D+16*dEXP(-10560/temperature)
            kd = 1d1**(dLOG10(0.36d0)/(1+(dLOG10(ka/kb) &
                    /(0.75-1.27*(dLOG10(0.36d0))))**2))
            qfor = (ka*kb)*kd/(ka+kb)
        CASE(15)  ! KMT15
            ka = 8.6D-29*SumMc*(temperature/3d2)**(-3.1d0)
            kb = 9.0D-12*(temperature/3d2)**(-0.85d0)
            kd = 1d1**(dLOG10(0.48d0)/(1+(dLOG10(ka/kb) &
                    /(0.75-1.27*(dLOG10(0.48d0))))**2))
            qfor = (ka*kb)*kd/(ka+kb)
        CASE(16)  ! KMT16
            ka = 8D-27*SumMc*(temperature/3d2)**(-3.5d0)
            kb = 3.0D-11*(temperature/3d2)**(-1d0)
            kd = 1d1**(dLOG10(0.5d0)/(1+(dLOG10(ka/kb) &
                    /(0.75-1.27*(dLOG10(0.5d0))))**2))
            qfor = (ka*kb)*kd/(ka+kb)
        CASE(17)  ! KMT17
            ka = 5.0D-30*SumMc*(temperature/3d2)**(-1.5d0)
            kb = 1.0D-12
            kd = 1d1**(dLOG10(0.17*dEXP(-51/temperature) &
                +dEXP(-temperature/204))/(1+(dLOG10(ka/kb) &
                /(0.75-1.27*(dLOG10(0.17*dEXP(-51/temperature) &
                +dEXP(-temperature/204))))**2)))
            qfor = (ka*kb)*kd/(ka+kb)
        CASE(18)  ! KMT18
            ka = 9.5D-39*SumMc*2d-1*dEXP(5270/temperature) !O2
            kb = 1+7.5D-29*SumMc*2d-1*dEXP(5610/temperature) !O2
            qfor = ka/kb

    !C------------------------------------------------------------------------
    !C  Formation and decomposition of PAN
    !C  IUPAC 2001
    !C  Extracted from mcm_3-3-1_unix/mcm_3-3-1_fortran_complete.txt
    !C------------------------------------------------------------------------
        CASE (21) ! KFPAN 
            ka = 3.28D-28*SumMc*(temperature/3d2)**(-6.87d0)
            kb = 1.125D-11*(temperature/3d2)**(-1.105d0)
            kd = 10**(dLOG10(0.30d0)/(1+(dLOG10(ka/kb) &
                    /(0.75-1.27*(dLOG10(0.30d0))))**2))
            qfor = (ka*kb)*kd/(ka+kb)
        CASE (22) ! KBPAN
            ka = 1.10D-05*SumMc*dEXP(-101d2/temperature)
            kb = 1.90D17*dEXP(-141d2/temperature)
            kd = 10**(dLOG10(0.30d0)/(1+(dLOG10(ka/kb) &
                    /(0.75-1.27*(dLOG10(0.30d0))))**2))
            qfor = (ka*kb)*kd/(ka+kb)
        CASE (23) ! KBPPN
            ka = 1.7D-03*dEXP(-1128d1/temperature)*SumMc
            kb = 8.3D+16*dEXP(-1394d1/temperature)
            kd = 10**(dLOG10(0.36d0)/(1+(dLOG10(ka/kb) &
                    /(0.75-1.27*(dLOG10(0.36d0))))**2))
            qfor = (ka*kb)*kd/(ka+kb)
        
        !CASE (30) ! KRO2 => move to 94
        !    qfor = 1.26D-12 * RO2
            
        CASE DEFAULT
          print*, 'Error: MCM Complex Rate Coefficients id non known: ', &
                    label,ind,ire,iex
          STOP
        END SELECT

    CASE DEFAULT
      print*, '--error-- in ssh_mcm_rate. Type unknown: ', &
                label,ire,iex
      STOP 

  END SELECT
  
  ! add a ratio
  kinetic_rate(ire) = qfor *  Arrhenius(ire,1)
  
end subroutine ssh_mcm_rate


subroutine ssh_genoa_spec(ire, iex)

!C------------------------------------------------------------------------
!C
!C     -- DESCRIPTION: Additional kientic rate compuration used in GENOA
!C     or GENOA related mechanism
!C
!C     with Arrihenis(ire,1) the ratio
!C
!C
!C     WZZ 01/12/2024
!C
!C------------------------------------------------------------------------

    IMPLICIT NONE 

    ! INPUT                                                                 
    INTEGER, INTENT(in) :: ire, iex
    integer :: label
    double precision :: ka, kb, qfor

  ! Get label from extra_coeff
  label = INT(extra_coeff(iex, 2)) ! 2nd level label
  SELECT CASE(label)

    CASE (1) ! 3.8D-13*dexp(780/TEMP)*(1-1/(1+498*dexp(-1160/TEMP))) - MCM
        qfor = 3.8D-13*dexp(78d1/temperature)*(1d0-1d0/(1d0+498d0* &
              dexp(-116d1/temperature)))
    CASE (2) ! 3.8D-13*dexp(780/TEMP)*(1/(1+498*dexp(-1160/TEMP)))
        qfor = 3.8D-13*dexp(78d1/temperature)*(1d0/(1d0+498d0* &
              dexp(-116d1/temperature)))
    CASE (3) ! 1.03D-13*math.dexp(365/TEMP)*(1-7.18*dexp(-885/TEMP))
        qfor = 1.03D-13*dexp(365d0/temperature)*(1d0-7.18d0* &
              dexp(-885d0/temperature))
    CASE (4) ! 8.8D-12*dexp(-1320/TEMP) + 1.7D-14*dexp(423/TEMP)
        qfor = 8.8D-12*dexp(-132d1/temperature) + 1.7D-14*dexp(423d0/temperature)
    CASE (5) ! 5.00E-12*O2*3.2*(1-dexp(-550/TEMP)) - need to go with TB O2
        qfor = 5.00E-12*3.2d0*(1d0-dexp(-550d0/temperature))
    CASE (6) ! 1.03D-13*0.5*math.dexp(365/TEMP)*(1-7.18*dexp(-885/TEMP))
        qfor = 1.03D-13*0.5d0*dexp(365d0/temperature)*(1d0-7.18d0* &
              dexp(-885d0/temperature))
    CASE (7) ! 2.20E+10*EXP(-8174/TEMP)*EXP(1.00E+8/TEMP@3)
        qfor = 2.20E+10*dexp(-8174d0/temperature)*dexp(1d0+8d0/ &
              temperature**3)
    CASE (8) ! 8.14E+9*EXP(-8591/TEMP)*EXP(1.00E+8/TEMP@3)
        qfor = 8.14E+9*dexp(-8591d0/temperature)*dexp(1d0+8d0/ &
              temperature**3)

    CASE (9) ! rewrite from ssh_IRDICARB 
!C     Write kinetic rate corresponding to irreversible pathway for MGLY condensation.
!C     depending on RH     
!C     See SI of Lannuque et al. (2023) for more informations
!C     Author: VICTOR LANNUQUE, 2022.
        ka = 6.112d2 * dexp(1.767d1 * (temperature - 273.15d0) &
             /(temperature - 29.65d0))
        kb = humidity * pressure / ((0.62197d0 * (1.d0-humidity)+humidity) * ka ) * 100.
        qfor = 1d-9*kb**3 - 1d-7*kb**2 + 3.0d-7*kb + 0.0003
    CASE DEFAULT
          print*, '--error-- in ssh_genoa_spec. Type unknown: ', &
                    label, ire, iex
          STOP
    END SELECT

  ! Add a ratio
  kinetic_rate(ire) = qfor *  Arrhenius(ire, 1)

end subroutine ssh_genoa_spec

subroutine ssh_spack_spec(ire, iex, label)

!C------------------------------------------------------------------------
!C
!C     -- DESCRIPTION: define the reaction rates when
!C                     the keyword EXTRA is used in the reactions list.
!C
!C------------------------------------------------------------------------
  
    IMPLICIT NONE 

    ! INPUT                                                                 
    INTEGER, INTENT(in) :: ire, iex, label
    
    integer :: ind
    double precision :: ka, kb, kd, qfor
    double precision :: psat, facteur, xlw
    double precision :: masmol, cstar, aw, cbar, awc, denom, kwon !! Wall

    if (humidity == -999.d0) then
       call ssh_compute_sh(relative_humidity, temperature, &
            pressure, humidity)
    endif
    
    xlw = humidity

    ! Label for different mechanisms
    ! 10: cb05, 20: racm2, 30: melchior2 ... 
    ind = INT(extra_coeff(iex, 2)) ! 2nd level label
    
    if (label .eq. 10) then ! CB05 - rewrite from ssh_WSPEC_CB0590
      SELECT CASE(ind)
        CASE (1)
           qfor = SumMc * 6.0d-34 * (temperature/3.d2) ** (-2.4d0)
           qfor = qfor * SumMc * 0.2d0
!C     YS 17/11/2008 values given by NASA/JPL 2003
        CASE (2)
            qfor = 2.3d-13 * dexp(600.0d0 / temperature) &
                + 1.7d-33* SumMc * dexp(1000.0d0 / temperature)
        CASE (3)
            qfor = 3.22d-34 * dexp(2800.0d0 / temperature) &
                 + 2.38d-54 * SumMc * dexp(3200.0d0 / temperature)
            qfor = qfor * YlH2O
                 
        CASE (4)
            ka = 2.4d-14 * dexp(460.0d0 / temperature)
            kb = 2.7d-17 * dexp(2199.0d0 / temperature)
            kd = 6.5d-34 * dexp(1335.0d0 / temperature) * SumMc
            qfor = ka + kd / (1d0 + kd / kb)
!C     YS 17/11/2008 values given by NASA/JPL 2003
        CASE (5)
            qfor = 1.44d-13 * (1.0d0 + 2.381d-20 * 8.0d-1 * SumMc)
!C     YS 19/11/2008 value given by IUPAC 2005
        CASE (6)
           ! NO2 + OH --> HNO3
           ka = 2.0d-30 * (temperature / 3.d2)**(-3.d0)
           kb = 2.5d-11 * (temperature / 3.d2)**(0.d0)
           qfor = (ka * SumMc / (1.0d0 + ka * SumMc / &
                kb)) * 0.6d0 ** (1.0d0 / (1.0d0 + &
                (dlog10(ka * SumMc / kb))**2))
           qfor = qfor * 0.885d0
        CASE (7)
            qfor = 1.8d-39 * YlH2O * YlH2O
!C     YS 26/11/2008 value given by IUPAC 2005
        CASE DEFAULT
          print*, '--error-- in ssh_spack_spec. CB05 Type unknown: ', &
                    label,ire,iex
          STOP
      END SELECT
    elseif (label .eq. 20) then ! RACM2 - rewrite from ssh_WSPEC_RACM90
     SELECT CASE(ind)
        CASE (1)
           qfor = sumMc * 5.74d-34 * (temperature/3.d2) ** (-2.6d0)
           qfor = qfor * SumMc * 0.2d0
!C     BS 05/02/2003 values given by RACM
        CASE (2)
          qfor = 2.2d-13 * dexp(6d2 / temperature) &
               + 1.9d-33* SumMc * dexp(980.0d0/ temperature)
        CASE (3)
          qfor = 3.08d-34 * dexp(28d2 / temperature) &
               + 2.59d-54 * SumMc * dexp(3180.0d0 / temperature)
          qfor = qfor * YlH2O
!C      MODIF BS 06/06/2003 on the basis of CMAQ
        CASE (4)
          ka = 2.4d-14 * dexp(460d0 / temperature)
          kb = 2.7d-17 * dexp(2199d0 / temperature)
          kd = 6.5d-34 * dexp(1335d0 / temperature) * SumMc
          qfor = ka + kd / (1d0 + kd / kb)
        CASE (5)
          ! CO +   HO        ->      HO2
          qfor = 1.44d-13 * (1.0d0 + 8.0d-1 * SumMc / 4.0d19)
        CASE (6)
          ! NO + HO2 --> HNO3
          ka = 3.43d-12 * dexp(270.0d0 / temperature)
          kb = (530.0d0 / temperature) + 4.8d-6 * pressure - 1.73
          qfor = ka * kb / 100.d0
        CASE (7)
          qfor = 2.0d-39 * YlH2O * YlH2O
        CASE (8)
          ! ACT +  HO        ->      ACTP
          qfor = 1.39d-13 + 3.72d-11 * dexp(-2.044d3 / temperature)
        CASE (9)
          ! N2O5            ->      NO2   +      NO3
           qfor = dexp(  0.6117554523749536D+02 &
                - (  0.1100000000000000D+05 )/temperature)
           ka =  0.2200000000000000D-29* (temperature / 3.d2) &
                **(- ( 0.4400000000000000D+01))
           kb =  0.1400000000000000D-11* (temperature / 3.d2) &
                **(- ( 0.7000000000000000D+00))
           kd = (ka * SumMc / ( 1.0d0 + ka * SumMc / &
                kb)) * 0.6d0 ** (1.0d0 / (1.0d0 + &
                (dlog10(ka * SumMc / kb))**2.0d0))
           qfor = kd * qfor
        CASE (10)
           ! HNO4            ->      HO2    +     NO2
           qfor =  dexp(  0.6142746008608852D+02 &
                - (  0.1090000000000000D+05 )/temperature)
           ka =  0.2000000000000000D-30* (temperature / 3.d2) &
                **(- ( 0.3400000000000000D+01))
           kb =  0.2900000000000000D-11* (temperature / 3.d2) &
                **(- ( 0.1100000000000000D+01))
           kd = (ka * SumMc / ( 1.0d0 + ka * SumMc / &
                kb)) * 0.6d0 ** (1.0d0 / (1.0d0 + &
                (dlog10(ka * SumMc / kb))**2))
           qfor = kd * qfor
        CASE (11)
           ! PAN             ->      ACO3   +     NO2
           qfor =  dexp(  0.6462080260895155D+02 &
                - (  0.1395400000000000D+05 )/temperature)
           ka =  0.9700000000000000D-28* (temperature / 3.d2) &
                **(- ( 0.5600000000000000D+01))
           kb =  0.9300000000000000D-11* (temperature / 3.d2) &
                **(- ( 0.1500000000000000D+01))
           kd = (ka * SumMc / ( 1.0d0 + ka * SumMc / &
                kb)) * 0.6d0 ** (1.0d0 / (1.0d0 + &
                (dlog10(ka * SumMc / kb))**2))
           qfor = kd * qfor
        CASE (12)
           ! DK3000 -> IRDK3000
           Psat = 611.2d0 * dexp(17.67d0 * (temperature - 273.15d0) &
                /(temperature - 29.65d0))
           facteur = xlw*Pressure / ((0.62197d0 * (1.d0-xlw)+xlw) &
                * Psat ) * 100.
           qfor = 1d-9*facteur**3 - 1d-7*facteur**2 &
                + 3.0d-7*facteur + 0.0003 
        CASE (13)           
           ! TOL    -> WTOL
           Psat =  0.16000D+02/760.*dexp( 0.40000D+02*1000. &
                / 8.314*(1./298.-1./temperature))
           Masmol =  0.92000D+02 
           cstar = Masmol*1d+6*Psat/(8.205d-5*temperature)
           aw = 10**(-0.27440D+01)*cstar**(-0.14066D+01)
           cbar = (8000*8.314*temperature/(3.14*Masmol))**0.5
           awc = (aw*cbar/4.)
           denom = awc / (( 0.42800D-02* 0.50000D-05)**0.5)
           kwon =  0.33330D+02 * awc / (1.+1.5708*denom)
           qfor = kwon
        CASE (14)
           ! NO2 + OH --> HNO3
           ka = 1.8d-30 * (temperature / 3.d2)**(-3.d0)
           kb = 2.8d-11 * (temperature / 3.d2)**(0.d0)
           qfor = (ka * SumMc / (1.0d0 + ka * SumMc / &
                kb)) * 0.6d0 ** (1.0d0 / (1.0d0 + &
                (dlog10(ka * SumMc / kb))**2))
           qfor = qfor * 0.868d0
           
        CASE DEFAULT
          print*, '--error-- in ssh_spack_spec. RACM2 Type unknown: ', &
                    label,ire,iex
          STOP
       END SELECT
    elseif (label .eq. 30) then ! MELCHIOR2 - rewrite from ssh_WSPEC_MELCHIOR2
       SELECT CASE(ind)
       CASE (1)
          ! NO2 + OH --> HNO3
          ka = 3.4d-30 * (3.d2 / temperature)**(3.d2) * SumMc
          kb = ka / (4.77d-11 * (3.d2 / temperature)**1.4)
          qfor = (ka / (1.0d0 + kb)) * 0.3 ** &
               (1.0d0 / (1.0d0 + &
               ((dlog10(kb) - 0.12) / 1.2)**2))
       CASE (2)
          ! N2O5 -> 2. HNO3
          qfor = 2.0d-39 * YlH2O * YlH2O
       CASE (3)
          ! NO2 -> HONO + NO2
          ka  = 15.0d0
          kb = 0.00033d0
          qfor = 0.5d0 * kb / ka
          
       END SELECT
    end if
  
  ! add a ratio
  kinetic_rate(ire) = qfor *  Arrhenius(ire,1)
   
end subroutine ssh_spack_spec


subroutine ssh_hetero_rate(i, reaction_ind, hetero_ind, &
     ns, nbin_aer, temp, press, &
     wetdiam, granulo, &
     dsf_aero, wmol, lwcmin, ktot)
  
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This routine computes the kinetic rates for the heterogeous reactions.
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     wmol: molecular weight
!     ICLD: 0: no, 1: vsrm, 2: simple
!     LWCmin: lwc_cloud_threshold
!     temp: temperature
!     press: pressure
!     dllwc: liquid_water_content in kg/kg (QCLOUD from WRF)
!     dsf_aero: diameter
!     WetDiam: wet diameter      
!
!     -- INPUT/OUTPUT VARIABLES
!
!     -- OUTPUT VARIABLES
!
!     ktot: heterogeneous reaction rate
!
!------------------------------------------------------------------------

      INCLUDE 'CONST.INC'
      INCLUDE 'CONST_A.INC'
      INCLUDE 'hetrxn.inc'

      INTEGER Ns,Nbin_aer

      integer :: reaction_ind, hetero_ind
      INTEGER i,js
      DOUBLE PRECISION Ndroplet

      DOUBLE PRECISION temp,press
      DOUBLE PRECISION vi
      DOUBLE PRECISION ktot,A
      DOUBLE PRECISION DIFF
      DOUBLE PRECISION avdiammeter

      DOUBLE PRECISION WetDiam(Nbin_aer)
      DOUBLE PRECISION WetDiammeter(Nbin_aer)
      DOUBLE PRECISION granulo(Nbin_aer)
      DOUBLE PRECISION dsf_aero(Nbin_aer)

      DOUBLE PRECISION lwctmp, dllwc 

      INTEGER ICLOUD
      INTEGER ICLD

      DOUBLE PRECISION sigm_tmp,parm_tmp,wmol_tmp
      DOUBLE PRECISION Wmol(Ns),LWCmin

      double precision dactiv, avdiam

      DOUBLE PRECISION :: lwca
      
     
      !     WET DIAMETERS OF THE TWO SECTIONS
      
!     If bulk, diameter = avdiam
      parameter (avdiam = 20)   ! in \mu m   
!     
!     ACTIVATION DIAMETER (Dry)
!     
      parameter (dactiv = 0.2D0) ! in \mu m

      if (i > 4) then
         write(*,*) "Index", i, "is not accepted for the heterogenous reactions."
         stop
      endif

      
!     Constants.
      avdiammeter = 1.d-6 * avdiam
      
      sigm_tmp = SIGM_NO2
      parm_tmp = PARM_NO2
      ICLOUD=0
      
      !C     Cloud liquid water content (g/m^3)
      !  Converts from kg / kg to g / m^3.
      lwctmp=DLLWC*1000.d0* &
           press/101325.d0*28.97d0/Pr/ &
           temp
      
      IF (ICLD.GE.1.AND.(lwctmp.GE.LWCmin)) THEN
         ICLOUD=1
      ENDIF      
      
      !     REACTION PROBABILITIES
      ktot = 0.d0

      wmol_tmp = Wmol(hetero_ind + 1)
      call SSH_COMPUTE_GAS_DIFFUSIVITY(TEMP,PRESS, &
           sigm_tmp, wmol_tmp, parm_tmp, DIFF)

      call SSH_COMPUTE_QUADRATIC_MEAN_VELOCITY(TEMP, &
           wmol_tmp,vi)

      if (gamma(i) .eq. 0.d0) then
         ktot = 0.d0
      else
         do js=1,Nbin_aer
            IF ((ICLOUD.NE.1).OR.(dsf_aero(js).lt.dactiv*1.d6)) then
               WetDiammeter(js) = 1.d-6 * WetDiam(js)

               !     CALCULATE k FOR EACH SECTION
               
               A = PI*(WetDiammeter(js)**2.d0)*granulo(js) !surface * nb_aero

               !     CALCULATE k FOR AEROSOL DISTRIBUTION

               IF(Gamma(i).ne.0.d0) then
                  ktot = ktot &
                       + (((WetDiammeter(js) / (2.d0 * DIFF)) &
                       + 4.d0 / (vi * Gamma(i)))**(-1.d0)) * A ! rate constant for each rxn
               ENDIF
            ENDIF
         enddo
         
         !     If we are in a cloud, then reactions on droplet surface for N2O5.
         IF ((ICLOUD.EQ.1).and.(i.eq.4)) THEN
            Ndroplet = lwctmp &
                 / (RHOwater * 1.d3 & !kg.m-3 -> g.m-3
                 * pi / 6.d0 * avdiammeter**3.d0)
            A = PI * avdiammeter**2.d0 * Ndroplet !surface * nb_droplet
            IF(Gamma(i).ne.0.d0) then
               ktot = ktot + ((avdiammeter / (2.d0 * DIFF) + &
                    4.d0 / (vi*Gamma(i)))**(-1.d0))*A ! rate constant for each rxn
            ENDIF
         ENDIF
      endif

      
end subroutine ssh_hetero_rate






subroutine ssh_dicarb_rate(i, reaction_ind, irdi_ind, &
!                           nbin_aer, temp, press, icld, &
                           nbin_aer, temp, press, &
                           wetdiam, granulo, &
!                           dsf_aero, lwcmin, dllwc, OHgasmass,ktot)
                           OHgasmass,ktot)
  
!------------------------------------------------------------------------
! 
!     -- DESCRIPTION
! 
!     This routine computes the kinetic rates for the heterogeous reactions
!     leading to irreversible condensation of dicarbonyles.
! 
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     wmol: molecular weight
!     ICLD: 0: no, 1: vsrm, 2: simple
!     LWCmin: lwc_cloud_threshold
!     temp: temperature
!     press: pressure
!     dllwc: liquid_water_content in kg/kg (QCLOUD from WRF)
!     dsf_aero: diameter
!     WetDiam: wet diameter      
!
!     -- INPUT/OUTPUT VARIABLES
!
!     -- OUTPUT VARIABLES
!
!     ktot: heterogeneous reaction rate
!
!------------------------------------------------------------------------

      INCLUDE 'CONST.INC'
      INCLUDE 'CONST_A.INC'

      INTEGER Nbin_aer

      integer :: reaction_ind, irdi_ind
      INTEGER i,js
!      DOUBLE PRECISION Ndroplet

      DOUBLE PRECISION temp,press
      DOUBLE PRECISION vi
      DOUBLE PRECISION ktot
!      DOUBLE PRECISION avdiammeter

      DOUBLE PRECISION WetDiam(Nbin_aer)
      DOUBLE PRECISION WetDiammeter(Nbin_aer)
      DOUBLE PRECISION granulo(Nbin_aer)
!      DOUBLE PRECISION dsf_aero(Nbin_aer)

!      DOUBLE PRECISION lwctmp, dllwc 
!      INTEGER ICLOUD
!      INTEGER ICLD

      DOUBLE PRECISION wmol_tmp, OHgasmass
      
!      DOUBLE PRECISION LWCmin
!      double precision dactiv, avdiam
!      DOUBLE PRECISION :: lwca
      
      ! data for heterogeneous reactivity of dicarbonyls (GLY and MGLY)
      DOUBLE PRECISION     Rcst, Daq, alpha,ve,OHaq,OHgas,kaq,kl,HeffOH
      DOUBLE PRECISION     Rpart,Asurf,qp,cothqp,Heff,gammad,Navogadro
      DOUBLE PRECISION     Heff298,deltaHenry,Heff298OH,deltaHenryOH
      
     
      !     WET DIAMETERS OF THE TWO SECTIONS
      
!     If bulk, diameter = avdiam
!      parameter (avdiam = 20)   ! in micro m   
!     
!     ACTIVATION DIAMETER (Dry)   
!      parameter (dactiv = 0.2D0) ! in micro m

      if (irdi_ind > 2) then
         write(*,*) "Index", i, "is not accepted for the dicarb. heterogenous reactions."
         stop
      endif

      
!     Constants.
!      avdiammeter = 1.d-6 * avdiam
!      ICLOUD=0
!      
!      !C     Cloud liquid water content (g/m^3)
!      !  Converts from kg / kg to g / m^3.
!      lwctmp = DLLWC *1000.d0 *press /101325.d0 *28.97d0 /Pr /temp
!      
!      IF (ICLD.GE.1.AND.(lwctmp.GE.LWCmin)) THEN
!         ICLOUD=1
!      ENDIF

! PARAMETERS

Rcst = 8.207d-5        !universal gas constant in m3 atm mol-1 K-1
Navogadro = 6.02214d23 !mol-1
Daq = 1.0d-9           !aqueous-phase diffusion coefficient (m2 s−1)
!The aqueous-phase diffusion coefficient used for both GLY and MGLY was Daq = 10−9 m2 s−1. 
!Daq does not vary much for small species, and this value is typical for small organics (Bird et al., 2006). 
![Curry et al. 2018 10.5194/acp-18-9823-2018]
alpha = 0.02d0           !mass accommodation coefficient ()
!The mass accommodation coefficient used was α = 0.02. This value of α is an estimate based on the assumption 
!that α values for GLYX and MGLY are similar to that of formaldehyde uptake to water (Jayne et al., 1992).
![Curry et al. 2018 10.5194/acp-18-9823-2018]

!Henry constant of dicarb. adapted from (Sander 2015 10.5194/acp-15-4399-2015) to match with Hu et al. 2022 observations.
if (irdi_ind .eq. 1) then
  Heff298 = 360000.d0 ; deltaHenry = 4500.d0 ; wmol_tmp = 58.d0        !GLYOXAL
elseif (irdi_ind .eq. 2) then
  Heff298 = 42000.d0  ; deltaHenry = 4500.d0 ; wmol_tmp = 72.d0        !METHYLGLYOXAL
else
  print*, "error: dicarb. properties missing in ssh_dicarb_rate of ModuleChemkinetic.f90" ; stop
endif
Heff  =  Heff298 * dexp( deltaHenry * ((1.d0/temp)-(1.d0/298.d0)) )

!Henry constant of OH (Sander 2015 10.5194/acp-15-4399-2015)
Heff298OH = 29.d0
deltaHenryOH = 4300.d0 !4300.d0
HeffOH = Heff298OH*dexp(deltaHenryOH*((1.d0/temp)-(1.d0/298.d0)))

!gas-phase thermal velocity of glyoxal/methylglyoxal (m s−1)
!ve = dsqrt(8000.d0*Rcst*101300.d0*temp/(pi*wmol_tmp))! with Mw (g mol-1)
call ssh_compute_quadratic_mean_velocity(temp,wmol_tmp,ve)


OHgas = OHgasmass/1.d12/17.d0*Navogadro
!OHgas =  1.5d6 !typical mean value in atmosphere (lelieveld et al. 2016)

OHaq = OHgas * Rcst * temp * HeffOH / Navogadro * 1.0d6 !OH conc in aq phase according to Henry's law (mol l-1)
!OHaq=3.0d-12 !typical mean value in atmosphere (Hermann et al. 2010)
kaq = 1.1d9  !aqueous reaction rate of MGLYaq or GLYaq with OHaq (l mol-1 s-1)
kl = kaq * OHaq !first-order aqueous loss rate(s−1)

ktot=0.d0 !effective uptake rate (s-1)
!calculs en lien avec la taille des aerosols:
do js=1,Nbin_aer
!  if ((ICLOUD.NE.1).OR.(dsf_aero(js).lt.dactiv*1.d6)) then
    if (WetDiam(js).ne.0.0) then
      Rpart = 1.d-6 * WetDiam(js) / 2.d0         !particle radius (m)
      Asurf = 4.d0*pi*(Rpart**2.d0) * granulo(js)   !aerosol surface area density (m2 m−3): surface * nb_aero
      qp = Rpart / dsqrt(Daq/kl)                 !parameter for measuring in-particle diffusion limitations ()

      !gammad = uptake coefficient ()
      cothqp = cosh(qp)/sinh(qp)
      
      gammad = 1.d0/( 1.d0/alpha + (ve/(4.d0*Rcst*temp*Heff*dsqrt(kaq*Daq)))/(cothqp - 1.d0/qp) )
      
      ktot = ktot + (ve*gammad*Asurf)/4.d0 

    endif
!  endif
enddo
         
!if (ICLOUD.EQ.1) then !     If we are in a cloud.
!  Ndroplet = lwctmp / (RHOwater * 1.d3 * pi / 6.d0 * avdiammeter**3.d0)
!  Rpart = avdiammeter / 2.d0              !particle radius (m)
!  Asurf = 4.d0*pi*(Rpart**2.d0) * Ndroplet   !aerosol surface area density (m2 m−3): surface * nb_aero
!  qp = Rpart / dsqrt(Daq/kl)              !parameter for measuring in-particle diffusion limitations ()
!  
!  !gammad = uptake coefficient ()
!  cothqp = cosh(qp)/sinh(qp)
!  gammad = 1.d0/( 1.d0/alpha + (ve/(4.d0*Rcst*temp*Heff*dsqrt(kaq*Daq)))/(cothqp - 1.d0/qp) )
!  
!  ktot = ktot + (ve*gammad*Asurf)/4.d0 
! 
!endif


   
end subroutine ssh_dicarb_rate











subroutine ssh_compute_granulo(fixed_density_aer, IDENS, &
     nbin_aer, ns_aer, dlconc_aer, mass_density_aer, &
     DLnumconc_aer, ncomp_aer, bin_bound_aer, dsf)

!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This routine computes the aerosol diameters for the heterogeous reactions.
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     
!     fixed_density_aer: aerosol density which is set in namelist.
!     IDENS: 0 for fixed density, 1 for varing density
!     nbin_aer: number of aerosol bins
!     ns_aer: number of aerosol species
!     dlconc_aer: aerosol mass concentrations in ug/m3
!     mass_density_aer: density of each aerosol species 
!     DLnumconc_aer: aerosol number concentrations in #/m3
!     ncomp_aer: number of aerosol composition   
!     bin_bound_aer: bin bound of the aerosol size bin
!  
!     -- INPUT/OUTPUT VARIABLES
!
!     -- OUTPUT VARIABLES
!
!     dsf: diameter of aerosol bins
!
!------------------------------------------------------------------------

  
      IMPLICIT NONE

      INCLUDE 'paraero.inc'
      INCLUDE 'CONST_A.INC'
      
      integer :: jb, jsp, nbin_aer, ns_aer, ncomp_aer
      integer :: IDENS

      double precision :: rhoa, fixed_density_aer
      double precision :: conc_tot
      double precision :: rho_dry(nbin_aer)
      DOUBLE PRECISION :: DLconc_aer(nbin_aer,ns_aer)
      double precision :: DLnumconc_aer(nbin_aer)
      DOUBLE PRECISION :: MSF(nbin_aer)
      DOUBLE PRECISION :: DSF(nbin_aer)
      DOUBLE PRECISION :: DBF((nbin_aer/ncomp_aer)+1)
      DOUBLE PRECISION :: bin_bound_aer((nbin_aer/ncomp_aer) + 1)
      INTEGER :: idx_bs(nbin_aer)
      double precision :: mass_density_aer(ns_aer)
      
!C     Aerosol density converted in microg / microm^3.
      RHOA = fixed_density_aer * 1.D-09

!C     Aerosol discretization converted in microm.
      DO Jb=1,(nbin_aer/ncomp_aer)+1
         DBF(Jb) = bin_bound_aer(Jb) * 1.D06
      ENDDO

!C     relations between bin idx and size idx
      DO Jb=1,nbin_aer
         idx_bs(Jb)=(Jb-1)/ncomp_aer+1
      ENDDO
      
      
!C     Compute aerosol density
      rho_dry = RHOA
      IF (IDENS.EQ.1) THEN   ! for varying density
         DO Jb=1,nbin_aer
            CALL SSH_COMPUTE_DENSITY(nbin_aer,ns_aer, ns_aer, TINYM, &
                 DLconc_aer, &
                 mass_density_aer,Jb,rho_dry(Jb)) 
         ENDDO
      ENDIF
      DO Jb = 1, nbin_aer
         conc_tot = 0.d0
         DO Jsp = 1, Ns_aer
            conc_tot = conc_tot + DLconc_aer(Jb,Jsp)
         ENDDO

         !C     Compute mass and diameter of each section
         IF (DLnumconc_aer(Jb) .GT. 0.d0) THEN
            MSF(Jb) = conc_tot/DLnumconc_aer(Jb)
         ELSE
            MSF(Jb) = 0.d0
         ENDIF

         if ((DLnumconc_aer(Jb).GT. TINYN .or.  &
              conc_tot.GT.TINYM) &
              .AND. IDENS .EQ. 1) then
            DSF(Jb) = (MSF(Jb)/cst_PI6/rho_dry(Jb))**cst_FRAC3
         else
            DSF(Jb) = DSQRT(DBF(idx_bs(Jb))* DBF(idx_bs(Jb)+1)) !sz
         endif

         if (DSF(Jb) .LT. DBF(idx_bs(Jb)) .or. &
              DSF(Jb) .GT. DBF(idx_bs(Jb)+1)) THEN
            DSF(Jb) =  DSQRT(DBF(idx_bs(Jb)) * DBF(idx_bs(Jb)+1))
         endif

      ENDDO

  
end subroutine ssh_compute_granulo


end module mod_sshchemkinetic

