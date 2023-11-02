module mod_sshchemkinetic

  use aInitialization, only : attenuation, n_reaction, n_gas, &
                          humidity, temperature, pressure, &
                          photo_rcn, TB_rcn, fall_rcn, extra_rcn, &
                          index_RCT, index_PDT, &
                          fall_coeff, extra_coeff, &
                          photo_ratio, Arrhenius, &
                          ratio_PDT, kinetic_rate, chem_prod, chem_loss, &
                          drv_knt_rate, rcn_rate, gas_yield, species_name, &
                          concentration_gas_all, SumMc, YlH2O
  
  implicit none
  
contains

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
!     n_reaction: reaction number.
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

  integer i,j,k,ntot
  
  drv_knt_rate=0.d0
  ntot = size(drv_knt_rate)
  
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k)
  do i=1, ntot
    j = index_RCT(i,1) ! ircn
    k = index_RCT(i,2) ! isps
    drv_knt_rate(i) = kinetic_rate(j)*gas_yield(k)
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
 
  integer i,j,k,ntot

! Chemical loss terms.
  chem_loss=0.d0
  ntot = size(index_RCT,1)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,k)
  do i=1, ntot
    !j = index_RCT(i,1) ! ircn
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

  integer i,j,k,ntot

! Chemical production terms.
  chem_prod=0.d0
  ntot = size(index_PDT,1)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k)
  do i=1, ntot
    j = index_PDT(i,1) ! index in reaction list
    k = index_PDT(i,2) ! index in species list
    if (k /= 0) then ! remove 'NOTHING'
        chem_prod(k) = chem_prod(k) + rcn_rate(j) * ratio_PDT(i)
    endif
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
!     n_reaction: reaction number.
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

  integer i,j,k,ntot 

  rcn_rate = kinetic_rate  ! init
  ntot = size(index_RCT,1) ! number of reactants
  
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k)
  do i=1, ntot
    j = index_RCT(i,1) ! ircn
    k = index_RCT(i,2) ! index in species list
    rcn_rate(j) = rcn_rate(j) * gas_yield(k)
  enddo
!$OMP END PARALLEL DO

END subroutine ssh_rates

! =================================================================

subroutine ssh_kinetic(azi,RO2)
     
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
!     n_reaction: reaction number.
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

  integer :: i,j,k,s,tag
  double precision,INTENT(IN) :: RO2,azi

  !integer,INTENT(IN) :: photo_rcn(:,:), TB_rcn(:,:)
  !integer,INTENT(IN) :: fall_rcn(:), extra_rcn(:)

  ! for photolysis
  integer, parameter :: nsza = 11! sza size
  double precision, parameter :: szas(11) = (/0d0,1d1,2d1,3d1,4d1, &
                                      5d1,6d1,7d1,7.8d1,8.6d1,9d1/)
  double precision :: photo

  ! Arrhenius ! k = C1 * T**C2 * exp(-C3/T)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
  do i=1, n_reaction ! reactions
    kinetic_rate(i) = Arrhenius(i,1)*(temperature**Arrhenius(i,2)) &
                      * dexp(-Arrhenius(i,3)/temperature)
  enddo
!$OMP END PARALLEL DO 

  ! TBs(6) = ["RO2", "O2 ", "H2O", "M  ", "N2 ", "H2 "]
  do i=1,size(TB_rcn,1) ! no.TB reactions
  
    j = TB_rcn(i,1) ! reaction index
    k = TB_rcn(i,2) ! TB index

    select case (k)
        case(1) !RO2
            kinetic_rate(j) = kinetic_rate(j) * RO2
        case(2) ! O2
            kinetic_rate(j) = kinetic_rate(j) * SumMc * 0.2d0
        case(3) ! H2O
            kinetic_rate(j) = kinetic_rate(j) * YlH2O
        case(4) ! M
            kinetic_rate(j) = kinetic_rate(j) * SumMc
        case(5) ! N2
            kinetic_rate(j) = kinetic_rate(j) * SumMc * 0.8d0
        case(6) ! H2
            kinetic_rate(j) = kinetic_rate(j) * SumMc * 5.8d-7
    end select
  enddo

  ! Photolysis
  if (azi .lt. 9d1) then ! with photolysis

    do i=1,size(photo_rcn,1)
        j = photo_rcn(i,1) ! reaction index
        k = photo_rcn(i,2) ! pholysis index
        
        tag = 0 ! tag for computing photolysis rate

        do s=1,nsza-1 ! get photolysis rate
            if (azi.ge.szas(s) .and. azi.lt.szas(s+1)) then ! find zone
                photo = photo_ratio(k,s,4)
                photo = photo_ratio(k,s,3) + (azi-szas(s))*photo
                photo = photo_ratio(k,s,2) + (azi-szas(s))*photo
                photo = photo_ratio(k,s,1) + (azi-szas(s))*photo
                tag = 1
                exit
            endif
        enddo
        
        if (tag.eq.0) then
            print*, "Not find photolysis info: ",j,k,azi
            stop 
        endif

       ! Cloud attenuation.
       kinetic_rate(j) = kinetic_rate(j)*max(photo*attenuation, 0.d0) ! add a ratio
    enddo

  else ! no photolysis

    do i=1,size(photo_rcn,1) ! no.TB reactions
        j = photo_rcn(i,1) ! reaction index
        kinetic_rate(j) = 0.d0 ! set to zero
    enddo

  endif

  ! falloff reactions
  do i=1,size(fall_rcn)
  
    j = fall_rcn(i) ! reaction index

    call akkfo(j,i)

  enddo
  
  ! extra reactions
  do i=1,size(extra_rcn)
  
    j = extra_rcn(i) ! reaction index

    call akkextra4(j,i)

  enddo
END subroutine ssh_kinetic 
! =================================================================

SUBROUTINE akkfo(ire,ifo)

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

     xki = Arrhenius(ire,1)*((temperature/300.)**Arrhenius(ire,2)) &
    &      *dexp(-Arrhenius(ire,3)/temperature)                                   
     xkom = fall_coeff(ifo,1)*((temperature/300.)**fall_coeff(ifo,2)) &
    &     * dexp(-fall_coeff(ifo,3)/temperature)*SumMc                               
      
     if (xki .ne. 0.d0) then
       kratio = xkom/xki
     else
       print*,ire,ifo,"MCk: xki is zero",fall_coeff(ifo,:)
       stop
     endif 
      
     factor = 1./(1.+(dlog10(kratio))**2.) 
                                                                        
     IF ((fall_coeff(ifo,4)/=0.d0) .AND. (fall_coeff(ifo,5)==0.d0)) THEN 
       qfor = (xkom/(1.+kratio)) * (fall_coeff(ifo,4)**factor) 
                                                                        
! The following IF use a different equation for MCM rates constant      
     ELSE IF ((fall_coeff(ifo,4)/=0.d0) .AND. (fall_coeff(ifo,5)==1.d0)) THEN 
       factor2 = 1./(1.+(dlog10(kratio)/(0.75-1.27*dlog10(fall_coeff(ifo,4))))**2.)                                              
       qfor = (xkom/(1.+kratio)) * (fall_coeff(ifo,4)**factor2) 
                                                                        
     ELSE IF ((fall_coeff(ifo,4)==0.d0) .AND. (fall_coeff(ifo,5)==204.) .AND. &
    &          (fall_coeff(ifo,6)==0.17) .AND. (fall_coeff(ifo,7)==51.)) THEN      
       fcent = dexp(-temperature/fall_coeff(ifo,5)) + fall_coeff(ifo,6) &
    &        + dexp(-fall_coeff(ifo,7)/temperature)                                  
       factor2 = 1./(1.+(dlog10(kratio)/(0.75-1.27*dlog10(fcent)))**2.) 
       qfor = (xkom/(1.+kratio)) * (fcent**factor2) 
     ELSE
       print*, fall_coeff(ifo,:)
       print*, "MCk: unexpected setup of the auxiliary fall_coeff in falloff rxn",ire,ifo
       STOP 
     ENDIF 

    kinetic_rate(ire) = qfor * fall_coeff(ifo,8) ! ADD A RATIO
    
    !qln = log(qfor) 
                                                                        
END SUBROUTINE akkfo                                  

subroutine akkextra4(ire, iex)

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

    IMPLICIT NONE 
    
    integer :: label
    double precision :: kd, xk1, xk2, xk3, qfor
    double precision :: water_dimer

    ! INPUT                                                                 
    INTEGER, INTENT(in) :: ire, iex

! CASE SELECTOR FOR LABEL
! -----------------------------
  label=NINT(extra_coeff(iex,1))
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
      kd = 4.7856E-4*exp(1851/temperature-5.10485E-3*temperature) 
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

      xk2=extra_coeff(iex,2)*(temperature**extra_coeff(iex,3))* &
 &       dEXP(-extra_coeff(iex,4)/temperature)
      xk3=extra_coeff(iex,5)*(temperature**extra_coeff(iex,6))* &
 &       dEXP(-extra_coeff(iex,7)/temperature)
      qfor = qfor + xk3*SumMc / (1.+ ((xk3*SumMc)/xk2) )
! Add: ISOM
    CASE (10) 
      qfor = qfor*(extra_coeff(iex,2)*temperature**4 +  &
 &      extra_coeff(iex,3)*temperature**3 + &
 &      extra_coeff(iex,4)*temperature**2 + &
 &      extra_coeff(iex,5)*temperature + extra_coeff(iex,6))
  
!********************************************************               
! identificateur inconnu                                *               
!********************************************************               
! unidentified label
    CASE DEFAULT
      print*, '--error-- in akkextra. Type non-known: ', &
                label,ire,iex
      STOP 
  END SELECT
  
  ! add a ratio
  kinetic_rate(ire) = qfor * extra_coeff(iex,8)
  
end subroutine akkextra4

end module mod_sshchemkinetic

