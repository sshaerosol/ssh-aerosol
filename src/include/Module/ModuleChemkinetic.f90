module mod_sshchemkinetic

  use aInitialization, only : attenuation, n_reaction, n_gas, &
                          humidity, temperature, pressure, &
                          photo_rcn, TB_rcn, fall_rcn, extra_rcn, &
                          Nsps_rcn, index_RCT,index_PDT, &
                          fall_coeff, photo_ratio, Arrhenius, &
                          ratio_PDT, kinetic_rate, chem_prod, chem_loss, &
                          drv_knt_rate, rcn_rate, gas_yield, species_name, &
                          concentration_gas_all
  
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

  integer i,j,k,p,ind
  
  drv_knt_rate=0.d0
  
  k = 0
  do i=1, n_reaction ! reactions
    p = Nsps_rcn(i,1) ! no.reactants
    do j=1, p
        k = k + 1 ! traverse products
        ind = index_RCT(k) ! index in species list
        drv_knt_rate(k) = kinetic_rate(i)*gas_yield(ind)
    enddo
  enddo
  
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
 
  integer i,j,k,p,ind

! Chemical loss terms.
  chem_loss=0.d0
  
  k = 0
  do i=1, n_reaction ! no.reactions
    p = Nsps_rcn(i,1) ! no.reactants
    do j=1, p
        k = k + 1 ! traverse products
        ind = index_RCT(k) ! index in species list
        chem_loss(ind) = chem_loss(ind) + drv_knt_rate(k)
    enddo
  enddo
  
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

  integer i,j,k,ind,p

! Chemical production terms.
  chem_prod=0.d0

  k = 0
  do i=1, n_reaction ! reactions
    p = Nsps_rcn(i,2) ! no.products, default = 1

    do j=1, p
        k = k + 1 ! traverse products
        ind = index_PDT(k) ! index in species list
        if (ind /= 0) then ! remove nothing
            chem_prod(ind) = chem_prod(ind) + rcn_rate(i) * ratio_PDT(k)
        endif
    enddo

  enddo

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

  integer i,j,k,p,ind 

  rcn_rate=0.d0
  
  k = 0 ! index
  do i=1, n_reaction ! reactions
    
    p = Nsps_rcn(i,1) ! no.reactants
    rcn_rate(i) = kinetic_rate(i)
    
    do j=1, p
        k = k + 1 ! traverse reactants
        ind = index_RCT(k) ! index in species list
        rcn_rate(i) = rcn_rate(i) * gas_yield(ind)

    enddo
  enddo
  
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
!     TEMP: temperature ([K]).
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
  !integer,INTENT(IN) :: fall_rcn(:,:), extra_rcn(:,:)

  ! for photolysis
  integer, parameter :: nsza = 11! sza size
  double precision, parameter :: szas(11) = (/0d0,1d1,2d1,3d1,4d1, &
                                      5d1,6d1,7d1,7.8d1,8.6d1,9d1/)
  double precision :: SumM, cosX, secX, YlH2O, photo, temp
  
! Compute third body.
! Conversion = Avogadro*1d-6/Perfect gas constant.
! PRESS in Pascal, SUMM in molecules/cm3, TEMP in Kelvin

  temp = temperature
  SumM = pressure * 7.243D16 / temp

! Number of water molecules computed from the massic fraction
! (absolute humidity)
 
  YlH2O = 29.d0*SumM*humidity/(18.d0+11.d0*humidity)

  kinetic_rate = 0.d0 ! init

  ! Arrhenius
  do i=1, n_reaction ! reactions
    ! k = C1T**C2exp(-C3/T)
    kinetic_rate(i) = Arrhenius(i,1) * temp**Arrhenius(i,2)
    if (Arrhenius(i,3).ne.0.d0) then
      kinetic_rate(i) = kinetic_rate(i) * dexp(-1.d0*temp/Arrhenius(i,3))
    end if
  enddo
  
  ! TBs(6) = ["RO2", "O2 ", "H2O", "M  ", "N2 ", "H2 "]
  do i=1,size(TB_rcn,1) ! no.TB reactions
  
    j = TB_rcn(i,1) ! reaction index
    k = TB_rcn(i,2) ! TB index

    select case (k)
        case(1) !RO2
            kinetic_rate(j) = kinetic_rate(j) * RO2
        case(2) ! O2
            kinetic_rate(j) = kinetic_rate(j) * SumM * 0.2d0
        case(3) ! H2O
            kinetic_rate(j) = kinetic_rate(j) * YlH2O
        case(4) ! M
            kinetic_rate(j) = kinetic_rate(j) * SumM
        case(5) ! N2
            kinetic_rate(j) = kinetic_rate(j) * SumM * 0.8d0
        case(6) ! H2
            kinetic_rate(j) = kinetic_rate(j) * SumM * 5.8d-7
    end select
  enddo

  ! Photolysis
  if (azi .lt. 9d1) then ! with photolysis

    do i=1,size(photo_rcn,1) ! no.TB reactions
        j = photo_rcn(i,1) ! reaction index
        k = photo_rcn(i,2) ! pholysis index
        
        tag = 0 ! tag for computing photolysis rate

        do s=1,nsza-1 ! get photolysis rate
            if (azi.ge.szas(s) .and. azi.lt.szas(s+1)) then ! find zone
                photo = photo_ratio(k,s,4)
                photo = photo_ratio(k,s,3) + (azi-szas(s)*photo)
                photo = photo_ratio(k,s,2) + (azi-szas(s)*photo)
                photo = photo_ratio(k,s,1) + (azi-szas(s)*photo)
                tag = 1
                exit
            endif
        enddo
        
        if (tag.eq.0) then
            print*, "Not find photolysis info: ",j,k,azi
            stop 
        endif
        
        kinetic_rate(j) = photo
    enddo

    !     Cloud attenuation.
    kinetic_rate(j) = kinetic_rate(j) * attenuation

  else ! no photolysis

    do i=1,size(photo_rcn,1) ! no.TB reactions
        j = photo_rcn(i,1) ! reaction index
        kinetic_rate(j) = 0.d0 ! set to zero
    enddo

  endif

  ! other reactions
  
END subroutine ssh_kinetic 
! =================================================================

end module mod_sshchemkinetic

