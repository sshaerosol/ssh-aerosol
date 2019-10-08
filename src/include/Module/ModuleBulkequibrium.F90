!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------
!!
!!     -- DESCRIPTION
!!    This module contains methods to solve particle condensation/evaporation
!!    based on the gas/aerosol bulk equilibrium method.
!!-----------------------------------------------------------------------
Module iBulkequibrium
  use aInitialization
  use cThermodynamics
  use kSOAP
  implicit none
contains

  subroutine bulkequi_org(nesp_eq, &
       lwc, watorg, ionic, proton, liquid)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine solves gas/aerosol bulk equilibrium for organic species
!     based on H2O method.
!     equilibrium will be established between all size bin and bulk gas
!     mass flux will be redistributed into each cell based on rates
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!------------------------------------------------------------------------   
    implicit none
    
    integer::j,jesp,iter,s,k
    integer ::nesp_eq!number of species at equilibirum    
    double precision::dq(N_aerosol),qext(N_aerosol),qextold(N_aerosol)
    double precision::qaero(N_aerosol),qgas(N_aerosol)
    double precision:: Kelvin_effect(N_size,N_aerosol)
    double precision::ce_kernal_coef_tot(N_aerosol)
    double precision organion, watorg, proton, lwc, d_ms
    double precision rhop_tmp,emw_tmp,wet_diam
    integer :: eq_species(nesp_eq)! species compute with equilibrium
    double precision:: total_ms(N_aerosol)

    !**** SOAP ****
    double precision :: liquid(12), ionic
!!    DOUBLE PRECISION :: q(N_size*(1+N_aerosol)+N_aerosol)
    double precision qaerona,qaerocl,qgascl


!!     ******zero init
    do jesp=1,(N_aerosol-1)
	dq(jesp)=0.D0
	qgas(jesp)=0.d0
	qext(jesp)=0.D0
	qextold(jesp)=0.D0
	qaero(jesp)=0.d0
	ce_kernal_coef_tot(jesp)=0.D0!TINYA
	total_ms(jesp)=0.d0
    end do
    Kelvin_effect=1.000001
  !     ****** allocate equilibrium species list      

    do s=1,nesp_aec
      eq_species(s)=aec_species(s)
    enddo

  !compute local equi

  !compute only c/e coefficients
    do jesp=1,(N_aerosol-1)
      if (aerosol_species_interact(jesp).GT.0) then
	emw_tmp = molecular_weight_aer(jesp) * 1.D-6 ! g/mol
		  
	do j = 1,N_size	!FOR OGANIC
 	 wet_diam=wet_diameter(j)
	 rhop_tmp = rho_wet_cell (j) * 1.d9 !1400 ! kg/m3
 	 call COMPUTE_KELVIN_COEFFICIENT(&
		  Temperature,&          ! temperature (Kelvin)
		  emw_tmp,&       ! ext mol weight (g.mol-1)
		  surface_tension(jesp),&   ! surface tension (N.m-1) from INC
		  wet_diam,&         ! wet aero diameter (µm)
		  rhop_tmp,&      ! aerosol density (kg.m-3)
		  Kelvin_effect(j,jesp) )   ! kelvin effect coef (adim)
	  call COMPUTE_CONDENSATION_TRANSFER_RATE(&
		diffusion_coef(jesp), &! diffusion coef (m2.s-1)
		quadratic_speed(jesp),& ! quadratic mean speed (m.s-1)
		accomodation_coefficient(jesp),& ! accomadation coef (adim)
		wet_diameter(j),   & ! wet aero diameter (Âµm)
		ce_kernal_coef(j,jesp) ) ! c/e kernel coef (m3.s-1)
	  if(Kelvin_effect(j,jesp).lt.1.d0) Kelvin_effect(j,jesp)=1.000001
	    ce_kernal_coef_tot(jesp)= ce_kernal_coef_tot(jesp)&! compute total ce_kernal_coef coef
	                +ce_kernal_coef(j,jesp)*concentration_number(j)!&
	  !              *(1.d0/(Kelvin_effect(j,jesp)-1.d0))
	  !KS: ce_kernal_coef_tot(jesp)= ce_kernal_coef_tot(jesp)&! compute total ce_kernal_coef coef
	  ! 	      +ce_kernal_coef(j,jesp)*concentration_number(j) 
	enddo
      endif
    enddo

    do s=1,(N_aerosol_layers-1)
      jesp=List_species(s)
      do j=1,N_size
	qaero(jesp)=qaero(jesp)+concentration_mass(j,jesp)
      enddo
      qgas(jesp)=concentration_gas(jesp)!initial gas for H2O
      total_ms(jesp)=qaero(jesp)+qgas(jesp)
      qextold(jesp)=qaero(jesp)
    end do
    jesp=isorropia_species(2)
    qgas(jesp)=0.d0
    qgas(EH2O)=0.0
    organion = 0.D0
    watorg = 0.D0
    proton = 0.D0
    ionic = 0.D0
    lwc = 0.D0

#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
    qaerona = qaero(ENa)
    qaerocl = qaero(ECl)
    qgascl = qgas(ECl)
    qaero(ENa) = 0.D0
    qaero(ECl) = 0.D0
    qgas(ECl) = 0.D0
#endif

    call isoropia_drv(N_aerosol,&
         qaero,qgas,organion, watorg, ionic, proton, lwc, Relative_Humidity, Temperature, &
         liquid)

    if (ISOAPDYN.eq.0) then
       call soap_eq(watorg, lwc, Relative_Humidity, ionic, proton, &
            Temperature, qaero, qgas, liquid)
    endif

#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
    qaero(ENa) = qaerona
    qaero(ECl) = qaerocl
    qgas(ECl) = qgascl
#endif
    
    qgas(EH2O)=0.0
    
      !for organic
    do s=1,nesp_eq
      jesp=eq_species(s)
      if(aerosol_species_interact(jesp).GT.0) then
	qext(jesp)=qaero(jesp)
	concentration_gas(jesp)=qgas(jesp)!new qgas is used in N_aerosol bin
	if(qext(jesp).gt.0.d0) then
	  dq(jesp)=qext(jesp)-qextold(jesp)! compute delta aero conc
	else
	  dq(jesp)=-qextold(jesp)
	endif
      endif
    enddo

!     ******redistribute on each cell according to Rates
    !call bulkequi_redistribution_anck(concentration_number,concentration_mass,&
   ! nesp_eq,eq_species,N_size,dq,ce_kernal_coef,ce_kernal_coef_tot,Kelvin_effect)
    call bulkequi_redistribution(concentration_number,concentration_mass,&
    nesp_eq,eq_species,N_size,dq,ce_kernal_coef,ce_kernal_coef_tot)

  end subroutine bulkequi_org

  ! subroutine bulkequi_inorg(nesp_eq, lwc, ionic, proton, liquid, lwc_nsize)
  subroutine bulkequi_inorg(nesp_eq, lwc, ionic, proton, liquid)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine solves gas/aerosol bulk equilibrium for inorganic species
!     based on ISORROPIA method.
!     equilibrium will be established between all size bin < ICUT and bulk gas
!     mass flux will be redistributed into each cell based on rates
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!------------------------------------------------------------------------   
    implicit none
    integer::j,jesp,s
    integer ::nesp_eq!number of species at equilibirum    
    double precision::dq(N_aerosol),qext(N_aerosol),qextold(N_aerosol)
    double precision::qaero(N_aerosol),qgas(N_aerosol)
    double precision::ce_kernal_coef_tot(N_aerosol)
    double precision::aatoteq,qgasa,qgasi
    double precision:: Kelvin_effect(N_size,N_aerosol)
    double precision organion, watorg, proton, lwc, d_ms
    double precision rhop_tmp,emw_tmp,wet_diam    
    integer :: eq_species(nesp_eq)! species compute with equilibrium
    double precision:: total_ms(N_aerosol)
    double precision :: liquid(12), ionic
      
!!     ******zero init
    do jesp=1, N_aerosol !nesp_isorropia
      dq(jesp)=0.D0
      qgas(jesp)=0.d0
      qext(jesp)=0.D0
      qextold(jesp)=0.D0
      qaero(jesp)=0.d0
      ce_kernal_coef_tot(jesp)=0.D0
      total_ms(jesp)=0.d0
    end do
    Kelvin_effect=1.000001
  !     ****** if sulfate computed dynamically avoid it
  !     ****** allocate equilibrium species list
    do jesp=1,nesp_isorropia
      eq_species(jesp)=isorropia_species(jesp)
    enddo

  !compute only c/e coefficients

    ce_kernal_coef=0.d0

    do s=1, nesp_isorropia
       jesp=isorropia_species(s)
       ce_kernal_coef_tot(jesp) = 0.d0

      if (aerosol_species_interact(jesp).GT.0) then
	    ! compute total ce_kernal_coef coef (s-1)	
!	if(jesp.ne.isorropia_species(2)) then ! need to compute sulfate in case of sulfate equilibrium 
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
         IF (jesp.NE.ECl) THEN
#endif
	do j = 1,ICUT	!FOR INOGANIC
	  rhop_tmp = rho_wet_cell (j) * 1.d9 !1400 ! kg/m3
	  emw_tmp = molecular_weight_aer(jesp) * 1.D-6 ! g/mol

 	  wet_diam=wet_diameter(j)

        ! if (wet_diam .lt. 1.d-3) then
        !   write(*,*) "bulkequi_inorg: too small wet_diameter",wet_diam
         ! endif 
	  call COMPUTE_KELVIN_COEFFICIENT(&
		    Temperature,&          ! temperature (Kelvin)
		    emw_tmp,&       ! ext mol weight (g.mol-1)
		    surface_tension(jesp),&   ! surface tension (N.m-1) from INC
		    wet_diam,&         ! wet aero diameter (µm)
		    rhop_tmp,&      ! aerosol density (kg.m-3)
		    Kelvin_effect(j,jesp) )   ! kelvin effect coef (adim)
	  call COMPUTE_CONDENSATION_TRANSFER_RATE(&
		diffusion_coef(jesp), &! diffusion coef (m2.s-1)
		quadratic_speed(jesp),& ! quadratic mean speed (m.s-1)
		accomodation_coefficient(jesp),& ! accomadation coef (adim)
		wet_diameter(j),   & ! wet aero diameter (Âµm)
		ce_kernal_coef(j,jesp) ) ! c/e kernel coef (m3.s-1)
	    if(Kelvin_effect(j,jesp).lt.1.d0) Kelvin_effect(j,jesp)=1.000001
	    ce_kernal_coef_tot(jesp)= ce_kernal_coef_tot(jesp)&! compute total ce_kernal_coef coef
			+ce_kernal_coef(j,jesp)*concentration_number(j) !&
!!			*(1.d0/(Kelvin_effect(j,jesp)-1.d0)) ! KS
	enddo
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
        ENDIF
#endif
!	endif
      endif
    enddo

    if ((sulfate_computation.eq.0).AND.(with_nucl.EQ.0)) then
       jesp=isorropia_species(2)
       aatoteq = 0.d0
      ce_kernal_coef_tot(jesp) = 0.d0
      do j=1,N_size
	ce_kernal_coef_tot(jesp)= ce_kernal_coef_tot(jesp)&
	    +ce_kernal_coef(j,jesp)*concentration_number(j)
      end do
      qgasa=0.d0
      qgasi=concentration_gas(jesp)
    else
  !     ****** if sulfate computed dynamically avoid it
      jesp=isorropia_species(2)
      qgasi=concentration_gas(jesp)
    endif

    !compute total mass of each species
    do s=1,nesp_isorropia!inorganic
      jesp=isorropia_species(s)
      do j=1,ICUT
	qaero(jesp)=qaero(jesp)+concentration_mass(j,jesp)
      enddo
      qgas(jesp)=concentration_gas(jesp)
      total_ms(jesp)=qaero(jesp)+qgas(jesp)
      qextold(jesp)=qaero(jesp)
    end do

    if ((sulfate_computation.eq.0).AND.(with_nucl.EQ.0)) then
      !compute apparent gas concentration of sulfate
!     i.e. the SO4 gas conc actually
!     seen by equilibrium aerosols
      jesp=isorropia_species(2)
      do j=1,ICUT
	aatoteq=aatoteq+ce_kernal_coef(j,jesp)*concentration_number(j)
      enddo
      if(ce_kernal_coef_tot(jesp).gt.0.d0) then
	qgas(jesp)=concentration_gas(jesp)*aatoteq/ce_kernal_coef_tot(jesp)
      endif
      qgasa=qgas(jesp)
   else
      jesp=isorropia_species(2)
      qgas(jesp)=0.d0
      qgasa = 0.d0
    endif

#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
      qaero(ENa) = 0.D0
      qaero(ECl) = 0.D0
      qgas(ECl) = 0.D0
#endif

    organion= 0.D0
    watorg = 0.D0
    proton= 0.D0
    ionic = 0.D0
    lwc= 0.0

  !eqilibirum for inorganic
    call isoropia_drv(N_aerosol,&
         qaero,qgas,organion, watorg, ionic, proton, lwc, Relative_Humidity, Temperature, &
         liquid)

    qaero(N_aerosol)=lwc
    qgas(EH2O)=0.0

  !     ******redistribute on each cell according to Rates

    if ((sulfate_computation.eq.0).AND.(with_nucl.EQ.0)) then
      jesp=isorropia_species(2)
      ce_kernal_coef_tot(jesp)=aatoteq!for later redistribution
    endif
    do s=1, nesp_isorropia
      jesp=eq_species(s)
      if(aerosol_species_interact(jesp).GT.0) then      
	qext(jesp)=qaero(jesp)
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
       if(jesp.ne.ECl) then
#endif
	concentration_gas(jesp)=qgas(jesp)!new qgas is used in N_aerosol bin
	if(qext(jesp).gt.0.d0) then
	  dq(jesp)=qext(jesp)-qextold(jesp)! compute delta aero conc
        else
   	  dq(jesp)=-qextold(jesp)! compute delta aero conc
	endif
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
 	endif
#endif	
      endif
   enddo
   
  ! give back initial SO4 gas conc
  ! minus that consumed by equi bins
    if ((sulfate_computation.eq.0).AND.(with_nucl.EQ.0)) then
      jesp=isorropia_species(2)
      concentration_gas(jesp)=qgasi-qgasa
   else
      jesp=isorropia_species(2)
      concentration_gas(jesp)=qgasi
   endif

!    call bulkequi_redistribution_anck(concentration_number,concentration_mass,&
!    nesp_isorropia,eq_species,N_size,dq,ce_kernal_coef,ce_kernal_coef_tot,Kelvin_effect)

    call bulkequi_redistribution(concentration_number,concentration_mass,&
      nesp_isorropia,eq_species,ICUT,dq,ce_kernal_coef,ce_kernal_coef_tot)

  end subroutine bulkequi_inorg


  
  subroutine bulkequi_redistribution(c_number,c_mass,nesp_eq,eq_species,end_bin,dq,AAi,ce_kernal_coef_tot)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine redistribute bulk mass variations into each bins
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     c_number: aerosol number concentration(#/m^3)
!     c_mass: aerosol mass concentration(µg/m^3)
!     nesp_eq: number of species at equilibirum
!     eq_species: the list of species pointer at equilibirum
!     end_bin: marks the number of bins concerned during the redistribution
!     dq: bulk mass variations
!     AAi: c/e kernel coefficient          ([m3.s-1]).
!     ce_kernal_coef_tot: sum of c/e kernel coefficient          ([m3.s-1]).
!------------------------------------------------------------------------
    implicit none
    integer::j,s,jesp,end_bin,iclip
    integer::nesp_eq
    integer::eq_species(nesp_eq)
    double precision::totaer,temp_mass,totaa
    double precision::ce_kernal_coef_tot(N_aerosol)
    double precision::dq(N_aerosol)
    double precision::AAi(N_size,N_aerosol)
    double precision::frac(N_size,N_aerosol)
    double precision::c_number(N_size)
    double precision::c_mass(N_size,N_aerosol_layers)
    double precision::frac_bin(N_sizebin,N_aerosol)
    double precision::dm_bin(N_sizebin,N_aerosol)
    integer::iclip_bin(N_sizebin,N_aerosol)
    integer::k,ibegin
 
    frac_bin=0.d0
    dm_bin=0.d0
    iclip_bin=0
    if ((sulfate_computation.eq.0).AND.(with_nucl.EQ.0)) then
       ibegin = -1
    else
       ibegin = isorropia_species(2)
    endif
   
    do s=1, nesp_eq
      jesp=eq_species(s)
      if (aerosol_species_interact(jesp).GT.0) then

       if(jesp.NE.ibegin) then
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
       IF (jesp.NE.ECl .and. jesp.ne.ENa) THEN
#endif
         iclip=0
         do j=1,end_bin!judgment
	  if(ce_kernal_coef_tot(jesp).gt.0.d0) then
  	    frac(j,jesp)= AAi(j,jesp)*c_number(j)/ce_kernal_coef_tot(jesp)
	    temp_mass=c_mass(j,jesp)+dq(jesp)*frac(j,jesp)
	    if(temp_mass.lt.0.d0) iclip=1!case of over evaporation
          endif
         enddo
         if(iclip.eq.1) then !over evaporate
  	  totaer=0.d0
	  do j=1,end_bin
	    totaer=totaer+c_mass(j,jesp)
	  enddo
	  if(totaer.gt.0.d0) then
	    do j=1,end_bin
	      frac(j,jesp)=c_mass(j,jesp)/totaer
	      c_mass(j,jesp)=c_mass(j,jesp)+dq(jesp)*frac(j,jesp)
	      if(dq(jesp)*frac(j,jesp).ne.0.d0) then
	        k=concentration_index(j,1)
	        frac_bin(k,jesp)=frac_bin(k,jesp)+frac(j,jesp)
	        dm_bin(k,jesp)=dm_bin(k,jesp)+dq(jesp)*frac(j,jesp)
	        iclip_bin(k,jesp)=iclip
             endif
	    enddo
	  endif
        else!normal case
	  do j=1,end_bin
	   if(ce_kernal_coef_tot(jesp).gt.0.d0) then
	    frac(j,jesp)= AAi(j,jesp)*c_number(j)/ce_kernal_coef_tot(jesp)
	    c_mass(j,jesp)=c_mass(j,jesp)+dq(jesp)*frac(j,jesp)
	    if(dq(jesp)*frac(j,jesp).ne.0.d0) then
	      k=concentration_index(j,1)
	      frac_bin(k,jesp)=frac_bin(k,jesp)+frac(j,jesp)
	      dm_bin(k,jesp)=dm_bin(k,jesp)+dq(jesp)*frac(j,jesp)
	      iclip_bin(k,jesp)=iclip
	    endif
           endif
	  enddo
        endif
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
       ENDIF
#endif
     endif
    endif
   enddo

  end subroutine bulkequi_redistribution

  subroutine bulkequi_redistribution_anck(c_number,c_mass,nesp_eq,eq_species,end_bin,dq,AAi,ce_kernal_coef_tot,ck_ef)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine redistribute bulk mass variations into each bins
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     c_number: aerosol number concentration(#/m^3)
!     c_mass: aerosol mass concentration(µg/m^3)
!     nesp_eq: number of species at equilibirum
!     eq_species: the list of species pointer at equilibirum
!     end_bin: marks the number of bins concerned during the redistribution
!     dq: bulk mass variations
!     AAi: c/e kernel coefficient          ([m3.s-1]).
!     ce_kernal_coef_tot: sum of c/e kernel coefficient          ([m3.s-1]).
!------------------------------------------------------------------------   
    implicit none
    integer::j,s,jesp,end_bin,iclip
    integer::nesp_eq
    integer::eq_species(nesp_eq)
    double precision::totaer,temp_mass
    double precision::ce_kernal_coef_tot(N_aerosol)
    double precision::dq(N_aerosol)
    double precision::AAi(N_size,N_aerosol)
    double precision::frac(N_size,N_aerosol)
    double precision::c_number(N_size)
    double precision::c_mass(N_size,N_aerosol_layers)
    double precision::ck_ef(N_size,N_aerosol)
    double precision::frac_bin(N_sizebin,N_aerosol)
    double precision::dm_bin(N_sizebin,N_aerosol)
    integer::iclip_bin(N_sizebin,N_aerosol)
    integer::k,ibegin
  
    frac_bin=0.d0
    dm_bin=0.d0
    iclip_bin=0
    if ((sulfate_computation.eq.0).AND.(with_nucl.EQ.0)) then
       ibegin = -1
    else
       ibegin = isorropia_species(2)
    endif
    
    do s=1, nesp_eq
      jesp=eq_species(s)
      if (aerosol_species_interact(jesp).GT.0) then

       if(jesp.NE.ibegin) then
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
      IF (jesp.NE.ECl) THEN
      IF (jesp.NE.ENa) THEN
#endif
         iclip=0
         do j=1,end_bin!judgment
	  if(ce_kernal_coef_tot(jesp).gt.0.d0) then
  	    frac(j,jesp)= AAi(j,jesp)*c_number(j)*(1.d0/(ck_ef(j,jesp)-1.d0))/ &
                 ce_kernal_coef_tot(jesp)
	    temp_mass=c_mass(j,jesp)+dq(jesp)*frac(j,jesp)
	    if(temp_mass.lt.0.d0) iclip=1!case of over evaporation
          endif
         enddo
         if(iclip.eq.1) then !over evaporate
  	  totaer=0.d0
	  do j=1,end_bin
	    totaer=totaer+c_mass(j,jesp)
	  enddo
	  if(totaer.gt.0.d0) then
	    do j=1,end_bin
             if(totaer.GT.0.0) then
	      frac(j,jesp)=c_mass(j,jesp)/totaer
	      c_mass(j,jesp)=c_mass(j,jesp)+dq(jesp)*frac(j,jesp)
	      if(dq(jesp)*frac(j,jesp).ne.0.d0) then
	        k=concentration_index(j,1)
	        frac_bin(k,jesp)=frac_bin(k,jesp)+frac(j,jesp)
	        dm_bin(k,jesp)=dm_bin(k,jesp)+dq(jesp)*frac(j,jesp)
	        iclip_bin(k,jesp)=iclip
	      endif
             endif
	    enddo
	  endif
        else!normal case
	  do j=1,end_bin
	   if(ce_kernal_coef_tot(jesp).gt.0.d0) then
	    frac(j,jesp)= AAi(j,jesp)*c_number(j)*(1.d0/(ck_ef(j,jesp)-1.d0))/ &
                   ce_kernal_coef_tot(jesp)
	    c_mass(j,jesp)=c_mass(j,jesp)+dq(jesp)*frac(j,jesp)
	    if(dq(jesp)*frac(j,jesp).ne.0.d0) then
	      k=concentration_index(j,1)
	      frac_bin(k,jesp)=frac_bin(k,jesp)+frac(j,jesp)
	      dm_bin(k,jesp)=dm_bin(k,jesp)+dq(jesp)*frac(j,jesp)
	      iclip_bin(k,jesp)=iclip
	    endif
           endif
	  enddo
        endif
        endif
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
      ENDIF
      ENDIF
#endif
      endif
    enddo

  end subroutine bulkequi_redistribution_anck

end Module iBulkequibrium
  
  
