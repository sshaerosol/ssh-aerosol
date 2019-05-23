!!-----------------------------------------------------------------------
!!     Copyright (C) 2012-2018, ENPC - EDF R&D - INERIS
!!     Author(s): Shupeng Zhu, Edouard Debry, Karine Sartelet
!!
!!     This file is part of the Size Composition Resolved Aerosol Model (SCRAM), a
!!     component of the SSH-aerosol model.
!!
!!     SSH-aerosol is a free software; you can redistribute it and/or modify
!!     it under the terms of the GNU General Public License as published
!!     by the Free Software Foundation; either version 2 of the License,
!!     or (at your option) any later version.
!!
!!     SSH-aerosol is distributed in the hope that it will be useful, but
!!     WITHOUT ANY WARRANTY; without even the implied warranty of
!!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
!!     General Public License for more details.
!!
!!-----------------------------------------------------------------------

Module cThermodynamics
  use aInitialization
  implicit none
contains
  subroutine compute_wet_mass_diameter(start_bin,end_bin,c_mass,c_number,c_inti, &
      wet_m,wet_d,wet_v)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine computes particle wet diameter for each
!     bin between start_bin and end_bin as well as their water content
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     start_bin: index of the first bin need the computation
!     end_bin: index of the last bin need the computation
!     c_number: aerosol number concentration(#/m^3)
!     c_mass: aerosol mass concentration(µg/m^3)
!
!     -- OUTPU VARIABLES
!
!     c_inti: aerosol internal species mass concentration(µg/m^3)
!     wet_m: particle wet mass
!     wet_d: particle wet diameter
!     wet_v: particle wet volume
!     dry_d: particle dry diameter
!------------------------------------------------------------------------
    implicit none
    integer :: jesp,s,j,k,i,start_bin,end_bin
    double precision :: c_number(N_size)
    double precision :: c_mass(N_size,N_aerosol)
    double precision :: c_inti(N_size,N_inside_aer)
    double precision :: wet_m(N_size)!wet mass
    double precision :: wet_d(N_size)!wet diameter
    double precision :: wet_v(N_size)!wet volume
    double precision :: dry_d(N_size)!dry diameter
    double precision :: qext(N_aerosol)
    double precision :: qinti(N_inside_aer)
    double precision :: lwc,qti,vad
    double precision :: aero(5)
    double precision :: rhoaer!aerosol density based on the internal composition
    !print*,"wet_diam_estimation",wet_diam_estimation
    double precision :: rhoaer_dry

!      OPEN(UNIT=10,FILE="test-report.txt")
    do j=start_bin,end_bin
      !initialization
      do jesp=1,N_aerosol
	qext(jesp)=c_mass(j,jesp)
      enddo
      c_mass(j,EH2O)=0.d0
      do jesp=1,N_inside_aer
	qinti(jesp)=0.d0
      enddo
!     ******total dry mass
      qti=0.D0
      do s=1,(N_aerosol-1)
	jesp=List_species(s)
	qti=qti+qext(jesp)     !total dry mass µg.m-3
      end do
      
      ! if (qti.gt.TINYM*N_aerosol) then
      ! The threshold 0.0 for the minimum aerosol concentration 
      ! should be avoided because it leads to a too small particle diameter.
      ! The threshold is set to 1.d-6 (YK).
      if (qti.gt.1.d-6) then 
         if(wet_diam_estimation.eq.0) then !with isorropia
	    do i=1,nesp_isorropia
               jesp=isorropia_species(i)
               aero(i)=qext(jesp)
            enddo
	    call calculatewater(aero,qinti,lwc)
! 	    do i=1,nesp_isorropia
! 	      jesp=isorropia_species(i)
! 	      qext(jesp)=aero(i)
! 	      !if(qext(jesp).gt.1.d3) print*,"huge_aer",j,jesp,qext(jesp),lwc
! 	    enddo
	    qext(EH2O)=lwc
	    call VOLAERO(N_aerosol,qext,qinti,rhoaer,rhoaer_dry) !! YK
	    if(rhoaer.gt.0.d0) then
	      rho_wet_cell(j)=rhoaer
	    else
	      rhoaer=density_aer_bin(j)
	    endif
            do jesp=1,N_inside_aer
               c_inti(j,jesp)=qinti(jesp)
            enddo

            if(rhoaer.gt.0.d0.and.c_number(j).gt.TINYN) then
               ! vad=qti/rhoaer!qti total dry mass !! YK (rhoaer is aerosol density including water)
               vad=qti/rhoaer_dry!qti total dry mass
               wet_v(j)=vad+qext(EH2O)/rhoaer!: wet volume aerosol concentration (µm3/m3).
               ! aerosol diameter
               ! qn cannot be zero, checked in eqpart routine
               dry_d(j)=(vad/c_number(j)/cst_pi6)**cst_FRAC3 ! dry aerosol dimaeter µm
               k=concentration_index(j, 1)
!	    if(dry_d(j).lt.diam_bound(k)) then
!	      dry_d(j)=size_diam_av(k)
!	      c_number(j)=vad/(cst_pi6*dry_d(j)**3.0)
!	    endif
               wet_d(j)=(wet_v(j)/c_number(j)/cst_pi6)**cst_FRAC3 ! wet aerosol diameter µm
               wet_d(j)=DMAX1(wet_d(j),dry_d(j))!wet diameter is always larger than dry diameter
               c_mass(j,EH2O)=qext(EH2O)
               wet_m(j)=(qti+qext(EH2O))/c_number(j) ! single wet mass (µg)
            else
               c_mass(j,EH2O)=qext(EH2O)
               k=concentration_index(j, 1)
               wet_d(j)=size_diam_av(k)
               wet_m(j)=size_mass_av(k)
               !	    if(qext(EH2O).gt.1.d0) print*,"Mod_Thy nub_mas abnormal",qti,&
               !	      c_number(j),qext(EH2O)
            endif

	 elseif(wet_diam_estimation.eq.1) then!with fastdiam
            ! if(rhoaer.gt.0.d0.and.c_number(j).gt.TINYN) then !! YR: rhoaer is not initialized.
           if(c_number(j).gt.TINYN) then
              rhoaer=density_aer_bin(j)
              vad=qti/rhoaer ! dry aerosol volume (µm3.m-3)
              dry_d(j)=(vad/c_number(j)/cst_pi6)**cst_FRAC3 ! µm, cannot be less than tinyn

              CALL gerber_wet_diameter(Relative_humidity, &! relative humidity (0<adim<1)
                   Temperature,     &     ! temperature (Kelvin)
                   dry_d(j),   &   ! dry aero diam (m)
                   wet_d(j))      ! wet aero diam (m)
              wet_d(j)=DMAX1(wet_d(j),dry_d(j))!wet diameter is always larger than dry diameter
              wet_v(j)=c_number(j)*cst_pi6*wet_d(j)*wet_d(j)*wet_d(j)!: wet volume aerosol concentration (µm3/m3).
              qext(EH2O)=DMAX1(wet_v(j)-vad,0.D0)*mass_density(EH2O) ! water conc (µg.m-3)
              c_mass(j,EH2O)=qext(EH2O)
              wet_m(j)=(qti+qext(EH2O))/c_number(j) ! single wet mass (µg)
              if (wet_d(j) .gt. 1.D100) then
                 print *, j," fastdiam, after gerber",&
                      " wet diam", wet_d(j), "dry diameter", dry_d(j),&
                      "TEMP", Temperature, "RH", Relative_humidity
              endif
              ! if (wet_d(j) .lt. 1.D-3) then
              !    print *, j," fastdiam, after gerber: ",&
              !         "wet diam", wet_d(j), "dry diameter", dry_d(j),&
              !         "TEMP", Temperature, "RH", Relative_humidity,&
              !         "dry conc", qti, "water conc", qext(EH2O),"number",c_number(j)
              ! endif
	   else
	      c_mass(j,EH2O)=qext(EH2O)
	      k=concentration_index(j, 1)
	      wet_d(j)=size_diam_av(k)
	      wet_m(j)=size_mass_av(k)
           endif
	 endif
      else
	! if too few aerosols or too few mass
	! we set variables of given bins as
	! its initial fixed ones,
	! thus avoiding zero values
	k=concentration_index(j, 1)
	wet_d(j)=size_diam_av(k)
	wet_m(j)=size_mass_av(k)
      endif
    enddo

  end subroutine compute_wet_mass_diameter

  subroutine calculatewater(aero,qinti,lwc)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine computes particle water content
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     aero: aerosol mass concentration(µg/m^3)
!
!     -- OUTPU VARIABLES
!
!     qinti: aerosol internal species mass concentration(µg/m^3)
!     lwc: particle water content
!------------------------------------------------------------------------
    implicit none
    integer jesp,j
    double precision aero(5)
    double precision organion, watorg, proton
    double precision lwc,ionic

    double precision wi(5),w(5),gas(3),cntrl(2), other(6)
    double precision liquid(N_liquid),solid(N_solid)
    double precision qinti(N_inside_aer)

    organion = 0.D0
    watorg = 0.D0
    proton = 0.D0

    cntrl(1) = 1.D0!reverse mode
    cntrl(2) = 1.D0

    gas(1)=0.d0
    gas(2)=0.d0
    gas(3)=0.d0

!     conversion unit for isorropia needed in mol.m-3
    do j=1,nesp_isorropia
	jesp=isorropia_species(j)
	if(molecular_weight_aer(jesp).gt.0.d0) then
	  wi(j)= aero(j)& ! µg.m-3
	      /molecular_weight_aer(jesp)!&  ! µg.mol-1
	endif
    end do

!     call isorropia fortran routine
    call ISOROPIA(wi, Relative_Humidity, Temperature, cntrl, w, gas,&
	  liquid, solid, other, organion, watorg)
!     clipping to tinym

!     Aqueous phase total liquid water content and pH (proton) concentration
    lwc = liquid(IH2O) * imw(IH2O) ! microg.m-3
    ionic = other(5)
    proton = liquid(IH) * imw(IH) !* gammaH  ! microg.m-3 but equivalent to micromol.m-3

    if (gas(1).lt.0.d0) gas(1)=tinym
    if (gas(2).lt.0.d0) gas(2)=tinym
    if (gas(3).lt.0.d0) gas(3)=tinym

!     Aqueous phase total liquid water content and pH (proton) concentration
    do jesp=IH,IOH
      qinti(jesp)= DMAX1(liquid(jesp),0.D0)*imw(jesp)   ! moles to µg MOLAR WEIGHT
    end do
			! solid inorg aerosol
    do jesp=SNaNO3,SLC
      qinti(jesp)= DMAX1(solid(jesp-12),0.D0)&
	  *smw(jesp)        ! moles to µg !molecular_weight_solid(jesp)
    end do
			      ! liquid water content
    !lwc= qinti(IH2O)+qinti(IOH)*1.05882352941D0 ! mwh2o/mwioh
  end subroutine calculatewater

  subroutine update_wet_diameter(start_bin,end_bin,c_mass,c_inti,c_number,wet_m,&
      wet_d,wet_v,dry_d)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine computes particle wet diameter for each
!     bin between start_bin and end_bin based on known water content
!     It also update the total water and pH of the system
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     start_bin: index of the first bin need the computation
!     end_bin: index of the last bin need the computation
!     c_number: aerosol number concentration(#/m^3)
!     c_mass: aerosol mass concentration(µg/m^3)
!     c_inti: aerosol internal species mass concentration(µg/m^3)
!
!     -- OUTPUT VARIABLES
!
!     wet_m: particle wet mass
!     wet_d: particle wet diameter
!     wet_v: particle wet volume
!     dry_d: particle dry diameter
!------------------------------------------------------------------------
    implicit none
    integer :: jesp,s,j,k,start_bin,end_bin
    double precision ::c_mass(N_size,N_aerosol)
    double precision ::c_number(N_size)
    double precision :: c_inti(N_size,N_inside_aer)
    double precision :: wet_m(N_size)!wet mass
    double precision :: wet_d(N_size)!wet diameter
    double precision :: wet_v(N_size)!wet volume
    double precision :: dry_d(N_size)!dry diameter
    double precision :: qext(N_aerosol)
    double precision :: qinti(N_inside_aer)
    double precision :: qti,vad,rhoaer
    double precision :: rhoaer_dry ! aerosol density without water (YK)

    total_water=0.d0
    total_IH=0.d0
    do j=start_bin,end_bin
      ! rho_wet_cell(j)=0.d0 ! YK
      rho_wet_cell(j)=fixed_density
      total_water=total_water+c_mass(j,EH2O)
      total_IH=total_IH+c_inti(j,IH)
      qti=0.D0
      do s=1,(N_aerosol-1)
	jesp=List_species(s)
	qti=qti+c_mass(j,jesp)     ! µg.m-3
      end do
      if (c_number(j).gt.TINYN.and.qti.gt.TINYM ) then
	do jesp=1,N_aerosol
	  qext(jesp)=c_mass(j,jesp)
	enddo
	do jesp=1,N_inside_aer
	  qinti(jesp)=c_inti(j,jesp)
	enddo
	call VOLAERO(N_aerosol,qext,qinti,rhoaer,rhoaer_dry)

	  if(rhoaer.gt.0.d0) then
	    rho_wet_cell(j)=rhoaer
	  else
	    rhoaer=density_aer_bin(j)
 	  endif
	if(rhoaer.gt.0.d0) then
	  ! vad=qti/rhoaer!qti total dry mass
	  vad=qti/rhoaer_dry!qti total dry mass YK: aerosol density without water.
	  wet_v(j)=vad+c_mass(j,EH2O)/rhoaer!: wet volume aerosol concentration (µm3/m3).
	  ! aerosol diameter
	  ! qn cannot be zero, checked in eqpart routine
	  dry_d(j)=(vad/c_number(j)/cst_pi6)**cst_FRAC3 ! dry aerosol dimaeter µm
	  wet_d(j)=(wet_v(j)/c_number(j)/cst_pi6)**cst_FRAC3 ! wet aerosol diameter µm
	  wet_d(j)=DMAX1(wet_d(j),dry_d(j))!wet diameter is always larger than dry diameter
	  wet_m(j)=(qti+c_mass(j,EH2O))/c_number(j) ! single wet mass (µg)
	endif
      else
	! if too few aerosols or too few mass
	! we set variables of given bins as
	! its initial fixed ones,
	! thus avoiding zero values
	k=concentration_index(j, 1)
	wet_d(j)=size_diam_av(k)
	wet_m(j)=size_mass_av(k)
      endif
    enddo
    if(total_IH*total_water.gt.0.d0) total_PH=-log10((total_IH/1.D6)/ (total_water/1.d9))

  end subroutine update_wet_diameter

  subroutine update_wet_diameter_liquid(start_bin,end_bin,c_mass,c_number,wet_m,&
      wet_d,wet_v,dry_d)

!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine computes particle wet diameter for each
!     bin between start_bin and end_bin based on known water content
!     It also update the total water and pH of the system
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     start_bin: index of the first bin need the computation
!     end_bin: index of the last bin need the computation
!     c_number: aerosol number concentration(#/m^3)
!     c_mass: aerosol mass concentration(µg/m^3)
!     c_inti: aerosol internal species mass concentration(µg/m^3) from isoroppia
!     lwcorg = organic liquid water content
!
!     -- OUTPUT VARIABLES
!
!     wet_m: particle wet mass
!     wet_d: particle wet diameter
!     wet_v: particle wet volume
!     dry_d: particle dry diameter
!------------------------------------------------------------------------
    implicit none
    integer :: jesp,s,j,k,start_bin,end_bin

    double precision ::c_mass(N_size,N_aerosol)
    double precision ::c_number(N_size)
    double precision :: c_inti(N_size,N_inside_aer)
    double precision :: wet_m(N_size)!wet mass
    double precision :: wet_d(N_size)!wet diameter
    double precision :: wet_v(N_size)!wet volume
    double precision :: dry_d(N_size)!dry diameter
    double precision :: qext(N_aerosol)
    double precision :: qinti(N_inside_aer)
    double precision :: qti,vad,rhoaer

    total_water=0.d0
    do j=start_bin,end_bin
      ! rho_wet_cell(j)=0.d0 ! YK
      rho_wet_cell(j)=fixed_density
      total_water=total_water+c_mass(j,EH2O)
      qti=0.D0
      do s=1,(N_aerosol-1)
	jesp=List_species(s)
	qti=qti+c_mass(j,jesp)     ! µg.m-3
      end do
      if (c_number(j).gt.TINYN.and.qti.gt.TINYM ) then
	do jesp=1,N_aerosol
	  qext(jesp)=c_mass(j,jesp)
	enddo

      ! call compute_density(N_sizebin,N_aerosol,EH2O,TINYM,c_mass,&
      !        mass_density_aer,j,rhoaer)
      call compute_density(N_size,N_aerosol,EH2O,TINYM,c_mass,&
             mass_density,j,rhoaer)

	  if(rhoaer.gt.0.d0) then
	    rho_wet_cell(j)=rhoaer
	  else
	    rhoaer=density_aer_bin(j)
	    rho_wet_cell(j)=rhoaer
 	  endif
	if(rhoaer.gt.0.d0) then
	  vad=qti/rhoaer!qti total dry mass
	  wet_v(j)=vad+c_mass(j,EH2O)/rhoaer!: wet volume aerosol concentration (µm3/m3).
	  ! aerosol diameter
	  ! qn cannot be zero, checked in eqpart routine
	  dry_d(j)=(vad/c_number(j)/cst_pi6)**cst_FRAC3 ! dry aerosol dimaeter µm
	  wet_d(j)=(wet_v(j)/c_number(j)/cst_pi6)**cst_FRAC3 ! wet aerosol diameter µm
	  wet_d(j)=DMAX1(wet_d(j),dry_d(j))!wet diameter is always larger than dry diameter
	  wet_m(j)=(qti+c_mass(j,EH2O))/c_number(j) ! single wet mass (µg)
	endif
      else
	! if too few aerosols or too few mass
	! we set variables of given bins as
	! its initial fixed ones,
	! thus avoiding zero values
	k=concentration_index(j, 1)
	wet_d(j)=size_diam_av(k)
	wet_m(j)=size_mass_av(k)
      endif
    enddo
    if(total_IH*total_water.gt.0.d0) total_PH=-log10((total_IH/1.D6)/ (total_water/1.d9))
  end subroutine update_wet_diameter_liquid

  subroutine VOLAERO(nesp_aer,qext,qinti,rhoaer, rhoaer_dry)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine computes the dry and wet aerosol volumes and the
!     aerosol density according to the internal composition.
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     nesp_aer: number of aerosol species
!     qext: aerosol mass concentration(µg/m^3)
!     qinti: aerosol internal species mass concentration(µg/m^3)
!
!     -- OUTPU VARIABLES
!
!     rhoaer: aerosol density according to the internal composition.
!------------------------------------------------------------------------
    implicit none
    integer nesp_aer
    double precision qext(nesp_aer),qinti(N_inside_aer),rhoaer

    integer jesp,s
    double precision vid,viw,vis,vil,vod,vad,vaw,sumint
    double precision rhoaer_dry, sum_dry !! YK

!!!     ******inorganic volume
    vid=0.D0
    vis=0.D0
    vil=0.D0
    viw=0.D0
    vod=0.D0
    sumint=0.d0
    sum_dry = 0.d0
    !print*,mass_density
!!     Mineral Dust
    vis = vis + qext(EMD)/mass_density(EMD)
    sumint = sumint + qext(EMD)

!!     Black Carbon
    vis = vis + qext(EBC)/mass_density(EBC)
    sumint = sumint + qext(EBC)

!!     Inorganic
    ! compute solid aerosol volume
    do jesp=SNaNO3,SLC
      vis=vis+qinti(jesp)/SMD(jesp)
    end do
    ! compute liquid aerosol volume
    ! sodium volume
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
     vil=vil + qext(ENa)/mass_density(ENa)
#endif
! 
#ifndef WITHOUT_NACL_IN_THERMODYNAMICS
     vil=vil + qinti(INa)/mass_density(ENa)
#endif
! 
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
     vil = vil + qext(ECl)/mass_density(ECl) ! HCl volume
#endif
			      ! ammonium volume
    vil=vil+( qinti(INH4)&
	+qinti(INH4)&
	*0.944444444444D0&    ! mwnh3/mwinh4
	)/mass_density(ENH4)          ! µg.µm-3

			      ! nitric acid volume
    vil=vil+( qinti(IHNO3)&
	+qinti(INO3)&
	*1.01612903226D0&     ! mwhno3/mwino3
	)/mass_density(ENO3)          ! µg.µm-3

			      ! chlorhydric acid volume
#ifndef WITHOUT_NACL_IN_THERMODYNAMICS
    vil=vil+( qinti(IHCl)&
	+qinti(ICl)&
	*1.02816901408D0&     ! mwhcl/mwicl
	)/mass_density(ECl)           ! µg.µm-3
#endif
			      ! sulfuric acid volume
    vil=vil+( qinti(IHSO4)&
	*1.01030927835D0&     ! mwh2so4/mwihso4
	+qinti(ISO4)&
	*1.02083333333D0&     ! mwh2so4/mwiso4
	)/mass_density(ESO4)          ! µg.µm-3

    vid=vil+vis               ! dry inorg vol µm3.m-3

			      ! water volume
    !print*,'qext(EH2O)',qext(EH2O),mass_density(EH2O)
    viw=qext(EH2O)/mass_density(EH2O)

    ! total inorg internal mass
    do jesp=1,N_inside_aer
      sumint=sumint+qinti(jesp)
    enddo
    ! correction for water
    sumint=sumint-qinti(IH2O)-qinti(IOH)+qext(EH2O)

#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
! 			      ! correction when no NaCl in internal composition
     sumint = sumint - qinti(INa) + qext(ENa)
     sumint = sumint - qinti(IHCl) - qinti(ICl)*1.02816901408D0&
 	+ qext(ECl)
#endif

!!     dry organic volume and total internal mass
    do s=1,nesp_aec
      jesp=aec_species(s)
      vod=vod+qext(jesp)/mass_density(jesp)
      sumint=sumint+qext(jesp)
    enddo

    do s=1,nesp_pankow
      jesp=pankow_species(s)
      vod=vod+qext(jesp)/mass_density(jesp)
      sumint=sumint+qext(jesp)
    enddo

    do s=1,nesp_pom
      jesp=poa_species(s)
      vod=vod+qext(jesp)/mass_density(jesp)
      sumint=sumint+qext(jesp)
    enddo
!!     ******tot dry and wet vol, µm3.m-3
    vad=vod+vid
    vaw=vad+viw

!!     Notice that the density is based on the internal composition
!!     and as such bigger than the minimal density and less that the maximal
!!     one
    if ((vaw.gt.0.d0).AND.(sumint.gt.0.D0)) then
	rhoaer=sumint/vaw
        rhoaer_dry = (sumint - qext(EH2O)) / vad
    else
	rhoaer=fixed_density
        rhoaer_dry = fixed_density
    endif
  end subroutine VOLAERO

   subroutine EQINORG(nesp_aer,qext,qinti,surface_equilibrium_conc,lwc,ionic,proton,aerliq)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine computes the equilibrium between inorganic aerosols
!     and gas-phase (reverse mode).
!     It calls ISORROPIA by Nenes et al.
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     nesp_aer: number of aerosol species
!     qext: aerosol mass concentration(µg/m^3)
!     qinti: aerosol internal species mass concentration(µg/m^3)
!
!     -- OUTPU VARIABLES
!
!     surface_equilibrium_conc: surface equilibrium concentration of aerosol species
!------------------------------------------------------------------------
    implicit none
    double precision GREAT
    parameter (GREAT=100.D0)

    integer nesp_aer
    double precision qext(nesp_aer),qinti(N_inside_aer)
    double precision surface_equilibrium_conc(nesp_aer)
    double precision lwc,ionic,proton

    integer jesp,j
    double precision wi(nesp_isorropia),w(nesp_isorropia)
    double precision aerliq(N_liquid),aersld(N_solid)
    double precision gas(3),cntrl(2)
    double precision other(6),coef
    !      CHARACTER*15 scase
    double precision organion,watorg

    organion = 0.D0
    watorg = 0.D0
!     Inputs  for Isorropia
    gas(1)=0.d0
    gas(2)=0.d0
    gas(3)=0.d0
    cntrl(1)=1.D0             ! reverse mode
    cntrl(2)=1.D0            ! metastable option
      ! convert µg to moles
!     conversion unit for isorropia needed in mol.m-3
    do j=1,nesp_isorropia
      jesp=isorropia_species(j)
      wi(j)= qext(jesp)& ! µg.m-3
	  /molecular_weight_aer(jesp)!&  ! µg.mol-1
    end do

#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
      wi(1) = 0.D0              !Do not consider sea salt in isoropia
      wi(5) = 0.D0
#endif

    call ISOROPIA(wi,Relative_Humidity,Temperature,cntrl,w,gas,&
	aerliq,aersld,other,organion,watorg)
!     clipping to tinym
    if (gas(1).lt.0.d0) gas(1)=tinym
    if (gas(2).lt.0.d0) gas(2)=tinym
    if (gas(3).lt.0.d0) gas(3)=tinym
!     Aqueous phase total liquid water content and pH (proton) concentration
    lwc = aerliq(IH2O) * imw(IH2O) ! microg.m-3 

    ionic = other(5)
    proton = aerliq(IH) * imw(IH) !* gammaH  ! microg.m-3 but equivalent to micromol.m-3


!     Outputs isorropia
			      ! sulfate surf conc always 0. µg.m-3
    surface_equilibrium_conc(ESO4)=0.D0

			      ! convert moles.m-3 to µg.m-3
    surface_equilibrium_conc(ENH4)=gas(1)*molecular_weight_aer(ENH4)
    ! write(*,*) "NH4 in eqinorg after isorropia",gas(1), &
    !      molecular_weight_aer(ENH4), qext(ENH4),wi(3)
    surface_equilibrium_conc(ENO3)=gas(2)*molecular_weight_aer(ENO3)
    surface_equilibrium_conc(ECl) =gas(3)*molecular_weight_aer(ECl)
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
    surface_equilibrium_conc(ECl) = 0.d0 
#endif
    
			      ! liquid inorg aerosol
    do jesp=IH,IOH
      qinti(jesp)= DMAX1(aerliq(jesp),0.D0)*imw(jesp)  ! moles to µg

    end do
			      ! liquid water content
    qext(EH2O)= qinti(IH2O)&
	+qinti(IOH)*1.05882352941D0 ! mwh2o/mwioh
			! solid inorg aerosol
    do jesp=SNaNO3,SLC
      qinti(jesp)= DMAX1(aersld(jesp-12),0.D0)&
	  *smw(jesp)        ! moles to µg  
     !molecular_weight_solid(jesp)
    end do

   end subroutine EQINORG

   subroutine KLIMIT(q,c_gas,k,ce_kernal_coef)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine limits the condensation/evaporation rates for
!     aerosol and gases in order to avoid clippings.
!     Two kinds of limitation are performed:
!     - the 1st is aerosol clipping : as it may reduce evaporation then
!     enlarge condensation;, it is done before condensation limitation.
!     - the 2nd is gas clipping : in practice it may reduce, per species,
!     aerosol condensation only in bins that leads to gas clipping.
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     q: aerosol mass concentration(µg/m^3)
!     c_gas: aerosol gas concentration(µg/m^3)
!     ce_kernal_coef: c/e kernel coefficient          ([m3.s-1]).
!
!     -- OUTPU VARIABLES
!
!     k: particle mass derivation(µg/m^3/s)
!------------------------------------------------------------------------
    implicit none
    double precision:: q(N_size,N_aerosol)!1th order mass concentration
    double precision:: k(N_size,N_aerosol)
    double precision:: c_gas(N_aerosol)
    double precision:: ce_kernal_coef(N_size,N_aerosol)

    integer:: jesp,j,s
    double precision:: ksum,klim,ktlim
    double precision:: qsum,ce_kernal_coef_tot
    double precision:: frac,qnew

!     ****** prevent aerosol clipping
    do s=1,(N_aerosol-1)
      jesp=List_species(s)
      do j=(ICUT+1),N_size
	  ! only when evaporation
	if (k(j,jesp).lt.0.D0) then
	  ! if q(j,s) is <=TINYM or =0
	  ! then k should be >=0, but
	  ! due to bad matrix inversion
	  ! this case may occur
	  if (q(j,jesp).lt.TINYM) k(j,jesp)=0.D0
	  ! test clipping in other cases
	  qnew=q(j,jesp)+k(j,jesp)
	  if (qnew.lt.0.D0) then
	    ! we are sure that q>=tinym
	    ! otherwise k would be = 0
	    ! from previous case
	    k(j,jesp)=(TINYM-q(j,jesp)) !/timestep_splitting ! <=0 µg.m-3.s-1
	    ! we force q to be a
	    ! 'little' more than TINYM
	    k(j,jesp)=0.99D0*k(j,jesp)
	  endif
	endif
      enddo
    enddo
!     ****** prevent gas clipping
    do s=1,(N_aerosol-1)
      jesp=List_species(s)
      ! compute total mass rate per species
      ksum=0.D0
      do j=ICUT+1,N_size
	  ksum=ksum+k(j,jesp)     !µg.m-3.s-1
      enddo
      ! we perform limiting in
      ! case of condensation only
      if(ksum.gt.0.D0) then
	! this is the total lumped mass
	! to perserve from clipping
	qsum=c_gas(jesp)
	do j=1,ICUT
	    qsum=qsum+q(j,jesp)  ! µg.m-3
	enddo
	! test if clipping occurs
	! then perform the limitation
	if (ksum.gt.qsum) then
	    ! sum of ce_kernal_coef(*) c/e coefficient
	  ce_kernal_coef_tot=0.D0
	  do j= 1, N_size
	    ce_kernal_coef_tot=ce_kernal_coef_tot+ce_kernal_coef(j,jesp)
	  enddo
	 ! tot max rate, µg.m-3.s-1
	  ktlim=qsum       !/timestep_splitting
	  do j=ICUT+1,N_size
	    ! fraction, adim
            if (ce_kernal_coef_tot .ne. 0.d0) then !! YK
               frac=ce_kernal_coef(j,jesp)/ce_kernal_coef_tot ! ce_kernal_coef_tot != 0
            else
               frac = 0.d0
            endif !! YK
	    ! we allow only a given fraction
	    ! of ktlim to condense on given bin
	    klim=ktlim*frac
	    ! apply the limit
	    ! only if necessary
	    if (k(j,jesp).gt.klim) k(j,jesp)=klim
	  enddo
	endif
      endif
    enddo

  end subroutine KLIMIT


  subroutine HPLFLIM(alfa,qih,N_size_loc,init_bulk_gas,&
   ce_kernal_coef_i,Kelvin_effect,surface_equilibrium_conc,ce_kernel)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine computes the flux limitation for the
!     condensation/evaporation flux. The algorithm is based
!     on the limitation of the aerosol acidity rate.
!
!     The details may be found in the PhD Work of Edouard Debry,
!     Chapter 10 (section 10.1.5) or in the reference:
!     Pilinis et al: MADM, a new multicomponent aerosol dynamic model
!     Aerosol Science and Technology 32, 482:502, 2000.
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     alfa: percentage of H+ allowed to c/e(0.1)
!     qih: int H+ conc (µg)
!     N_size_loc: size of vectors following below
!     init_bulk_gas: bulk gas conc (µg.m-3)
!     surface_equilibrium_conc: surface equilibrium concentration of aerosol species(µg/m^3)
!     Kelvin_effect:kelvin coef (adim)
!     ce_kernal_coef_i: c/e kernel coefficient          ([m3.s-1]).
!
!     -- OUTPU VARIABLES
!
!     ce_kernel: modified c/e kernel (µg/m^3)
!------------------------------------------------------------------------
    implicit none
    integer N_size_loc
    double precision init_bulk_gas(N_size_loc),ce_kernal_coef_i(N_size_loc)
    double precision Kelvin_effect(N_size_loc),surface_equilibrium_conc(N_size_loc)
    double precision ce_kernel(N_size_loc),qih,alfa

    integer jesp,s
    double precision maa(N_size_loc),mkercd(N_size_loc)
    double precision cfa,cfb,cfc,cc
    double precision mih,melec,mlim,q

!     Compute mol fluxes

    do s=1,nesp_isorropia
      jesp=isorropia_species(s)
      ! maa(*) in m3.mol.s-1.µg-1
      maa(jesp)= ce_kernal_coef_i(jesp)&   ! m3.s-1
	  /molecular_weight_aer(jesp)        ! µg.mol-1
      ! mkercd(*) in mol.s-1
      mkercd(jesp)= ce_kernel(jesp)& ! µg.s-1
	  /molecular_weight_aer(jesp)        ! µg.mol-1
    end do

!     H+ limitation

    mih=qih/imw(IH)           ! mol of H+ in aerosol


			      ! maximum of mih variation tolerated
    mlim=mih*alfa             ! mol.s-1 alfa=(0.1)

			      ! electroneutrality  ! mol.s-1
    melec= 2.D0*mkercd(ESO4)+mkercd(ENO3)&
	+mkercd(ECl)-mkercd(ENH4)

    ! correction factor default value
    cc=0.D0
    ! correction calculation
    if (DABS(melec).gt.mlim) then
      ! we give to mlim the sign of melec
      mlim=DSIGN(mlim,melec)! returns the value of mlim with the sign of melec
      !to judge + or -
      ! cfa,cfb,cfc are coefficients of
      ! 2nd order eq : cfa*cc^2+cfb*cc+cfc=0
      ! satisfied by the correction factor cc
      cfa= maa(ENO3)*surface_equilibrium_conc(ENO3)*Kelvin_effect(ENO3)&
	  +maa(ECl)*surface_equilibrium_conc(ECl)*Kelvin_effect(ECl)
      cfb=mlim-2.D0*maa(ESO4)*init_bulk_gas(ESO4)& ! mol.s-1
	  -maa(ENO3)*init_bulk_gas(ENO3)&
	  -maa(ECl) *init_bulk_gas(ECl)&
	  +maa(ENH4)*init_bulk_gas(ENH4)
      cfc=-maa(ENH4)*surface_equilibrium_conc(ENH4)*Kelvin_effect(ENH4)
      ! one can note cfa>=0 and cfc<=0
      ! then there always exist a + but
      ! possibly zero root
      ! root computation
      if (cfa.gt.0.D0) then
	if (cfb*cfb-4.D0*cfa*cfc.lt.0.D0) then
	    WRITE(6,*)'(hplflim.f): sqrt(<0)'
	    WRITE(6,*)cfb,cfa,cfc,mlim,'Time',current_sub_time
	    WRITE(6,*)'ASO4',ce_kernal_coef_i(ESO4),'ANH4',ce_kernal_coef_i(ENH4),'ANO3',&
	    ce_kernal_coef_i(ENO3),'ACl',ce_kernal_coef_i(ECl)
	    WRITE(6,*)'SO4',mkercd(ESO4),'NH4',mkercd(ENH4),&
	    'NO3',mkercd(ENO3),'Cl',mkercd(ECl)
	    WRITE(6,*)'SO4',maa(ESO4),init_bulk_gas(ESO4)
	    WRITE(6,*)'NH4',maa(ENH4),init_bulk_gas(ENH4),surface_equilibrium_conc(ENH4)
	    WRITE(6,*)'NO3',maa(ENO3),init_bulk_gas(ENO3),surface_equilibrium_conc(ENO3)
	    WRITE(6,*)'Cl',maa(ECl),init_bulk_gas(ECl),surface_equilibrium_conc(ECl)
	    STOP
	endif
	q=-5.D-01*(cfb+DSIGN(1.D0,cfb)&
	    *DSQRT(cfb*cfb-4.D0*cfa*cfc))
	cc=DMAX1(q/cfa,cfc/q) ! we select the + root
      else
	if (cfb.ne.0.D0) cc=-cfc/cfb
      endif
    endif

!     A correction is done if only upper calculation
!     root has changed cc to a strictly positive
!     value, otherwise it is considered as non stiff cases

    if (cc.gt.0.D0) then
	surface_equilibrium_conc(ENH4)=surface_equilibrium_conc(ENH4)/cc
	surface_equilibrium_conc(ENO3)=surface_equilibrium_conc(ENO3)*cc
	surface_equilibrium_conc(ECl) =surface_equilibrium_conc(ECl)*cc
	do jesp=ENH4,ECl
	  ce_kernel(jesp)= ce_kernal_coef_i(jesp)*(init_bulk_gas(jesp)-surface_equilibrium_conc(jesp)*Kelvin_effect(jesp))
	end do
    endif

  end subroutine HPLFLIM


  ! YK 
  subroutine equi_const(TEMP, XK3, XK4, XK6, XK8, XK9, XK10)
    implicit none
    double precision TEMP, T0, T0T, COEF
    double precision XK3, XK4, XK6, XK8, XK9, XK10

      XK3  = 1.971D6   ! HCL(g)           <==> H(aq)     + CL(aq)
      XK4  = 2.511e6   ! HNO3(g)          <==> H(aq)     + NO3(aq) ! ISORR
      XK6  = 1.086D-16 ! NH4CL(s)         <==> NH3(g)    + HCL(g)
      XK8  = 37.661D0  ! NACL(s)          <==> NA(aq)    + CL(aq)
      XK9  = 11.977D0  ! NANO3(s)         <==> NA(aq)    + NO3(aq)
      XK10 = 5.746D-17 ! NH4NO3(s)        <==> NH3(g)    + HNO3(g) ! ISORR

      IF (INT(TEMP) .NE. 298) THEN   ! FOR T != 298K or 298.15K
         T0  = 298.15D0
         T0T = T0/TEMP
         COEF= 1.0+LOG(T0T)-T0T

         XK3 = XK3 *EXP( 30.20*(T0T-1.0) + 19.910*COEF)
         XK4 = XK4 *EXP( 29.17*(T0T-1.0) + 16.830*COEF) !ISORR
         XK6 = XK6 *EXP(-71.00*(T0T-1.0) +  2.400*COEF)
         XK8 = XK8 *EXP( -1.56*(T0T-1.0) + 16.900*COEF)
         XK9 = XK9 *EXP( -8.22*(T0T-1.0) + 16.010*COEF)
         XK10= XK10*EXP(-74.38*(T0T-1.0) +  6.120*COEF) ! ISORR
      ENDIF

  end subroutine equi_const
  
  subroutine DRYIN(Temp,qinti,N_size_loc,init_bulk_gas,ce_kernal_coef_i,&
   Kelvin_effect,surface_equilibrium_conc,ce_kernel)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine computes the aerosol surface gas-phase
!     concentration (through equilibrium) for dry aerosols.
!
!     The algorithms are detailed in Chapter 10 (section 10.1.6) of
!     the PhD work of Edouard Debry.
!     See also the reference:
!     Pilinis et al: MADM, a new multicomponent aerosol dynamic model
!     Aerosol Science and Technology 32, 482:502, 2000.
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     qih: int H+ conc (µg)
!     N_size_loc: size of vectors following below
!     init_bulk_gas: bulk gas conc (µg.m-3)
!     surface_equilibrium_conc: surface equilibrium concentration of aerosol species(µg/m^3)
!     Kelvin_effect:kelvin coef (adim)
!     ce_kernal_coef_i: c/e kernel coefficient          ([m3.s-1]).
!
!     -- OUTPU VARIABLES
!
!     ce_kernel: modified c/e kernel (µg/m^3)
!------------------------------------------------------------------------
! ISORROPIA commons needed by this routine,
! directly taken from isrpia.inc.
    implicit none
    double precision XK3,XK4,XK6,XK8,XK9,XK10
    ! COMMON /EQUK/ XK1,XK2,XK3,XK4,XK5,XK6,XK7,XK8,&
    !     	  XK9,XK10,XK11,XK12,XK13,XK14,&
    !     	  XKW,XK21,XK22,XK31,XK32,XK41,XK42 ! YK
!!!$OMP THREADPRIVATE(/EQUK/)
    integer N_size_loc
    double precision qinti(N_inside_aer),Kelvin_effect(N_size_loc)
    double precision ce_kernal_coef_i(N_size_loc),init_bulk_gas(N_size_loc)
    double precision surface_equilibrium_conc(N_size_loc),ce_kernel(N_size_loc)
    double precision Temp
    LOGICAL leq1,leq2,leq3,lr47,lr56
    integer icase,jesp,s
    double precision rgas1,maa(N_size_loc)
    double precision rk1,rk2,rk3
    double precision sat1,sat2,sat3a,sat3b
    double precision mkercd(N_size_loc),msat
    double precision cfa,cfb,cfc

!     Initialization:
!     1 stands for nh4no3 equilibrium
!     2 stands for nh4cl equilibrium
!     3 stands for nacl/nano3 equilibrium
!     47 for nacl and nh4cl reactions
!     56 for nano3 and nh4no3 reactions

    leq1=.false.
    leq2=.false.
    leq3=.false.
    lr47=.false.
    lr56=.false.

    rgas1=ATM/(RGAS*Temp)

    ! YK
    call equi_const(Temp, XK3, XK4, XK6, XK8, XK9, XK10)

    rk1= XK10*rgas1*rgas1*molecular_weight_aer(ENH4)&
    *molecular_weight_aer(ENO3)*Kelvin_effect(ENH4)*Kelvin_effect(ENO3)
			      ! rk1 in (µg.m-3)2

    rk2=XK6*rgas1*rgas1*molecular_weight_aer(ENH4)&
    *molecular_weight_aer(ECl)*Kelvin_effect(ENH4)*Kelvin_effect(ECl)
			      ! rk2 in (µg.m-3)2

    rk3=XK4*XK8/(XK3*XK9) ! adim

    sat1=init_bulk_gas(ENH4)*init_bulk_gas(ENO3) ! (µg.m-3)2
    sat2=init_bulk_gas(ENH4)*init_bulk_gas(ECl) ! (µg.m-3)2

    sat3a=init_bulk_gas(ECl)          ! µg.m-3
    sat3b=rk3*init_bulk_gas(ENO3)     ! µg.m-3

    do s=2,nesp_isorropia
      jesp=isorropia_species(s)
			    ! maa(*) in m3.mol.s-1.µg-1
      maa(jesp)= ce_kernal_coef_i(jesp)&   ! m3.s-1
	  /molecular_weight_aer(jesp)        ! µg.mol-1

			    ! mkercd(*) in mol.s-1
      mkercd(jesp)= ce_kernel(jesp)& ! µg.s-1
	  /molecular_weight_aer(jesp)        ! µg.mol-1
    end do

    msat=2.D0*mkercd(ESO4)-mkercd(ENH4) ! mol.s-1

!     Determine which reaction is active

    if (qinti(SNH4NO3).gt.0.D0.OR.sat1.gt.rk1) then
      leq1=.true.
    endif

    if (qinti(SNH4Cl).gt.0.D0.OR.sat2.gt.rk2) then
      leq2=.true.
    endif

    if (qinti(SNaNO3).gt.0.D0.AND.qinti(SNaCl).gt.0.D0) then
      leq3=.true.
    endif

    if (qinti(SNaNO3).gt.0.D0.AND. sat3a.gt.sat3b) then
      leq3=.true.
    endif

    if (qinti(SNaCl).gt.0.D0.AND.sat3a.lt.sat3b) then
      leq3=.true.
    endif

    if (init_bulk_gas(ESO4).gt.0.D0) then
      if (qinti(SNaCl).gt.0.D0.OR.qinti(SNH4Cl).gt.0.D0) then
	lr47=.true.
      endif

      if (qinti(SNaNO3).gt.0.D0.OR.qinti(SNH4NO3).gt.0.D0) then
	lr56=.true.
      endif
    endif

    if (leq1.AND.leq2.AND.leq3) then
      PRINT *,'Warning from dryin.f: << solid : leq123 >>'
    endif

!     Determine which case is relevant

    if (leq2.AND.leq3) then
	icase=1                ! R2 and R3 active
    elseif (leq1.AND.leq2) then
	icase=2                ! R1 and R2 active
    elseif (leq1.AND.leq3) then
	icase=3                ! R1 and R3 active
    elseif (leq1) then
	icase=4                ! only R1 active
    elseif (leq2) then
	icase=5                ! only R2 active
    elseif (leq3) then
	icase=6                ! only R3 active
    else                      ! no active equilibrium
      if (lr47) then
	icase=7             ! R4 or R7 active
      elseif (lr56) then
	icase=8             ! R5 or R6 active
      else                   ! nothing active

	if (msat.lt.0.D0) then
	    icase=9          ! enough nh3 to neutralize so4
	else
	    icase=10         ! not enough nh3 to neutralize so4
			    ! in this case aerosol become acidic
	endif
      endif
    endif
!     Solve each case
			      ! icase 3 is not physical but used
			      ! to determine the real icase
    if (icase.eq.3) then
      cfa=( maa(ENO3)*Kelvin_effect(ENO3)+maa(ECl)*Kelvin_effect(ECl)*rk3 )

      cfb=2.D0*mkercd(ESO4)+mkercd(ENO3)+mkercd(ECl)-mkercd(ENH4)

      cfc=rk1*maa(ENH4)*Kelvin_effect(ENH4)

      if (cfb*cfb+4.D0*cfa*cfc.le.0.D0) then
	WRITE(6,*)'(dryin.f): (1) sqrt(<0) '
	STOP
      endif
      surface_equilibrium_conc(ENO3)= (cfb+DSQRT(cfb*cfb+4.D0*cfa*cfc)) /(2.D0*cfa)
      !surface_equilibrium_conc(ECl)=rk3*surface_equilibrium_conc(ENO3)
      !surface_equilibrium_conc(ENH4)=rk1/surface_equilibrium_conc(ENO3)
      ce_kernel(ENO3)= ce_kernal_coef_i(ENO3) *( init_bulk_gas(ENO3)&
      -surface_equilibrium_conc(ENO3)*Kelvin_effect(ENO3) )
      ! test if no3 condenses
      if (ce_kernel(ENO3).gt.0.D0) then
	icase=2             ! if nh4no3 forms then real case=2
      else
	icase=1             ! if nacl forms then real case=1
      endif
    endif

			      ! other cases
    if (icase.eq.1) then
      cfa=( maa(ENO3)*Kelvin_effect(ENO3)+maa(ECl)*Kelvin_effect(ECl)*rk3 )
      cfb=2.D0*mkercd(ESO4)+mkercd(ENO3)+mkercd(ECl)-mkercd(ENH4)
      cfc=rk2/rk3*maa(ENH4)*Kelvin_effect(ENH4)
      if (cfb*cfb+4.D0*cfa*cfc.le.0.D0) then
	WRITE(6,*)'(dryin.f): (2) sqrt(<0) '
	STOP
      endif
      surface_equilibrium_conc(ENO3)= (cfb+DSQRT(cfb*cfb+4.D0*cfa*cfc))/(2.D0*cfa)
      surface_equilibrium_conc(ECl)=rk3*surface_equilibrium_conc(ENO3)
      surface_equilibrium_conc(ENH4)=rk2/rk3/surface_equilibrium_conc(ENO3)
    elseif (icase.eq.2) then
      cfa=( maa(ENO3)*Kelvin_effect(ENO3)+maa(ECl)*Kelvin_effect(ECl)*rk2/rk1 )
      cfb=2.D0*mkercd(ESO4)+mkercd(ENO3)+mkercd(ECl)-mkercd(ENH4)
      cfc=rk1*maa(ENH4)*Kelvin_effect(ENH4)
      if (cfb*cfb+4.D0*cfa*cfc.le.0.D0) then
	WRITE(6,*)'(dryin.f): (3) sqrt(<0) '
	STOP
      endif
      surface_equilibrium_conc(ENO3)= (cfb+DSQRT(cfb*cfb+4.D0*cfa*cfc))/(2.D0*cfa)
      surface_equilibrium_conc(ENH4)=rk1/surface_equilibrium_conc(ENO3)
      surface_equilibrium_conc(ECl)=rk2/rk1*surface_equilibrium_conc(ENO3)
    elseif (icase.eq.4) then
      cfa=maa(ENO3)*Kelvin_effect(ENO3)
      cfb= 2.D0*mkercd(ESO4)+mkercd(ENO3)-mkercd(ENH4)
      cfc=rk1*maa(ENH4)*Kelvin_effect(ENH4)
      if (cfb*cfb+4.D0*cfa*cfc.le.0.D0) then
	WRITE(6,*)'(dryin.f): (4) sqrt(<0) '
	STOP
      endif
      surface_equilibrium_conc(ENO3)= (cfb+DSQRT(cfb*cfb+4.D0*cfa*cfc))/(2.D0*cfa)
      surface_equilibrium_conc(ENH4)=rk1/surface_equilibrium_conc(ENO3)
      surface_equilibrium_conc(ECl)=init_bulk_gas(ECl)/Kelvin_effect(ECl)
    elseif (icase.eq.5) then
      cfa=maa(ENH4)*Kelvin_effect(ENH4)
      cfb= 2.D0*mkercd(ESO4)+mkercd(ECl)-mkercd(ENH4)
      cfc=maa(ECl)*rk2*Kelvin_effect(ECl)
      if (cfb*cfb+4.D0*cfa*cfc.le.0.D0) then
	WRITE(6,*)'(dryin.f): (5) sqrt(<0) '
	STOP
      endif
      surface_equilibrium_conc(ENH4)= (cfb+DSQRT(cfb*cfb+4.D0*cfa*cfc))/(2.D0*cfa)
      surface_equilibrium_conc(ECl)=rk2/surface_equilibrium_conc(ENH4)
      surface_equilibrium_conc(ENO3)=init_bulk_gas(ENO3)/Kelvin_effect(ENO3)
    elseif (icase.eq.6) then
      cfa=( maa(ENO3)*Kelvin_effect(ENO3)+rk3*maa(ECl)*Kelvin_effect(ECl) )
      cfb=2.D0*mkercd(ESO4)+mkercd(ECl)+mkercd(ENO3)
      surface_equilibrium_conc(ENO3)=cfb/cfa
      surface_equilibrium_conc(ECl)=rk3*surface_equilibrium_conc(ENO3)
      surface_equilibrium_conc(ENH4)=init_bulk_gas(ECl)/Kelvin_effect(ECl)
    elseif (icase.eq.7) then
      surface_equilibrium_conc(ENH4)=init_bulk_gas(ENH4)/Kelvin_effect(ENH4)
      surface_equilibrium_conc(ENO3)=init_bulk_gas(ENO3)/Kelvin_effect(ENO3)

      ce_kernel(ENH4)=0.D0
      ce_kernel(ENO3)=0.D0

      ce_kernel(ECl)=-molecular_weight_aer(ECl)*2.D0*mkercd(ESO4)
      surface_equilibrium_conc(ECl)= ( init_bulk_gas(ECl)-ce_kernel(ECl)&
      /ce_kernal_coef_i(ECl) )/Kelvin_effect(ECl)
    elseif (icase.eq.8) then
      surface_equilibrium_conc(ENH4)=init_bulk_gas(ENH4)/Kelvin_effect(ENH4)
      surface_equilibrium_conc(ECl)=init_bulk_gas(ECl)/Kelvin_effect(ECl)
      ce_kernel(ENH4)=0.D0
      ce_kernel(ECl)=0.D0
      ce_kernel(ENO3)=-molecular_weight_aer(ECl)*2.D0 *mkercd(ESO4)
      surface_equilibrium_conc(ENO3)= ( init_bulk_gas(ENO3)-ce_kernel(ENO3)&
      /ce_kernal_coef_i(ENO3) )/Kelvin_effect(ENO3)
    elseif (icase.eq.9) then
      surface_equilibrium_conc(ENO3)=init_bulk_gas(ENO3)/Kelvin_effect(ENO3)
      surface_equilibrium_conc(ECl)=init_bulk_gas(ECl)/Kelvin_effect(ECl)
      ce_kernel(ENO3)=0.D0
      ce_kernel(ECl)=0.D0
      ce_kernel(ENH4)= molecular_weight_aer(ENH4)*2.D0 *mkercd(ESO4)
      surface_equilibrium_conc(ENH4)= ( init_bulk_gas(ENH4)-ce_kernel(ENH4)&
      /ce_kernal_coef_i(ENH4) )/Kelvin_effect(ENH4)
    elseif (icase.eq.10) then
      surface_equilibrium_conc(ENO3)=init_bulk_gas(ENO3)/Kelvin_effect(ENO3)
      surface_equilibrium_conc(ECl)=init_bulk_gas(ECl)/Kelvin_effect(ECl)
      ce_kernel(ENO3)=0.D0
      ce_kernel(ECl)=0.D0
!    no more electroneutrality in this case
    endif
!     Giving out the kernel for case <=6
    if (icase.LE.6) then
      do s=3,nesp_isorropia
	jesp=isorropia_species(s)
	ce_kernel(jesp)= ce_kernal_coef_i(jesp)*( init_bulk_gas(jesp)-&
	surface_equilibrium_conc(jesp)*Kelvin_effect(jesp) )
      end do
    endif
  end subroutine DRYIN

  subroutine POA_DRV(nesp_aer,aero, gas, &
   soa_part_coef, Temp, vaporization_enthalpy)

!------------------------------------------------------------------------
!
!     -- DES!RIPTION
!
!     This subroutine computes the equilibrium between gas and particle
!     phase using an absorption partioning model (Poa, 1994a, 1994b)
!     for organic species which are not managed by AEC model (Pun et al 2001).
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     Nesp_aer: Number of species
!     Temperature: temperature ([Kelvin]).
!     soa_part_coef: partition coefficient ([m^3.\mu g^-1]).
!     vaporization_enthalpy: vaporization enthalpy
!
!     -- INPUT/OUTPUT VARIABLES
!
!     AERO: aerosol bulk concentration ([\mu g.m^-3]).
!     GAS: gas concentration ([\mu g.m^-3]).
!
!     -- OUTPUT VARIABLES
!
!------------------------------------------------------------------------

    IMPLICIT NONE

    integer nesp_aer
    double precision aero(nesp_aer), gas(nesp_aer), soa_part_coef(nesp_aer)
    double precision vaporization_enthalpy(nesp_aer)
    double precision ctot(nesp_pom), caer(nesp_pom)
    double precision cgas(nesp_pom), kpart2(nesp_pom)
    double precision emw2(nesp_pom), Temp
    integer i,j
    double precision totom
    double precision paom,soam

!     Fill concentration vectors
    paom=0.0
    soam=0.0
    do i = 1,nesp_pom
      j = poa_species(i)
      caer(i) = aero(j)
      cgas(i) = gas(j)
      emw2(i) = molecular_weight_aer(j)
      paom=paom+caer(i)
      !print*,aero(j),gas(j),molecular_weight_aer(j),paom
    enddo

    if(paom<0.1) then
      paom=0.1
    endif

! Set secondary organic mass quantity in microg/m3
    do i = N_hydrophilic+1,nesp_aec ! only dry AEC species
      j = aec_species(i)
      soam = soam + aero(j)
    enddo

    do i = 1,nesp_pankow
      j = pankow_species(i)
      soam = soam + aero(j)
    enddo

! aero in microg/m3
    do i = 1,nesp_pom
      j = poa_species(i)
      kpart2(i) = soa_part_coef(j)
    enddo

    do i = 1,nesp_pom
      ctot(i) = caer(i) + cgas(i)
    enddo

    do j=1,NITER_POA
      totom = paom+soam
      paom=0
      do i=1,nesp_pom
	if(kpart2(i)*totom.ne.-1.d0) then
	  caer(i) = ctot(i)*kpart2(i)*totom/(1+kpart2(i)*totom)
	endif
	paom=paom+caer(i)
      enddo
    enddo
    do j = 1,nesp_pom
      cgas(j) = ctot(j) - caer(j)
    enddo

    do i = 1,nesp_pom
      j = poa_species(i)
      aero(j) = caer(i)
      gas(j) = cgas(i)
    enddo

  end subroutine POA_DRV

  subroutine ISOROPIA_DRV(nesp_aer,&
          aero, gas, organion, watorg, ionic, proton, lwc, rh, Temp, &
          liquid)

!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This subroutine computes the equilibrium between inorganic aerosols
!     and gas-phase (forward mode), taking in account organic liquid
!     water content and organic ions.
!     It calls ISORROPIA by Nenes et al.
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     ORGANION: organic ions ([\mu mol.m^-3]).
!     WATORG: organic liquid water content ([\mu g.m^-3]).
!     Relative_Humidity: relative humidity 0< <1 ([]).
!     Temperature: temperature ([Kelvin]).
!
!     -- INPUT/OUTPUT VARIABLES
!
!     AERO: aerosol bulk concentration ([\mu g.m^-3]).
!     GAS: gas concentration ([\mu g.m^-3]).
!
!     -- OUTPUT VARIABLES
!
!     PROTON: hydronium ion concentration ([\mu g.m^-3]).
!     LWC: total liquid water content ([\mu g.m^-3]).
!------------------------------------------------------------------------

    IMPLICIT NONE

    integer nesp_aer
    double precision aero(nesp_aer), gas(nesp_aer)
    double precision organion, watorg, proton
    double precision lwc, rh, Temp
    double precision wi(N_inorganic),w(N_inorganic),gas2(3),cntrl(2), other(6)
    double precision liquid(N_liquid),solid(N_solid)
    double precision organion2, watorg2, ionic ,gammaH
    integer i,idx !,j
!     mol neg charge in mol.m-3 */
    organion2 = organion * 1.D-6
    liquid=0.d0
    solid=0.d0
    w=0.d0
    gas2=0.d0
    other=0.d0
!     organic water content converted from
!     microg/m3 (aec output) to kg/m3 (isorropia input)
    watorg2 = watorg * 1.D-9
    cntrl(1) = 0.D0
    cntrl(2) = MTSBL
!     concentration in microg.m-3
    gas(ENa) = 0.D0

    do i=1,nesp_isorropia
	idx = isorropia_species(i)
	wi(i) = aero(idx) + gas(idx)
    enddo
!     conversion unit for isorropia needed in mol.m-3
    do i=1,nesp_isorropia
	idx = isorropia_species(i)
	wi(i) = wi(i) / molecular_weight_aer(idx) ! microg.m-3 / microg.mol-1 = mol.m-3
    enddo
    
    call ISOROPIA(wi, rh, Temp, cntrl, w, gas2,&
	liquid, solid, other, organion2, watorg2)

!     Aqueous phase total liquid water content and pH (proton) concentration
    lwc = liquid(IH2O) * imw(IH2O) ! microg.m-3 
    ionic = other(5)
!    gammaH = 10**(-0.511 * (298.0/Temperature)**1.5 * sqrt(ionic)/(1+sqrt(ionic)))
!! YK gammaH will be calculated later in SOAP

    proton = liquid(IH) * imw(IH) !* gammaH  ! microg.m-3 but equivalent to micromol.m-3
!    print*,'lwc=',lwc,'proton=',proton,'SO4',w(2),'NH4',w(3) - gas2(1),&
!	'NO3',w(4) - gas2(2),'Cl',w(5) - gas2(3)

    gas(ESO4) = 0.D0
    do i=3,nesp_isorropia
	idx = isorropia_species(i)
	gas(idx) = gas2(i-2) * molecular_weight_aer(idx)
    enddo

    aero(ESO4) = w(2) * molecular_weight_aer(ESO4)
    do i=3,nesp_isorropia
	idx = isorropia_species(i)
	aero(idx) = (w(i) - gas2(i-2)) * molecular_weight_aer(idx)
    enddo

  end subroutine ISOROPIA_DRV

  subroutine PANKOW_DRV(nesp_aer,aero, gas, soa_part_coef)

!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This subroutine computes the equilibrium between gas and particle
!     phase using an absorption partioning model (Pankow, 1994a, 1994b)
!     for organic species which are not managed by AEC model (Pun et al 2001).
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     soa_part_coef: partition coefficient ([m^3.\mu g^-1]).
!
!     -- INPUT/OUTPUT VARIABLES
!
!     AERO: aerosol bulk concentration ([\mu g.m^-3]).
!     GAS: gas concentration ([\mu g.m^-3]).
!
!     -- OUTPUT VARIABLES
!------------------------------------------------------------------------
    IMPLICIT NONE

    integer nesp_aer
    double precision aero(nesp_aer), gas(nesp_aer), soa_part_coef(nesp_aer)
    double precision ctot(nesp_pankow), caer(nesp_pankow)
    double precision cgas(nesp_pankow), kpart2(nesp_pankow)
    double precision emw2(nesp_pankow)
    integer i,j
    double precision totmol, totmol2
    double precision a, b, c, deter, q, paom

!     Fill concentration vectors
    do i = 1,nesp_pankow
	j = pankow_species(i)
	caer(i) = aero(j)
	cgas(i) = gas(j)
	emw2(i) = molecular_weight_aer(j)
    enddo

! Set primary organic molar quantity in mol/m3
    paom=0.0
    do i = 1,nesp_pom
	j = poa_species(i)
	paom=paom+aero(j)/molecular_weight_aer(j)
    enddo
    do i = N_hydrophilic+1,nesp_aec ! only dry AEC species
	j = aec_species(i)
	paom = paom + aero(j) / molecular_weight_aer(j)
    enddo

! aero in microg/m3, molecular_weight_aer in microg/mol
    do i = 1,nesp_pankow
	j = pankow_species(i)
	kpart2(i) = soa_part_coef(j)
    enddo

    do i = 1,nesp_pankow
	ctot(i) = caer(i) + cgas(i)
    enddo

    do j=1,NITER_PKW
	totmol = paom
	do i=1,nesp_pankow
	  totmol = totmol + caer(i) / emw2(i)
	enddo
	do i=1,nesp_pankow
	  totmol2 = totmol - caer(i) / emw2(i)
	  a = 1.D0 / emw2(i)
	  b = (1.D0 / kpart2(i) - ctot(i)) / emw2(i) + totmol2
	  c = - ctot(i) * totmol2
	  deter = b * b - 4.D0 * a * c
	  if (deter.lt.0.d0) STOP 'pankow_drv.f: deter < 0'
	  q= - 0.5D0 * ( b + DSIGN(1.D0,b) * DSQRT(deter))
	  caer(i) = DMAX1(q / a, c / q)
	enddo
    enddo
    do j = 1,nesp_pankow
	cgas(j) = ctot(j) - caer(j)
    enddo
    do i = 1,nesp_pankow
	j = pankow_species(i)
	aero(j) = caer(i)
	gas(j) = cgas(i)
    enddo

  end subroutine PANKOW_DRV
! imw  molecular_weight_inside
! smw  molecular_weight_solid
! solmd  mass_density_solid
End module cThermodynamics
