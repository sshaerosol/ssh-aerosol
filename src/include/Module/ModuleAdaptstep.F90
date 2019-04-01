!!-----------------------------------------------------------------------
!!     Copyright (C) 20012-2018, ENPC - EDF R&D - INERIS
!!     Author(s): Shupeng Zhu, Karine Sartelet, Edouard Debry
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
!!-----------------------------------------------------------------------
!!
!!     -- DESCRIPTION
!!    This module contains numerical solves for time integration
!!-----------------------------------------------------------------------
Module jAdaptstep
  use aInitialization
  use bCoefficientRepartition
  use cThermodynamics
  use dPhysicalbalance
  use eRedistribution
  use fCondensation
  use gCoagulation
  use hCongregation
  use iBulkequibrium
  use kSOAP
   implicit none
contains
  subroutine aerodyn(start_time,end_time)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This is the main subroutine for solving the GDE for aerosol
!     with a size-composition-resolved method in a given grid cell.
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     start_time: input timestep ([s]).
!     end_time: output timestep ([s]).
!
!     -- OUTPUT VARIABLES
!------------------------------------------------------------------------
    implicit none
    Integer jesp,j,i,j1,j2
    Integer solver
    Double precision emw_tmp,tmp
    doUBLE PRECISION start_time,end_time
    double precision :: timestep_coag,timestep_cond    

    !**** SOAP ****
    double precision :: inorg_total, inorg_bin(N_size) 
    integer :: neq !! YK
    DOUBLE PRECISION :: q(N_size*(1+N_aerosol)+N_aerosol)
    integer :: iq(N_aerosol, N_size) 
    ! integer :: ig(N_aerosol)
    double precision :: ionic, lwc, proton
!    double precision :: lwc_Nsize(N_size),ionic_Nsize(N_size)
   ! double precision :: ionic_Nsize(N_size)
    double precision :: lwcorg_Nsize(N_size)
    ! double precision :: proton_Nsize(N_size),liquid_Nsize(12,N_size)
    double precision :: liquid(12),rhoaer,lwcorg
    double precision :: qext(N_aerosol),surface_equilibrium_conc_tmp(N_aerosol)
    double precision :: qinti_tmp(N_inside_aer)
    integer k,s

    ionic_nsize = 0.d0
    lwc_nsize = 0.d0

    lwcorg_nsize = 0.d0
   
! Initialize the density of aerosols

    do j1 = 1, N_size
      ! call compute_density(N_sizebin,N_aerosol,EH2O,TINYM,concentration_mass,&
      !        mass_density_aer,j1,rhoaer) !! YK
      call compute_density(N_size,N_aerosol,EH2O,TINYM,concentration_mass,&
             mass_density,j1,rhoaer)
      rho_wet_cell(j1) = rhoaer
      if(rho_wet_cell(j1).LT.0.1d-6) rho_wet_cell(j1)=density_aer_bin(j1)
    enddo

    call update_wet_diameter_liquid(1,N_size,concentration_mass, concentration_number,wet_mass,wet_diameter,wet_volume,cell_diam_av)

    call mass_conservation(concentration_mass,concentration_number,concentration_gas, total_mass)

    IF (with_coag.EQ.1) THEN
       do j1 = 1, N_size
          do j2 = 1, N_size
             call compute_bidisperse_coagulation_kernel(Temperature,air_free_mean_path,&
                  wet_diameter(j1),wet_diameter(j2),&
                  wet_mass(j1),wet_mass(j2), kernel_coagulation(j1,j2))
          enddo
       enddo
    ENDIF

!!     Update thermo and microphysical parameters.
    call COMPUTE_RELATIVE_HUMIDITY(humidity,Temperature,Pressure,Relative_Humidity)
    Relative_Humidity = DMIN1(DMAX1(Relative_Humidity, Threshold_RH_inf), Threshold_RH_sup)
! air free mean path

    if (with_coag.EQ.1) then
	call COMPUTE_AIR_FREE_MEAN_PATH(Temperature,Pressure,&
	air_free_mean_path,viscosity)
    endif

    if (with_cond.EQ.1) then

      quadratic_speed=0.D0
      diffusion_coef=0.D0
      soa_sat_conc=0.D0
      soa_part_coef=0.D0

	do jesp=1,N_aerosol
	  if (aerosol_species_interact(jesp).GT.0) then
	      emw_tmp = molecular_weight_aer(jesp) * 1.D-6 ! g/mol
	      call COMPUTE_GAS_DIFFUSIVITY(Temperature,Pressure,&
					  molecular_diameter(jesp),emw_tmp,&
					  collision_factor_aer(jesp),diffusion_coef(jesp) ) ! gas diff coef in air

	      call COMPUTE_QUADRATIC_MEAN_VELOCITY(Temperature,&
		  emw_tmp, quadratic_speed(jesp) ) ! gas quad mean speed in air
	  endIF
	enddo

	do i=1,nesp_pankow
	  jesp=pankow_species(i)
	  emw_tmp = molecular_diameter(jesp) * 1.D-6 ! g/mol
	  call COMPUTE_SATURATION_CONCENTRATION(Temperature,&
	      emw_tmp, vaporization_enthalpy(jesp), saturation_pressure(jesp), soa_sat_conc(jesp))
	enddo

	do i=1,nesp_pom
	  jesp=poa_species(i)
	  emw_tmp = molecular_diameter(jesp) * 1.D-6 ! g/mol
	  call COMPUTE_SATURATION_CONCENTRATION(Temperature,&
	      emw_tmp, vaporization_enthalpy(jesp), saturation_pressure(jesp), soa_sat_conc(jesp))
	enddo
	tmp = 1.D0/RGAS * (1.D0 / Temperature - 1.D0 / 298.D0)!RGAS       perfect gas constant

	do i=1,nesp_pankow
	  jesp=pankow_species(i)
	  soa_part_coef(jesp) = Temperature / 298.D0& ! temperature dependency of partition coefficient
	      * partition_coefficient(jesp) * DEXP(vaporization_enthalpy(jesp) * tmp)
	enddo
	tmp = 1.D0/RGAS * (1.D0 / Temperature - 1.D0 / 300.D0) !Reference at 300K for prymary organic aerosol
	do i=1,nesp_pom
	  jesp=poa_species(i)
	  soa_part_coef(jesp) = Temperature / 300.D0 &! temperature dependency of partition coefficient
	      * partition_coefficient(jesp) * DEXP(vaporization_enthalpy(jesp) * tmp)
	enddo
    endIF
    
    !call mass_conservation(concentration_mass,concentration_number,concentration_gas, total_mass)

    ! write(*,*) "before initstep in ModuleAdaptstep",concentration_gas(9)
    
    initial_time_splitting = 0.D0
    final_time = end_time-start_time
    timestep_splitting=0.D0
    current_sub_time=0.D0
    final_sub_time=0.D0

    timestep_coag = 0.d0
    timestep_cond = 0.d0

!!     **************************************************
!!     SOLVE GAS/AEROSOL GDE LAGRANGIAN EQS
!!     **************************************************
    do while (initial_time_splitting.lt.final_time)
      ! Solve coagulation, condensation/evaporation and nucleation
      ! Compute the characteristic time step of each physical processes.

       ! write(*,*) "timestep_coag and _cond", timestep_coag, timestep_cond
       call initstep(concentration_mass,concentration_number,&
            concentration_gas,&
            timestep_coag,timestep_cond, timestep_splitting, final_time)
       ! write(*,*) "after initstep"

!       call check_nan_inf(6)
      ! Solve with the slowest process (coagulation)
      if (with_coag.eq.1) then
	current_sub_time=initial_time_splitting
	final_sub_time=current_sub_time+timestep_splitting
	tag_coag=1
	tag_cond=0
	tag_nucl=0
	sub_timestep_splitting=timestep_coag
	solver=1            ! only etr for coagulation

!	call  processaero(solver,ionic_Nsize,proton_Nsize,liquid_Nsize)
	call  processaero(solver)

!  	call mass_conservation(concentration_mass,concentration_number,&
!  	concentration_gas, total_mass)
      endif

      ! write(*,*) "after processaero_coag in ModuleAdaptstep",concentration_gas(9)

      if(N_fracmax.gt.1) then
	call redistribution_fraction()!fraction redistribution
      endif

      ! Solve the fastest process (condensation/evaporation and nucleation).
      if (with_cond+with_nucl.gt.0) then
	current_sub_time   = initial_time_splitting
	final_sub_time  = current_sub_time + timestep_splitting
	tag_coag   = 0
	tag_cond   = with_cond
	tag_nucl   = with_nucl
	solver=dynamic_solver
	sub_timestep_splitting=timestep_cond
!	call  processaero(solver,ionic_Nsize,proton_Nsize,liquid_Nsize)
	call  processaero(solver)

        ! write(*,*) "after processaero_cond in ModuleAdaptstep",concentration_gas(9)
        ! write(*,*) "lwc_nsize after processaero cond", lwc_nsize

      endif

      if (with_cond+with_nucl+with_coag.eq.0) then
	final_sub_time  = initial_time_splitting + timestep_splitting!pure emission
      endif

!       ! ! Check mass conservation
!       call mass_conservation(concentration_mass,concentration_number,concentration_gas, total_mass) !! YK removed ?

      initial_time_splitting = final_sub_time

    enddo

    !do C/E equilibrium after each emission
    if (with_cond.gt.0) then

      if(ICUT.ge.1) then

        ! call  bulkequi_inorg(nesp_isorropia,& 
        !        lwc, ionic, proton, liquid,lwc_nsize)!equlibrium for inorganic

        call  bulkequi_inorg(nesp_isorropia,& 
               lwc, ionic, proton, liquid)!equlibrium for inorganic

        ! call redistribution_lwc(lwc,ionic,proton,liquid, &
        !        ionic_Nsize,proton_Nsize,liquid_Nsize)

        call redistribution_lwc(lwc,ionic,proton,liquid)

!      else
!       if inorganics are computed dynamically but organics are at equilibrium,
!       then lwc, proton, ionic and liquid are estimated by isorropia in bulkequi_org      
!        lwc = 0.d0
!       proton = 0.d0
!       ionic = 0.d0
!       liquid = 0.d0
!       do j1 = 1, N_size
!          lwc = lwc + lwc_Nsize(j1)
!          proton = proton + proton_Nsize(j1)
!          ionic = ionic + ionic_Nsize(j1)
!          liquid = liquid + liquid_Nsize(j1)
!       enddo
!       ionic = ionic/N_size
!       liquid = liquid/N_size
      endif

      if (ISOAPDYN.eq.0) then

! ******** equilibrium SOA even if inorganic aerosols are estimated dynamically


         call  bulkequi_org(nesp_pom+nesp_pankow+nesp_aec, &
              lwc,lwcorg,ionic,proton,liquid)!equilibrium for organic


         call redistribution_lwcorg(lwcorg,lwcorg_Nsize)

      else 
  	do jesp=1,N_aerosol
           qext(jesp) = 0.d0
           surface_equilibrium_conc_tmp(jesp) = 0.d0
           do j=1,N_size
              qext(jesp) = qext(jesp) + concentration_mass(j,jesp)
           enddo
        enddo
        call EQINORG(N_aerosol,qext,qinti_tmp,surface_equilibrium_conc_tmp,lwc,ionic,proton,liquid)
        ! call redistribution_lwc(lwc,ionic,proton,liquid, &
        !        lwc_Nsize,ionic_Nsize,proton_Nsize,liquid_Nsize) ! YK
        call redistribution_lwc(lwc,ionic,proton,liquid)
! *** SOA are dynamically partitioned even if inorganic aerosols are estimated by equilibrium.
!
         ! call SOAP_DYN(Relative_Humidity,&
         !      ionic, proton, lwc,lwcorg,&
         !      Temperature, timestep_cond,&
         !      cell_diam_av, neq, q, iq, liquid,&
         !      ionic_Nsize,proton_Nsize,liquid_Nsize)

         call SOAP_DYN(Relative_Humidity,&
              ionic, proton, lwc,lwcorg,&
              Temperature, timestep_cond,&
              cell_diam_av, neq, q, iq, liquid)

      endif

!update wet diameter at the end of current c/e cycle
       call update_wet_diameter_liquid(1,N_size,concentration_mass, &
            concentration_number,wet_mass,wet_diameter,wet_volume,cell_diam_av)
 

      if(N_fracmax.gt.1) then
	call redistribution_fraction()!fraction redistribution
      endif

        call redistribution_size(redistribution_method)!size redistribution

    endif

!    write(*,*) "after redistribution in ModuleAdaptstep",concentration_gas(9)


! IS THIS UPDATE WET DIAMETER REQUIRED OR NOT?
!update wet diameter at the end of current c/e cycle
!	call update_wet_diameter_liquid(1,N_size,concentration_mass, concentration_number,wet_mass,wet_diameter,wet_volume,cell_diam_av)
    call mass_conservation(concentration_mass,concentration_number,&
    concentration_gas, total_mass)

!    write(*,*) "after mass_conservation in ModuleAdaptstep",concentration_gas(9)


  end subroutine aerodyn

!  subroutine processaero(solver,ionic_Nsize,proton_Nsize,liquid_Nsize)
  subroutine processaero(solver)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine solves different aerosol process based on
!     the chosen numerical solver.
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     solver: type of chosen numerical solver
!
!     -- OUTPUT VARIABLES
!------------------------------------------------------------------------
    implicit none
    
    integer:: solver,itrat
    double precision:: time_t,time_p
    integer :: j
!    double precision:: lwc_Nsize(N_size),ionic_Nsize(N_size)
!    double precision:: ionic_Nsize(N_size)
!    double precision:: proton_Nsize(N_size),liquid_Nsize(12,N_size)
!    double precision :: dqdt(N_size, N_aerosol)

    time_p=0.d0
    time_t=final_sub_time-current_sub_time
    itrat=0
!    ******Dynamic time loop
    do while ( current_sub_time .lt. final_sub_time )

     itrat=itrat+1
     ! write(*,*) "time ", current_sub_time, final_sub_time
      if (tag_cond.eq.1) then
                   !     Compute gas mass conservation.
	call mass_conservation(concentration_mass,concentration_number,&
            concentration_gas, total_mass)
         !     solve sulfate dynamically
        if (sulfate_computation.eq.1) then
           if (ICUT.lt.N_size) then
 	      call SULFDYN(concentration_mass_tmp,concentration_mass,concentration_number_tmp,&
	         concentration_number,concentration_gas,dqdt,sub_timestep_splitting)
	   else
	      call SULFDYN(concentration_mass_tmp,concentration_mass,concentration_number_tmp,&
	         concentration_number,concentration_gas,dqdt,time_t)	  
	   endif
        endif
      endif

      ! write(*,*) "after conden processaero ", solver

      if(ICUT.eq.N_size.AND.tag_nucl.eq.0) then
	  sub_timestep_splitting=final_sub_time-current_sub_time
	  current_sub_time=final_sub_time
   
      elseif (solver.eq.0) then
			  ! euler solver used for coagulation
	! call Euler_solver(concentration_mass,concentration_number,concentration_gas,dqdt,& 
        !              proton_Nsize,ionic_Nsize,liquid_Nsize)
	call Euler_solver(concentration_mass,concentration_number,concentration_gas,dqdt)

      elseif (solver.eq.1) then
			    ! etr dynamic solver
	! call Etr_solver(concentration_mass_tmp,concentration_mass,&
	!       concentration_number_tmp,concentration_number,concentration_gas,dqdt,&
        !       proton_Nsize,ionic_Nsize,liquid_Nsize)

         ! write(*,*) "call etr_solver"
	call Etr_solver(concentration_mass_tmp,concentration_mass,&
	      concentration_number_tmp,concentration_number,concentration_gas,dqdt)
      elseif (solver.eq.2) then
			    ! ros2 dynamic solver used for condensation/evaporation of inorganics
	call Ros2_solver(concentration_mass_tmp,concentration_mass,&
	      concentration_number_tmp,concentration_number,&
              concentration_gas,dqdt)
	! call Ros2_solver(concentration_mass_tmp,concentration_mass,&
	!       concentration_number_tmp,concentration_number,&
        !       concentration_gas,dqdt,&
        !       proton_Nsize,ionic_Nsize,liquid_Nsize)
      endif

      ! write(*,*) "after coag processaero "

	call mass_conservation(concentration_mass,concentration_number,&
	  concentration_gas, total_mass)

      time_p=time_p+sub_timestep_splitting

!!     check the need to call ADAPTIME (not for the last timestep)
      If (current_sub_time.lt.final_sub_time.and.ICUT.ne.N_size) then!.and.current_time.gt.(final_time*1.d-2)
	    call adaptime(concentration_mass_tmp,concentration_mass,concentration_number_tmp,&
		concentration_number,sub_timestep_splitting)
      endif
    end do

  end subroutine processaero

  subroutine initstep(c_mass,c_number,c_gas,time_coag,time_cond,&
       time_splitting,t_total)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine performs the initialization of the timestep
!     for the integration of the GDE. The criterion is related to
!     the timescales of the aerosol processes.
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     c_mass: aerosol mass concentration(µg/m^3)
!     c_number: aerosol number concentration(#/m^3)
!     c_gas: aerosol gas phase concentration(µg/m^3)
!     t_total: total time limitation in current loop(s)
!
!     -- OUTPUT VARIABLES
!
!     time_coag: timescales of coagulation process(s)
!     time_cond: timescales of condensation process(s)
!     time_splitting: splitting time step
!------------------------------------------------------------------------
    implicit none

    integer:: j,jesp,s
    double precision:: c_mass(N_size,N_aerosol)
    double precision:: dqdt1(N_size,N_aerosol)
    double precision:: c_number(N_size)
    double precision:: dndt1(N_size)
    double precision:: c_gas(N_aerosol)
    double precision:: time_coag
    double precision:: time_cond
    double precision:: time_splitting
    double precision:: tmp,tscale
    double precision:: t_total!total time of the simunlation
!    double precision:: lwc_Nsize(N_size),ionic_Nsize(N_size)
!    double precision:: ionic_Nsize(N_size)
!    double precision:: proton_Nsize(N_size),liquid_Nsize(12,N_size)
    integer::ICUT_tmp
    double precision :: ce_kernal_coef(N_size, N_aerosol)

    ICUT_tmp=ICUT!avoid the influence of hybrid method

!     Coagulation time step.

    time_coag=t_total
    time_cond=t_total

    if (with_coag.eq.1) then
      tag_coag=1
      tag_cond=0
      tag_nucl=0
      call fgde(c_mass,c_number,c_gas,dqdt1,dndt1,ce_kernal_coef)
!      call fgde(c_mass,c_number,c_gas,dqdt1,dndt1,ce_kernal_coef,ionic_Nsize, &
!                proton_Nsize,liquid_Nsize)

      do j=1,N_size
	tmp=c_number(j)*dndt1(j)
	if (tmp.ne.0.D0) then
	  tscale=c_number(j)/DABS(dndt1(j))
	  time_coag=DMIN1(time_coag,tscale)
	endif
	do s= 1, (N_aerosol-1)
	  jesp=List_species(s)
	  tmp=c_mass(j,jesp)*dqdt1(j,jesp)
	  if (tmp.ne.0.D0.and.c_mass(j,jesp).gt.TINYM) then
	    tscale=c_mass(j,jesp)/DABS(dqdt1(j,jesp))
	    time_coag=DMIN1(time_coag,tscale)
	  endif
	enddo
      end do
    endif

!    write(*,*) "initstep step 1"
!     Cond/evap time step.
    if (with_cond.eq.1.and.ICUT.lt.N_size) then
      ICUT=0
      tag_coag=0
      tag_cond=1
      tag_nucl=with_nucl
      ! call fgde(c_mass,c_number,c_gas,dqdt1,dndt1,ce_kernal_coef,ionic_Nsize, &
      !           proton_Nsize,liquid_Nsize)
      call fgde(c_mass,c_number,c_gas,dqdt1,dndt1,ce_kernal_coef)

   
!      write(*,*) "initstep step 1-1"
      do j=ICUT+1,N_size
	do s= 1, (N_aerosol-1)!s=G1,G2
	  jesp=List_species(s)	  
	  tmp=c_mass(j,jesp)*dqdt1(j,jesp)
	  if (DABS(dqdt1(j,jesp)).gt.0.d0.and.c_mass(j,jesp).gt.TINYM) then
	    tscale=c_mass(j,jesp)/DABS(dqdt1(j,jesp))
	    time_cond=DMIN1(time_cond,tscale)
	    if (time_cond.lt.DTAEROMIN) then
	      time_cond=DTAEROMIN
	    endif
	  endif
	enddo
      end do
    endif

 !   write(*,*) "initstep step 2", time_coag

 !     Minimal time step. chose the bigger one
    if(time_coag.gt.time_cond) then
      time_splitting=time_coag
    else
      time_splitting=time_cond
    endif

    time_splitting=DMIN1(time_splitting,t_total-initial_time_splitting)


    ! if only one process then set
    ! splitting step to whole timestep
    time_splitting = DMAX1(time_splitting,DTAEROMIN)
    time_cond = DMIN1(time_splitting,DMAX1(time_cond,DTAEROMIN))
    time_coag = DMIN1(time_splitting,DMAX1(time_coag,DTAEROMIN))
    if(with_coag.eq.0.OR.with_cond.eq.0)  time_splitting=t_total-initial_time_splitting
    if(ICUT.eq.N_size) then
      time_splitting=t_total-initial_time_splitting
      time_cond=time_splitting
    endif

    ICUT=ICUT_tmp

  end subroutine initstep

  subroutine adaptime(q1,q2,n1,n2,T_dt)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine computes N_aerosol time step for time
!     integration, based on the difference between the
!     first and second order evaluations
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     q1: aerosol mass concentration of first order evaluation(µg/m^3)
!     n1: aerosol number concentration of first order evaluation(#/m^3)
!     q2: aerosol mass concentration of second order evaluation(µg/m^3)
!     n2: aerosol number concentration of second order evaluation(#/m^3)
!
!     -- OUTPUT VARIABLES
!
!     T_dt: dynamic time step for N_aerosol integration
!------------------------------------------------------------------------
    implicit none

    integer j,jesp,s
    double precision:: q1(N_size,N_aerosol)!tmp mass concentration
    double precision:: q2(N_size,N_aerosol)!tmp mass concentration
    double precision:: n1(N_size)!1th order number concentration
    double precision:: n2(N_size)!2d order number concentration
    double precision tmp,n2err,R
    double precision T_dt
!!!     ******zero init
    n2err=0.D0
!!!     ******local error estimation
    do j=1,N_size
      if(n2(j).gt.TINYN) then
	tmp=(n2(j)-n1(j))/(n2(j)+TINYN)
	n2err=n2err+tmp*tmp
      endif
    end do
    do j=1,N_size
      do s= 1, (N_aerosol-1)
	jesp=List_species(s)
	if(q2(j,jesp).gt.TINYM) then
	  tmp=(q2(j,jesp)-q1(j,jesp))/(q2(j,jesp))
	  n2err=n2err+tmp*tmp
	endif
      enddo
    enddo
    n2err=DSQRT(n2err)
!!!     ******compute new time step
    ! first we constrain norm2 error
    ! in order to prevent division by zero
    ! and to keep new time step between
    ! DTMIN and DTMAX defined in time.inc
    R = (1.D+02)/(1.D-05)
    tmp=R*R
    !EPSER= 0.01D0  relative error_precision
    n2err=DMIN1(n2err,EPSER*tmp)
    n2err=DMAX1(EPSER/tmp,n2err)
    ! formula to compute new time step
    T_dt=T_dt*DSQRT(EPSER/n2err)
    T_dt = DMIN1( (final_sub_time-current_sub_time) , DMAX1(DTAEROMIN, T_dt) )
  end subroutine adaptime


  ! subroutine Etr_solver(q1,q2,n1,n2,c_gas,dqdt,proton_Nsize,ionic_Nsize,liquid_Nsize)
  subroutine Etr_solver(q1,q2,n1,n2,c_gas,dqdt)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine solves a system of Ordinary Differential Equations
!     provided by the GDE for aerosols with the Explicit Trapezoidal Rule
!     algorithm (ETR).
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     c_gas: aerosol gas phase concentration(µg/m^3)
!
!     -- INPUT/OUTPUT VARIABLES
!
!     q2: aerosol mass concentration of second order evaluation(µg/m^3)
!     n2: aerosol number concentration of second order evaluation(#/m^3)
!
!     -- OUTPUT VARIABLES
!
!     q1: aerosol mass concentration of first order evaluation(µg/m^3)
!     n1: aerosol number concentration of first order evaluation(#/m^3)
!     dqdt: mass derivation(µg/m^3/s)
!
!------------------------------------------------------------------------
    implicit none

    integer::s,j,jesp!s is the species index of ESO4
    double precision:: dqdt(N_size,N_aerosol)
    double precision:: dq1dt(N_size,N_aerosol)
    double precision:: dq2dt(N_size,N_aerosol)
    double precision:: dn1dt(N_size)
    double precision:: dn2dt(N_size)
    double precision:: c_gas(N_aerosol)
    double precision:: c_gas_t(N_aerosol)
    double precision:: t_mass(N_aerosol)
    double precision:: dtetr,tmp
    double precision:: n1(N_size)!1th order number concentration
    double precision:: n2(N_size)!2d order number concentration
    double precision:: q1(N_size,N_aerosol)!1th order mass concentration
    double precision:: q2(N_size,N_aerosol)!2d order mass concentration
!    double precision:: lwc_Nsize(N_size),ionic_Nsize(N_size)
!    double precision:: ionic_Nsize(N_size)
!    double precision:: proton_Nsize(N_size),liquid_Nsize(12,N_size)
!    double precision :: ce_kernal_coef(N_size,N_aerosol)

    !for condensation or coagulation
    ! call fgde(q2,n2,c_gas,dq1dt,dn1dt,ce_kernal_coef,lwc_Nsize,ionic_Nsize, &
    !           proton_Nsize, liquid_Nsize)

    ! call fgde(q2,n2,c_gas,dq1dt,dn1dt,ce_kernal_coef,ionic_Nsize, &
    !           proton_Nsize, liquid_Nsize)

    ! write(*,*) "before fgde in etr_solver"
    call fgde(q2,n2,c_gas,dq1dt,dn1dt,ce_kernal_coef)
    ! write(*,*) "after fgde in etr_solver"
    !     First step
    do j=1,N_size
      n1(j)=n2(j)+sub_timestep_splitting*dn1dt(j)
      if(n1(j).lt.0.d0) n1(j)=0.d0
      do s=1,(N_aerosol-1)
	jesp=List_species(s)
	t_mass(jesp)=total_mass(jesp)
	q1(j,jesp)=q2(j,jesp)+sub_timestep_splitting*dq1dt(j,jesp)
	if(q1(j,jesp).lt.0.d0) q1(j,jesp)=0.d0
      enddo
    enddo

    ! write(*,*) " after 1st step"

    !     Second step
    ! call fgde(q1,n1,c_gas_t,dq2dt,dn2dt,ce_kernal_coef,lwc_Nsize,ionic_Nsize, &
    !           proton_Nsize, liquid_Nsize)
    ! call fgde(q1,n1,c_gas_t,dq2dt,dn2dt,ce_kernal_coef,ionic_Nsize, &
    !           proton_Nsize, liquid_Nsize)
    call fgde(q1,n1,c_gas_t,dq2dt,dn2dt,ce_kernal_coef)

    ! write(*,*) " after 2nd step fgde, sub_timestep_splitting", sub_timestep_splitting

    dtetr=sub_timestep_splitting*5.0D-01
    current_sub_time = current_sub_time + sub_timestep_splitting
    do j=1,N_size
      tmp=dtetr*(dn1dt(j)+dn2dt(j))
      if((n2(j)+tmp).lt.0.d0) then
	n2(j)=n1(j)
      else
	n2(j)=n2(j)+tmp
      endif
      ! write(*,*) j, n2(j), dn1dt(j), dn2dt(j)
      n_grow_coag=n_grow_coag+tmp
      do s=1,(N_aerosol-1)
	jesp=List_species(s)
        ! write(*,*) dtetr, dq1dt(j,jesp),dq2dt(j,jesp)
	tmp=dtetr*(dq1dt(j,jesp)+dq2dt(j,jesp))
	dqdt(j,jesp)=tmp/sub_timestep_splitting
        ! write(*,*) j, jesp !, dqdt(j,jesp)
	if((q2(j,jesp)+tmp).lt.0.d0) then
	  q2(j,jesp)=q1(j,jesp)
	else
	  q2(j,jesp)=q2(j,jesp)+tmp
	endif
      enddo
    enddo

    ! write(*,*) " end of etr_solver"
  end subroutine Etr_solver

  ! subroutine Euler_solver(c_mass,c_number,c_gas,dqdt,lwc_Nsize,proton_Nsize,ionic_Nsize,liquid_Nsize)
  ! subroutine Euler_solver(c_mass,c_number,c_gas,dqdt,proton_Nsize,ionic_Nsize,liquid_Nsize)
  subroutine Euler_solver(c_mass,c_number,c_gas,dqdt)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine solves a system of Ordinary Differential Equations
!     provided by the GDE for aerosols with Euler algorithm.
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     c_gas: aerosol gas phase concentration(µg/m^3)
!
!     -- INPUT/OUTPUT VARIABLES
!
!     c_mass: aerosol mass concentration (µg/m^3)
!     c_number: aerosol number concentration (#/m^3)
!
!     -- OUTPUT VARIABLES
!
!     dqdt: mass derivation(µg/m^3/s)
!
!------------------------------------------------------------------------
    implicit none

    integer::s,j,jesp
    double precision:: dqdt(N_size,N_aerosol)
    double precision:: dndt(N_size)
    double precision:: c_number(N_size)!number concentration
    double precision:: c_mass(N_size,N_aerosol)!micg/m^-3
    double precision:: c_gas(N_aerosol)!micg/m^-3
    double precision:: tmp
    ! double precision:: lwc_Nsize(N_size),ionic_Nsize(N_size)
!    double precision:: ionic_Nsize(N_size)
!    double precision:: proton_Nsize(N_size),liquid_Nsize(12,N_size)
    double precision :: ce_kernal_coef(N_size,N_aerosol)

    !for condensation
    ! call fgde(c_mass,c_number,c_gas,dqdt,dndt,ce_kernal_coef,lwc_Nsize,ionic_Nsize,&
    !           proton_Nsize,liquid_Nsize)
    ! call fgde(c_mass,c_number,c_gas,dqdt,dndt,ce_kernal_coef,ionic_Nsize,&
    !           proton_Nsize,liquid_Nsize)
    call fgde(c_mass,c_number,c_gas,dqdt,dndt,ce_kernal_coef)


    do j=1,N_size
      c_number(j)=c_number(j)+dndt(j)*sub_timestep_splitting
      do s=1,(N_aerosol-1)
	jesp=List_species(s)
	tmp=sub_timestep_splitting*dqdt(j,jesp)
	c_mass(j,jesp)=c_mass(j,jesp)+tmp
      enddo
    enddo

    current_sub_time = current_sub_time + sub_timestep_splitting
  end subroutine Euler_solver

  ! subroutine Ros2_solver(q1,q2,n1,n2,c_gas,dqdt,lwc_Nsize,proton_Nsize,ionic_Nsize,liquid_Nsize)
  ! subroutine Ros2_solver(q1,q2,n1,n2,c_gas,dqdt,proton_Nsize,ionic_Nsize,liquid_Nsize)
  subroutine Ros2_solver(q1,q2,n1,n2,c_gas,dqdt)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutines solves the GDE with the scond-order Rosenbrock
!     algorithm (ROS2).
!
!     The canonic form of ROS2 is:
!     (I-Gamma*Jn*h) * k1 = fn*h
!     (I-Gamma*J1*h) * k2 = f1*h - Gamma*(Jn+J1)*h*k1
!     c(n+1) = c(n) + (k1+k2)/2
!
!     where c is the variable, f the time derivative and J an
!     approximation of the Jacobian matrix.
!
!     To spare jacobian matrix computation
!     the same approximation can be used at both states
!     (only one LU factorisation):
!     (I-Gamma*Jn*h) * k1 = fn*h
!     (I-Gamma*Jn*h) * k2 = f1*h - 2*Gamma*Jn*h*k1
!     c(n+1) = c(n) + (k1+k2)/2
!
!     A cheaper implementation is furthermore possible
!     (avoiding a matrix-vector product)
!     (I-Gamma*Jn*h) * k1 = fn*h
!     (I-Gamma*Jn*h) * k2 = f1*h - 2*k1
!     c(n+1) = c(n) + (3* k1+k2)/2
!
!     The use of a poor approximation is another way to spare
!     computation. Here only the diagonal part is estimated
!     WITH PHYSICAL ASSUMPTIONS attempting to ensure positivity
!     if f < 0 , it is assumed zero-valued in zero (decreasing may stop)
!     f > 0 , it is assumed zero-valued for total available amount
!     of the considered specie (condensation may stop)
!     Both of these assumptions implies negative values for
!     jacobian diagonal part.
!     With these approximations, the use of canonic form is no longer
!     computational expensive.
!
!     ! Notations.
!     in Verwer (1)   , k is variation rate of variable c
!     in this software, k is first derivative, and then variation
!     (step) for the variable q
!
!     Reference:
!     J.G. Verwer, E.J. Spee, J.G. Blom, W.H. Hundsdorfer
!     "A second order Rosenbrock method applied to photochemical
!     dispersion problems" Report MAS-R9717 August 1997
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     c_gas: aerosol gas phase concentration(µg/m^3)
!
!     -- INPUT/OUTPUT VARIABLES
!
!     q2: aerosol mass concentration of second order evaluation(µg/m^3)
!     n2: aerosol number concentration of second order evaluation(#/m^3)
!
!     -- OUTPUT VARIABLES
!
!     q1: aerosol mass concentration of first order evaluation(µg/m^3)
!     n1: aerosol number concentration of first order evaluation(#/m^3)
!     dqdt: mass derivation(µg/m^3/s)
!
!------------------------------------------------------------------------
    implicit none

    integer::j,jesp,s
    double precision:: dqdt(N_size,N_aerosol)
    double precision:: k1(N_size,N_aerosol)
    double precision:: k2(N_size,N_aerosol)
    double precision:: dn1dt(N_size)
    double precision:: dn2dt(N_size)
    double precision:: Jdn(N_size,N_aerosol)
    double precision:: Jd1(N_size,N_aerosol)
    double precision:: c_gas(N_aerosol)
    double precision:: c_gas_t(N_aerosol)
    double precision:: tmp
    double precision:: q1(N_size,N_aerosol)!1th order mass concentration
    double precision:: q2(N_size,N_aerosol)!2d order mass concentration
    double precision:: n1(N_size)!1th order number concentration
    double precision:: n2(N_size)!2d order number concentration
    double precision Gamma,temp1,temp2
    parameter ( Gamma= 1.7071D0)
    ! double precision:: lwc_Nsize(N_size),ionic_Nsize(N_size)
!    double precision:: ionic_Nsize(N_size)
!    double precision:: proton_Nsize(N_size),liquid_Nsize(12,N_size)
!    double precision :: ce_kernal_coef(N_size,N_aerosol)

  !for condensation
    ! call fgde(q2,n2,c_gas,k1,dn1dt,ce_kernal_coef,ionic_Nsize, &
    !           proton_Nsize,liquid_Nsize)!compute first order derivative
    call fgde(q2,n2,c_gas,k1,dn1dt,ce_kernal_coef)!compute first order derivative


!    write(*,*) "ros2 step 1"
!     Every dynamical variable protected against vanishing
    do j = 1 , N_size
      do s= 1, nesp_isorropia!(N_aerosol-1)
         jesp=List_species(s)
         Jdn(j,jesp) = 0.d0
         if ( q2(j,jesp).gt.1.d-26 .and. k1(j,jesp).lt.0.d0 ) then
	    Jdn(j,jesp) = k1(j,jesp) / q2(j,jesp)!mass variation percentage of each bin
	    !for equilibrium bin k1(j,jesp)=0.d0
         endif
      enddo
    enddo

    ! write(*,*) "ros2 step 2"
!     Condensation limited by available amount (temp1) in both gas phase
!     and small bins (assumed in equilibrium)
!     sum of condensation rates (temp2) is distributed between
!     dynamical bins
    do s= 1, nesp_isorropia!(N_aerosol-1)
       jesp=List_species(s)
       if (aerosol_species_interact(jesp).gt.0) then
	  temp1 = c_gas(jesp)
	  if(jesp.ge.E1.and.jesp.le.E2) then!Gas phase
             do j = 1, ICUT !small bins
                temp1 = temp1 + q2(j,jesp)! the sum of gas and equilibrium bin mass
             enddo
          endif

	  temp2 = 0.d0
	  do j = (ICUT+1), N_size!for dynamic bins
	    if (k1(j,jesp).gt.0.d0) then
	      temp2 = temp2 + k1(j,jesp)! the sum of dynamic bin mass variation
	    endif
	  enddo

	  do j = (ICUT+1), N_size
	    if (temp1-q2(j,jesp).gt.1.d-26 .AND. k1(j,jesp).gt.0.D0) then!in case of condensation
	      Jdn(j,jesp) = -temp2 / (temp1-q2(j,jesp))
	    endif
	  enddo
     endif
    enddo

    ! write(*,*) "ros2 step 3"
!     ******2. first evaluation

    do j = 1 , N_size
      dn1dt(j)=dn1dt(j)*sub_timestep_splitting!dndt always 0
      do s= 1, nesp_isorropia!(N_aerosol-1)
	jesp=List_species(s)
	k1(j,jesp) = k1(j,jesp) * sub_timestep_splitting / ( 1.d0 - Gamma * Jdn(j,jesp) * sub_timestep_splitting )
      enddo
    enddo

    call KLIMIT(q2,c_gas,k1,ce_kernal_coef)

    do j = 1 , N_size
      n1(j) = DMAX1 ( 0.01d0*n2(j) , (n2(j)+dn1dt(j)) )
      do s= 1, nesp_isorropia!(N_aerosol-1)
	jesp=List_species(s)
	q1(j,jesp) = DMAX1 ( 0.01d0*q2(j,jesp) , (q2(j,jesp)+k1(j,jesp)) )
      enddo
    enddo

!     ******3. compute the derivative k2 in (q+k1 ; Tin2+sub_timestep_splitting)
!     ******&  the approximation of jacobian diagonal

    current_sub_time = current_sub_time + sub_timestep_splitting
    ! call fgde(q1,n1,c_gas_t,k2,dn2dt,ce_kernal_coef,lwc_Nsize,ionic_Nsize, &
    !           proton_Nsize,liquid_Nsize)

    ! call fgde(q1,n1,c_gas_t,k2,dn2dt,ce_kernal_coef,ionic_Nsize, &
    !           proton_Nsize,liquid_Nsize)

    call fgde(q1,n1,c_gas_t,k2,dn2dt,ce_kernal_coef)

    do j = 1 , N_size
      do s= 1, nesp_isorropia!(N_aerosol-1)
	jesp=List_species(s)
        if (isnan(k2(j,jesp))) then
           write(*,*) j,s,list_species(s),jd1(j,jesp),k2(j,jesp),q1(j,jesp)
           stop
        endif

	if (IsNaN(k2(j,jesp))) then 
           k2(j,jesp) = 0.d0 !NAN
        endif

        Jd1(j,jesp) = 0.d0
	if ( q1(j,jesp).gt.1.d-26 .and. k2(j,jesp).lt.0.d0 ) then
	   Jd1(j,jesp) = k2(j,jesp) / q1(j,jesp)
        endif 
      enddo
    enddo

    do s= 1, nesp_isorropia!(N_aerosol-1)
      jesp=List_species(s)
      if (aerosol_species_interact(jesp).gt.0) then
	temp1 = c_gas(jesp)
 	do j = 1, ICUT
 	  temp1 = temp1 + q1(j,jesp)
 	enddo

	temp2 = 0.d0
	do j = (ICUT+1), N_size
	  if (k2(j,jesp).gt.0.d0) temp2 = temp2 + k2(j,jesp)!Temperature 2 sum of drivation
	enddo
	do j = (ICUT+1), N_size
	  if ((temp1-q1(j,jesp)).gt.1.d-26 .and. k2(j,jesp).gt.0.d0) then
	    Jd1(j,jesp) = -temp2  / (temp1-q1(j,jesp))
	  endif
	enddo
      endif
    enddo

!     ******4. second evaluation

    do j = 1 , N_size
      dn2dt(j)=dn2dt(j)*sub_timestep_splitting
      do s= 1, nesp_isorropia!(N_aerosol-1)
	jesp=List_species(s)
	k2(j,jesp) = (k2(j,jesp)*sub_timestep_splitting- Gamma*sub_timestep_splitting&
	    *(Jdn(j,jesp)+Jd1(j,jesp))*k1(j,jesp))&
	    / (1.d0- Gamma*Jd1(j,jesp) *sub_timestep_splitting  )
      enddo
    enddo

    call KLIMIT(q2,c_gas,k2,ce_kernal_coef)

    do j = 1 , N_size
      tmp=n2(j)
      n2(j)=DMAX1 ( 0.01d0*n2(j) , (n2(j)+0.5d0*(dn1dt(j)+dn2dt(j))) )
      n_grow_nucl=n_grow_nucl+(n2(j)-tmp)
      do s = 1, nesp_isorropia
	  jesp=isorropia_species(s)
	  tmp=q2(j,jesp)
	  q2(j,jesp)=DMAX1(0.01d0*q2(j,jesp),(q2(j,jesp)+0.5d0*(k1(j,jesp)+k2(j,jesp))) )
	  tmp=q2(j,jesp)-tmp
	  dqdt(j,jesp)=tmp/sub_timestep_splitting
	  m_grow_cond=m_grow_cond+tmp
      enddo
    enddo

  end subroutine Ros2_solver
  
end Module jAdaptstep
