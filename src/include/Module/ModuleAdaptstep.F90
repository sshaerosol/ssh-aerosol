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
  subroutine aerodyn(start_time,delta_t)
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
    integer :: neq
    DOUBLE PRECISION :: q(N_size*(1+N_aerosol)+N_aerosol)
    integer :: iq(N_aerosol, N_size) 
    double precision :: ionic, lwc, proton
    double precision :: lwcorg_Nsize(N_size)
    double precision :: liquid(12),rhoaer,lwcorg
    double precision :: qext(N_aerosol),surface_equilibrium_conc_tmp(N_aerosol)
    double precision :: qinti_tmp(N_inside_aer)
    double precision :: delta_t
    integer :: k,s,b,nb_iter

    double precision :: timestep_splitting,sub_timestep_splitting
    double precision :: initial_time_splitting,current_sub_time,final_sub_time

    lwcorg_nsize = 0.d0

    ! Initialize the density of aerosols
    if (with_fixed_density.ne.1) then

       do j1 = 1, N_size
          call compute_density(N_size,N_aerosol,EH2O,TINYM,concentration_mass,&
               mass_density,j1,rhoaer)
          rho_wet_cell(j1) = rhoaer
          if(rho_wet_cell(j1).LT.0.1d-6) rho_wet_cell(j1)=density_aer_bin(j1)
       enddo
    end if

    call update_wet_diameter_liquid(1,N_size,concentration_mass, concentration_number,wet_mass,wet_diameter,wet_volume,cell_diam_av)

    call mass_conservation(concentration_mass,concentration_number,concentration_gas, total_mass)

    if (with_coag.EQ.1) then
       call COMPUTE_AIR_FREE_MEAN_PATH(Temperature,Pressure,&
            air_free_mean_path,viscosity)
       do j1 = 1, N_size
          do j2 = 1, N_size
             call compute_bidisperse_coagulation_kernel(Temperature,air_free_mean_path,&
                  wet_diameter(j1),wet_diameter(j2),&
                  wet_mass(j1),wet_mass(j2), kernel_coagulation(j1,j2))
          enddo
       enddo
    endif

    !   Update thermo and microphysical parameters.

    if (with_cond.EQ.1) then
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
    endif

    initial_time_splitting = 0.D0
    end_time = start_time + delta_t
    timestep_splitting=0.D0
    current_sub_time=0.D0
    final_sub_time=0.D0

    timestep_coag = 0.d0
    timestep_cond = 0.d0

    !!     **************************************************
    !!     SOLVE GAS/AEROSOL GDE LAGRANGIAN EQS
    !!     **************************************************

    nb_iter = 0
    do while (initial_time_splitting.lt.delta_t)
       if(splitting == 0) then ! Initial time step when processes are splitted
          if(nb_iter == 0) then
             ! Solve coagulation, condensation/evaporation and nucleation
             ! Compute the characteristic time step of each physical process.
             call initstep(concentration_mass,concentration_number,&
                  concentration_gas,timestep_coag,timestep_cond, &
                  initial_time_splitting,timestep_splitting,delta_t)
             nb_iter = 1
          else
             timestep_splitting = DMAX1(timestep_coag,timestep_cond)
             timestep_splitting = DMIN1(timestep_splitting,delta_t-current_sub_time)
          endif
       else  
          call initstep_coupled(concentration_mass,concentration_number,&
               concentration_gas,sub_timestep_splitting, &
               initial_time_splitting,timestep_splitting,delta_t)
       endif

       if(splitting == 0) then ! Split processes
          ! Solve with the slowest process (coagulation)
          if (with_coag.eq.1) then
             current_sub_time = initial_time_splitting
             final_sub_time = current_sub_time+timestep_splitting
             tag_coag = 1
             tag_cond = 0
             tag_nucl = 0
             solver=1            ! only etr for coagulation
             call  processaero(solver,current_sub_time,timestep_coag,final_sub_time,splitting)
          endif
          if(N_fracmax.gt.1) then
             call redistribution_fraction()!fraction redistribution
          endif

          ! Solve the fastest process (condensation/evaporation and nucleation).
          if (with_cond+with_nucl.gt.0) then
             current_sub_time = initial_time_splitting
             final_sub_time = current_sub_time + timestep_splitting
             tag_coag = 0
             tag_cond = with_cond
             tag_nucl = with_nucl
             solver=dynamic_solver
             call  processaero(solver,current_sub_time,timestep_cond,final_sub_time,splitting)
          endif

       else
          current_sub_time=initial_time_splitting
          final_sub_time = current_sub_time+timestep_splitting
          tag_coag = with_coag
          tag_cond = with_cond
          tag_nucl = with_nucl
          solver = 1            ! only etr for coagulation
          call  processaero(solver,current_sub_time,sub_timestep_splitting,final_sub_time,splitting)
       endif

       if (with_cond+with_nucl+with_coag.eq.0) then
          final_sub_time = initial_time_splitting + timestep_splitting !pure emission
       endif

       initial_time_splitting = final_sub_time

    enddo

    !do C/E equilibrium after each emission
    if (with_cond.gt.0) then

       if(ICUT.ge.1) then

          call  bulkequi_inorg(nesp_isorropia,& 
               lwc, ionic, proton, liquid) !equlibrium for inorganic

          call redistribution_lwc(lwc,ionic,proton,liquid)

       endif

       if (ISOAPDYN.eq.0) then

          ! ******** equilibrium SOA even if inorganic aerosols are estimated dynamically

          call  bulkequi_org(nesp_eq_org,lwc,lwcorg,ionic,proton,liquid)!equilibrium for organic

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
          call redistribution_lwc(lwc,ionic,proton,liquid)

          ! *** SOA are dynamically partitioned even if inorganic aerosols are estimated by equilibrium.
          !
          call SOAP_DYN(Relative_Humidity,&
               ionic, proton, lwc,lwcorg,&
               Temperature, delta_t,&
               cell_diam_av, neq, q, iq, liquid)

       endif

       call update_wet_diameter_liquid(1,N_size,concentration_mass, &
            concentration_number,wet_mass,wet_diameter,wet_volume,cell_diam_av)

       if(N_fracmax.gt.1 .and. redistribution_method.ne.0) then
          call redistribution_fraction()!fraction redistribution
       endif

       if (redistribution_method.ne.0) call redistribution_size(redistribution_method)!size redistribution

    endif

    call update_wet_diameter_liquid(1,N_size,concentration_mass,concentration_number,&
         wet_mass,wet_diameter,wet_volume,cell_diam_av) 

    ! AVOID mass conservation after redistribution
    ! call mass_conservation(concentration_mass,concentration_number,&
    !   concentration_gas, total_mass)

  end subroutine aerodyn

  subroutine processaero(solver,current_sub_time,sub_timestep_splitting,&
       final_sub_time,splitting)
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

    integer:: solver,splitting
    double precision:: time_t,time_p,time_step_sulf
    double precision:: final_sub_time,current_sub_time,sub_timestep_splitting
    integer :: i,j

    time_p=0.d0
    time_t=final_sub_time-current_sub_time
    !    ******Dynamic time loop
    do while ( current_sub_time .lt. final_sub_time )
       call update_wet_diameter_liquid(1,N_size,concentration_mass, &
            concentration_number,wet_mass,wet_diameter,wet_volume,cell_diam_av)
       time_step_sulf = 0.0

       if (tag_cond.eq.1) then
          !!if (tag_emis .ne. 0) call emission(sub_timestep_splitting)
          !     Compute gas mass conservation.
          call mass_conservation(concentration_mass,concentration_number,&
               concentration_gas, total_mass)
          !     solve sulfate dynamically
          if (sulfate_computation.eq.1) then
             call SULFDYN(concentration_mass_tmp,concentration_mass,concentration_number_tmp,&
                  concentration_number,concentration_gas,dqdt,sub_timestep_splitting,time_step_sulf)
          endif
       endif

       if(ICUT.eq.N_size.AND.tag_nucl.eq.0) then
          if ((sulfate_computation.eq.1).AND.(time_step_sulf>0)) then !Case full equilibrium but sulfate computed dynamically
             time_p=time_p+sub_timestep_splitting
             current_sub_time = current_sub_time + sub_timestep_splitting
             sub_timestep_splitting = DMIN1(time_step_sulf,final_sub_time-current_sub_time)
          else ! case full equilibrium and sulfate computed in full equilibrium
             sub_timestep_splitting=final_sub_time-current_sub_time
             current_sub_time=final_sub_time
          endif

       elseif (solver.eq.0) then
          ! euler solver used for coagulation
          call Euler_solver(concentration_mass,concentration_number,concentration_gas,dqdt,&
               current_sub_time,sub_timestep_splitting)

       elseif (solver.eq.1) then
          ! etr dynamic solver
          call Etr_solver(concentration_mass_tmp,concentration_mass,&
               concentration_number_tmp,concentration_number,concentration_gas,dqdt,&
               current_sub_time,sub_timestep_splitting)
       elseif (solver.eq.2) then
          ! ros2 dynamic solver used for condensation/evaporation of inorganics
          call Ros2_solver(concentration_mass_tmp,concentration_mass,&
               concentration_number_tmp,concentration_number,&
               concentration_gas,dqdt,&
               current_sub_time,sub_timestep_splitting)
       endif

       call mass_conservation(concentration_mass,concentration_number,&
            concentration_gas, total_mass)

       time_p=time_p+sub_timestep_splitting

       If (current_sub_time.le.final_sub_time.and.ICUT.ne.N_size) then
          call adaptime(concentration_mass_tmp,concentration_mass,concentration_number_tmp,&
               concentration_number,sub_timestep_splitting,time_step_sulf,current_sub_time,&
               final_sub_time)
          ! Need to redistribute onto fixed grid if nucleation is solved with 
          ! condensation/evaporation or if processes are not splitted
          if((tag_nucl.EQ.1).OR.(splitting.EQ.1)) then 
             call update_wet_diameter_liquid(1,N_size,concentration_mass, &
                  concentration_number,wet_mass,wet_diameter,wet_volume,cell_diam_av)
             if (redistribution_method.ne.0) call redistribution_size(redistribution_method)!size redistribution
             if(N_fracmax.gt.1) then
                call redistribution_fraction()!fraction redistribution
             endif
          endif
       endif
    end do

  end subroutine processaero

  subroutine initstep_coupled(c_mass,c_number,c_gas,sub_time_splitting,&
       initial_time_splitting,time_splitting,t_total)
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
    !     c_mass: aerosol mass concentration(\B5g/m^3)
    !     c_number: aerosol number concentration(#/m^3)
    !     c_gas: aerosol gas phase concentration(\B5g/m^3)
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
    double precision:: sub_time_splitting
    double precision:: time_splitting,initial_time_splitting
    double precision:: tmp,tscale
    double precision:: t_total!total time of the simunlation
    integer::ICUT_tmp
    double precision :: rate,ce_kernal_coef(N_size, N_aerosol)

    ICUT_tmp = ICUT!avoid the influence of hybrid method

    !     Coagulation time step.

    sub_time_splitting=t_total-initial_time_splitting

    tag_coag = with_coag
    tag_cond = with_cond
    tag_nucl = with_nucl
    call fgde(c_mass,c_number,c_gas,dqdt1,dndt1,ce_kernal_coef)

    do j=1,N_size
       tmp=c_number(j)*dndt1(j)
       if (tmp.ne.0.D0) then
          tscale = 0.1d0 * c_number(j)/DABS(dndt1(j))
          sub_time_splitting = DMIN1(sub_time_splitting,tscale)
       endif
       do s= 1, (N_aerosol-1)
          jesp=List_species(s)
          tmp=c_mass(j,jesp)*dqdt1(j,jesp)
          if (tmp.ne.0.D0.and.c_mass(j,jesp).gt.TINYM) then
             tscale = 0.1d0 * c_mass(j,jesp)/DABS(dqdt1(j,jesp))
             sub_time_splitting = DMIN1(sub_time_splitting,tscale)
          endif
       enddo
    end do
    sub_time_splitting = DMIN1(sub_time_splitting,t_total-initial_time_splitting)

    if (sub_time_splitting.lt.DTAEROMIN) then
       sub_time_splitting=DTAEROMIN
    endif
    sub_time_splitting = DMIN1(sub_time_splitting,t_total-initial_time_splitting)

    ! if only one process then set
    ! splitting step to whole timestep
    ! or if processes are solved coupled
    time_splitting = t_total-initial_time_splitting

  end subroutine initstep_coupled

  subroutine initstep(c_mass,c_number,c_gas,time_coag,time_cond,&
       initial_time_splitting,time_splitting,t_total)
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
    !     c_mass: aerosol mass concentration(\B5g/m^3)
    !     c_number: aerosol number concentration(#/m^3)
    !     c_gas: aerosol gas phase concentration(\B5g/m^3)
    !     t_total: total time limitation in current loop(s)
    !
    !     -- OUTPUT VARIABLES
    !
    !     time_coag: timescales of coagulation process(s)
    !     time_cond: timescales of condensation process(s)
    !     time_splitting: splitting time step
    !------------------------------------------------------------------------
    implicit none

    integer:: j,jesp,s,splitting
    double precision:: c_mass(N_size,N_aerosol)
    double precision:: dqdt1(N_size,N_aerosol)
    double precision:: c_number(N_size)
    double precision:: dndt1(N_size)
    double precision:: c_gas(N_aerosol)
    double precision:: time_coag
    double precision:: time_cond
    double precision:: time_splitting,initial_time_splitting
    double precision:: tmp,tscale
    double precision:: t_total!total time of the simunlation
    integer::ICUT_tmp
    double precision :: rate,ce_kernal_coef(N_size, N_aerosol)

    ICUT_tmp = ICUT!avoid the influence of hybrid method

    !     Coagulation time step.

    time_coag=t_total-initial_time_splitting
    time_cond=t_total-initial_time_splitting

    if (with_coag.eq.1) then
       tag_coag=1
       tag_cond=0
       tag_nucl=0

       call fgde(c_mass,c_number,c_gas,dqdt1,dndt1,ce_kernal_coef)

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
       time_coag = DMIN1(time_coag,t_total-initial_time_splitting)
    endif

    !     Cond/evap time step.
    if (with_cond.eq.1.OR.with_nucl.eq.1) then
       ICUT=0
       tag_coag=0
       tag_cond=with_cond
       tag_nucl=with_nucl
       call fgde(c_mass,c_number,c_gas,dqdt1,dndt1,ce_kernal_coef)
       if (sulfate_computation.eq.1) then
          do j = 1,N_size! Reassigned distribution by mass of each species
             call compute_condensation_transfer_rate(diffusion_coef(ESO4), &
                  quadratic_speed(ESO4), accomodation_coefficient(ESO4), &
                  wet_diameter(j), rate)
             dqdt1(j,ESO4) =  dqdt1(j,ESO4) + rate * c_number(j)*c_gas(ESO4)
          enddo
       endif
       do j=1,N_size  !Loop from 1 in case of nucleation - if bins at equilibrium then dqdt1 = 0 from fgde
	  if (DABS(dndt1(j)).gt.0.d0.and.c_number(j).gt.TINYN) then
             tscale=c_number(j)/DABS(dndt1(j))
             time_cond=DMIN1(time_cond,tscale)
	  endif
	  do s= 1, (N_aerosol-1)
             jesp=List_species(s)  
             tmp=c_mass(j,jesp)*dqdt1(j,jesp)
             if (DABS(dqdt1(j,jesp)).gt.0.d0.and.c_mass(j,jesp).gt.TINYM) then
                tscale=c_mass(j,jesp)/DABS(dqdt1(j,jesp))
                time_cond=DMIN1(time_cond,tscale)
             endif
	  enddo
       end do

       if (time_cond.lt.DTAEROMIN) then
          time_cond=DTAEROMIN
       endif
       time_cond = DMIN1(time_cond,t_total-initial_time_splitting)
    endif

    if(time_coag.gt.time_cond) then
       time_splitting=time_coag
    else
       time_splitting=time_cond
    endif

    time_splitting=DMIN1(time_splitting,t_total-initial_time_splitting)
    time_splitting = DMAX1(time_splitting,DTAEROMIN)
    time_cond = DMIN1(time_splitting,DMAX1(time_cond,DTAEROMIN))
    time_coag = DMIN1(time_splitting,DMAX1(time_coag,DTAEROMIN))

    ! if only one process then set
    ! splitting step to whole timestep
    ! or if processes are solved coupled
    if(with_coag.eq.0.OR.with_cond.eq.0)  then
       time_splitting = t_total-initial_time_splitting
    endif
    if(ICUT.eq.N_size) then
       time_splitting=t_total-initial_time_splitting
       time_cond=time_splitting
    endif

    ICUT=ICUT_tmp

  end subroutine initstep

  subroutine adaptime(q1,q2,n1,n2,T_dt,time_step_sulf,current_sub_time,&
       final_sub_time)
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
    !     q1: aerosol mass concentration of first order evaluation(\B5g/m^3)
    !     n1: aerosol number concentration of first order evaluation(#/m^3)
    !     q2: aerosol mass concentration of second order evaluation(\B5g/m^3)
    !     n2: aerosol number concentration of second order evaluation(#/m^3)
    !
    !     -- OUTPUT VARIABLES
    !
    !     T_dt: dynamic time step for N_aerosol integration
    !------------------------------------------------------------------------
    implicit none

    integer j,jesp,s
    double precision :: q1(N_size,N_aerosol)!tmp mass concentration
    double precision :: q2(N_size,N_aerosol)!tmp mass concentration
    double precision :: n1(N_size)!1th order number concentration
    double precision :: n2(N_size)!2d order number concentration
    double precision :: tmp,n2err,R
    double precision :: T_dt,time_step_sulf
    double precision ::  current_sub_time,final_sub_time
    !     ******zero init
    n2err=0.d0
    !     ******local error estimation
    do j=1,N_size
       if(n2(j).gt.TINYN) then
          tmp=(n2(j)-n1(j))/(n2(j)+TINYN)
          n2err=n2err+tmp*tmp
       endif
    end do
    do j=1,N_size
       do s= 1, (N_aerosol-1) ! Do not consider water for time step
          jesp=List_species(s)
          if(q2(j,jesp).gt.TINYM) then
             tmp=(q2(j,jesp)-q1(j,jesp))/(q2(j,jesp))
             n2err=n2err+tmp*tmp
          endif
       enddo
    enddo
    n2err=DSQRT(n2err)

!!!     ******compute new time step
    ! formula to compute new time step
    if(n2err.NE.0.d0) then
       T_dt=T_dt*DSQRT(EPSER/n2err)
    else
       T_dt=final_sub_time-current_sub_time
    endif
    if(time_step_sulf > 0) T_dt = DMIN1(T_dt,time_step_sulf)
    T_dt = DMAX1(DTAEROMIN, T_dt)
    if(final_sub_time-current_sub_time.GT.1.d-20) then ! Estimate time step for next splitting iteration
       T_dt = DMIN1( (final_sub_time-current_sub_time), T_dt)
    endif

  end subroutine adaptime

  subroutine Etr_solver(q1,q2,n1,n2,c_gas,dqdt,current_sub_time,sub_timestep_splitting)
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
    !     c_gas: aerosol gas phase concentration(\B5g/m^3)
    !
    !     -- INPUT/OUTPUT VARIABLES
    !
    !     q2: aerosol mass concentration of second order evaluation(\B5g/m^3)
    !     n2: aerosol number concentration of second order evaluation(#/m^3)
    !
    !     -- OUTPUT VARIABLES
    !
    !     q1: aerosol mass concentration of first order evaluation(\B5g/m^3)
    !     n1: aerosol number concentration of first order evaluation(#/m^3)
    !     dqdt: mass derivation(\B5g/m^3/s)
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
    double precision:: dtetr,tmp,current_sub_time,sub_timestep_splitting
    double precision:: n1(N_size)!1th order number concentration
    double precision:: n2(N_size)!2d order number concentration
    double precision:: q1(N_size,N_aerosol)!1th order mass concentration
    double precision:: q2(N_size,N_aerosol)!2d order mass concentration

    !for condensation or coagulation

    call fgde(q2,n2,c_gas,dq1dt,dn1dt,ce_kernal_coef)
    !     First step
    do j=1,N_size
       if(n2(j)+dn1dt(j)*sub_timestep_splitting.GE.TINYN) then
          n1(j)=n2(j)+sub_timestep_splitting*dn1dt(j)
       else
          n1(j) = 0.0
          ! write(*,*) 'pb ETR n1 < 0',n2(j),sub_timestep_splitting*dn1dt(j),sub_timestep_splitting
       endif
       do s=1,(N_aerosol) 
          jesp=List_species(s)
          t_mass(jesp)=total_mass(jesp)
          if(q2(j,jesp)+dq1dt(j,jesp)*sub_timestep_splitting.GE.TINYM) then
             q1(j,jesp)=q2(j,jesp)+sub_timestep_splitting*dq1dt(j,jesp)
          else
             q1(j,jesp) = 0.0
             !if (dq1dt(j,jesp).GT.1.d-10) write(*,*) 'pb ETR q1 < 0',q2(j,jesp),sub_timestep_splitting*dq1dt(j,jesp),sub_timestep_splitting
          endif
       enddo
    enddo

    !     Second step
    call fgde(q1,n1,c_gas_t,dq2dt,dn2dt,ce_kernal_coef)

    dtetr=sub_timestep_splitting*5.0D-01
    current_sub_time = current_sub_time + sub_timestep_splitting
    do j=1,N_size
       tmp=dtetr*(dn1dt(j)+dn2dt(j))
       if(n2(j)+tmp.GE.TINYN) then
          n2(j)=n2(j)+tmp
       else
          n2(j) = 0.0
          ! write(*,*) 'pb ETR n2 < 0',n2(j),sub_timestep_splitting*dn2dt(j),sub_timestep_splitting
       endif
       do s=1,(N_aerosol)
          jesp=List_species(s)
          tmp=dtetr*(dq1dt(j,jesp)+dq2dt(j,jesp))
          if(tmp+q2(j,jesp).GE.TINYM) then
             dqdt(j,jesp)=tmp/sub_timestep_splitting
             q2(j,jesp)=q2(j,jesp)+tmp
          else
             q2(j,jesp) = 0.0
             !   write(*,*) 'pb ETR q2 < 0',q2(j,jesp),sub_timestep_splitting*dq2dt(j,jesp),sub_timestep_splitting
          endif
       enddo
    enddo

  end subroutine Etr_solver

  subroutine Euler_solver(c_mass,c_number,c_gas,dqdt,current_sub_time,sub_timestep_splitting)
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
    !     c_gas: aerosol gas phase concentration(\B5g/m^3)
    !
    !     -- INPUT/OUTPUT VARIABLES
    !
    !     c_mass: aerosol mass concentration (\B5g/m^3)
    !     c_number: aerosol number concentration (#/m^3)
    !
    !     -- OUTPUT VARIABLES
    !
    !     dqdt: mass derivation(\B5g/m^3/s)
    !
    !------------------------------------------------------------------------
    implicit none

    integer::s,j,jesp
    double precision:: dqdt(N_size,N_aerosol)
    double precision:: dndt(N_size)
    double precision:: c_number(N_size)!number concentration
    double precision:: c_mass(N_size,N_aerosol)!micg/m^-3
    double precision:: c_gas(N_aerosol)!micg/m^-3
    double precision:: tmp,current_sub_time,sub_timestep_splitting
    double precision :: ce_kernal_coef(N_size,N_aerosol)

    !for condensation
    call fgde(c_mass,c_number,c_gas,dqdt,dndt,ce_kernal_coef)

    do j=1,N_size
       if(c_number(j)+dndt(j)*sub_timestep_splitting.GE.0.d0) then
          c_number(j) = c_number(j)+dndt(j)*sub_timestep_splitting
       else
          c_number(j)= 0.d0
       endif
       do s=1,(N_aerosol) 
          jesp=List_species(s)
          tmp=sub_timestep_splitting*dqdt(j,jesp)
          if(tmp + c_mass(j,jesp).GE.0.d0) then
             c_mass(j,jesp)=c_mass(j,jesp)+tmp
          endif
       enddo
    enddo

    current_sub_time = current_sub_time + sub_timestep_splitting
  end subroutine Euler_solver

  subroutine Ros2_solver(q1,q2,n1,n2,c_gas,dqdt,current_sub_time,sub_timestep_splitting)
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
    !     c_gas: aerosol gas phase concentration(\B5g/m^3)
    !
    !     -- INPUT/OUTPUT VARIABLES
    !
    !     q2: aerosol mass concentration of second order evaluation(\B5g/m^3)
    !     n2: aerosol number concentration of second order evaluation(#/m^3)
    !
    !     -- OUTPUT VARIABLES
    !
    !     q1: aerosol mass concentration of first order evaluation(\B5g/m^3)
    !     n1: aerosol number concentration of first order evaluation(#/m^3)
    !     dqdt: mass derivation(\B5g/m^3/s)
    !
    !------------------------------------------------------------------------
    implicit none

    integer::j,jesp,s
    double precision:: dqdt(N_size,N_aerosol)
    double precision:: dq1dt(N_size,N_aerosol)
    double precision:: dq2dt(N_size,N_aerosol)
    double precision:: dn1dt(N_size)
    double precision:: dn2dt(N_size)
    double precision:: Jdq1(N_size,N_aerosol)
    double precision:: Jdq2(N_size,N_aerosol)
    double precision:: Jdn1(N_size)
    double precision:: Jdn2(N_size)
    double precision:: c_gas(N_aerosol)
    double precision:: c_gas_t(N_aerosol)
    double precision:: tmp
    double precision:: q1(N_size,N_aerosol)!1th order mass concentration
    double precision:: q2(N_size,N_aerosol)!2d order mass concentration
    double precision:: n1(N_size)!1th order number concentration
    double precision:: n2(N_size)!2d order number concentration
    double precision:: Gamma,temp1,temp2
    double precision:: current_sub_time,sub_timestep_splitting
    parameter ( Gamma= 1.7071D0)

    !for condensation
    call fgde(q2,n2,c_gas,dq1dt,dn1dt,ce_kernal_coef)!compute first order derivative

    !     Every dynamical variable protected against vanishing
    do j = 1 , N_size
       do s= 1, N_aerosol !!nesp_isorropia!(N_aerosol-1)
          jesp=List_species(s)
          Jdq1(j,jesp) = 0.d0
          if ( q2(j,jesp).gt.TINYM .and. dq1dt(j,jesp).lt.0.d0 ) then ! case of evaporation
             Jdq1(j,jesp) = dq1dt(j,jesp) / q2(j,jesp)!mass variation percentage of each bin
             !for equilibrium bin dq1dt(j,jesp)=0.d0
          endif
          Jdn1(j) = 0.d0
          if ( n2(j).gt.TINYN)then
             Jdn1(j) = dn1dt(j) / (n2(j)+TINYN)
          endif
       enddo
    enddo

    !     Condensation limited by available amount (temp1) in both gas phase
    !     (and small bins (assumed in equilibrium))
    !     sum of condensation rates (temp2) is distributed between
    !     dynamical bins
    do s= 1, N_aerosol !! nesp_isorropia!(N_aerosol-1)
       jesp=List_species(s)
       if (aerosol_species_interact(jesp).gt.0) then
	  temp1 = c_gas(jesp)
   !if(jesp.ge.E1.and.jesp.le.E2) then!Gas phase
   !   do j = 1, ICUT !small bins
   !      temp1 = temp1 + q2(j,jesp)! the sum of gas and equilibrium bin mass
   !   enddo
   !endif

	  temp2 = 0.d0
   !	  do j = (ICUT+1), N_size!for dynamic bins !! Pb in case of nucleation => dq1dt(1) non pris en compte
	  do j = 1, N_size
             if (dq1dt(j,jesp).gt.0.d0) then
                temp2 = temp2 + dq1dt(j,jesp)! the sum of dynamic bin mass variation
             endif
	  enddo

   !	  do j = (ICUT+1), N_size
	  do j = 1, N_size
             if (temp1-q2(j,jesp).gt.TINYM .AND. dq1dt(j,jesp).gt.0.D0) then!in case of condensation
                Jdq1(j,jesp) = -temp2 / (temp1-q2(j,jesp))
             endif
	  enddo
       endif
    enddo

    !     ******2. first evaluation

    do j = 1 , N_size
       dn1dt(j)=dn1dt(j)*sub_timestep_splitting/ ( 1.d0 - Gamma * Jdn1(j) * sub_timestep_splitting )
       do s= 1,N_aerosol !! nesp_isorropia!(N_aerosol-1)
          jesp=List_species(s)
          dq1dt(j,jesp) = dq1dt(j,jesp) * sub_timestep_splitting / ( 1.d0 - Gamma * Jdq1(j,jesp) * sub_timestep_splitting )
       enddo
    enddo

    call KLIMIT(q2,c_gas,dq1dt,ce_kernal_coef)

    do j = 1 , N_size
       n1(j) = DMAX1 ( 1.d-12*n2(j) , (n2(j)+dn1dt(j)) )
       do s= 1, N_aerosol!! nesp_isorropia!(N_aerosol-1)
          jesp=List_species(s)
          q1(j,jesp) = DMAX1 (1.d-12*q2(j,jesp) , (q2(j,jesp)+dq1dt(j,jesp)) )
       enddo
    enddo

    !     ******3. compute the derivative dq2dt in (q+dq1dt ; Tin2+sub_timestep_splitting)
    !     ******&  the approximation of jacobian diagonal

    current_sub_time = current_sub_time + sub_timestep_splitting

    call fgde(q1,n1,c_gas_t,dq2dt,dn2dt,ce_kernal_coef)

    do j = 1 , N_size
       do s= 1, N_aerosol !! nesp_isorropia!(N_aerosol-1)
          jesp=List_species(s)
          if (isnan(dq2dt(j,jesp))) then
             write(*,*) j,s,list_species(s),jdq1(j,jesp),dq2dt(j,jesp),q1(j,jesp)
             stop
          endif

          if (IsNaN(dq2dt(j,jesp))) then 
             dq2dt(j,jesp) = 0.d0 !NAN
          endif

          Jdq2(j,jesp) = 0.d0
          if ( q1(j,jesp).gt.TINYM .and. dq2dt(j,jesp).lt.0.d0 ) then
             Jdq2(j,jesp) = dq2dt(j,jesp) / q1(j,jesp)
          endif
          Jdn2(j) = 0.d0
          if ( n1(j).gt.TINYN) &
               Jdn2(j) = dn2dt(j) / (n1(j)+TINYN)
       enddo
    enddo

    do s= 1, N_aerosol !!nesp_isorropia!(N_aerosol-1)
       jesp=List_species(s)
       if (aerosol_species_interact(jesp).gt.0) then
          temp1 = c_gas(jesp)
          !          do j = 1, ICUT
          !             temp1 = temp1 + q1(j,jesp)
          !          enddo

          temp2 = 0.d0
          !          do j = (ICUT+1), N_size
          do j = 1, N_size
             if (dq2dt(j,jesp).gt.0.d0) temp2 = temp2 + dq2dt(j,jesp)
          enddo
          do j = 1, N_size
             if ((temp1-q1(j,jesp)).gt.TINYM .and. dq2dt(j,jesp).gt.0.d0) then
                Jdq2(j,jesp) = -temp2  / (temp1-q1(j,jesp))
             endif
          enddo
       endif
    enddo

    !     ******4. second evaluation

    do j = 1 , N_size
       dn2dt(j)=(dn2dt(j)*sub_timestep_splitting- Gamma*sub_timestep_splitting&
	    *(Jdn1(j)+Jdn2(j))*dn1dt(j))&
	    / (1.d0- Gamma*Jdn2(j) *sub_timestep_splitting  )
       do s= 1, N_aerosol !!nesp_isorropia!(N_aerosol-1)
          jesp=List_species(s)
          dq2dt(j,jesp) = (dq2dt(j,jesp)*sub_timestep_splitting- Gamma*sub_timestep_splitting&
               *(Jdq1(j,jesp)+Jdq2(j,jesp))*dq1dt(j,jesp))&
               / (1.d0- Gamma*Jdq2(j,jesp) *sub_timestep_splitting  )
       enddo
    enddo

    call KLIMIT(q2,c_gas,dq2dt,ce_kernal_coef)

    do j = 1 , N_size
       tmp=n2(j)
       n2(j)=DMAX1 ( 1.d-12*n2(j) , (n2(j)+0.5d0*(dn1dt(j)+dn2dt(j))) )
       do s = 1, N_aerosol !! nesp_isorropia
          jesp=List_species(s)
	  tmp=q2(j,jesp)
	  q2(j,jesp)=DMAX1(1.d-12*q2(j,jesp),(q2(j,jesp)+0.5d0*(dq1dt(j,jesp)+dq2dt(j,jesp))) )
	  tmp=q2(j,jesp)-tmp
	  dqdt(j,jesp)=tmp/sub_timestep_splitting
       enddo
    enddo

  end subroutine Ros2_solver

end Module jAdaptstep
