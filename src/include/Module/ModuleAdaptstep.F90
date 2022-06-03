!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
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
  subroutine ssh_aerodyn(start_time,delta_t)
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
    Double precision emw_tmp
    doUBLE PRECISION start_time,end_time
    double precision :: timestep_coag,timestep_cond    

    !**** SOAP ****
    double precision :: inorg_total, inorg_bin(N_size) 
    integer :: neq
    double precision :: ionic, lwc, proton
    double precision :: lwcorg_Nsize(N_size)
    double precision :: liquid(12),rhoaer,lwcorg
    double precision :: qext(N_aerosol),surface_equilibrium_conc_tmp(N_aerosol)
    double precision :: qinti_tmp(N_inside_aer)
    double precision :: delta_t
    integer :: k,s,b,nb_iter

    double precision :: timestep_splitting,sub_timestep_splitting
    double precision :: initial_time_splitting,current_sub_time,final_sub_time
    double precision :: qaero(N_aerosol)
    double precision :: watorg,organion,qgas(N_aerosol)
    double precision :: dry_density,wet_density,dry_mass,wmass,dry_to_wet

    lwcorg_nsize = 0.d0

    frac_oligo=0.d0
  
    iter_eqconc(:) = 0
    ratio_water(:)=0.d0
    ratio_eqconc(:,:)=0.d0   
    iter_water(:) = 0    

    ! Initialize the density of aerosols
    if (with_fixed_density.eq.0) then

       do j1 = 1, N_size
          call ssh_compute_density(N_size,N_aerosol_layers,EH2O_layers,TINYM,concentration_mass,&
               mass_density_layers,j1,rhoaer)
          rho_wet_cell(j1) = rhoaer
          if(rho_wet_cell(j1).LT.0.1d-6) rho_wet_cell(j1)=density_aer_bin(j1)
       enddo
    else if (with_fixed_density.eq.2) then

       do j1 = 1, N_size
          do jesp=1,N_aerosol
             qext(jesp)=0.d0
          enddo
          do jesp=1,N_aerosol_layers
             s = List_species(jesp)
             qext(s) = qext(s) + concentration_mass(j1,jesp)
          enddo

          call ssh_get_nonlinear_density(qext,dry_density,wet_density,dry_mass,wmass,dry_to_wet)
          rho_wet_cell(j1) = wet_density
          if(rho_wet_cell(j1).LT.0.1d-6) rho_wet_cell(j1)=density_aer_bin(j1)
       enddo
    end if  

    call ssh_update_wet_diameter_liquid(N_size,concentration_mass, concentration_number,&
         wet_mass,wet_diameter,wet_volume,cell_diam_av)

    !call ssh_mass_conservation(concentration_mass,concentration_number,concentration_gas, total_mass)

    if (with_coag.EQ.1) then
       call ssh_COMPUTE_AIR_FREE_MEAN_PATH(Temperature,Pressure,&
            air_free_mean_path,viscosity)
       do j1 = 1, N_size
          do j2 = 1, N_size
             call ssh_compute_bidisperse_coagulation_kernel(Temperature,air_free_mean_path,&
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
             call ssh_COMPUTE_GAS_DIFFUSIVITY(Temperature,Pressure,&
                  molecular_diameter(jesp),emw_tmp,&
                  collision_factor_aer(jesp),diffusion_coef(jesp) ) ! gas diff coef in air

             call ssh_COMPUTE_QUADRATIC_MEAN_VELOCITY(Temperature,&
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

    ! set cut-off diameter if cond/evap is considered and tag_ciut >= 1
    if (with_cond.eq.1.and.tag_icut.ge.1) then
	set_icut = 0 !set tag to recompute icut each time when enter aerodyn()
	icut = 0     !initialise icut to aero
    endif

    do while (initial_time_splitting.lt.delta_t)
       if(splitting == 0) then ! Initial time step when processes are splitted
          if(nb_iter == 0) then
             ! Solve coagulation, condensation/evaporation and nucleation
             ! Compute the characteristic time step of each physical process.
             !call initstep(concentration_mass,concentration_number,&
             !     concentration_gas,timestep_coag,timestep_cond, &
             !     initial_time_splitting,timestep_splitting,delta_t)
             nb_iter = 1
             timestep_coag = DTAEROMIN 
             timestep_cond = DTAEROMIN 
             timestep_splitting = DTAEROMIN 
             timestep_splitting = DMIN1(timestep_splitting,delta_t-initial_time_splitting) 
             if(with_coag.eq.0.OR.with_cond.eq.0)  then 
                timestep_splitting = delta_t-initial_time_splitting 
             endif 
          else
             timestep_coag = DMIN1(timestep_coag,delta_t-current_sub_time)
             timestep_cond = DMIN1(timestep_cond,delta_t-current_sub_time)
             timestep_splitting = DMAX1(timestep_coag,timestep_cond)
             timestep_splitting = DMIN1(timestep_splitting,delta_t-current_sub_time)
          endif
       else  
          !call initstep_coupled(concentration_mass,concentration_number,&
          !     concentration_gas,sub_timestep_splitting, &
          !     initial_time_splitting,timestep_splitting,delta_t)
          sub_timestep_splitting = DTAEROMIN
          sub_timestep_splitting = DMIN1(sub_timestep_splitting,delta_t-initial_time_splitting)
          timestep_splitting = delta_t-initial_time_splitting
       endif

       if(splitting == 0) then ! Split processes
          ! Solve with the slowest process (coagulation)
          if ((with_coag.eq.1).AND.(nlayer.eq.1)) then
             current_sub_time = initial_time_splitting
             final_sub_time = current_sub_time+timestep_splitting
             tag_coag = 1
             tag_cond = 0
             tag_nucl = 0
             solver=1            ! only etr for coagulation
             call ssh_processaero(solver,current_sub_time,timestep_coag,final_sub_time,splitting)
          endif
          if(N_fracmax.gt.1) then
             call ssh_redistribution_fraction()!fraction redistribution
          endif

          ! Solve the fastest process (condensation/evaporation and nucleation).
          if (with_cond+with_nucl.gt.0) then
             current_sub_time = initial_time_splitting
             final_sub_time = current_sub_time + timestep_splitting
             tag_coag = 0
             tag_cond = with_cond
             tag_nucl = with_nucl
             solver=dynamic_solver
             call ssh_processaero(solver,current_sub_time,timestep_cond,final_sub_time,splitting)
          endif

       else
          current_sub_time=initial_time_splitting
          final_sub_time = current_sub_time+timestep_splitting
          tag_coag = with_coag
          tag_cond = with_cond
          tag_nucl = with_nucl
          solver = 1            ! only etr for coagulation
          call ssh_processaero(solver,current_sub_time,sub_timestep_splitting,final_sub_time,splitting)
       endif

       if (with_cond+with_nucl+with_coag.eq.0) then
          final_sub_time = initial_time_splitting + timestep_splitting !pure emission
       endif

       initial_time_splitting = final_sub_time

    enddo
   
    !do C/E equilibrium after each emission
    if (with_cond.gt.0) then
       if(ICUT.ge.1.and.soap_inorg==0) then ! use bulk equilibrium method for inorganics in section <= ICUT

          call ssh_bulkequi_inorg(nesp_isorropia,& 
               lwc, ionic, proton, liquid) !equlibrium for inorganic
          !call ssh_mass_conservation(concentration_mass,concentration_number,&
          !                       concentration_gas, total_mass)

          call ssh_redistribution_lwc(lwc,ionic,proton,liquid,0,ICUT)

       endif

       if (ISOAPDYN.eq.0.or.soap_inorg==1) then !.and.nlayer.eq.1) then

          if (ISOAPDYN.eq.1.and.soap_inorg==1) then
             soap_inorg_loc=-1               !If ISOAPDYN=1 and soap_inorg=1, I just want to compute inorganics in this part
          else
             soap_inorg_loc=soap_inorg
          endif
          
          ! ******** equilibrium SOA even if inorganic aerosols are estimated dynamically
          if (ICUT_org>0) then

             call  ssh_bulkequi_org(nesp_aec,lwc,lwcorg,ionic,proton,liquid,delta_t)!equilibrium for organic
             !call ssh_mass_conservation(concentration_mass,concentration_number,&
             !     concentration_gas, total_mass)

             if (soap_inorg==0) then
                call ssh_redistribution_lwcorg(lwcorg,lwcorg_Nsize)
             else
                call ssh_redistribution_lwc(lwc,ionic,proton,liquid,0,ICUT_org)
                call ssh_redistribution_lwcorg(lwcorg,lwcorg_Nsize)
             endif
          endif

          if (ICUT_org+1<=N_size) then
             call SSH_SOAP_DYN_ICUT(Relative_Humidity,&
                  ionic, proton, lwc,lwcorg,&
                  Temperature, delta_t,&
                  cell_diam_av(ICUT_org+1:N_size), neq, liquid)
          endif
       endif
       
       if (ISOAPDYN.eq.1) then
          if (soap_inorg==1) then
             soap_inorg_loc=0      !only compute organic
          endif

          organion=0.d0
          watorg=0.d0
          proton=0.d0
          ionic=0.d0
          lwc=0.d0
          qaero=0.d0
          do jesp=1,N_aerosol
             do j=1,N_size
                qaero(jesp) = qaero(jesp) + concentration_mass(j,jesp)
             enddo
             qgas(jesp)=concentration_gas(jesp)
          enddo
          
          jesp=isorropia_species(2)
          qgas(jesp)=0.d0
          qgas(EH2O)=0.0
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
          qaero(ENa) = 0.d0
          qaero(ECl) = 0.d0
          qgas(ECl) = 0.d0
#endif
          call ssh_isoropia_drv(N_aerosol,&
               qaero,qgas,organion, watorg, ionic, proton, lwc, Relative_Humidity, Temperature, &
               liquid)
          call ssh_redistribution_lwc(lwc,ionic,proton,liquid,1,N_size)

          ! *** SOA are dynamically partitioned even if inorganic aerosols are estimated by equilibrium.
          !
          call ssh_SOAP_DYN(Relative_Humidity,&
               ionic, proton, lwc,lwcorg,&
               Temperature, delta_t,&
               cell_diam_av, neq, liquid)
       endif

       call ssh_update_wet_diameter_liquid(N_size,concentration_mass, &
            concentration_number,wet_mass,wet_diameter,wet_volume,cell_diam_av)

       if(N_fracmax.gt.1 ) then !.and. redistribution_method.ne.0) then
          call ssh_redistribution_fraction()!fraction redistribution
       endif
       
       call ssh_compute_average_diameter() ! Update average_diame
       if (redistribution_method.ne.0) call ssh_redistribution_size(redistribution_method)!size redistribution

    endif

    call ssh_update_wet_diameter_liquid(N_size,concentration_mass,concentration_number,&
         wet_mass,wet_diameter,wet_volume,cell_diam_av) 

    ! AVOID mass conservation after redistribution
    ! call mass_conservation(concentration_mass,concentration_number,&
    !   concentration_gas, total_mass)

  end subroutine ssh_aerodyn

  subroutine ssh_processaero(solver,current_sub_time,sub_timestep_splitting,&
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
    double precision:: time_t,time_step_sulf
    double precision:: final_sub_time,current_sub_time,sub_timestep_splitting
    ! update ICUT
    integer:: j,s,jesp,icut_tmp
    double precision:: tau(N_size),cond_time(N_size,3),tmp,tmp1
    !integer :: tag_cond_save

    cond_time = 0.0
    
    ! set icut=0 in case of processing coagulation under split numerical scheme
    if (splitting.eq.0.and.tag_coag.eq.1.and.tag_cond.eq.0) then
	icut_tmp = icut
	icut = 0 !set icut to zero when solve coagulation splitted
    endif

    time_t=final_sub_time-current_sub_time
    !    ******Dynamic time loop
    !tag_cond_save=tag_cond
    !if (soap_inorg==1) then
    !   tag_cond=0
    !endif
    do while ( current_sub_time .lt. final_sub_time )
       time_step_sulf = 0.0
       !if (tag_cond.eq.1) then
          !!if (tag_emis .ne. 0) call ssh_emission(sub_timestep_splitting)
          !     Compute gas mass conservation.
       !   call ssh_mass_conservation(concentration_mass,concentration_number,&
       !        concentration_gas, total_mass)
       !endif

       if (solver.eq.0) then
          ! euler solver used for coagulation
          call ssh_Euler_solver(concentration_mass,concentration_number,concentration_gas,dqdt,&
               current_sub_time,sub_timestep_splitting)

       elseif (solver.eq.1) then
          ! etr dynamic solver
          call ssh_Etr_solver(concentration_mass_tmp,concentration_mass,&
               concentration_number_tmp,concentration_number,concentration_gas,dqdt,&
               current_sub_time,sub_timestep_splitting,cond_time)
       elseif (solver.eq.2) then
          ! ros2 dynamic solver used for condensation/evaporation of inorganics
          call ssh_Ros2_solver(concentration_mass_tmp,concentration_mass,&
               concentration_number_tmp,concentration_number,&
               concentration_gas,dqdt,&
               current_sub_time,sub_timestep_splitting)
       endif
       
       !call ssh_mass_conservation(concentration_mass,concentration_number,&
       !           concentration_gas, total_mass)

       if (current_sub_time.le.final_sub_time) then

	  ! compute criteria for icut computation if need
	  if ((tag_cond .eq. 1).and.(set_icut.eq.0)) then

	     ! c/e characteristic timescale criteria
	     if (tag_icut.eq.1) then
		!init tau to a huge value
	        tau = 1d99
		! compute tau
		do j=1,N_size
		   do s=1,3 ! only for inorganics ENH4, ENO3, ECL
                      if(cond_time(j,s).gt.TINYM) tau(j) = dmin1(tau(j),cond_time(j,s)) !sectional c/e characteriatic timescales
		   enddo
		enddo
		! set sections with no timescale back to 0
	        do j=1,N_size
		  if (tau(j) .eq. 1d99) tau(j)=0.d0
	        enddo

	      !weighted QSSA criteria
	      elseif (tag_icut.eq.3) then
		!init tau
	        tau = 0.d0
		do j=1,N_size
                   tmp = 0.d0 !weighted QSSA
                   tmp1 =0.d0 !total mass
	           do s=1,3
		  	jesp=cond_time_index(s)
                        if(cond_time(j,s).gt.TINYM.and.concentration_mass(j,jesp).GT.TINYM) then
			   tmp = tmp + cond_time(j,s)*concentration_mass(j,jesp)
			   tmp1  = tmp1 + concentration_mass(j,jesp)
                        endif
		   enddo
                  if (tmp1.GT.TINYM) tau(j) = tmp/tmp1
		enddo
	     endif

	  endif


             call ssh_adaptime(concentration_mass_tmp,concentration_mass,concentration_number_tmp,&
                  concentration_number,sub_timestep_splitting,time_step_sulf,current_sub_time,&
                  final_sub_time,tau)
             if((tag_nucl.EQ.1).OR.((tag_cond.EQ.1).AND.(tag_coag.EQ.1))) then 
             ! Need to redistribute onto fixed grid if nucleation is solved with 
             ! condensation/evaporation or if processes are not splitted
               call ssh_update_wet_diameter_liquid(N_size,concentration_mass, &
                       concentration_number,wet_mass,wet_diameter,wet_volume,cell_diam_av)
               if(N_fracmax.gt.1) then
                    call ssh_redistribution_fraction()!fraction redistribution
               endif
                 
               call ssh_compute_average_diameter() ! Update average_diame
               if (redistribution_method.ne.0) call ssh_redistribution_size(redistribution_method)!size redistribution    
             endif     
       endif
    end do

    ! restore icut in case of processing coagulation under split numerical scheme
    if (tag_coag.eq.1.and.tag_cond.eq.0.and.splitting.eq.0) icut = icut_tmp
    !if (soap_inorg==1) then
    !   tag_cond=tag_cond_save
    !endif
    
  end subroutine ssh_processaero

  subroutine ssh_initstep_coupled(c_mass,c_number,c_gas,sub_time_splitting,&
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
    double precision:: c_mass(N_size,N_aerosol_layers)
    double precision:: dqdt1(N_size,N_aerosol_layers)
    double precision:: c_number(N_size),qH2O(N_size)
    double precision:: dndt1(N_size)
    double precision:: c_gas(N_aerosol),cond_time(N_size,3)
    double precision:: sub_time_splitting
    double precision:: time_splitting,initial_time_splitting
    double precision:: tmp,tscale
    double precision:: t_total!total time of the simunlation
    double precision :: rate,ce_kernal_coef(N_size, N_aerosol)

    !     Coagulation time step.

    sub_time_splitting=t_total-initial_time_splitting
    
    tag_coag = with_coag
    tag_cond = with_cond
    tag_nucl = with_nucl
    call ssh_fgde(c_mass,c_number,c_gas,dqdt1,dndt1,ce_kernal_coef,qH2O,cond_time,0)

    do j=1,N_size
       tmp=c_number(j)*dndt1(j)
       if (tmp.ne.0.D0) then
          tscale = 0.1d0 * c_number(j)/DABS(dndt1(j))
          sub_time_splitting = DMIN1(sub_time_splitting,tscale)
       endif
       do jesp= 1, N_aerosol_layers
          if(jesp.NE.EH2O_layers) then
             tmp=c_mass(j,jesp)*dqdt1(j,jesp)
             if (tmp.ne.0.D0.and.c_mass(j,jesp).gt.TINYM) then
                tscale = 0.1d0 * c_mass(j,jesp)/DABS(dqdt1(j,jesp))
                sub_time_splitting = DMIN1(sub_time_splitting,tscale)
             endif
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

  end subroutine ssh_initstep_coupled

  subroutine ssh_initstep(c_mass,c_number,c_gas,time_coag,time_cond,&
       initial_time_splitting,time_splitting,t_total)
    !------------------------------------------------------------------------
    !
    !     -- DESCRIPTION
    !     This subroutine ssh_performs the initialization of the timestep
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
    double precision:: c_mass(N_size,N_aerosol_layers)
    double precision:: dqdt1(N_size,N_aerosol_layers)
    double precision:: c_number(N_size),qH2O(N_size)
    double precision:: dndt1(N_size)
    double precision:: c_gas(N_aerosol)
    double precision:: time_coag,cond_time(N_size,3)
    double precision:: time_cond
    double precision:: time_splitting,initial_time_splitting
    double precision:: tmp,tscale
    double precision:: t_total!total time of the simunlation
    double precision :: rate,ce_kernal_coef(N_size, N_aerosol)

    !     Coagulation time step.

    time_coag=t_total-initial_time_splitting
    time_cond=t_total-initial_time_splitting

    if ((with_coag.eq.1).AND.(nlayer.eq.1)) then
       tag_coag=1
       tag_cond=0
       tag_nucl=0

       call ssh_fgde(c_mass,c_number,c_gas,dqdt1,dndt1,ce_kernal_coef,qH2O,cond_time,0)

       do j=1,N_size
          tmp=c_number(j)*dndt1(j)
          if (tmp.ne.0.D0) then
             tscale=c_number(j)/DABS(dndt1(j))
             time_coag=DMIN1(time_coag,tscale)
          endif
          do jesp = 1,N_aerosol_layers
             if(jesp.NE.EH2O_layers) then
                tmp=c_mass(j,jesp)*dqdt1(j,jesp)
                if (tmp.ne.0.D0.and.c_mass(j,jesp).gt.TINYM) then
                   tscale=c_mass(j,jesp)/DABS(dqdt1(j,jesp))
                   time_coag=DMIN1(time_coag,tscale)
                endif
             endif
          enddo
       end do
       time_coag = DMIN1(time_coag,t_total-initial_time_splitting)
    endif

    !     Cond/evap time step.
    if (with_cond.eq.1.OR.with_nucl.eq.1) then
       tag_coag=0
       tag_cond=with_cond
       tag_nucl=with_nucl
       call ssh_fgde(c_mass,c_number,c_gas,dqdt1,dndt1,ce_kernal_coef,qH2O,cond_time,0)
       do j=1,N_size  !Loop from 1 in case of nucleation - if bins at equilibrium then dqdt1 = 0 from fgde
	  if (DABS(dndt1(j)).gt.0.d0.and.c_number(j).gt.TINYN) then
             tscale=c_number(j)/DABS(dndt1(j))
             time_cond=DMIN1(time_cond,tscale)
	  endif
          do jesp = 1,N_aerosol_layers
             if(jesp.NE.EH2O_layers) then
                tmp=c_mass(j,jesp)*dqdt1(j,jesp)
                if (DABS(dqdt1(j,jesp)).gt.0.d0.and.c_mass(j,jesp).gt.TINYM) then
                   tscale=c_mass(j,jesp)/DABS(dqdt1(j,jesp))
                   time_cond=DMIN1(time_cond,tscale)
                endif
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
    time_cond = DMIN1(time_splitting,DMAX1(time_cond,DTAEROMIN))
    time_coag = DMIN1(time_splitting,DMAX1(time_coag,DTAEROMIN))

    time_splitting=DMIN1(time_splitting,t_total-initial_time_splitting)
    time_splitting = DMAX1(time_splitting,DTAEROMIN)


    ! if only one process then set
    ! splitting step to whole timestep
    ! or if processes are solved coupled
    if(with_coag.eq.0.OR.with_cond.eq.0)  then
       time_splitting = t_total-initial_time_splitting
    endif
!    if((ICUT.eq.N_size).AND.(with_nucl.EQ.0)) then 
!       time_splitting=t_total-initial_time_splitting
!       time_cond=time_splitting
!    endif

  end subroutine ssh_initstep

  subroutine ssh_adaptime(q1,q2,n1,n2,T_dt,time_step_sulf,current_sub_time,&
       final_sub_time,tau)
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
    double precision :: q1(N_size,N_aerosol_layers)!tmp mass concentration
    double precision :: q2(N_size,N_aerosol_layers)!tmp mass concentration
    double precision :: n1(N_size)!1th order number concentration
    double precision :: n2(N_size)!2d order number concentration
    double precision :: tmp,n2err,R
    double precision :: T_dt,time_step_sulf
    double precision ::  current_sub_time,final_sub_time

    ! update ICUT
    integer k,f,ind,l,icut_tmp
    double precision :: tau(N_size),criteria
    double precision :: n2errg,c_gas1(N_aerosol),c_gas2(N_aerosol)

    ! ICUT computation 
    if (set_icut.eq.0.and.tag_cond.eq.1) then
	!init icut
	icut = 0
	!init criteria
	criteria = Cut_dim

	! characteristic timescale or weighted QSSA criteria
	if (tag_icut.eq.1.or.tag_icut.eq.3) then
		do k=1, N_size
			if (tau(k).lt.criteria.and.(k-1).eq.icut) icut = k
		enddo

	! ETR solver efficientcy criteria
	elseif (tag_icut .eq. 2) then 
	      tau = 0d0 !initialization
	      ! compute the error from large aerosol section to small aerosol section
  	      do j=1,N_size
	   	 ind = N_size - j + 1 ! section index
		 tmp = 0.0
		 do s= 1,N_aerosol_layers
			if(s.NE.EH2O_layers) then! Do not consider water for time step
			    if(q2(ind,s).ne.q1(ind,s).and.q2(ind,s).gt.TINYM) then
    			        tmp = tmp + ((q2(ind,s)-q1(ind,s))/q2(ind,s))**2
			    endif
			endif
		 enddo
		 ! sum errors
		 if (ind .ne. N_size) then
			tau(ind) = tau(ind+1)+tmp
		 else
			tau(ind) = tmp
		 endif
	      enddo

	      ! in case no mass in the last few bins
	      ind = N_size
	      do j = 1, N_size
		 k = N_size - j + 1
		 if (tau(k).eq.0.d0.and.ind.eq.k) ind = k-1
	      enddo

	      !compute ICUT
	      do k=1,N_size
		 tmp = (tau(k)/tau(ind))**0.25 !ETR criteria
		 !tmp = T_dt*DSQRT(EPSER)/(tau(k)**0.25) This is the time scale
		 if (tmp.gt.criteria.and.icut.eq.(k-1)) icut = k
	      enddo

	endif

	!error checking to see if the new icut can be accepted
	if (icut.gt.0) then

		! get gas-phase conc.
		call ssh_mass_conservation(q1,n1,c_gas1, total_mass)
		call ssh_mass_conservation(q2,n2,c_gas2, total_mass)
		! error for gas-phase NH3, HCL and HNO3
		n2errg=0.d0
		do jesp=1,3
		   s=cond_time_index(jesp)
		   if(c_gas2(s).gt.TINYM) then
		      	tmp=(c_gas2(s)-c_gas1(s))/c_gas2(s)
		      	n2errg=n2errg+tmp*tmp
		   endif
		end do

		! error for particle ENH4, ENO3 and ECL in bin<=ICUT
		n2err = 0.d0
		do j=1,N_size
		   if(concentration_index(j, 1) <= ICUT) then
			do s=1,3
		   		jesp=cond_time_index(s)
				if(q2(j,jesp).gt.TINYM) then
		                	tmp=(q2(j,jesp)-q1(j,jesp))/q2(j,jesp)
		                	n2err=n2err+tmp*tmp
				endif
		        enddo
		   endif
		 enddo

		! Do not accept icut if it leads to lower time step, that means larger variations, than those of the gas phase
		if(n2err.ge.n2errg) then
			icut = 0
		else
			set_icut=1 !set tag so in this time step do not compute icut again
		endif
	endif
	! set_icut=1 !do not recompute icut even icut=0
    endif

    !     ******zero init
    n2err=0.d0
    !     ******local error estimation
    do j=1,N_size
     !if(concentration_index(j, 1) > ICUT) then
       if(n2(j).gt.TINYN) then
          tmp=(n2(j)-n1(j))/(n2(j)+TINYN)
          n2err=n2err+tmp*tmp
       endif
     !endif
    end do
    do j=1,N_size
     !if(concentration_index(j, 1) > ICUT) then
       do jesp= 1,N_aerosol_layers ! Do not consider water for time step
          if(jesp.NE.EH2O_layers) then
            if(q2(j,jesp).gt.TINYM) then
               tmp=(q2(j,jesp)-q1(j,jesp))/(q2(j,jesp))
               n2err=n2err+tmp*tmp
            endif
          endif
       enddo
     !endif
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

  end subroutine ssh_adaptime

  subroutine ssh_Etr_solver(q1,q2,n1,n2,c_gas,dqdt,current_sub_time,sub_timestep_splitting,cond_time)
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
    double precision:: dqdt(N_size,N_aerosol_layers)
    double precision:: dq1dt(N_size,N_aerosol_layers)
    double precision:: dq2dt(N_size,N_aerosol_layers)
    double precision:: dn1dt(N_size)
    double precision:: dn2dt(N_size),qH2O(N_size)
    double precision:: c_gas(N_aerosol)
    double precision:: c_gas_t(N_aerosol),cond_time(N_size,3)
    double precision:: t_mass(N_aerosol)
    double precision:: dtetr,tmp,current_sub_time,sub_timestep_splitting
    double precision:: n1(N_size)!1th order number concentration
    double precision:: n2(N_size)!2d order number concentration
    double precision:: q1(N_size,N_aerosol_layers)!1st order mass concentration
    double precision:: q2(N_size,N_aerosol_layers)!2d order mass concentration

    qh2o = 0.0
    
    !for condensation or coagulation
   
    call ssh_fgde(q2,n2,c_gas,dq1dt,dn1dt,ce_kernal_coef,qH2O,cond_time,0)
    
    !     First step
    do j=1,N_size
       if(n2(j)+dn1dt(j)*sub_timestep_splitting.GE.TINYN) then
          n1(j)=n2(j)+sub_timestep_splitting*dn1dt(j)
       else
          n1(j) = 0.0
       endif
       do s=1,N_aerosol_layers
          jesp=List_species(s)
          t_mass(jesp)=total_mass(jesp)
          if(q2(j,s)+dq1dt(j,s)*sub_timestep_splitting.GE.TINYM) then
             q1(j,s)=q2(j,s)+sub_timestep_splitting*dq1dt(j,s)
          else
             q1(j,s) = 0.0
          endif
       enddo
!       q1(j,N_aerosol_layers) = qH2O(j) ! water is the last species
    enddo

    call ssh_mass_conservation(q1,n1,c_gas_t, total_mass)
    !     Second step
    call ssh_fgde(q1,n1,c_gas_t,dq2dt,dn2dt,ce_kernal_coef,qH2O,cond_time,0)
    ! KS: call fgde with the last argument = 1 if surface equilibrium concentrations
    ! are not recomputed save CPU time 
    !call ssh_fgde(q1,n1,c_gas_t,dq2dt,dn2dt,ce_kernal_coef,qH2O,cond_time,1)

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
       do s=1,N_aerosol_layers
          tmp=dtetr*(dq1dt(j,s)+dq2dt(j,s))
          if(tmp+q2(j,s).GE.TINYM) then
             dqdt(j,s)=tmp/sub_timestep_splitting
             q2(j,s)=q2(j,s)+tmp
          else
             q2(j,s) = 0.0
             !   write(*,*) 'pb ETR q2 < 0',q2(j,jesp),sub_timestep_splitting*dq2dt(j,jesp),sub_timestep_splitting
          endif
       enddo
!       q2(j,N_aerosol_layers) = qH2O(j)!water is the last species
    enddo

    call ssh_mass_conservation(q2,n2,c_gas, total_mass) !! YK

  end subroutine ssh_Etr_solver

  subroutine ssh_Euler_solver(c_mass,c_number,c_gas,dqdt,current_sub_time,sub_timestep_splitting)
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
    double precision:: dqdt(N_size,N_aerosol_layers)
    double precision:: dndt(N_size),qH2O(N_size),cond_time(N_size,3)
    double precision:: c_number(N_size)!number concentration
    double precision:: c_mass(N_size,N_aerosol_layers)!micg/m^-3
    double precision:: c_gas(N_aerosol)!micg/m^-3
    double precision:: tmp,current_sub_time,sub_timestep_splitting
    double precision :: ce_kernal_coef(N_size,N_aerosol)

    !for condensation
    call ssh_fgde(c_mass,c_number,c_gas,dqdt,dndt,ce_kernal_coef,qH2O,cond_time,0)

    do j=1,N_size
       if(c_number(j)+dndt(j)*sub_timestep_splitting.GE.0.d0) then
          c_number(j) = c_number(j)+dndt(j)*sub_timestep_splitting
       else
          c_number(j)= 0.d0
       endif
       do s=1,N_aerosol_layers
          tmp=sub_timestep_splitting*dqdt(j,s)
          if(tmp + c_mass(j,s).GE.0.d0) then
             c_mass(j,s)=c_mass(j,s)+tmp
          endif
       enddo
!       c_mass(j,N_aerosol_layers) = qH2O(j) !water is the last species
     enddo

    current_sub_time = current_sub_time + sub_timestep_splitting
  end subroutine ssh_Euler_solver

  subroutine ssh_Ros2_solver(q1,q2,n1,n2,c_gas,dqdt,current_sub_time,sub_timestep_splitting)
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
    double precision:: dqdt(N_size,N_aerosol_layers)
    double precision:: dq1dt(N_size,N_aerosol_layers)
    double precision:: dq2dt(N_size,N_aerosol_layers)
    double precision:: dn1dt(N_size)
    double precision:: dn2dt(N_size)
    double precision:: Jdq1(N_size,N_aerosol_layers)
    double precision:: Jdq2(N_size,N_aerosol_layers)
    double precision:: Jdn1(N_size)
    double precision:: Jdn2(N_size),qH2O(N_size)
    double precision:: c_gas(N_aerosol)
    double precision:: c_gas_t(N_aerosol)
    double precision:: tmp,cond_time(N_size,3)
    double precision:: q1(N_size,N_aerosol_layers)!1th order mass concentration
    double precision:: q2(N_size,N_aerosol_layers)!2d order mass concentration
    double precision:: n1(N_size)!1th order number concentration
    double precision:: n2(N_size)!2d order number concentration
    double precision:: Gamma,temp1,temp2
    double precision:: current_sub_time,sub_timestep_splitting
    parameter ( Gamma= 1.7071D0)

    !for condensation
    call ssh_fgde(q2,n2,c_gas,dq1dt,dn1dt,ce_kernal_coef,qH2O,cond_time,0)!compute first order derivative

    !     Every dynamical variable protected against vanishing
    do j = 1 , N_size
       do jesp= 1, N_aerosol_layers !!nesp_isorropia!(N_aerosol-1)
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
    do jesp= 1, N_aerosol_layers !! nesp_isorropia!(N_aerosol-1)
       if (aerosol_species_interact(List_species(jesp)).gt.0) then
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
       do jesp= 1,N_aerosol_layers !! nesp_isorropia!(N_aerosol-1)
          dq1dt(j,jesp) = dq1dt(j,jesp) * sub_timestep_splitting / ( 1.d0 - Gamma * Jdq1(j,jesp) * sub_timestep_splitting )
       enddo
    enddo

    call ssh_KLIMIT(q2,c_gas,dq1dt,ce_kernal_coef)

    do j = 1 , N_size
       n1(j) = DMAX1 ( 1.d-12*n2(j) , (n2(j)+dn1dt(j)) )
       do jesp= 1, N_aerosol_layers !! nesp_isorropia!(N_aerosol-1)
          q1(j,jesp) = DMAX1 (1.d-12*q2(j,jesp) , (q2(j,jesp)+dq1dt(j,jesp)) )
       enddo
    enddo
!    do j = 1 , N_size
!       q1(j,N_aerosol_layers) = qH2O(j) !water is the last species
!    enddo


    !     ******3. compute the derivative dq2dt in (q+dq1dt ; Tin2+sub_timestep_splitting)
    !     ******&  the approximation of jacobian diagonal

    current_sub_time = current_sub_time + sub_timestep_splitting

    call ssh_fgde(q1,n1,c_gas_t,dq2dt,dn2dt,ce_kernal_coef,qH2O,cond_time,1)

    do j = 1 , N_size
       do jesp= 1, N_aerosol_layers !! nesp_isorropia!(N_aerosol-1)
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

    do jesp= 1, N_aerosol_layers !!nesp_isorropia!(N_aerosol-1)
       s = List_species(jesp)
       if (aerosol_species_interact(s).gt.0) then
          temp1 = c_gas(s)
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
       do jesp= 1, N_aerosol_layers !!nesp_isorropia!(N_aerosol-1)
          dq2dt(j,jesp) = (dq2dt(j,jesp)*sub_timestep_splitting- Gamma*sub_timestep_splitting&
               *(Jdq1(j,jesp)+Jdq2(j,jesp))*dq1dt(j,jesp))&
               / (1.d0- Gamma*Jdq2(j,jesp) *sub_timestep_splitting  )
       enddo
    enddo

    call ssh_KLIMIT(q2,c_gas,dq2dt,ce_kernal_coef)

    do j = 1 , N_size
       tmp=n2(j)
       n2(j)=DMAX1 ( 1.d-12*n2(j) , (n2(j)+0.5d0*(dn1dt(j)+dn2dt(j))) )
       do jesp = 1, N_aerosol_layers !! nesp_isorropia
	  tmp=q2(j,jesp)
	  q2(j,jesp)=DMAX1(1.d-12*q2(j,jesp),(q2(j,jesp)+0.5d0*(dq1dt(j,jesp)+dq2dt(j,jesp))) )
	  tmp=q2(j,jesp)-tmp
	  dqdt(j,jesp)=tmp/sub_timestep_splitting
       enddo
    enddo
!    do j = 1 , N_size
!       q2(j,N_aerosol_layers) = qH2O(j)  !water is the last species
!    enddo

  end subroutine ssh_Ros2_solver

end Module jAdaptstep
