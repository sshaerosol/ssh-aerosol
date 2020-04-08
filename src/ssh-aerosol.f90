!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------

PROGRAM SSHaerosol

  use aInitialization
  use jAdaptstep
  use bCoefficientRepartition
  use mEmissions
  use netcdf
  use lDiscretization
  use Resultoutput
  use gCoagulation
  use mod_photolysis

  implicit none

  integer :: t, j, s,jesp,day
  character (len=40) :: namelist_ssh  ! Configuration file
  double precision, dimension(:), allocatable :: timer

  double precision :: t_since_update_photolysis, t0

  ! Initial time (seconds)
  call cpu_time(t0)

  ! Initialisation: discretization and distribution
  call getarg(1, namelist_ssh)

  ! Read the number of gas-phase species and chemical reactions
  call dimensions(N_gas, n_reaction, n_photolysis)

  call read_namelist(namelist_ssh)

  call read_inputs()                                
 
  call init_parameters()  
 
  call init_distributions()  

  call init_output_conc() 

  call save_report()

  call save_concentration() 

  if ((tag_chem .ne. 0).AND.(option_photolysis.eq.2)) then
    ! Allocate arrays for photolysis
    call allocate_photolysis()    
    ! Read photolysis rate for current day
    current_time = initial_time
    call init_photolysis() 
    call interpol_photolysis()
  endif
  ! **** simulation starts 
  t_since_update_photolysis = 0.d0

  ! Initialization is finished
  allocate(timer(nt+3))
  timer(1) = t0
  call cpu_time(t0)
  timer(2) = t0

  do t = 1, nt

     current_time = initial_time + (t - 1) * delta_t
     t_since_update_photolysis = t_since_update_photolysis +  delta_t

     if (ssh_standalone) write(*,*) "Performing iteration #" // trim(str(t)) // "/" // trim(str(nt))
     if (ssh_logger) write(logfile,*) "Performing iteration #" // trim(str(t)) // "/" // trim(str(nt))

     ! Read the photolysis rates.
     if ((tag_chem .ne. 0).AND.(option_photolysis.eq.2)) then
       if (t_since_update_photolysis >= time_update_photolysis) then
        call interpol_photolysis()
        t_since_update_photolysis = 0.d0
       endif
     endif

     ! Emissions
     if (tag_emis .ne. 0) call emission(delta_t)

     ! Gas-phase chemistry

     ! initial and final physical parameters are set same
     ! no volumetric emission
     ! 0 : vertical gas volumetric emission    1 : with number 
     ! 0 : not take into account cloud    0.d0 : air water content fracion sets to 0  

     if (tag_chem .ne. 0) then
       call chem(n_gas, n_reaction, n_photolysis, photolysis_reaction_index,&
          ns_source, source_index, conversionfactor, conversionfactorjacobian,&
          0, lwc_cloud_threshold, molecular_weight, &
          current_time, attenuation, &
          humidity, temperature,&
          pressure, source, &
          photolysis_rate, delta_t, attenuation,&
          humidity, temperature,&
          pressure, source, &
          photolysis_rate, longitude,&
          latitude, concentration_gas_all,&
          0, with_heterogeneous, n_aerosol, n_size, n_fracmax,&
          0.d0,&
          diam_bound, fixed_density, &
          wet_diameter, &
          heterogeneous_reaction_index, &
          concentration_mass,&
          with_adaptive, adaptive_time_step_tolerance,&
          min_adaptive_time_step, option_photolysis, ind_jbiper, ind_kbiper,&
          1, not(with_fixed_density), concentration_number, &
          mass_density)
      end if

	! re-calculate total_mass(N_aerosol) because mass change due to gas-phase chemistry  
    total_aero_mass = 0.d0
    total_mass = 0.d0
    do s = 1, N_aerosol_layers
       jesp = List_species(s)
       do j=1,N_size
         total_aero_mass(jesp) = total_aero_mass(jesp) + concentration_mass(j,s)
       enddo
    enddo
	! update mass conc. of aerosol precursors
	! concentration_gas_all(precursor_index) -> concentration_gas(n_aerosol)
    do s = 1, N_aerosol
       if (aerosol_species_interact(s) .gt. 0) then
          concentration_gas(s) = concentration_gas_all(aerosol_species_interact(s))
       end if
          total_mass(s) = total_mass(s) + concentration_gas(s) + total_aero_mass(s)
    end do

    ! Aerosol dynamic
    CALL AERODYN(current_time,delta_t)

	! update mass conc. of aerosol precursors
	! concentration_gas(n_aerosol) -> concentration_gas_all(precursor_index)
    do s = 1, N_aerosol
       if (aerosol_species_interact(s) .gt. 0) then
          concentration_gas_all(aerosol_species_interact(s)) = concentration_gas(s)
       end if
    end do

    call save_concentration()         ! Text or Binary format outout

    ! Time step is finished
    call cpu_time(t0)
    timer(t+2) = t0

  end do			! finsh simulation


  !call delete_empty_file() ! delete empty output files

  call free_allocated_memory()
  IF (with_coag.EQ.1) call DeallocateCoefficientRepartition()
 
    ! Desallocate arrays for photolysis
  if ((tag_chem .ne. 0).AND.(option_photolysis.eq.2)) then
    call deallocate_photolysis()    
  endif

  if (ssh_standalone) write(*,*) "============================================"
  if (ssh_standalone) write(*,*) "==== SSH-aerosol simulation completed  ====="
  if (ssh_standalone) write(*,*) "============================================"
  if (ssh_logger) write(logfile,*) "============================================"
  if (ssh_logger) write(logfile,*) "==== SSH-aerosol simulation completed  ====="
  if (ssh_logger) write(logfile,*) "============================================"

  ! Simulation is finished
  call cpu_time(t0)
  timer(nt+3) = t0

  ! Print various times
  if (ssh_standalone) then
    write(*,*) ""
    write(*,*) "Total simulation time in seconds : ", timer(nt+3) - timer(1)
    write(*,*) "Initialization time in seconds : ", timer(2)-timer(1)
    write(*,*) "Average time per time step in seconds : ", sum(timer(3:nt+2) - timer(2:nt+1))/dble(nt)
    write(*,*) "Maximal time per time step in seconds : ", maxval(timer(3:nt+2) - timer(2:nt+1))
    write(*,*) "Minimal time per time step in seconds : ", minval(timer(3:nt+2) - timer(2:nt+1))
  endif
  if (ssh_logger) then
    write(logfile,*) ""
    write(logfile,*) "Total simulation time in seconds : ", timer(nt+3) - timer(1)
    write(logfile,*) "Initialization time in seconds : ", timer(2)-timer(1)
    write(logfile,*) "Average time per time step in seconds : ", sum(timer(3:nt+2) - timer(2:nt+1))/dble(nt)
    write(logfile,*) "Maximal time per time step in seconds : ", maxval(timer(3:nt+2) - timer(2:nt+1))
    write(logfile,*) "Minimal time per time step in seconds : ", minval(timer(3:nt+2) - timer(2:nt+1))
  endif

  ! Free memory
  deallocate(timer)
  
end PROGRAM SSHaerosol
