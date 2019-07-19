!!-----------------------------------------------------------------------
!!     Copyright (C) 2012-2019, ENPC - EDF R&D - INERIS
!!
!!     This file is part of the SSH-aerosol model.
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
!!
!!     -- DESCRIPTION
!!    This module propose wrappers to interface SSH-aerosol with Code_Saturne
!!
!!-----------------------------------------------------------------------

module SSHSaturne

  implicit none

  contains

! =============================================================
!
! External code can force SSH-aerosol to:
!   read simulation settings from file
!   call initialize functions
!
! input : namelist.ssh file
!   maximum length of input set to 40 chars (cf read_namelist)
! =============================================================

    subroutine cs_initialize(input_namelist_file) bind(c, name='cs_sshaerosol_initialize_')

      use iso_c_binding
      use aInitialization, only : read_namelist, read_inputs, N_gas, n_reaction, n_photolysis
      use lDiscretization, only : init_parameters, init_distributions

      implicit none

      integer, parameter :: size_namelist_file = 40
      character(kind=c_char), intent(in) :: input_namelist_file(size_namelist_file)
      character(len=size_namelist_file) :: namelist_file

      ! N_gas = 93; N_reaction = 206; N_photolysis = 24 
      call dimensions(N_gas, n_reaction, n_photolysis)  

      ! Read SSH simulation settings file
      namelist_file = transfer(input_namelist_file(1:size_namelist_file), &
                               namelist_file)
      call read_namelist(namelist_file)

      ! Read inputs
      call read_inputs()

      ! Initialize parameters
      call init_parameters()

      ! Initialize distributions
      call init_distributions()

    end subroutine cs_initialize

! =============================================================
!
! External code can force SSH-aerosol to:
!   deallocate memory
!   close the log file
!
! =============================================================

    subroutine cs_finalize() bind(c, name='cs_sshaerosol_finalize_')

      use iso_c_binding
      use aInitialization, only : free_allocated_memory, with_coag, ssh_logger, logfile
      use bCoefficientRepartition, only : DeallocateCoefficientRepartition

      implicit none

      integer :: ierr

      ! Free memory
      call free_allocated_memory()
      if (with_coag.eq.1) call DeallocateCoefficientRepartition()

      ! Close log file if needed
      if (ssh_logger) then
        close(unit = logfile, iostat = ierr)
        if (ierr.ne.0) then
          write(*,*) "SSH-aerosol: error when closing log file."
          stop
        endif
      endif

    end subroutine cs_finalize

! =============================================================
!
! External code can set the SSH-aerosol time step
!
! input : time step in seconds
! =============================================================

    subroutine set_dt(time_step) bind(c, name='cs_set_sshaerosol_dt_')

      use iso_c_binding
      use aInitialization, only : delta_t, dt

      implicit none

      real(c_double), intent(in) :: time_step

      delta_t = time_step
      dt = time_step ! this is probably useless

    end subroutine set_dt

! =============================================================
!
! External code can get the SSH-aerosol time step
!
! return value : time step in seconds
! =============================================================

    function get_dt() bind(c, name='cs_get_sshaerosol_dt_')

      use iso_c_binding
      use aInitialization, only : delta_t

      implicit none

      real(c_double) :: get_dt

      get_dt = delta_t

    end function get_dt

! =============================================================
!
! External code can set the flag to declare SSH-aerosol is not running standalone
!
! input : true if running standalone (default), false otherwise
! =============================================================

    subroutine set_standalone(flag) bind(c, name='cs_set_sshaerosol_standalone_')

      use iso_c_binding
      use aInitialization, only : ssh_standalone

      implicit none

      logical(kind=c_bool), intent(in) :: flag
      
      ssh_standalone = flag

    end subroutine set_standalone

! =============================================================
!
! External code can get the flag to check if SSH-aerosol is running standalone
!
! return value : true if running standalone (default), false otherwise
! =============================================================

    function standalone() bind(c, name='cs_get_sshaerosol_standalone_')

      use iso_c_binding
      use aInitialization, only : ssh_standalone

      implicit none

      logical(kind=c_bool) :: standalone
      
      standalone = ssh_standalone

    end function standalone

! =============================================================
!
! External code can set the flag to decide if SSH-aerosol is logging informations
!
! Important: This subroutine must be called before cs_initialize
!
! input : true if logging to a file, false (default) otherwise
! =============================================================

    subroutine set_logger(flag) bind(c, name='cs_set_sshaerosol_logger_')

      use iso_c_binding
      use aInitialization, only : ssh_logger, ssh_logger_file, logfile

      implicit none

      logical(kind=c_bool), intent(in) :: flag
      integer :: ierr
      logical :: log_file_exists
      
      ssh_logger = flag

      ! Create / replace log file if needed
      if (ssh_logger) then
        inquire(file = trim(ssh_logger_file), exist = log_file_exists, iostat = ierr)
        if (ierr.ne.0) then
          write(*,*) "SSH-aerosol: error when inquiring log file."
          stop
        endif
        if (log_file_exists) then
          open(unit = logfile, file = trim(ssh_logger_file), status = "replace", iostat = ierr)
        else
          open(unit = logfile, file = trim(ssh_logger_file), status = "new", iostat = ierr)
        endif
        if (ierr.ne.0) then
          write(*,*) "SSH-aerosol: error when creating / replacing log file."
          stop
        endif
      endif


    end subroutine set_logger

! =============================================================
!
! External code can get the flag to check if SSH-aerosol is logging informations
!
! return value : true if logging to a file, false (default) otherwise
! =============================================================

    function logger() bind(c, name='cs_get_sshaerosol_logger_')

      use iso_c_binding
      use aInitialization, only : ssh_logger

      implicit none

      logical(kind=c_bool) :: logger
      
      logger = ssh_logger

    end function logger

! =============================================================
!
! External code can get the number of gas species
!
! return value : number of gas species
! =============================================================

    function cs_get_ngas() bind(c, name='cs_get_sshaerosol_ngas_')

      use iso_c_binding
      use aInitialization, only : N_gas

      implicit none

      integer(kind=c_int) :: cs_get_ngas
      
      cs_get_ngas = N_gas

    end function cs_get_ngas

! =============================================================
!
! External code can get the number of aerosols species
!
! return value : number of aerosols species
! =============================================================

    function cs_get_naero() bind(c, name='cs_get_sshaerosol_naero_')

      use iso_c_binding
      use aInitialization, only : N_aerosol

      implicit none

      integer(kind=c_int) :: cs_get_naero
      
      cs_get_naero = N_aerosol

    end function cs_get_naero

! =============================================================
!
! External code can get the number of aerosols size bins
!
! return value : number of aerosols size bins
! =============================================================

    function cs_get_nsizebin() bind(c, name='cs_get_sshaerosol_nsizebin_')

      use iso_c_binding
      use aInitialization, only : N_sizebin

      implicit none

      integer(kind=c_int) :: cs_get_nsizebin
      
      cs_get_nsizebin = N_sizebin

    end function cs_get_nsizebin

! =============================================================
!
! External code can set the gaseous concentrations
!
! input : array of concentrations in micrograms / m^3
! =============================================================

    subroutine cs_set_gas_concentration(array) bind(c, name='cs_set_sshaerosol_gas_concentration_')

      use iso_c_binding
      use aInitialization, only : N_gas, concentration_gas_all

      implicit none

      real(kind=c_double), intent(in), dimension(N_gas) :: array

      concentration_gas_all(:) = array(:)

    end subroutine cs_set_gas_concentration

! =============================================================
!
! External code can get the gaseous concentrations
!
! output : array of concentrations in micrograms / m^3
! =============================================================

    subroutine cs_get_gas_concentration(array) bind(c, name='cs_get_sshaerosol_gas_concentration_')

      use iso_c_binding
      use aInitialization, only : N_gas, concentration_gas_all

      implicit none

      real(kind=c_double), intent(out), dimension(N_gas) :: array
      
      array(:) = concentration_gas_all(:)

    end subroutine cs_get_gas_concentration

! =============================================================
!
! External code can set the aerosols concentrations
!
! input : array of concentrations in micrograms / m^3
! =============================================================

    subroutine cs_set_aero_concentration(array) bind(c, name='cs_set_sshaerosol_aero_concentration_')

      use iso_c_binding
      use aInitialization, only : N_aerosol, concentration_gas

      implicit none

      real(kind=c_double), intent(in), dimension(N_aerosol) :: array

      concentration_gas(:) = array(:)

    end subroutine cs_set_aero_concentration

! =============================================================
!
! External code can get the aerosols concentrations
!
! output : array of concentrations in micrograms / m^3
! =============================================================

    subroutine cs_get_aero_concentration(array) bind(c, name='cs_get_sshaerosol_aero_concentration_')

      use iso_c_binding
      use aInitialization, only : N_aerosol, concentration_gas

      implicit none

      real(kind=c_double), intent(out), dimension(N_aerosol) :: array
      
      array(:) = concentration_gas(:)

    end subroutine cs_get_aero_concentration

! =============================================================
!
! External code can call the chemistry scheme
!
! input : current time in seconds (GMT, computed from January 1st)
! =============================================================

    subroutine cs_call_ssh_chemistry(time) bind(c, name='cs_call_sshaerosol_chemistry_')

      use iso_c_binding
      use aInitialization

      implicit none

      real(kind=c_double), intent(in) :: time
      double precision :: current_time

      current_time = time

      call chem(n_gas, n_reaction, n_photolysis, photolysis_reaction_index,&
          ns_source, source_index, conversionfactor, conversionfactorjacobian,&
          0, lwc_cloud_threshold, molecular_weight, &
          current_time, attenuation, &
          humidity, temperature,&
          pressure, source, &
          photolysis, delta_t, attenuation,&
          humidity, temperature,&
          pressure, source, &
          photolysis, longitude,&
          latitude, concentration_gas_all,&
          0, with_heterogeneous, n_aerosol, n_size, n_fracmax,&
          0.d0,&
          diam_bound, fixed_density, &
          wet_diameter, &
          heterogeneous_reaction_index, &
          concentration_mass,&
          with_adaptive, adaptive_time_step_tolerance,&
          min_adaptive_time_step, with_photolysis, ind_jbiper, ind_kbiper,&
          1, not(with_fixed_density), concentration_number, &
          mass_density)

    end subroutine cs_call_ssh_chemistry

end module SSHSaturne
