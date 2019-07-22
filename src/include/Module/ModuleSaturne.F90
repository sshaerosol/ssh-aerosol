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
      use aInitialization, only : free_allocated_memory, with_coag, ssh_logger, close_logger
      use bCoefficientRepartition, only : DeallocateCoefficientRepartition

      implicit none

      integer :: ierr

      ! Free memory
      call free_allocated_memory()
      if (with_coag.eq.1) call DeallocateCoefficientRepartition()

      ! Close log file if needed
      if (ssh_logger) then
        call close_logger()
      endif

    end subroutine cs_finalize

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

    subroutine cs_set_logger(cflag) bind(c, name='cs_set_sshaerosol_logger_')

      use iso_c_binding
      use aInitialization, only : set_logger

      implicit none

      logical(kind=c_bool), intent(in) :: cflag
      logical :: flag

      flag = cflag
      call set_logger(flag)

    end subroutine cs_set_logger

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
! External code can set the SSH-aerosol time step
!
! input : time step in seconds
! =============================================================

    subroutine set_dt(val) bind(c, name='cs_set_sshaerosol_dt_')

      use iso_c_binding
      use aInitialization, only : delta_t, dt

      implicit none

      real(c_double), intent(in) :: val

      delta_t = val
      dt = val ! this is probably useless

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
! External code can set the SSH-aerosol initial time
!
! input : initial time in seconds (GMT, computed from January 1st)
! =============================================================

    subroutine set_initial_t(val) bind(c, name='cs_set_sshaerosol_initial_t_')

      use iso_c_binding
      use aInitialization, only : initial_time

      implicit none

      real(c_double), intent(in) :: val

      initial_time = val

    end subroutine set_initial_t

! =============================================================
!
! External code can get the SSH-aerosol initial time
!
! return value : initial time in seconds (GMT, computed from January 1st)
! =============================================================

    function get_initial_t() bind(c, name='cs_get_sshaerosol_initial_t_')

      use iso_c_binding
      use aInitialization, only : initial_time

      implicit none

      real(c_double) :: get_initial_t

      get_initial_t = initial_time

    end function get_initial_t

! =============================================================
!
! External code can force to update the specific humidity
!
! =============================================================

    subroutine cs_update_humidity() bind(c, name='cs_update_sshaerosol_humidity_')

      use iso_c_binding
      use aInitialization, only : temperature, humidity, pressure, pressure_sat, relative_humidity

      implicit none

      ! This is taken from the subroutine read_namelist
      ! TODO have a dedicated subroutine to avoid duplication of code
      pressure_sat = 611.2d0 * exp(17.67d0 * (temperature - 273.15d0) / (temperature - 29.65d0))

      ! This is taken from the subroutine read_namelist
      ! TODO have a dedicated subroutine to avoid duplication of code
      humidity = 1.d0/(Pressure/(pressure_sat *0.62197d0* Relative_Humidity)-1.d0)

    end subroutine cs_update_humidity

! =============================================================
!
! External code can get the specific humidity
!
! output : specific humidity in kg / kg
! =============================================================

    function cs_get_humidity() bind(c, name='cs_get_sshaerosol_humidity_')

      use iso_c_binding
      use aInitialization, only : humidity

      implicit none

      real(kind=c_double) :: cs_get_humidity
      
      cs_get_humidity = humidity

    end function cs_get_humidity

! =============================================================
!
! External code can set the relative humidity
!
! input : relative humidity in ??? TODO
! =============================================================

    subroutine cs_set_relhumidity(val) bind(c, name='cs_set_sshaerosol_relhumidity_')

      use iso_c_binding
      use aInitialization, only : relative_humidity

      implicit none

      include 'CONST_A.INC' ! needed for threshold_RH_inf and threshold_RH_sup

      real(kind=c_double) :: val

      ! This is taken from the subroutine read_namelist
      relative_humidity = min(max(val, threshold_RH_inf), threshold_RH_sup)

    end subroutine cs_set_relhumidity

! =============================================================
!
! External code can get the relative humidity
!
! output : relative humidity in ??? TODO
! =============================================================

    function cs_get_relhumidity() bind(c, name='cs_get_sshaerosol_relhumidity_')

      use iso_c_binding
      use aInitialization, only : relative_humidity

      implicit none

      real(kind=c_double) :: cs_get_relhumidity

      cs_get_relhumidity = relative_humidity

    end function cs_get_relhumidity

! =============================================================
!
! External code can set the temperature
!
! input : temperature in K
! =============================================================

    subroutine cs_set_temperature(val) bind(c, name='cs_set_sshaerosol_temperature_')

      use iso_c_binding
      use aInitialization, only : temperature

      implicit none

      real(kind=c_double), intent(in) :: val

      temperature = val

    end subroutine cs_set_temperature

! =============================================================
!
! External code can get the temperature
!
! output : temperature in K
! =============================================================

    function cs_get_temperature() bind(c, name='cs_get_sshaerosol_temperature_')

      use iso_c_binding
      use aInitialization, only : temperature

      implicit none

      real(kind=c_double) :: cs_get_temperature

      cs_get_temperature = temperature

    end function cs_get_temperature

! =============================================================
!
! External code can set the pressure
!
! input : pressure in Pa
! =============================================================

    subroutine cs_set_pressure(val) bind(c, name='cs_set_sshaerosol_pressure_')

      use iso_c_binding
      use aInitialization, only : pressure

      implicit none

      real(kind=c_double), intent(in) :: val

      pressure = val

    end subroutine cs_set_pressure

! =============================================================
!
! External code can get the pressure
!
! output : pressure in Pa
! =============================================================

    function cs_get_pressure() bind(c, name='cs_get_sshaerosol_pressure_')

      use iso_c_binding
      use aInitialization, only : pressure

      implicit none

      real(kind=c_double) :: cs_get_pressure

      cs_get_pressure = pressure

    end function cs_get_pressure

! =============================================================
!
! External code can set the pH
!
! input : pH in ??? TODO
! =============================================================

    subroutine cs_set_ph(val) bind(c, name='cs_set_sshaerosol_ph_')

      use iso_c_binding
      use aInitialization, only : ph

      implicit none

      real(kind=c_double), intent(in) :: val

      ph = val

    end subroutine cs_set_ph

! =============================================================
!
! External code can get the pH
!
! output : pH in ??? TODO
! =============================================================

    function cs_get_ph() bind(c, name='cs_get_sshaerosol_ph_')

      use iso_c_binding
      use aInitialization, only : ph

      implicit none

      real(kind=c_double) :: cs_get_ph

      cs_get_ph = ph

    end function cs_get_ph

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
! input : 2D array of concentrations in micrograms / m^3
! =============================================================

    subroutine cs_set_aero_concentration(array) bind(c, name='cs_set_sshaerosol_aero_concentration_')

      use iso_c_binding
      use aInitialization, only : N_size, N_aerosol, concentration_mass

      implicit none

      real(kind=c_double), intent(in), dimension(N_size, N_aerosol) :: array

      concentration_mass(:,:) = array(:,:)

    end subroutine cs_set_aero_concentration

! =============================================================
!
! External code can get the aerosols concentrations
!
! output : 2D array of concentrations in micrograms / m^3
! =============================================================

    subroutine cs_get_aero_concentration(array) bind(c, name='cs_get_sshaerosol_aero_concentration_')

      use iso_c_binding
      use aInitialization, only : N_size, N_aerosol, concentration_mass

      implicit none

      real(kind=c_double), intent(out), dimension(N_size, N_aerosol) :: array

      array(:,:) = concentration_mass(:,:)

    end subroutine cs_get_aero_concentration

! =============================================================
!
! External code can call the chemistry scheme
!
! input : current time in seconds (GMT, computed from January 1st)
! =============================================================

    subroutine cs_call_ssh_gaschemistry(time) bind(c, name='cs_call_sshaerosol_gaschemistry_')

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

    end subroutine cs_call_ssh_gaschemistry

end module SSHSaturne
