!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------
!!
!!     -- DESCRIPTION
!!    This module propose wrappers to interface SSH-aerosol with external tools
!!
!!-----------------------------------------------------------------------

module SSHaerosolAPI

  implicit none

  contains

! =============================================================
!
! External code can use a simplified initialization for SSH
!
! input : namelist.ssh file
!   maximum length of input set to 400 chars (cf read_namelist)
! =============================================================

    subroutine ssh_api_simple_initialize(input_namelist_file, ngas, &
         n_aerosol_layers, nsize) &
         bind(c, name='api_sshaerosol_simple_initialize_')

      use iso_c_binding

      implicit none

      integer, parameter :: size_namelist_file = 400
      character(kind=c_char), intent(in) :: input_namelist_file(size_namelist_file)
      integer(kind=c_int), intent(out) :: ngas, n_aerosol_layers, nsize
      logical(kind=c_bool), parameter :: ok=.true.
      logical(kind=c_bool), parameter :: nope=.false.
      real(c_double) :: time

      call ssh_set_standalone(nope)
      call ssh_api_set_logger(ok)
      call ssh_api_initialize(input_namelist_file)
      call ssh_api_init_distributions()
      ngas = api_get_ngas()
      n_aerosol_layers = api_get_n_aerosol_layers()
      nsize = api_get_nsize()
      call ssh_api_call_ssh_initoutput()
      call ssh_api_call_ssh_report()
      call ssh_api_call_ssh_output()
      call ssh_api_call_ssh_initphoto()
      time = get_initial_t()
      call ssh_set_current_t(time)

    end subroutine ssh_api_simple_initialize

! =============================================================
!
! External code can use a simplified time step for SSH
!
! input : time step
! =============================================================

    subroutine ssh_api_simple_dt(dt) bind(c, name='api_sshaerosol_simple_dt_')

      use iso_c_binding
      use aInitialization, only : logfile, ssh_standalone, ssh_logger

      implicit none

      real(c_double), intent(in) :: dt
      real(c_double) :: time

      if (.not.set_dt(dt)) then
        if (ssh_standalone) write(*,*) "Error. Decrease the time step."
        if (ssh_logger) write(logfile,*) "Error. Decrease the time step."
        stop
      endif
      time = get_current_t()
      call ssh_set_current_t(time + dt)
      call ssh_api_call_ssh_updatephoto()
      call ssh_api_call_ssh_emission()
      call ssh_api_call_ssh_gaschemistry()
      call ssh_api_call_ssh_aerochemistry()
      call ssh_api_call_ssh_output()

    end subroutine ssh_api_simple_dt

! =============================================================
!
! External code can force SSH-aerosol to:
!   read simulation settings from file
!   call initialize functions
!
! input : namelist.ssh file
!   maximum length of input set to 400 chars (cf read_namelist)
! =============================================================

    subroutine ssh_api_initialize(input_namelist_file) bind(c, name='api_sshaerosol_initialize_')

      use iso_c_binding
      use aInitialization, only : ssh_read_namelist, ssh_read_inputs, N_gas, n_reaction, n_photolysis
      use lDiscretization, only : ssh_init_parameters, ssh_init_distributions
      use mod_meteo, only : ssh_read_meteo

      implicit none

      integer, parameter :: size_namelist_file = 400
      character(kind=c_char), intent(in) :: input_namelist_file(size_namelist_file)
      character(len=size_namelist_file) :: namelist_file
      
      ! N_gas = 93; N_reaction = 206; N_photolysis = 24 
      call ssh_dimensions(N_gas, n_reaction, n_photolysis)  

      ! Read SSH simulation settings file
      namelist_file = transfer(input_namelist_file(1:size_namelist_file), &
                               namelist_file)
      call ssh_read_namelist(namelist_file)

      ! Read inputs
      call ssh_read_inputs()

      ! Read meteo
      call ssh_read_meteo()

      ! Initialize parameters
      call ssh_init_parameters()

      ! Initialize coag coefficients
      call ssh_api_init_coag()

    end subroutine ssh_api_initialize

! =============================================================
!
! External code can force SSH-aerosol to:
!   call ssh_init_distributions0D
!
! =============================================================

    subroutine ssh_api_init_coag() bind(c, name='api_sshaerosol_init_coag_')

      use iso_c_binding
      use lDiscretization, only : ssh_init_coag

      implicit none

      ! Initialize distributions
      call ssh_init_coag()

    end subroutine ssh_api_init_coag

    
! =============================================================
!
! External code can force SSH-aerosol to:
!   call ssh_init_distributions
!
! =============================================================

    subroutine ssh_api_init_distributions() bind(c, name='api_sshaerosol_init_distributions_')

      use iso_c_binding
      use lDiscretization, only : ssh_init_distributions

      implicit none

      ! Initialize distributions
      call ssh_init_distributions()

      
    end subroutine ssh_api_init_distributions

! =============================================================
!
! External code can get converted aerosol mass concentration
!
! =============================================================

    subroutine ssh_api_get_aero_concentration_conv(array) &
         bind(c, name='api_sshaerosol_get_aero_conv')

      use iso_c_binding
      use aInitialization, only : N_size, N_aerosol, &
           N_aerosol_layers, concentration_mass, list_species

      implicit none

      real(kind=c_double), intent(out), dimension(N_size, N_aerosol) :: array
      integer :: b, s, jesp
      
      array = 0.d0

      do b = 1, N_size
         do s = 1, N_aerosol_layers
            jesp = List_species(s)
            array(b, jesp) = array(b, jesp) + concentration_mass(b, s)
         enddo
      enddo
    
    end subroutine ssh_api_get_aero_concentration_conv

! =============================================================
!
! External code can set aerosol converted mass concentration
!
! =============================================================

    subroutine ssh_api_set_aero_concentration_layer(array) &
         bind(c, name='api_sshaerosol_set_aero_conv')

      use iso_c_binding
      use aInitialization, only : N_size, N_aerosol, N_nonorganics, nlayer, &
                                  concentration_mass, index_species, vlayer, &
                                  i_hydrophilic

      implicit none

      real(kind=c_double), intent(in), dimension(N_size, N_aerosol) :: array
      integer :: s, jesp, j, lay
      
      do j = 1, N_size
         do s = 1, N_aerosol             ! s : index of species         
            if (s .LE. N_nonorganics) then
               concentration_mass(j, s) = array(j, s) 
            else
               do lay = 1, nlayer
                  jesp = index_species(s, lay)
                  concentration_mass(j, jesp) = array(j, s) * Vlayer(lay)
               enddo
               if(i_hydrophilic==1) then
                  jesp = index_species(s, nlayer + i_hydrophilic)
                  concentration_mass(j, jesp) = 0.d0 ! Initialise the aqueous phase to zero
               endif
            endif
         end do
      end do
    
    end subroutine ssh_api_set_aero_concentration_layer

    
! =============================================================
!
! External code can force SSH-aerosol to intialize again
!   internal variables. This subroutine must be called after
!   updating thermodynamic variables or concentrations.
!   TODO most of the code here is duplicated from the
!   subroutine init_distributions.
!   TODO add a dedicated subroutine to avoid duplicate
!
! =============================================================

    subroutine ssh_api_reinitialize() bind(c, name='api_sshaerosol_init_again_')

      use iso_c_binding
      use aInitialization
      use dPhysicalbalance
      use cThermodynamics

      implicit none

      integer :: s, jesp, i
      double precision :: tmp

      ! Initialize the concentrations of aerosol with gas-phase percursors
      total_aero_mass = 0.0
      total_mass = 0.0
      do s = 1, N_aerosol_layers
        jesp=List_species(s)
        do i=1,N_size
          total_aero_mass(jesp)=total_aero_mass(jesp)+concentration_mass(i,s)
        enddo
      enddo
      do s = 1, N_aerosol
        if (aerosol_species_interact(s) .gt. 0) then
          concentration_gas(s) = concentration_gas_all(aerosol_species_interact(s)) !µg/m3
        end if
        total_mass(s)=total_mass(s) + total_aero_mass(s) + concentration_gas(s)
      end do

    end subroutine ssh_api_reinitialize

! =============================================================
!
! External code can force SSH-aerosol to:
!   deallocate memory
!   close the log file
!
! =============================================================

    subroutine ssh_api_finalize() bind(c, name='api_sshaerosol_finalize_')

      use iso_c_binding
      use aInitialization, only : ssh_free_allocated_memory, with_coag, ssh_logger, ssh_close_logger, tag_chem, option_photolysis
      use bCoefficientRepartition, only : ssh_DeallocateCoefficientRepartition
      use mod_photolysis, only : ssh_deallocate_photolysis
      
      implicit none

      integer :: ierr

      ! Free memory
      call ssh_free_allocated_memory()
      if (with_coag.eq.1) call ssh_DeallocateCoefficientRepartition()

      if ((tag_chem .ne. 0).AND.(option_photolysis.eq.2)) then
         call ssh_deallocate_photolysis()
      endif
      
      ! Close log file if needed
      if (ssh_logger) then
        call ssh_close_logger()
      endif

    end subroutine ssh_api_finalize

! =============================================================
!
! External code can set the flag to declare SSH-aerosol is not running standalone
!
! input : true if running standalone (default), false otherwise
! =============================================================

    subroutine ssh_set_standalone(flag) bind(c, name='api_sshaerosol_set_standalone_')

      use iso_c_binding
      use aInitialization, only : ssh_standalone

      implicit none

      logical(kind=c_bool), intent(in) :: flag
      
      ssh_standalone = flag

    end subroutine ssh_set_standalone

! =============================================================
!
! External code can get the flag to check if SSH-aerosol is running standalone
!
! return value : true if running standalone (default), false otherwise
! =============================================================

    function standalone() bind(c, name='api_sshaerosol_get_standalone_')

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
! Important: This subroutine must be called before api_initialize
!
! input : true if logging to a file, false (default) otherwise
! =============================================================

    subroutine ssh_api_set_logger(cflag) bind(c, name='api_sshaerosol_set_logger_')

      use iso_c_binding
      use aInitialization, only : ssh_set_logger

      implicit none

      logical(kind=c_bool), intent(in) :: cflag
      logical :: flag

      flag = cflag
      call ssh_set_logger(flag)

    end subroutine ssh_api_set_logger

! =============================================================
!
! External code can get the flag to check if SSH-aerosol is logging informations
!
! return value : true if logging to a file, false (default) otherwise
! =============================================================

    function logger() bind(c, name='api_sshaerosol_get_logger_')

      use iso_c_binding
      use aInitialization, only : ssh_logger

      implicit none

      logical(kind=c_bool) :: logger
      
      logger = ssh_logger

    end function logger

! =============================================================
!
! External code can get the flag to check if SSH-aerosol is running with gas-phase chemistry.
!
! return value : true if running with gas-phase chemistry
! =============================================================

    function get_tag_chem() bind(c, name='api_sshaerosol_get_tag_chem_')

      use iso_c_binding
      use aInitialization, only : tag_chem

      implicit none

      logical(kind=c_bool) :: get_tag_chem

      if (tag_chem == 0) then
         get_tag_chem = .false.
      else
         get_tag_chem = .true.
      endif
      
    end function get_tag_chem

! =============================================================    

! =============================================================
!
! External code can set option with_init_num
!
! input : 1 or 0
! =============================================================

    subroutine ssh_api_set_with_init_num(val) &
      bind (c, name='api_sshaerosol_set_with_init_num_')

      use iso_c_binding
      use aInitialization, only : with_init_num

      implicit none

      integer(kind=c_int), intent(in) :: val

      with_init_num = val

    end subroutine ssh_api_set_with_init_num

    
    
! =============================================================
!
! External code can get the number of gas species
!
! return value : number of gas species
! =============================================================

    function api_get_ngas() bind(c, name='api_sshaerosol_get_ngas_')

      use iso_c_binding
      use aInitialization, only : N_gas

      implicit none

      integer(kind=c_int) :: api_get_ngas
      
      api_get_ngas = N_gas

    end function api_get_ngas

! =============================================================
!
! External code can get the number of aerosols species
!
! return value : number of aerosols species
! =============================================================

    function api_get_naero() bind(c, name='api_sshaerosol_get_naero_')

      use iso_c_binding
      use aInitialization, only : N_aerosol

      implicit none

      integer(kind=c_int) :: api_get_naero
      
      api_get_naero = N_aerosol

    end function api_get_naero

! =============================================================
!
! External code can get the number of aerosols layers
!
! return value : number of aerosols layers
! =============================================================

    function api_get_nlayer() bind(c, name='api_sshaerosol_get_nlayer_')

      use iso_c_binding
      use aInitialization, only : nlayer

      implicit none

      integer(kind=c_int) :: api_get_nlayer
      
      api_get_nlayer = nlayer

    end function api_get_nlayer    

    
! =============================================================
!
! External code can get the number of aerosols species in layers
!
! return value : number of aerosol species in layers
! =============================================================

    function api_get_n_aerosol_layers() bind(c, name='api_sshaerosol_get_n_aerosol_layers_')

      use iso_c_binding
      use aInitialization, only : N_aerosol_layers

      implicit none

      integer(kind=c_int) :: api_get_n_aerosol_layers

      api_get_n_aerosol_layers = N_aerosol_layers

    end function api_get_n_aerosol_layers

! =============================================================
!
! External code can get the number of aerosols
!
! return value : number of aerosols
! =============================================================

    function api_get_nsize() bind(c, name='api_sshaerosol_get_nsize_')

      use iso_c_binding
      use aInitialization, only : N_size

      implicit none

      integer(kind=c_int) :: api_get_nsize

      api_get_nsize = N_size

    end function api_get_nsize

! =============================================================
!
! External code can get the number of aerosols size bins
!
! return value : number of aerosols size bins
! =============================================================

    function api_get_nsizebin() bind(c, name='api_sshaerosol_get_nsizebin_')

      use iso_c_binding
      use aInitialization, only : N_sizebin

      implicit none

      integer(kind=c_int) :: api_get_nsizebin
      
      api_get_nsizebin = N_sizebin

    end function api_get_nsizebin


! =============================================================
!
! External code can get the number of inorganic aerosols species
!
! return value : number of inorganic aerosols species
! =============================================================

    function api_get_nesp_isorropia() bind(c, name='api_sshaerosol_get_nesp_isorropia_')

      use iso_c_binding
      use aInitialization, only : nesp_isorropia

      implicit none

      integer(kind=c_int) :: api_get_nesp_isorropia
      
      api_get_nesp_isorropia = nesp_isorropia

    end function api_get_nesp_isorropia

! =============================================================
!
! External code can get the number of organic aerosols species 
!
! return value : number of organic aerosols species
! =============================================================

    function api_get_nesp_aec() bind(c, name='api_sshaerosol_get_nesp_aec_')

      use iso_c_binding
      use aInitialization, only : nesp_aec

      implicit none

      integer(kind=c_int) :: api_get_nesp_aec
      
      api_get_nesp_aec = nesp_aec

    end function api_get_nesp_aec

    
! =============================================================
!
! External code can get the number of photolysis reactions
!
! return value : number of photolysis reactions
! =============================================================
   
    
    function api_get_nphotolysis() bind(c, name='api_sshaerosol_get_nphotolysis_')

      use iso_c_binding
      use aInitialization, only : N_photolysis

      implicit none

      integer(kind=c_int) :: api_get_nphotolysis
      
      api_get_nphotolysis = N_photolysis

    end function api_get_nphotolysis
    
! =============================================================
!
! External code can get the name of given aerosol species
!
! input : number of the aerosol species
! return value : name of the aerosol species
! =============================================================

    subroutine ssh_api_get_aero_name(aero_num, c_string) bind(c, name='api_sshaerosol_get_aero_name_')

      use iso_c_binding
      use aInitialization, only : N_aerosol, aerosol_species_name, logfile, ssh_standalone, ssh_logger

      implicit none

      ! Arguments

      integer(kind=c_int), intent(in) :: aero_num
      character(len = 1, kind = c_char), dimension(81), intent(out) :: c_string

      ! Local variables

      integer i, strlen

      ! Safety bounds check
      if (1.gt.aero_num .or. aero_num.gt.N_aerosol) then
        if (ssh_standalone) write(*,*) "Error: given aerosol number out of bounds."
        if (ssh_standalone) write(*,*) "Given : ", aero_num
        if (ssh_standalone) write(*,*) "Bounds: [1,",N_aerosol,"]"
        if (ssh_logger) write(logfile,*) "Error: given aerosol number out of bounds."
        if (ssh_logger) write(logfile,*) "Given : ", aero_num
        if (ssh_logger) write(logfile,*) "Bounds: [1,",N_aerosol,"]"
        stop
      endif

      ! Fill with aerosol species name
      strlen = min(80,len_trim(aerosol_species_name(aero_num)))
      do i = 1, strlen
        c_string(i) = aerosol_species_name(aero_num)(i:i)
      enddo

      ! C termination
      c_string(strlen + 1) = C_NULL_CHAR

    end subroutine ssh_api_get_aero_name

! =============================================================
!
! External code can set the SSH-aerosol time step
!
! input : time step in seconds
! return value : false if the time step is too small
! =============================================================

    function set_dt(val) bind(c, name='api_sshaerosol_set_dt_')

      use iso_c_binding
      use aInitialization, only : delta_t, dt, DTAEROMIN

      implicit none

      real(c_double), intent(in) :: val
      logical(kind=c_bool) :: set_dt

      if (val > DTAEROMIN) then
        set_dt = .true.
        delta_t = val
        dt = val ! this is probably useless
      else
        set_dt = .false.
        delta_t = DTAEROMIN
        dt = DTAEROMIN ! this is probably useless
     endif
     
    end function set_dt

! =============================================================
!
! External code can get the SSH-aerosol time step
!
! return value : time step in seconds
! =============================================================

    function get_dt() bind(c, name='api_sshaerosol_get_dt_')

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

    subroutine ssh_set_initial_t(val) bind(c, name='api_sshaerosol_set_initial_t_')

      use iso_c_binding
      use aInitialization, only : initial_time

      implicit none

      real(c_double), intent(in) :: val

      initial_time = val

    end subroutine ssh_set_initial_t

! =============================================================
!
! External code can get the SSH-aerosol initial time
!
! return value : initial time in seconds (GMT, computed from January 1st)
! =============================================================

    function get_initial_t() bind(c, name='api_sshaerosol_get_initial_t_')

      use iso_c_binding
      use aInitialization, only : initial_time

      implicit none

      real(c_double) :: get_initial_t

      get_initial_t = initial_time

    end function get_initial_t

! =============================================================
!
! External code can set the SSH-aerosol current time
!
! input : current time in seconds (GMT, computed from January 1st)
! =============================================================

    subroutine ssh_set_current_t(val) bind(c, name='api_sshaerosol_set_current_t_')

      use iso_c_binding
      use aInitialization, only : current_time

      implicit none

      real(c_double), intent(in) :: val

      current_time = val

    end subroutine ssh_set_current_t

! =============================================================
!
! External code can get the SSH-aerosol current time
!
! return value : current time in seconds (GMT, computed from January 1st)
! =============================================================

    function get_current_t() bind(c, name='api_sshaerosol_get_current_t_')

      use iso_c_binding
      use aInitialization, only : current_time

      implicit none

      real(c_double) :: get_current_t

      get_current_t = current_time

    end function get_current_t

! =============================================================
!
! External code can force to update the specific humidity
!
! =============================================================

    subroutine ssh_api_update_humidity() bind(c, name='api_sshaerosol_update_humidity_')

      use iso_c_binding
      use aInitialization, only : temperature, humidity, pressure, pressure_sat, relative_humidity

      implicit none

      call ssh_compute_psat_sh(Relative_Humidity, temperature, Pressure, pressure_sat, humidity)

    end subroutine ssh_api_update_humidity

! =============================================================
!
! External code can set specific humidity
!
! input : specific humidity in kg / kg
! =============================================================

    subroutine ssh_api_set_humidity(val) &
         bind(c, name='api_sshaerosol_set_humidity')

      use iso_c_binding
      use aInitialization, only : humidity

      implicit none

      real(kind=c_double), intent(in) :: val

      humidity = val

    end subroutine ssh_api_set_humidity
   

! =============================================================
!
! External code can set attenuation
!
! input : attenuation in adimensional value
! =============================================================

    subroutine ssh_api_set_attenuation(val) &
         bind(c, name='api_sshaerosol_set_attenuation')

      use iso_c_binding
      use aInitialization, only : attenuation

      implicit none

      real(kind=c_double), intent(in) :: val

      attenuation = val

    end subroutine ssh_api_set_attenuation

    
! =============================================================
!
! External code can get the specific humidity
!
! output : specific humidity in kg / kg
! =============================================================

    function api_get_humidity() bind(c, name='api_sshaerosol_get_humidity_')

      use iso_c_binding
      use aInitialization, only : humidity

      implicit none

      real(kind=c_double) :: api_get_humidity
      
      api_get_humidity = humidity

    end function api_get_humidity

! =============================================================
!
! External code can set the relative humidity
!
! input : relative humidity in ??? TODO
! =============================================================

    subroutine ssh_api_set_relhumidity(val) bind(c, name='api_sshaerosol_set_relhumidity_')

      use iso_c_binding
      use aInitialization, only : relative_humidity

      implicit none

      include 'CONST_A.INC' ! needed for threshold_RH_inf and threshold_RH_sup

      real(kind=c_double) :: val

      ! This is taken from the subroutine read_namelist
      relative_humidity = min(max(val, threshold_RH_inf), threshold_RH_sup)

    end subroutine ssh_api_set_relhumidity

! =============================================================
!
! External code can get the relative humidity
!
! output : relative humidity in ??? TODO
! =============================================================

    function api_get_relhumidity() bind(c, name='api_sshaerosol_get_relhumidity_')

      use iso_c_binding
      use aInitialization, only : relative_humidity

      implicit none

      real(kind=c_double) :: api_get_relhumidity

      api_get_relhumidity = relative_humidity

    end function api_get_relhumidity

! =============================================================
!
! External code can set the temperature
!
! input : temperature in K
! =============================================================

    subroutine ssh_api_set_temperature(val) bind(c, name='api_sshaerosol_set_temperature_')

      use iso_c_binding
      use aInitialization, only : temperature

      implicit none

      real(kind=c_double), intent(in) :: val

      temperature = val

    end subroutine ssh_api_set_temperature

! =============================================================
!
! External code can get the temperature
!
! output : temperature in K
! =============================================================

    function api_get_temperature() bind(c, name='api_sshaerosol_get_temperature_')

      use iso_c_binding
      use aInitialization, only : temperature

      implicit none

      real(kind=c_double) :: api_get_temperature

      api_get_temperature = temperature

    end function api_get_temperature

! =============================================================
!
! External code can set the pressure
!
! input : pressure in Pa
! =============================================================

    subroutine ssh_api_set_pressure(val) bind(c, name='api_sshaerosol_set_pressure_')

      use iso_c_binding
      use aInitialization, only : pressure

      implicit none

      real(kind=c_double), intent(in) :: val

      pressure = val

    end subroutine ssh_api_set_pressure

! =============================================================
!
! External code can get the pressure
!
! output : pressure in Pa
! =============================================================

    function api_get_pressure() bind(c, name='api_sshaerosol_get_pressure_')

      use iso_c_binding
      use aInitialization, only : pressure

      implicit none

      real(kind=c_double) :: api_get_pressure

      api_get_pressure = pressure

    end function api_get_pressure

! =============================================================
!
! External code can set the pH
!
! input : pH in ??? TODO
! =============================================================

    subroutine ssh_api_set_ph(val) bind(c, name='api_sshaerosol_set_ph_')

      use iso_c_binding
      use aInitialization, only : ph

      implicit none

      real(kind=c_double), intent(in) :: val

      ph = val

    end subroutine ssh_api_set_ph

! =============================================================
!
! External code can get the pH
!
! output : pH in ??? TODO
! =============================================================

    function api_get_ph() bind(c, name='api_sshaerosol_get_ph_')

      use iso_c_binding
      use aInitialization, only : ph

      implicit none

      real(kind=c_double) :: api_get_ph

      api_get_ph = ph

    end function api_get_ph


! =============================================================
!
! External code can get initial diameter
!
! output : array of diam_input in micro m
! =============================================================

    subroutine ssh_api_get_diam_input(diam) bind(c, name='api_sshaerosol_get_diam_input_')

      use iso_c_binding
      use aInitialization, only : diam_input, N_sizebin

      implicit none

      real(kind=c_double), dimension(N_sizebin + 1) :: diam

      diam = diam_input

    end subroutine ssh_api_get_diam_input
    
    
! =============================================================
!
! External code can set the gaseous concentrations
!
! input : array of concentrations in micrograms / m^3
! =============================================================

    subroutine ssh_api_set_gas_concentration(array) bind(c, name='api_sshaerosol_set_gas_')

      use iso_c_binding
      use aInitialization, only : N_gas, concentration_gas_all

      implicit none

      real(kind=c_double), intent(in), dimension(N_gas) :: array

      concentration_gas_all(:) = array(:)

    end subroutine ssh_api_set_gas_concentration

! =============================================================
!
! External code can get the gaseous concentrations
!
! output : array of concentrations in micrograms / m^3
! =============================================================

    subroutine ssh_api_get_gas_concentration(array) bind(c, name='api_sshaerosol_get_gas_')

      use iso_c_binding
      use aInitialization, only : N_gas, concentration_gas_all

      implicit none

      real(kind=c_double), intent(out), dimension(N_gas) :: array
      
      array(:) = concentration_gas_all(:)

    end subroutine ssh_api_get_gas_concentration

! =============================================================
!
! External code can set the initial aerosols concentrations
!
! input : 2D array of concentrations in micrograms / m^3
! =============================================================

    subroutine ssh_api_set_init_bin_mass(array) bind(c, name='api_sshaerosol_set_init_bin_mass_')

      use iso_c_binding
      use aInitialization, only : N_sizebin, N_aerosol, init_bin_mass

      implicit none

      real(kind=c_double), intent(in), dimension(N_sizebin, N_aerosol) :: array

      init_bin_mass(:,:) = array(:,:)

    end subroutine ssh_api_set_init_bin_mass

! =============================================================
!
! External code can set the initial aerosols number concentrations
!
! input : 2D array of concentrations in particles / m^3
! =============================================================

    subroutine ssh_api_set_init_bin_number(array) bind(c, name='api_sshaerosol_set_init_bin_number_')

      use iso_c_binding
      use aInitialization, only : N_sizebin, init_bin_number

      implicit none

      real(kind=c_double), intent(in), dimension(N_sizebin) :: array

      init_bin_number(:) = array(:)

    end subroutine ssh_api_set_init_bin_number
    
! =============================================================
!
! External code can set the aerosols concentrations
!
! input : 2D array of concentrations in micrograms / m^3
! =============================================================

    subroutine ssh_api_set_aero_concentration(array) bind(c, name='api_sshaerosol_set_aero_')

      use iso_c_binding
      use aInitialization, only : N_size, N_aerosol_layers, concentration_mass

      implicit none

      real(kind=c_double), intent(in), dimension(N_size, N_aerosol_layers) :: array

      concentration_mass(:,:) = array(:,:)

    end subroutine ssh_api_set_aero_concentration

! =============================================================
!
! External code can get the aerosols concentrations
!
! output : 2D array of concentrations in micrograms / m^3
! =============================================================

    subroutine ssh_api_get_aero_concentration(array) bind(c, name='api_sshaerosol_get_aero_')

      use iso_c_binding
      use aInitialization, only : N_size, N_aerosol_layers, concentration_mass

      implicit none

      real(kind=c_double), intent(out), dimension(N_size, N_aerosol_layers) :: array

      array(:,:) = concentration_mass(:,:)

    end subroutine ssh_api_get_aero_concentration

! =============================================================
!
! External code can set the aerosols numbers per size bin
!
! input : 1D array in particles / m^3
! =============================================================

    subroutine ssh_api_set_aero_number(array) bind(c, name='api_sshaerosol_set_aero_num_')

      use iso_c_binding
      use aInitialization, only : N_size, concentration_number

      implicit none

      real(kind=c_double), intent(in), dimension(N_size) :: array

      concentration_number(:) = array(:)

    end subroutine ssh_api_set_aero_number

! =============================================================
!
! External code can get the aerosols numbers per size bin
!
! output : 1D array in particles / m^3
! =============================================================

    subroutine ssh_api_get_aero_number(array) bind(c, name='api_sshaerosol_get_aero_num_')

      use iso_c_binding
      use aInitialization, only : N_size, concentration_number

      implicit none

      real(kind=c_double), intent(out), dimension(N_size) :: array

      array(:) = concentration_number(:)

    end subroutine ssh_api_get_aero_number

! =============================================================
!
! External code can call the emissions subroutine
!
! =============================================================

    subroutine ssh_api_call_ssh_emission() bind(c, name='api_sshaerosol_emission_')

      use iso_c_binding
      use aInitialization, only : tag_emis, delta_t
      use mEmissions, only : ssh_emission

      implicit none

      if (tag_emis .ne. 0) call ssh_emission(delta_t)

    end subroutine ssh_api_call_ssh_emission

! =============================================================
!
! External code can set the position in (lon, lat)
!
! input : longitude and latitude in degree
! =============================================================

    subroutine ssh_api_set_lonlat(val) bind(c, name='api_sshaerosol_set_lonlat_')

      use iso_c_binding
      use aInitialization, only : longitude, latitude

      implicit none

      real(kind=c_double), intent(in), dimension(2) :: val

      longitude = val(1)
      latitude = val(2)

    end subroutine ssh_api_set_lonlat


! =============================================================
!
! External code can get aerosol density for each size-composition section.
!
! output : array of density_aer_bin in micrograms / microm^3
! =============================================================

    subroutine ssh_api_get_density_aer_bin(array) bind(c, name='api_sshaerosol_get_density_aer_bin')

      use iso_c_binding
      use aInitialization, only : N_size, density_aer_bin

      implicit none

      real(kind=c_double), intent(out), dimension(N_size) :: array
      
      array(:) = density_aer_bin(:)

    end subroutine ssh_api_get_density_aer_bin


! =============================================================
!
! External code can get aerosol density for each size section.
!
! output : array of density_aer_bin in micrograms / microm^3
! =============================================================

    subroutine ssh_api_get_density_aer_size(array) bind(c, name='api_sshaerosol_get_density_aer_size')

      use iso_c_binding
      use aInitialization, only : N_sizebin, density_aer_size

      implicit none

      real(kind=c_double), intent(out), dimension(N_sizebin) :: array
      
      array(:) = density_aer_size(:)

    end subroutine ssh_api_get_density_aer_size


    
! =============================================================
!
! External code can get fixed density
!
! output : fixed density in microgram/microm3
! =============================================================

    function api_get_fixed_density() bind(c, name='api_sshaerosol_get_fixed_density_')

      use iso_c_binding
      use aInitialization, only : fixed_density

      implicit none

      real(kind=c_double) :: api_get_fixed_density

      api_get_fixed_density = fixed_density

    end function api_get_fixed_density

    
    
! =============================================================
!
! External code can call the chemistry scheme
!
! =============================================================

    subroutine ssh_api_call_ssh_gaschemistry() bind(c, name='api_sshaerosol_gaschemistry_')

      use iso_c_binding
      use aInitialization

      implicit none

      double precision :: current_time_api

      current_time_api = current_time

      if (tag_chem .ne. 0) then
        call ssh_chem(n_gas, n_reaction, n_photolysis, photolysis_reaction_index,&
          ns_source, source_index, conversionfactor, conversionfactorjacobian,&
          0, lwc_cloud_threshold, molecular_weight, &
          current_time_api, attenuation, &
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
      endif
      
    end subroutine ssh_api_call_ssh_gaschemistry

! =============================================================
!
! External code can call the aerosol scheme
!
! =============================================================

    subroutine ssh_api_call_ssh_aerochemistry() bind(c, name='api_sshaerosol_aerodyn_')

      use iso_c_binding
      use aInitialization
      use jAdaptstep, only : ssh_AERODYN

      implicit none

      integer :: s, j, jesp
      double precision :: current_time_api

      current_time_api = current_time

      ! re-calculate total_mass(N_aerosol) because mass change due to gas-phase chemistry
      total_aero_mass = 0.d0
      total_mass = 0.d0
      do s = 1, N_aerosol_layers
        jesp = List_species(s)
        do j = 1, N_size
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
      CALL SSH_AERODYN(current_time_api,delta_t)

      ! update mass conc. of aerosol precursors
      ! concentration_gas(n_aerosol) -> concentration_gas_all(precursor_index)
      do s = 1, N_aerosol
        if (aerosol_species_interact(s) .gt. 0) then
          concentration_gas_all(aerosol_species_interact(s)) = concentration_gas(s)
        end if
      end do
      
    end subroutine ssh_api_call_ssh_aerochemistry

! =============================================================
!
! External code can call the save_report subroutine
!         
! =============================================================

    subroutine ssh_api_call_ssh_report() bind(c, name='api_sshaerosol_report_')

      use iso_c_binding
      use Resultoutput

      implicit none

      call ssh_save_report()

    end subroutine ssh_api_call_ssh_report

! =============================================================
!
! External code can call the init_output_conc subroutine
!         
! =============================================================
          
    subroutine ssh_api_call_ssh_initoutput() bind(c, name='api_sshaerosol_initoutput_')
      
      use iso_c_binding
      use Resultoutput 

      implicit none

      call ssh_init_output_conc()

    end subroutine ssh_api_call_ssh_initoutput

! =============================================================
!
! External code can call the save_concentration subroutine
!         
! =============================================================
          
    subroutine ssh_api_call_ssh_output() bind(c, name='api_sshaerosol_output_')

      use iso_c_binding
      use Resultoutput

      implicit none

      call ssh_save_concentration()

    end subroutine ssh_api_call_ssh_output

! =============================================================
!
! External code can initialize the photolysis
!
! =============================================================

    subroutine ssh_api_call_ssh_initphoto() bind(c, name='api_sshaerosol_initphoto_')

      use iso_c_binding
      use aInitialization
      use mod_photolysis

      implicit none

      if ((tag_chem .ne. 0).AND.(option_photolysis.eq.2)) then
        ! Allocate arrays for photolysis
        call ssh_allocate_photolysis()
        ! Read photolysis rate for current day
        current_time = initial_time
        call ssh_init_photolysis()
        call ssh_interpol_photolysis()
      endif

    end subroutine ssh_api_call_ssh_initphoto

! =============================================================
!    
! External code can update the photolysis
!
! =============================================================

    subroutine ssh_api_call_ssh_updatephoto() bind(c, name='api_sshaerosol_updatephoto_')

      use iso_c_binding
      use aInitialization
      use mod_photolysis

      implicit none

      ! Read the photolysis rates.
      if ((tag_chem .ne. 0).AND.(option_photolysis.eq.2)) then
        call ssh_interpol_photolysis()
      endif

    end subroutine ssh_api_call_ssh_updatephoto

! =============================================================
!
! External code can set the photolysis rate
!
! input : array of photolysis_rate  
! =============================================================

    subroutine ssh_api_set_photolysis_rate(array) bind(c, &
         name='api_sshaerosol_set_photolysis_rate_')

      use iso_c_binding
      use aInitialization, only : n_photolysis, photolysis_rate

      implicit none

      real(kind=c_double), intent(in), dimension(n_photolysis) :: array

      photolysis_rate(:) = array(:)

    end subroutine ssh_api_set_photolysis_rate

    

! =============================================================
!
! External code can call ssh_compute_all_density subroutine
!         
! =============================================================

    subroutine ssh_api_call_ssh_compute_all_density() bind(c, name='api_sshaerosol_compute_all_density_')

      use iso_c_binding
      use dPhysicalbalance

      implicit none

      call ssh_compute_all_density()

    end subroutine ssh_api_call_ssh_compute_all_density

! =============================================================
!
! External code can set aerosol density
!
! input : µg/µm3
! =============================================================

    subroutine ssh_api_set_fixed_density(val) bind(c, name='api_sshaerosol_set_fixed_density_')

      use iso_c_binding
      use aInitialization, only : fixed_density

      implicit none

      real(kind=c_double), intent(in) :: val

      fixed_density = val

    end subroutine ssh_api_set_fixed_density

! =============================================================
!
! External code can get concentration_index
!
! output : array of concentration_index
! =============================================================

    subroutine ssh_api_get_concentration_index(array) bind(c, name='api_sshaerosol_get_concentration_index_')

      use iso_c_binding
      use aInitialization, only : N_size, concentration_index

      implicit none

      integer(kind=c_int), intent(out), dimension(N_size, 2) :: array
      
      array(:,:) = concentration_index(:,:)

    end subroutine ssh_api_get_concentration_index


! =============================================================
!
! External code can get index for aerosol species
!
! output : array of list_speices
! =============================================================

    subroutine ssh_api_get_list_species(array) bind(c, name='api_sshaerosol_get_list_species_')

      use iso_c_binding
      use aInitialization, only : N_aerosol_layers, list_species

      implicit none

      integer(kind=c_int), intent(out), dimension(N_aerosol_layers) :: array
      
      array(:) = list_species(:)

    end subroutine ssh_api_get_list_species

! =============================================================
!
! External code can get the value for section_pass
!
! return value : section_pass
! =============================================================

    function api_get_section_pass() bind(c, name='api_sshaerosol_get_section_pass_')

      use iso_c_binding
      use aInitialization, only : section_pass

      implicit none

      integer(kind=c_int) :: api_get_section_pass
      
      api_get_section_pass = section_pass

    end function api_get_section_pass    

! =============================================================
!
! External code can call ssh_compute_number
!         
! =============================================================

    subroutine ssh_api_call_ssh_compute_number() bind(c, name='api_sshaerosol_compute_number_')

      use iso_c_binding
      use dPhysicalbalance

      implicit none

      call ssh_compute_number
      
    end subroutine ssh_api_call_ssh_compute_number


    
end module SSHaerosolAPI
