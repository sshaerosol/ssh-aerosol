!
! Compile ssh-aerosol using ./compile --sharedlib=yes
!
! How to compile this example :
! gfortran -g -o f_example fapi.f90 -ldl ../src/libssh-aerosol.so -rdynamic
!
! How to run :
! LD_LIBRARY_PATH=../src/:$LD_LIBRARY_PATH ./f_example
!
! See https://rosettacode.org/wiki/Call_a_function_in_a_shared_library#Fortran
!

module so_module

  use iso_c_binding
  
  implicit none
  
  private
  
  public :: loadlib, get_func, closelib

  interface
    function dlopen(filename,mode) bind(c,name="dlopen")
      ! void *dlopen(const char *filename, int mode);
      use iso_c_binding
      implicit none
      type(c_ptr) :: dlopen
      character(c_char), intent(in) :: filename(*)
      integer(c_int), value :: mode
    end function

    function dlsym(handle,name) bind(c,name="dlsym")
      ! void *dlsym(void *handle, const char *name);
      use iso_c_binding
      implicit none
      type(c_funptr) :: dlsym
      type(c_ptr), value :: handle
      character(c_char), intent(in) :: name(*)
    end function

    function dlclose(handle) bind(c,name="dlclose")
      ! int dlclose(void *handle);
      use iso_c_binding
      implicit none
      integer(c_int) :: dlclose
      type(c_ptr), value :: handle
    end function
  end interface

  contains

  ! Try to load the given shared library
  function loadlib(libname)

    use iso_c_binding
    
    implicit none

    ! Arguments / output
    character(len=*), intent(in) :: libname
    type(c_ptr) :: loadlib

    ! Local variables
    integer(kind=c_int), parameter :: rtld_lazy=1
    
    loadlib = dlopen( trim(libname)//c_null_char, rtld_lazy )
    if( .not. c_associated(loadlib) ) then
      write(*,*) "Could not load the shared library "//trim(libname)
      stop
    else
      ! The line below can be safely removed
      write(*,*) "Shared library loaded"
    endif

  end function loadlib
  
  ! Look for given function inside given library
  function get_func(lib, funcname)

    use iso_c_binding

    implicit none

    ! Arguments / output
    type(c_ptr), intent(in) :: lib
    character(len=*), intent(in) :: funcname
    type(c_funptr) :: get_func

    get_func = dlsym( lib, trim(funcname)//c_null_char )
    if( .not. c_associated(get_func) ) then
      get_func = dlsym( lib, trim(funcname)//"_"//c_null_char )
    endif

    if( .not. c_associated(get_func) ) then
      write(*,*) "Could not find subroutine "//trim(funcname)
      stop
    else
      ! The line below can be safely removed
      write(*,*) "Subroutine "//trim(funcname)//" found"
    endif

  end function get_func
  
  ! Close the given shared library
  subroutine closelib(lib)

    use iso_c_binding

    implicit none

    ! Arguments
    type(c_ptr), intent(in) :: lib

    ! Local variables
    integer(kind=c_int) :: retval

    retval = dlclose(lib)
    if (retval.eq.0) then
      ! The line below can be safely removed
      write(*,*) "Shared library is closed"
    else
      write(*,*) "Error when closing the shared library"
      stop
    endif

  end subroutine closelib

end module so_module

module modsshaero

  use iso_c_binding
  use so_module

  implicit none

  ! C pointer for the shared library
  type(c_ptr) :: libssh
  ! C pointer for the functions
  type(c_funptr) :: procaddr

  ! This is uded to map C functions in Fortran
  !procedure(send_bool), pointer :: send_b
  !procedure(recv_bool), pointer :: recv_b
  !procedure(send_int), pointer :: send_i
  !procedure(recv_int), pointer :: recv_i
  !procedure(send_dble), pointer :: send_db
  !procedure(recv_dble), pointer :: recv_db
  !procedure(send_recv_array), pointer :: snd_rcv_db
  !procedure(send_recv_char), pointer :: snd_rcv_chr
  !procedure(just_call), pointer :: jst_cll

  ! Default is public, a few quantities are private
  private :: procaddr

  abstract interface
    ! Send / receive logical
    subroutine send_bool(flag)
      implicit none
      logical, intent(in) :: flag
    end subroutine send_bool
    function recv_bool()
      implicit none
      logical :: recv_bool
    end function recv_bool
    ! Send / receive int
    subroutine send_int(val)
      implicit none
      integer, intent(in) :: val
    end subroutine send_int
    function recv_int()
      implicit none
      integer :: recv_int
    end function recv_int
    ! Send / receive double precision
    subroutine send_dble(val)
      implicit none
      double precision, intent(in) :: val
    end subroutine send_dble
    function recv_dble()
      implicit none
      double precision :: recv_dble
    end function recv_dble
    ! Exchange double precision array
    subroutine send_recv_array(val)
      use iso_c_binding
      implicit none
      double precision, intent(inout) :: val
    end subroutine send_recv_array
    ! Exchange character array
    subroutine send_recv_char(val)
      implicit none
      character(len=*), intent(inout) :: val
    end subroutine send_recv_char
    ! Just call a subroutine
    subroutine just_call
      implicit none
    end subroutine just_call
  end interface

  contains

  !
  ! Read the given namelist
  ! Initialize SSH
  ! Return ngas, naerolayer, nsize
  !
  subroutine modsshaero_initialize(nml_file, ngas, naerolayer, nsize)

    implicit none

    ! Arguments
    character(len=*) :: nml_file
    integer, intent(out) :: ngas, naerolayer, nsize

    ! Local variables
    integer, parameter :: size_namelist_file = 400
    character(len=size_namelist_file) :: cfile
    procedure(send_bool), pointer :: send_b
    procedure(recv_int), pointer :: recv_i
    procedure(send_recv_char), pointer :: snd_rcv_chr
    procedure(just_call), pointer :: jst_cll

    ! The file should end with c_null_char
    cfile = trim(nml_file)//c_null_char

    ! Load the shared library
    libssh = loadlib("libssh-aerosol.so")

    ! Set the standalone flag to false
    procaddr = get_func(libssh, "api_sshaerosol_set_standalone_")
    call c_f_procpointer(procaddr, send_b)
    call send_b(.false.)

    ! Set the logger to true
    procaddr = get_func(libssh, "api_sshaerosol_set_logger_")
    call c_f_procpointer(procaddr, send_b)
    call send_b(.true.)

    ! Initialize with the namelist
    procaddr = get_func(libssh, "api_sshaerosol_initialize_")
    call c_f_procpointer(procaddr, snd_rcv_chr)
    call snd_rcv_chr(cfile)

    ! call Init_distribution
    procaddr = get_func(libssh, "api_sshaerosol_init_distributions_")
    call c_f_procpointer(procaddr, jst_cll)
    call jst_cll()

    ! Get number of gas species
    procaddr = get_func(libssh, "api_sshaerosol_get_ngas_")
    call c_f_procpointer(procaddr, recv_i)
    ngas = recv_i()

    ! Get number of aerosol layers
    procaddr = get_func(libssh, "api_sshaerosol_get_n_aerosol_layers_")
    call c_f_procpointer(procaddr, recv_i)
    naerolayer = recv_i()

    ! Get number of size bins
    procaddr = get_func(libssh, "api_sshaerosol_get_nsize_")
    call c_f_procpointer(procaddr, recv_i)
    nsize = recv_i()

  end subroutine modsshaero_initialize

  !
  ! Initialize photolysis
  !
  subroutine modsshaero_initphoto()

    implicit none

    ! Local variables
    procedure(just_call), pointer :: jst_cll

    procaddr = get_func(libssh, "api_sshaerosol_initphoto_")
    call c_f_procpointer(procaddr, jst_cll)
    call jst_cll()

  end subroutine modsshaero_initphoto

  !
  ! Get current time from SSH
  !
  subroutine modsshaero_get_current_t(t)

    implicit none

    ! Argument
    double precision :: t

    ! Local variable
    procedure(recv_dble), pointer :: recv_db

    procaddr = get_func(libssh, "api_sshaerosol_get_current_t_")
    call c_f_procpointer(procaddr, recv_db)
    t = recv_db()

  end subroutine modsshaero_get_current_t

  !
  ! Set current time in SSH
  !
  subroutine modsshaero_set_current_t(t)

    implicit none

    ! Argument
    double precision :: t

    ! Local variable
    procedure(send_dble), pointer :: send_db

    procaddr = get_func(libssh, "api_sshaerosol_set_current_t_")
    call c_f_procpointer(procaddr, send_db)
    call send_db(t)

  end subroutine modsshaero_set_current_t

  !
  ! Set time step in SSH
  !
  subroutine modsshaero_set_dt(dt)

    implicit none

    ! Argument
    double precision, intent(in) :: dt

    ! Local variable
    procedure(send_dble), pointer :: send_db

    ! Set time step
    procaddr = get_func(libssh, "api_sshaerosol_set_dt_")
    call c_f_procpointer(procaddr, send_db)
    call send_db(dt)

  end subroutine modsshaero_set_dt

  !
  ! Get gas concentration from SSH
  !
  subroutine modsshaero_get_gas(cgas)

    implicit none

    ! Argument
    double precision, dimension(*) :: cgas

    ! Local variable
    procedure(send_recv_array), pointer :: snd_rcv_db

    procaddr = get_func(libssh, "api_sshaerosol_get_gas_")
    call c_f_procpointer(procaddr, snd_rcv_db)
    call snd_rcv_db(cgas(1))

  end subroutine modsshaero_get_gas

  !
  ! Set gas concentration in SSH
  !
  subroutine modsshaero_set_gas(cgas)

    implicit none

    ! Argument
    double precision, dimension(*) :: cgas

    ! Local variable
    procedure(send_recv_array), pointer :: snd_rcv_db

    procaddr = get_func(libssh, "api_sshaerosol_set_gas_")
    call c_f_procpointer(procaddr, snd_rcv_db)
    call snd_rcv_db(cgas(1))

  end subroutine modsshaero_set_gas

  !
  ! Get aerosol concentration from SSH
  !
  subroutine modsshaero_get_aero(conc_mass_aero)

    implicit none

    ! Argument
    double precision, dimension(*) :: conc_mass_aero

    ! Local variable
    procedure(send_recv_array), pointer :: snd_rcv_db

    procaddr = get_func(libssh, "api_sshaerosol_get_aero_")
    call c_f_procpointer(procaddr, snd_rcv_db)
    call snd_rcv_db(conc_mass_aero(1))

  end subroutine modsshaero_get_aero

  !
  ! Set aerosol concentration in SSH
  !
  subroutine modsshaero_set_aero(conc_mass_aero)

    implicit none

    ! Argument
    double precision, dimension(*) :: conc_mass_aero

    ! Local variable
    procedure(send_recv_array), pointer :: snd_rcv_db

    procaddr = get_func(libssh, "api_sshaerosol_set_aero_")
    call c_f_procpointer(procaddr, snd_rcv_db)
    call snd_rcv_db(conc_mass_aero(1))

  end subroutine modsshaero_set_aero

  !
  ! Get aerosol number from SSH
  !
  subroutine modsshaero_get_naero(conc_numb_aero)

    implicit none

    ! Argument
    double precision, dimension(*) :: conc_numb_aero

    ! Local variable
    procedure(send_recv_array), pointer :: snd_rcv_db

    procaddr = get_func(libssh, "api_sshaerosol_get_aero_num_")
    call c_f_procpointer(procaddr, snd_rcv_db)
    call snd_rcv_db(conc_numb_aero(1))

  end subroutine modsshaero_get_naero

  !
  ! Set aerosol number in SSH
  !
  subroutine modsshaero_set_naero(conc_numb_aero)

    implicit none

    ! Argument
    double precision, dimension(*) :: conc_numb_aero

    ! Local variable
    procedure(send_recv_array), pointer :: snd_rcv_db

    procaddr = get_func(libssh, "api_sshaerosol_set_aero_num_")
    call c_f_procpointer(procaddr, snd_rcv_db)
    call snd_rcv_db(conc_numb_aero(1))

  end subroutine modsshaero_set_naero

  !
  ! Get temperature, pressure and relative humidity from SSH
  !
  subroutine modsshaero_get_t_p_rh(temp, pres, relh)

    implicit none

    ! Arguments
    double precision :: temp, pres, relh

    ! Local variable
    procedure(recv_dble), pointer :: recv_db

    ! Set temperature
    procaddr = get_func(libssh, "api_sshaerosol_get_temperature_")
    call c_f_procpointer(procaddr, recv_db)
    temp = recv_db()

    ! Set pressure
    procaddr = get_func(libssh, "api_sshaerosol_get_pressure_")
    call c_f_procpointer(procaddr, recv_db)
    pres = recv_db()

    ! Set relative humidity
    procaddr = get_func(libssh, "api_sshaerosol_get_relhumidity_")
    call c_f_procpointer(procaddr, recv_db)
    relh = recv_db()

  end subroutine modsshaero_get_t_p_rh

  !
  ! Set temperature, pressure and relative humidity in SSH
  !
  subroutine modsshaero_set_t_p_rh(temp, pres, relh)

    implicit none

    ! Arguments
    double precision, intent(in) :: temp, pres, relh

    ! Local variable
    procedure(send_dble), pointer :: send_db
    procedure(just_call), pointer :: jst_cll

    ! Set temperature
    procaddr = get_func(libssh, "api_sshaerosol_set_temperature_")
    call c_f_procpointer(procaddr, send_db)
    call send_db(temp)

    ! Set pressure
    procaddr = get_func(libssh, "api_sshaerosol_set_pressure_")
    call c_f_procpointer(procaddr, send_db)
    call send_db(pres)

    ! Set relative humidity
    procaddr = get_func(libssh, "api_sshaerosol_set_relhumidity_")
    call c_f_procpointer(procaddr, send_db)
    call send_db(relh)

    ! Update saturation pressure
    procaddr = get_func(libssh, "api_sshaerosol_update_humidity_")
    call c_f_procpointer(procaddr, jst_cll)
    call jst_cll()


  end subroutine modsshaero_set_t_p_rh

  !
  ! Impose given time step, temperature, pressure and relative humidity
  !
  subroutine modsshaero_update_dt_t_p_rh(dt, temp, pres, relh)

    implicit none

    ! Arguments
    double precision, intent(in) :: dt, temp, pres, relh

    ! Set time step
    call modsshaero_set_dt(dt)

    ! Set temperature, pressure and relative humidity
    call modsshaero_set_t_p_rh(temp, pres, relh)

  end subroutine modsshaero_update_dt_t_p_rh

  !
  ! Impose given time step, temperature, pressure, relative humidity and concentrations
  ! Perform one time step of aerosol dynamic
  !
  subroutine modsshaero_aerosol_only(duration, temp, pres, relh, gas, aero, numb)

    implicit none

    ! Arguments
    double precision, intent(in) :: duration, temp, pres, relh
    double precision, dimension(*) :: gas, aero, numb

    ! Local variable
    procedure(just_call), pointer :: jst_cll

    ! Update time step and thermodynamic parameters
    call modsshaero_update_dt_t_p_rh(duration, temp, pres, relh)

    ! Set gas concentration
    call modsshaero_set_gas(gas)

    ! Set aerosol concentration
    call modsshaero_set_aero(aero)

    ! Set aerosol number
    call modsshaero_set_naero(numb)

    ! Add emission
    call modsshaero_emission()

    ! Perform one time step of aerosol dynamic
    procaddr = get_func(libssh, "api_sshaerosol_aerodyn_")
    call c_f_procpointer(procaddr, jst_cll)
    call jst_cll()

    ! Get aerosol concentration
    call modsshaero_get_aero(aero)

    ! Get aerosol number
    call modsshaero_get_naero(numb)

  end subroutine modsshaero_aerosol_only

  !
  ! Impose given time step, temperature, pressure, relative humidity and concentrations
  ! Perform one time step of gas chemistry and aerosol dynamic
  !
  subroutine modsshaero_aerochem(duration, temp, pres, relh, gas, aero, numb)

    implicit none

    ! Arguments
    double precision, intent(in) :: duration, temp, pres, relh
    double precision, dimension(*) :: gas, aero, numb

    ! Local variable
    procedure(just_call), pointer :: jst_cll

    ! Update time step and thermodynamic parameters
    call modsshaero_update_dt_t_p_rh(duration, temp, pres, relh)

    ! Set gas concentration
    call modsshaero_set_gas(gas)

    ! Set aerosol concentration
    call modsshaero_set_aero(aero)

    ! Set aerosol number
    call modsshaero_set_naero(numb)

    ! Add emission
    call modsshaero_emission()

    ! Perform one time step of gas chemistry
    procaddr = get_func(libssh, "api_sshaerosol_gaschemistry_")
    call c_f_procpointer(procaddr, jst_cll)
    call jst_cll()
    
    ! Perform one time step of aerosol dynamic
    procaddr = get_func(libssh, "api_sshaerosol_aerodyn_")
    call c_f_procpointer(procaddr, jst_cll)
    call jst_cll()

    ! Get gas concentration
    call modsshaero_get_gas(gas)

    ! Get aerosol concentration
    call modsshaero_get_aero(aero)

    ! Get aerosol number
    call modsshaero_get_naero(numb)

  end subroutine modsshaero_aerochem

  !
  ! Initialize the built-in processing
  !
  subroutine modsshaero_preprocessing()

    implicit none

    ! Local variable
    procedure(just_call), pointer :: jst_cll

    procaddr = get_func(libssh, "api_sshaerosol_initoutput_")
    call c_f_procpointer(procaddr, jst_cll)
    call jst_cll()

    procaddr = get_func(libssh, "api_sshaerosol_report_")
    call c_f_procpointer(procaddr, jst_cll)
    call jst_cll()

    call modsshaero_processing()

  end subroutine modsshaero_preprocessing

  !
  ! Perform built-in processing
  !
  subroutine modsshaero_processing()

    implicit none

    ! Local variable
    procedure(just_call), pointer :: jst_cll

    procaddr = get_func(libssh, "api_sshaerosol_output_")
    call c_f_procpointer(procaddr, jst_cll)
    call jst_cll()

  end subroutine modsshaero_processing

  !
  ! Update photolysis
  !
  subroutine modsshaero_updatephoto()

    implicit none

    ! Local variable
    procedure(just_call), pointer :: jst_cll

    procaddr = get_func(libssh, "api_sshaerosol_updatephoto_")
    call c_f_procpointer(procaddr, jst_cll)
    call jst_cll()

  end subroutine modsshaero_updatephoto

  !
  ! Perform emission
  !
  subroutine modsshaero_emission()

    implicit none

    ! Local variable
    procedure(just_call), pointer :: jst_cll

    ! Apply emissions
    procaddr = get_func(libssh, "api_sshaerosol_emission_")
    call c_f_procpointer(procaddr, jst_cll)
    call jst_cll()

  end subroutine modsshaero_emission

  !
  ! Finalize SSH-Aerosol and close the shared library
  !
  subroutine modsshaero_finalize()

    implicit none

    ! Local variable
    procedure(just_call), pointer :: jst_cll

    ! Finalize SSH-Aerosol
    procaddr = get_func(libssh, "api_sshaerosol_finalize_")
    call c_f_procpointer(procaddr, jst_cll)
    call jst_cll()

    ! Close the shared library
    call closelib(libssh)

  end subroutine modsshaero_finalize

end module modsshaero

program main

  use iso_c_binding
  use modsshaero

  implicit none

  ! Local variables
  character(50) :: nmlst = "namelist.ssh"
  integer :: it, nt, ngas, naerolayer, nsize
  double precision :: t, dt, temp, pres, relh
  double precision, allocatable :: cgas(:), conc_mass_aero(:,:), conc_numb_aero(:)

  ! Initialize the shared library
  ! This should be performed at the beginning of the simulation
  ! In case of parallelism, each MPI process should call
  call modsshaero_initialize(nmlst, ngas, naerolayer, nsize)

  ! Built-in post-processing
  ! In case of parallelism, only one MPI process should call
  call modsshaero_preprocessing()

  ! Initialize photolysis
  call modsshaero_initphoto()

  ! Get current time from SSH
  call modsshaero_get_current_t(t)

  ! Get thermodynamic parameters or set them
  if (.true.) then
    call modsshaero_get_t_p_rh(temp, pres, relh)
  else
    temp = 310.d0
    pres = 1013.d2
    relh = 0.123456789d0
    call modsshaero_set_t_p_rh(temp, pres, relh)
  endif

  ! Print some information and allocate memory
  write(*,*) "Number of gas species, aerosol layers and size bins : ", ngas, naerolayer, nsize
  allocate(cgas(ngas)); cgas = 0.d0
  allocate(conc_mass_aero(nsize, naerolayer)); conc_mass_aero = 0.d0
  allocate(conc_numb_aero(nsize)); conc_numb_aero = 0.d0

  ! Example: get concentrations and numbers from SSH
  call modsshaero_get_gas(cgas)
  call modsshaero_get_aero(conc_mass_aero)
  call modsshaero_get_naero(conc_numb_aero)

  ! Example: set concentrations and numbers in SSH
  call modsshaero_set_gas(cgas)
  call modsshaero_set_aero(conc_mass_aero)
  call modsshaero_set_naero(conc_numb_aero)

  ! Example: perform 60 time steps of gas chemistry + aerosol dynamic
  nt = 60
  ! Time step set to 60 seconds
  dt = 60.d0
  do it = 1, nt

    ! Set current time in SSH
    if (it > 1) then
      t = t + dt    
      call modsshaero_set_current_t(t)
    endif

    ! Update photolysis
    call modsshaero_updatephoto()

    ! The subroutine will update in SSH the given
    !   time step
    !   temperature
    !   pressure
    !   relative humidity
    !   gas concentration
    !   aerosol concentration
    !   aerosol number
    ! After the call
    !   gas concentration, aerosol concentration and number contain updated values
    call modsshaero_aerochem(dt, temp, pres, relh, cgas, conc_mass_aero, conc_numb_aero)

    ! Built-in post-processing
    ! In case of parallelism, only one MPI process should call
    call modsshaero_processing()

  enddo

  ! This should be performed at the end of the simulation
  ! Finalize SSH-Aerosol and close the shared library
  ! In case of parallelism, each MPI process should call
  call modsshaero_finalize()

  ! Free memory
  deallocate(cgas, conc_mass_aero, conc_numb_aero)

end program main
