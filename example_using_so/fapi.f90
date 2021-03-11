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
      write(*,*) "Shared library is closed"
    else
      write(*,*) "Error when closing the shared library"
      stop
    endif

  end subroutine closelib

end module so_module

program main

  use iso_c_binding
  use so_module

  implicit none

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

  ! This is needed to map C functions in Fortran
  procedure(send_bool), pointer :: send_b
  procedure(recv_bool), pointer :: recv_b
  procedure(send_int), pointer :: send_i
  procedure(recv_int), pointer :: recv_i
  procedure(send_dble), pointer :: send_db
  procedure(recv_dble), pointer :: recv_db
  procedure(send_recv_array), pointer :: snd_rcv_db
  procedure(send_recv_char), pointer :: snd_rcv_chr
  procedure(just_call), pointer :: jst_cll

  ! C pointer for the shared library
  type(c_ptr) :: libssh
  
  ! C pointer for the functions
  type(c_funptr) :: procaddr

  ! Local variables
  character(50) :: nmlst = trim("namelist.ssh")//c_null_char
  logical :: bool
  integer :: ngas, nlayer, nsize
  integer :: it
  double precision :: time, dt
  double precision, dimension(:), allocatable :: cgas

  ! Load the shared library
  libssh = loadlib("../src/libssh-aerosol.so")

  ! Set the standalone flag to false
  procaddr = get_func(libssh, "api_sshaerosol_set_standalone_")
  call c_f_procpointer(procaddr, send_b)
  call send_b(.false.)

  ! Read the standalone flag, check if the value is correct
  procaddr = get_func(libssh, "api_sshaerosol_get_standalone_")
  call c_f_procpointer(procaddr, recv_b)
  bool = recv_b()
  if (bool) then
    write(*,*) "API NOT OK"
  else
    write(*,*) "API OK"  
  endif

  ! Set the logger to true
  procaddr = get_func(libssh, "api_sshaerosol_set_logger_")
  call c_f_procpointer(procaddr, send_b)
  call send_b(.true.)

  ! Initialize with the namelist
  procaddr = get_func(libssh, "api_sshaerosol_initialize_")
  call c_f_procpointer(procaddr, snd_rcv_chr)
  call snd_rcv_chr(nmlst)

  ! Get number of gas species
  procaddr = get_func(libssh, "api_sshaerosol_get_ngas_")
  call c_f_procpointer(procaddr, recv_i)
  ngas = recv_i()

  ! Get number of aerosol layers
  procaddr = get_func(libssh, "api_sshaerosol_get_n_aerosol_layers_")
  call c_f_procpointer(procaddr, recv_i)
  nlayer = recv_i()

  ! Get number of size bins
  procaddr = get_func(libssh, "api_sshaerosol_get_nsize_")
  call c_f_procpointer(procaddr, recv_i)
  nsize = recv_i()

  ! Write some information and allocate memory
  write(*,*) "Number of gas species, aerosol layers and size bins : ", ngas, nlayer, nsize
  allocate(cgas(ngas))

  ! Read gas concentration
  write(*,*) "Gas concentration : ", cgas
  procaddr = get_func(libssh, "api_sshaerosol_get_gas_")
  call c_f_procpointer(procaddr, snd_rcv_db)
  call snd_rcv_db(cgas(1))
  write(*,*) "Gas concentration : ", cgas

  ! Read initial time and time step
  procaddr = get_func(libssh, "api_sshaerosol_get_initial_t_")
  call c_f_procpointer(procaddr, recv_db)
  time = recv_db()
  procaddr = get_func(libssh, "api_sshaerosol_get_dt_")
  call c_f_procpointer(procaddr, recv_db)
  dt = recv_db()

  ! Init_distributions
  procaddr = get_func(libssh, "api_sshaerosol_init_distributions_")
  call c_f_procpointer(procaddr, jst_cll)
  call jst_cll()

  ! Built-in reporting
  procaddr = get_func(libssh, "api_sshaerosol_initoutput_")
  call c_f_procpointer(procaddr, jst_cll)
  call jst_cll()
  procaddr = get_func(libssh, "api_sshaerosol_report_")
  call c_f_procpointer(procaddr, jst_cll)
  call jst_cll()
  procaddr = get_func(libssh, "api_sshaerosol_output_")
  call c_f_procpointer(procaddr, jst_cll)
  call jst_cll()

  ! Initialize photochemistry
  procaddr = get_func(libssh, "api_sshaerosol_initphoto_")
  call c_f_procpointer(procaddr, jst_cll)
  call jst_cll()

  ! time loop, 60 time steps
  do it = 1, 60

    ! update time in SSH
    if (it > 1) then
      time = time + dt
      procaddr = get_func(libssh, "api_sshaerosol_set_current_t_")
      call c_f_procpointer(procaddr, send_db)
      call send_db(time)
    endif

    ! Update photolysis
    procaddr = get_func(libssh, "api_sshaerosol_updatephoto_")
    call c_f_procpointer(procaddr, jst_cll)
    call jst_cll()

    ! Emission
    procaddr = get_func(libssh, "api_sshaerosol_emission_")
    call c_f_procpointer(procaddr, jst_cll)
    call jst_cll()

    ! Gas chemistry
    procaddr = get_func(libssh, "api_sshaerosol_gaschemistry_")
    call c_f_procpointer(procaddr, jst_cll)
    call jst_cll()

    ! Aerosol dynamic
    procaddr = get_func(libssh, "api_sshaerosol_aerodyn_")
    call c_f_procpointer(procaddr, jst_cll)
    call jst_cll()

    ! Built-in reporting
    procaddr = get_func(libssh, "api_sshaerosol_output_")
    call c_f_procpointer(procaddr, jst_cll)
    call jst_cll()

  enddo

  ! finalize
  procaddr = get_func(libssh, "api_sshaerosol_finalize_")
  call c_f_procpointer(procaddr, jst_cll)
  call jst_cll()

  ! Close the shared library
  call closelib(libssh)

  ! Free memory
  deallocate(cgas)

end program main
