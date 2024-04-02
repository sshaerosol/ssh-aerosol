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
  use mod_meteo
  use mod_sshchem, only: ssh_chem_twostep,ssh_chem
  use mod_sshchemkinetic, only: ssh_gck_compute_gas_phase_water
  
  implicit none

  integer :: t, j, s, jesp, day, tag_genoa_in  
  character (len=400) :: namelist_ssh  ! Configuration file
  double precision, dimension(:), allocatable :: timer

  double precision :: t_since_update_photolysis, t0

  ! genoa char to read index for initial set
  character (len=20) :: ivoc0
  ! genoa timestep for aerosol_dynamic used when keep_gp = 1
  double precision :: delta_t2
  double precision :: emw_tmp

  ! Initial time (seconds)
  call cpu_time(t0)
  
  ! Initialize input args used for genoa v3
  initID = "-" ! initial set id
  chemID = "-" ! chem id
  chemID2= "-" ! chemID2 = chemID/chemID
  resID  = "-" ! result id
  tag_genoa_in = 0 ! normal mode

  ! Initialisation: discretization and distribution
  if (iargc() == 0) then
    write(*,*) "usage: ssh-aerosol namelist.input"
    stop
  else
    call getarg(1, namelist_ssh) ! Read namelist file
    !!! Inputs for gneoa START !!!
    ! No.inputs can be: 2, 3, 5 - read up to 5
    if (iargc() .ge. 2) then ! Read flag and settings for GENOA reductions ! 
      ! Read the #2 input as tag_genoa
      call getarg(2, ivoc0)
      read(ivoc0, *, iostat=s) tag_genoa_in
      if (s /= 0) then
        print*, "Error: Input genoa mode is not an integer: ", trim(ivoc0)
        stop
      endif
      ! Update tag_genoa: 0 = no genoa, 1 = fast mode, or 2 = complete mode
      if (tag_genoa_in .eq. 1 .or. tag_genoa_in .eq. 2) then
        tag_genoa = tag_genoa_in
        print*, "Under GENOA mode!", tag_genoa
      else if (tag_genoa_in .ne. 0) then
        print*, "Error: Input unknown tag_genoa mode: ", tag_genoa_in
        stop
      endif

      ! Read #3 input for init id
      if (iargc() .ge. 3) then
        initID = trim(adjustl(ivoc0))
        print*,"Read init ID: ", trim(initID)
        ! Read #4 and #5 inputs for chem id & res id
        if (iargc() .ge. 5) then
          call getarg(3,ivoc0)
          chemID = trim(adjustl(ivoc0))
          chemID2 = trim(chemID)//"/"//trim(chemID)
          call getarg(4,ivoc0)
          resID = trim(adjustl(ivoc0))
          print*, "Read chem & result IDs:",trim(chemID)," ",trim(resID) 
        endif
      endif
    endif
    !!! Inputs for genoa END !!!
  end if

  ! Read the number of gas-phase species and chemical reactions
  !call ssh_dimensions(N_gas, n_reaction, n_photolysis)

  call ssh_read_namelist(namelist_ssh)

  call ssh_read_inputs(0)                                

  call ssh_read_meteo()
  
  ! Read the meteorological data before output
  temperature = temperature_array(1)
  pressure = pressure_array(1)
  humidity = humidity_array(1)
  relative_humidity = relative_humidity_array(1) 

  call ssh_init_parameters()
  
  call ssh_init_coag()  

  call ssh_init_distributions()  

  ! Initialize output
  if (tag_genoa.eq.1.or.output_type.eq.0)  then  ! genoa fast

     call ssh_init_output_genoa()

     call ssh_save_concentration_genoa()

  else ! default settings

     call ssh_init_output()
     
     call ssh_save_report()
     
     call ssh_save_concentration()

  endif

  ! save total soas in concs.txt (tag_genoa = 1/2)
  if (tag_genoa.gt.0) call ssh_save_total_soa_genoa()
  
  if ((tag_chem .ne. 0).AND.(option_photolysis.eq.2)) then
    ! Allocate arrays for photolysis
    call ssh_allocate_photolysis()    
    ! Read photolysis rate for current day
    current_time = initial_time
    call ssh_init_photolysis() 
    call ssh_interpol_photolysis()
  endif
  
  ! **** simulation starts 
  t_since_update_photolysis = 0.d0

  ! Initialization is finished
  allocate(timer(nt+3))
  timer(1) = t0
  call cpu_time(t0)
  timer(2) = t0

  if (with_cond.EQ.1.and.kwall_gas>0.d0) then
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

  do t = 1, nt

     current_time = initial_time + (t - 1) * delta_t
     t_since_update_photolysis = t_since_update_photolysis +  delta_t
     delta_t2 = delta_t

     if (ssh_standalone) write(*,*) "Performing iteration #" // trim(str(t)) // "/" // trim(str(nt))
     if (ssh_logger) write(logfile,*) "Performing iteration #" // trim(str(t)) // "/" // trim(str(nt))

     ! Read the photolysis rates.
     if ((tag_chem .ne. 0).AND.(option_photolysis.eq.2)) then
       if (t_since_update_photolysis >= time_update_photolysis) then
        call ssh_interpol_photolysis()
        t_since_update_photolysis = 0.d0
       endif
     endif

     ! Emissions
     if (tag_emis .ne. 0) call ssh_emission(delta_t)

     ! Read the meteorological data.
     temperature = temperature_array(t)
     pressure = pressure_array(t)
     humidity = humidity_array(t)
     relative_humidity = relative_humidity_array(t) 
     
     ! Compute sumM using ideal gas law
     ! Compute third body.
     ! Conversion = Avogadro*1d-6/Perfect gas constant.
     SumMc = pressure * 7.243D16 / temperature
     
     ! Number of water molecules computed from the massic fraction
! (absolute humidity)
     YlH2O = 29.d0*SumMc*humidity/(18.d0+11.d0*humidity)
     !call compute_gas_phase_water(temperature,relative_humidity,YlH2O)

     ! Gas-phase chemistry

     ! initial and final physical parameters are set same
     ! no volumetric emission
     ! 0 : vertical gas volumetric emission    1 : with number 
     ! 0 : not take into account cloud    0.d0 : air water content fracion sets to 0  

     if (tag_chem .ne. 0) then

       if (tag_twostep .ne. 1) then
           call ssh_chem()
       else
         if (t==1.and.keep_gp.eq.1) then
         
            call ssh_chem_twostep(current_time,0.001*delta_t,t) ! current time, timestep

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
            CALL SSH_AERODYN(current_time,0.001*delta_t)

            ! update mass conc. of aerosol precursors
            ! concentration_gas(n_aerosol) -> concentration_gas_all(precursor_index)
            do s = 1, N_aerosol
               if (aerosol_species_interact(s) .gt. 0) then
                  concentration_gas_all(aerosol_species_interact(s)) = concentration_gas(s)
               end if
            end do
            
            delta_t2=0.999*delta_t
            
         endif
        
         ! solve chemistry with the two-step time numerical solver if tag_twostep .eq. 1
         call ssh_chem_twostep(current_time,delta_t2,t)
         
      end if
    end if ! finish chem

    ! set cst_aero(n_species,n_size,n_step) if need. N_sizebin is assumed to be N_size. Only for internal mixing.
    if (ncst_aero .gt. 0) then
      do s = 1, ncst_aero
        jesp = cst_aero_index(s)
        do j=1,N_size
          concentration_mass(j,jesp) = cst_aero(s,j,t)
        enddo
      enddo
    endif

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
    CALL SSH_AERODYN(current_time,delta_t2)

    ! update mass conc. of aerosol precursors
    ! concentration_gas(n_aerosol) -> concentration_gas_all(precursor_index)
    do s = 1, N_aerosol
       if (aerosol_species_interact(s) .gt. 0) then
          concentration_gas_all(aerosol_species_interact(s)) = concentration_gas(s)
       end if
    end do

    ! Save outputs
    if (tag_genoa.eq.1.or.output_type.eq.0) then ! genoa fast
        call ssh_save_concentration_genoa() ! only output certain species
    else  ! default
        call ssh_save_concentration()
    endif

   ! save concentrations in concs.txt
   if (tag_genoa.gt.0) call ssh_save_total_soa_genoa()

    ! Time step is finished
    call cpu_time(t0)
    timer(t+2) = t0

  end do ! finsh simulation

  ! Compute errors - genoa
  if (nref_file.gt.0) call ssh_compute_error_genoa()

  ! Write outputs
  if (tag_genoa.ne.1) call ssh_write_output()  !Creation of .txt .bin or .nc output files

  if (tag_genoa.ge.1.or.output_type.eq.0) call ssh_write_output_genoa()

  call ssh_free_allocated_memory()
  IF (with_coag.EQ.1) call ssh_DeallocateCoefficientRepartition()
 
  ! Desallocate arrays for photolysis
  if ((tag_chem .ne. 0).AND.(option_photolysis.eq.2)) then
    call ssh_deallocate_photolysis()    
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
  if (allocated(timer)) deallocate(timer)
end PROGRAM SSHaerosol
