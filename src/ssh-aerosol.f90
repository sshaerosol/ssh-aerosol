
PROGRAM SSHaerosol


  use aInitialization
  use jAdaptstep
  use bCoefficientRepartition
  use mEmissions
  use netcdf
  use lDiscretization
  use Resultoutput
  use gCoagulation

  implicit none

  integer :: t, j, k, s, tag_test = 0
  character (len=40) :: namelist_ssh  ! Configuration file
  character(len=10)::vchar ! Current time 
  double precision :: current_time, ttmassaero = 0.d0, ttmass = 0.d0, totsulf = 0.d0

  ! Initialisation: discretization and distribution  
  call getarg(1, namelist_ssh) 

  ! N_gas = 93; N_reaction = 206; N_photolysis = 24 
  call dimensions(N_gas, n_reaction, n_photolysis)  

  call read_namelist(namelist_ssh)                  

  call read_inputs()                                

  call init_parameters()  

  call init_distributions()  

  call init_output_conc() 
  ! mass/number/total number/total mass conc. for all species over each time step

  call init_output_vl()
  ! to test values evolution : cell_diam_av, wet_diameter

  call save_report()

  call save_values() 

  call save_concentration() 


  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! simulation starts ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  do t = 1, nt
     current_time = initial_time + (t - 1) * delta_t

     write(*,*)

     write(*,*) "Performing iteration #" // trim(str(t)) // "/" // trim(str(nt))

     write(vchar, '(f8.0)') current_time

     write(*,*) "Current time: " // trim(adjustl(vchar))
     


     ! Emissions
     if (tag_emis .ne. 0) call emission(current_time, delta_t)


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
	print*, 'called chem()' 
      end if

	! re-calculate total_mass(N_aerosol) because mass change due to gas-phase chemistry  
    total_aero_mass = 0.d0
    total_mass = 0.d0
    do s = 1, N_aerosol
       do j=1,N_size
         total_aero_mass(s) = total_aero_mass(s) + concentration_mass(j,s)
       enddo
	! update mass conc. of aerosol precursors
	! concentration_gas_all(precursor_index) -> concentration_gas(n_aerosol)
       if (aerosol_species_interact(s) .gt. 0) then
          concentration_gas(s) = concentration_gas_all(aerosol_species_interact(s))
       end if
          total_mass(s) = total_mass(s) + concentration_gas(s) + total_aero_mass(s)
    end do


    ! Aerosol dynamic
    CALL AERODYN(current_time,delta_t)
	print*,'called aerodyn()'

	! update mass conc. of aerosol precursors
	! concentration_gas(n_aerosol) -> concentration_gas_all(precursor_index)
    do s = 1, N_aerosol
       if (aerosol_species_interact(s) .gt. 0) then
          concentration_gas_all(aerosol_species_interact(s)) = concentration_gas(s)
       end if
    end do



    call save_concentration()         ! Text or Binary format outout
    call save_values()
  end do			! finsh simulation


  OPEN(UNIT=10,FILE=trim(output_directory) // "/" // "report.txt",status='old', position = "append")
	write(unit=10,FMT=*)
	write(unit=10,FMT=*), 'AFTER SIMULATION', nt, 'steps :'
	write(unit=10,FMT=*), '====== after concentration_number : ======'
	write(unit=10,FMT=*), concentration_number
	write(unit=10,FMT=*)
	write(unit=10,FMT=*), '===== after diam_bounds : ====='
	write(unit=10,FMT=*), diam_bound
	write(unit=10,FMT=*)
	write(unit=10,FMT=*), '===== after cell_diam_av : ====='
	write(unit=10,FMT=*), cell_diam_av 
	write(unit=10,FMT=*)
	write(unit=10,FMT=*), '===== after wet_volume : ====='
	write(unit=10,FMT=*), wet_volume
	write(unit=10,FMT=*)
	write(unit=10,FMT=*), '===== after concentration_mass : ====='
	write(unit=10,FMT=*), aerosol_species_name(4)
	write(unit=10,FMT=*), (concentration_mass(j,4), j = 1, N_size)
	write(unit=10,FMT=*)
	write(unit=10,FMT=*), aerosol_species_name(34)
	write(unit=10,FMT=*), (concentration_mass(j,34), j = 1, N_size)


  CLOSE(10)
  call free_allocated_memory()

  IF (with_coag.EQ.1) call DeallocateCoefficientRepartition()
  

  write(*,*) "============================================"
  write(*,*) "==== SSH-aerosol simulation completed  ====="
  write(*,*) "============================================"
  
end PROGRAM SSHaerosol
