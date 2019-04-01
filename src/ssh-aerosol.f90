!!-----------------------------------------------------------------------
PROGRAM SSHaerosol

  use lDiscretization
  use aInitialization
  use jAdaptstep
  use Resultoutput
  use bCoefficientRepartition
  use mEmissions
!  use h2o
  ! use Coagulation
  ! use Congregation
  use netcdf
  
  implicit none
  
  integer ::solver, n
  character (len=40) :: configuration_file!initial configure file
  character(len=10)::vchar
  double precision :: current_time,equilibrium_time
  double precision :: timestep_coag,timestep_cond
  double precision :: cpu_t1,cpu_t2
  double precision :: temp_mass,temp_mass_2
  double precision :: total_ms
  integer :: j, b, f, g, k, id, jesp, jt, t
  double precision :: tfchem_aer, tschem_aer, dtchem_aer
  integer :: ns_source, ncycle
  integer, dimension(:), allocatable :: source_index
  double precision, dimension(:), allocatable :: conversionfactor
  double precision, dimension(:,:), allocatable :: conversionfactorjacobian
  double precision, dimension(:), allocatable :: source
  double precision :: attenuation, liquid_water
  integer :: heterogeneous_reaction_index(4)
  integer :: ind_jbiper, ind_kbiper

  call cpu_time(cpu_t2)

  ! Initialisation: discretization and distribution  
  call getarg(1, configuration_file)!obtain the path of configuration file

  ! Read the number of gas-phase species, chemical reactions and photolysis
  call dimensions(N_gas, n_reaction, n_photolysis)

!  call read_pointer()

  call read_configuration(configuration_file)

  call read_pointer()

  call read_parameter()

  call init_distribution()

  call init_output()

  ! Compute the number of iteration, nt
  nt = int(final_time / delta_t)
  write(*,*) "Number of iterations:", nt

  ns_source = 1
  allocate(source_index(ns_source))
  source_index = [1]
  allocate(source(ns_source))
  source = 0.d0

  allocate(conversionfactor(n_gas))
  allocate(conversionfactorjacobian(n_gas, n_gas))

  do j = 1, n_gas
     conversionfactor(j) = Navog * 1.d-12 / molecular_weight(j)
     do jesp = 1, n_gas
        conversionfactorjacobian(j, jesp) = molecular_weight(j) / &
             molecular_weight(jesp)
     enddo
  enddo

  ncycle = 1
  with_heterogeneous = 1
  attenuation = 0.d0
  liquid_water = 0.d0
  with_adaptive = 1
  adaptive_time_step_tolerance = 0.001d0
  min_adaptive_time_step = 10.d0

  ! This index is modified by adding one later
  ! See ispeclost in Chemistry/common/hetrxn.f
  heterogeneous_reaction_index(1)= 86-1 ! HO2
  heterogeneous_reaction_index(2)= 85-1 ! NO2
  heterogeneous_reaction_index(3)= 78-1 ! NO3
  heterogeneous_reaction_index(4)= 23-1 ! N2O5

  opt_photolysis = 2
  ind_jbiper = 14 ! aerosol species index
  ind_kbiper = 24 ! pholysis index
  
  do t = 1, nt

     current_time = initial_time + (t - 1) * delta_t

     write(*,*) "Performing iteration #" // trim(str(t)) // "/" // trim(str(nt))
     write(vchar, '(f8.0)') current_time
     write(*,*) "Current time: " // trim(adjustl(vchar))
     
     call emission(current_time, delta_t)

     write(*,*) "Concentration NO2 after emission", concentration_gas_all(94)

     ! Gas-phase chemistry
     call chem(n_gas, n_reaction, n_photolysis, photolysis_reaction_index,&
          ns_source, source_index, conversionfactor, conversionfactorjacobian,&
          0, lwc_cloud_threshold, molecular_weight, &
          current_time, 0.d0, &
          humidity, temperature,&
          pressure, source, &
          photolysis, delta_t, attenuation,&
          humidity, temperature,&
          pressure, source, &
          photolysis, ncycle, longitude,&
          latitude, concentration_gas_all,&
          0, with_heterogeneous, n_aerosol, n_size, n_fracmax,&
          liquid_water,&
          bin_bound, fixed_density, &
          wet_diameter, &
          heterogeneous_reaction_index, &
          concentration_mass,&
          with_adaptive, adaptive_time_step_tolerance,&
          min_adaptive_time_step, opt_photolysis, ind_jbiper, ind_kbiper,&
          with_number, not(with_fixed_density), concentration_number, &
          mass_density) 

     write(*,*) "Concentration NO2 after gas-phase chem", concentration_gas_all(94)

     ! Aerosol chemistry
     ncycle_aer = 1
     DO Jt=1,ncycle_aer
        
        tschem_aer=initial_time+(Jt-1)*delta_t/ncycle_aer
        tfchem_aer=tschem_aer+delta_t/ncycle_aer
        dtchem_aer=tfchem_aer-tschem_aer
        
        !initialize gas
        DO jesp=1,N_gas
           if(IsNaN(concentration_gas_all(jesp)*0.d0)) then
              print*,"Error of infinity/NaN initial",concentration_gas_all(jesp)
              stop
           endif
           IF (concentration_gas_all(jesp).LT.0.D0) THEN
              concentration_gas_all(jesp) = 0.D0
           ENDIF
        ENDDO
        
        !     Number concentration: loop on bins
        IF (with_number.EQ.1) THEN
           DO j=1, N_size
              if(IsNaN(DLnum_conc_aer(j)*0.d0)) then
                 print*,"Error of infinity/NaN initial",DLnum_conc_aer(j)
                 stop
              endif
              concentration_number(j) = DLnum_conc_aer(j)
           ENDDO
        ENDIF
         
        total_ms=0.d0
        total_aero_mass=0.d0
        DO j=1,N_size
           DO jesp=1,N_aerosol-1
              if(IsNaN(concentration_mass(j,jesp)*0.d0)) then
                 print*,"Error of infinity/NaN initial",concentration_mass(j,jesp)
                 stop
              endif
              total_aero_mass(jesp)=total_aero_mass(jesp)+concentration_mass(j,jesp)
              total_ms=total_ms+concentration_mass(j,jesp)
           ENDDO
        ENDDO

     
        DO jesp=1,N_aerosol
           IF(aerosol_species_interact(jesp).GT.0) THEN
              concentration_gas(jesp)=concentration_gas_all(aerosol_species_interact(jesp))
           ENDIF
           total_mass(jesp)=total_aero_mass(jesp)+concentration_gas(jesp)
           total_mass_old(jesp)= total_mass(jesp)
        ENDDO
        
     
        density_aer_bin=fixed_density
        density_aer_size=fixed_density

        if(with_fixed_density.ne.1) then! modified
           !fixed density become overall averaged density
           call compute_all_density()
        endif
        
        call check_nan_inf(1)

!      !!     Aerosol discretization converted into microm.
!      do b= 1, (N_sizebin+1)
! !        diam_bound(b)=bin_bound_aer(b)* 1.D06
! !        write(*,*) b, diameter(b), diam_bound(b)
!         if(b.eq.N_sizebin+1) then
!            mass_bound(b)= density_aer_size(b-1)* cst_pi6 * diam_bound(b)**3
!         else
!            mass_bound(b)= density_aer_size(b)* cst_pi6 * diam_bound(b)**3
!         endif
!         log_bound(b)=dlog10(mass_bound(b))
!      end do

     ! do b =1,N_sizebin
     !    size_diam_av(b)=dsqrt(diam_bound(b)*diam_bound(b+1))
     !    size_log_av(b)=(log_bound(b) + log_bound(b+1))*5.D-01
     !    size_sect(b)=log_bound(b+1) - log_bound(b)
     !    ! YK need to check if size_sect is ok.
     !   write(*,*) "size_sect in main",b,size_sect(b), log_bound(b+1), log_bound(b)
     ! enddo


        IF (with_number.NE.1) THEN
           call compute_number()
        ENDIF

        call compute_average_bin_diameter()
     
        do j =1,N_size
           b=concentration_index(j,1)
           cell_log_av(j)=size_log_av(b)
        enddo

        call compute_average_diameter()
     
        call compute_wet_mass_diameter(1,N_size,concentration_mass,concentration_number,&
          concentration_inti,wet_mass,wet_diameter,wet_volume)
        
         
        CALL AERODYN(tschem_aer,tfchem_aer)


     ! DO jesp=1,N_aerosol
     !    IF(aerosol_species_interact(jesp).GT.0) THEN
     !       concentration_gas(jesp)=concentration_gas_all(aerosol_species_interact(jesp))
     !    ENDIF
     !    total_mass(jesp)=total_aero_mass(jesp)+concentration_gas(jesp)
     ! ENDDO

        DO jesp=1,N_aerosol-1
           total_ms=total_mass(jesp)-total_mass_old(jesp)
           if((total_ms*total_ms).gt.(TINYN)) then
              print*,"warning! mass not conserved - species:",jesp,&
                   total_ms,total_aero_mass(jesp),concentration_gas(jesp)
           endif
           if(concentration_gas(jesp).lt.0.d0) then
              print*,"warning! gas<0 species:",jesp,concentration_gas(jesp)
              concentration_gas(jesp)=0.d0
           endif
           IF(aerosol_species_interact(jesp).GT.0) THEN
              concentration_gas_all(aerosol_species_interact(jesp))=concentration_gas(jesp)
           ENDIF
        ENDDO

     write(*,*) "Concentration NO2 after aerodyn", concentration_gas_all(94)

        
        DO j=1,N_size
           DO jesp=1,N_aerosol
              if(IsNaN(concentration_mass(j,jesp)*0.d0)) then
                 print*,"Error of infinity/NaN end",concentration_mass(j,jesp)
                 stop
              endif
           ENDDO
        ENDDO
        
        !!     Number concentration
        IF(with_number.EQ.1) THEN
           DO j=1,N_size
              DLnum_conc_aer(j) = concentration_number(j)
           ENDDO
        ENDIF
         
        DO j=1,N_size
           if(IsNaN(DLnum_conc_aer(j)*0.d0)) then
              print*,"Error of infinity/NaN end",DLnum_conc_aer(j)
              stop
           endif
        ENDDO
        
     ENDDO

     ! ! Text format outout
     ! call write_result()

     ! Binary format output
     call save_result()


  end do

  CALL free_allocated_memory()

  IF (with_coag.EQ.1) THEN
     call DeallocateCoefficientRepartition()
  ENDIF

  deallocate(source_index)
  deallocate(source)
  deallocate(conversionfactor)
  deallocate(conversionfactorjacobian)

  write(*,*) "============================================"
  write(*,*) "==== SSH-aerosol simulation completed  ====="
  write(*,*) "============================================"
  
end PROGRAM SSHaerosol
