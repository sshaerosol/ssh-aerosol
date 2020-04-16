!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------

module SSHaerosol
  implicit none

contains

  subroutine initialisation_modssh(namelist_ssh,iemis,naero,nspecies,nbins,nlayers,Vlayers)

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

    character (len=40) :: namelist_ssh  ! Configuration file
    integer naero,nspecies,nbins,iemis,nlayers
    double precision,allocatable,dimension(:) :: Vlayers
    
    ssh_standalone = .false.

    call dimensions(N_gas, n_reaction, n_photolysis)  

    call read_namelist(namelist_ssh)                  

    tag_emis=iemis
    i_compute_repart = 1
    i_write_repart = 0

    call read_inputs()

    call init_parameters()

    call init_distributions()

    naero=N_aerosol_layers
    nspecies=N_aerosol
    nbins=N_sizebin
    nlayers=nlayer
    allocate(Vlayers(nlayer))
    Vlayers=Vlayer

    !call free_allocated_memory()
    !IF (with_coag.EQ.1) call DeallocateCoefficientRepartition()
  end subroutine initialisation_modssh

!!$  subroutine launch_ssh_chemonly(namelist_ssh,ngas,nbins,lat,lon,duration,temp,pres,relh,atte,gas,gas_conc,aero_conc,number)
!!$    use aInitialization
!!$    use jAdaptstep
!!$    use bCoefficientRepartition
!!$    use mEmissions
!!$    use netcdf
!!$    use lDiscretization
!!$    use Resultoutput
!!$    use gCoagulation
!!$    use mod_photolysis
!!$
!!$    implicit none
!!$
!!$    integer :: t, j, k, s, js,jesp,day,nbins,nspecies,ngas
!!$    character (len=40) :: namelist_ssh  ! Configuration file
!!$    double precision :: duration,temp,pres,relh,atte
!!$    double precision :: ttmassaero = 0.d0, ttmass = 0.d0, totsulf = 0.d0
!!$    double precision,dimension(nbins) :: number
!!$    double precision,dimension(nspecies) :: gas_conc
!!$    double precision,dimension(nbins,nspecies) :: aero_conc
!!$
!!$    double precision :: t_since_update_photolysis                                                                                    
!!$
!!$    ! Read the number of gas-phase species and chemical reactions
!!$    call dimensions(N_gas, n_reaction, n_photolysis)  
!!$
!!$    call read_namelist(namelist_ssh)
!!$
!!$    if (.not.allocated(init_bin_number)) call read_namelist(namelist_ssh)   
!!$
!!$    tag_emis=0
!!$    i_compute_repart = 1
!!$    i_write_repart = 0
!!$
!!$    if (.not.allocated(molecular_weight)) call read_inputs() 
!!$    Temperature=temp
!!$    Pressure=pres   
!!$    Relative_Humidity=relh
!!$    initial_time=0.d0
!!$    final_time=duration
!!$    nt = int((final_time-initial_time) / delta_t)
!!$    pressure_sat = 611.2d0* dexp(17.67d0* (temperature - 273.15d0) / (temperature - 29.65d0))
!!$    Relative_Humidity = DMIN1(DMAX1(Relative_Humidity, Threshold_RH_inf), Threshold_RH_sup)
!!$    humidity =  1.d0/(Pressure/(pressure_sat *0.62197d0* Relative_Humidity)-1.d0)
!!$    attenuation=atte
!!$    concentration_gas_all=0.d0
!!$    init_bin_mass=0.d0
!!$    init_mass=0.d0                    
!!$
!!$    call init_parameters()  
!!$
!!$    call init_distributions()  
!!$
!!$    if ((tag_chem .ne. 0).AND.(option_photolysis.eq.2)) then
!!$       ! Allocate arrays for photolysis
!!$       call allocate_photolysis()    
!!$       ! Read photolysis rate for current day
!!$       current_time = initial_time
!!$       call init_photolysis() 
!!$       call interpol_photolysis()
!!$    endif
!!$    ! **** simulation starts 
!!$    t_since_update_photolysis = 0.d0
!!$
!!$    do t = 1, nt
!!$
!!$       current_time = initial_time + (t - 1) * delta_t
!!$       t_since_update_photolysis = t_since_update_photolysis +  delta_t
!!$
!!$       if (ssh_standalone) write(*,*) "Performing iteration #" // trim(str(t)) // "/" // trim(str(nt))
!!$       if (ssh_logger) write(logfile,*) "Performing iteration #" // trim(str(t)) // "/" // trim(str(nt))
!!$
!!$       ! Read the photolysis rates.
!!$       if ((tag_chem .ne. 0).AND.(option_photolysis.eq.2)) then
!!$          if (t_since_update_photolysis >= time_update_photolysis) then
!!$             call interpol_photolysis()
!!$             t_since_update_photolysis = 0.d0
!!$          endif
!!$       endif
!!$
!!$       ! Emissions
!!$       if (tag_emis .ne. 0) call emission(delta_t)
!!$
!!$       ! Gas-phase chemistry
!!$
!!$       ! initial and final physical parameters are set same
!!$       ! no volumetric emission
!!$       ! 0 : vertical gas volumetric emission    1 : with number 
!!$       ! 0 : not take into account cloud    0.d0 : air water content fracion sets to 0  
!!$
!!$       if (tag_chem .ne. 0) then
!!$          call chem(n_gas, n_reaction, n_photolysis, photolysis_reaction_index,&
!!$               ns_source, source_index, conversionfactor, conversionfactorjacobian,&
!!$               0, lwc_cloud_threshold, molecular_weight, &
!!$               current_time, attenuation, &
!!$               humidity, temperature,&
!!$               pressure, source, &
!!$               photolysis_rate, delta_t, attenuation,&
!!$               humidity, temperature,&
!!$               pressure, source, &
!!$               photolysis_rate, longitude,&
!!$               latitude, concentration_gas_all,&
!!$               0, with_heterogeneous, n_aerosol, n_size, n_fracmax,&
!!$               0.d0,&
!!$               diam_bound, fixed_density, &
!!$               wet_diameter, &
!!$               heterogeneous_reaction_index, &
!!$               concentration_mass,&
!!$               with_adaptive, adaptive_time_step_tolerance,&
!!$               min_adaptive_time_step, option_photolysis, ind_jbiper, ind_kbiper,&
!!$               1, not(with_fixed_density), concentration_number, &
!!$               mass_density)
!!$       end if
!!$
!!$       ! re-calculate total_mass(N_aerosol) because mass change due to gas-phase chemistry  
!!$       total_aero_mass = 0.d0
!!$       total_mass = 0.d0
!!$       do s = 1, N_aerosol_layers
!!$          jesp = List_species(s)
!!$          do j=1,N_size
!!$             total_aero_mass(jesp) = total_aero_mass(jesp) + concentration_mass(j,s)
!!$          enddo
!!$       enddo
!!$       ! update mass conc. of aerosol precursors
!!$       ! concentration_gas_all(precursor_index) -> concentration_gas(n_aerosol)
!!$       do s = 1, N_aerosol
!!$          if (aerosol_species_interact(s) .gt. 0) then
!!$             concentration_gas(s) = concentration_gas_all(aerosol_species_interact(s))
!!$          end if
!!$          total_mass(s) = total_mass(s) + concentration_gas(s) + total_aero_mass(s)
!!$       end do
!!$    end do			! finsh simulation
!!$
!!$  end subroutine launch_ssh_chemonly

  subroutine launch_ssh_aerosolonly(namelist_ssh,iemis,naero,nbins,nspecies,duration,temp,pres,relh,gas_conc,aero_conc,number)    
    ! Computes aerosol formation
    !    With internal mixing
    !    Inviscid organic aerosol (only one layer) in that case naero=nspecies=number of aerosol species
    !    Viscous organic aerosol, 1 species for each organic aerosol layers

    !iemis=1 ; to use emissions. For test cases only

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

    integer :: t, j, k, s, js,jesp,day,nbins,nspecies,iemis,naero
    character (len=40) :: namelist_ssh  ! Configuration file
    double precision :: duration,temp,pres,relh
    double precision :: ttmassaero = 0.d0, ttmass = 0.d0, totsulf = 0.d0
    double precision,dimension(nbins) :: number
    double precision,dimension(nspecies) :: gas_conc
    double precision,dimension(nbins,naero) :: aero_conc

    double precision :: t_since_update_photolysis

    ! Initialisation: discretization and distribution  
    !call getarg(1, namelist_ssh)

    ssh_standalone = .false.

    ! Read the number of gas-phase species and chemical reactions    
    call dimensions(N_gas, n_reaction, n_photolysis)    

    if (.not.allocated(init_bin_number)) call read_namelist(namelist_ssh)    

    tag_emis=iemis
    i_compute_repart = 1
    i_write_repart = 0

    if (.not.allocated(molecular_weight)) call read_inputs()    

    if (N_aerosol_layers.ne.naero) then
       print*,"Number of aerosol species (including layers) should be equal to : ",N_aerosol_layers
       stop
    endif
    
    if (N_aerosol.ne.nspecies) then
       print*,"Number of aerosol species should be equal to : ",N_aerosol
       stop
    endif

    if (N_sizebin.ne.nbins) then
       print*,"Number of bins species should be equal to : ",N_sizebin
       stop
    endif

    Temperature=temp
    Pressure=pres   
    Relative_Humidity=relh
    initial_time=0.d0
    final_time=duration
    nt = int((final_time-initial_time) / delta_t)
    pressure_sat = 611.2d0* dexp(17.67d0* (temperature - 273.15d0) / (temperature - 29.65d0))
    Relative_Humidity = DMIN1(DMAX1(Relative_Humidity, Threshold_RH_inf), Threshold_RH_sup)
    humidity =  1.d0/(Pressure/(pressure_sat *0.62197d0* Relative_Humidity)-1.d0)
    concentration_gas_all=0.d0
    init_bin_mass=0.d0
    init_mass=0.d0

    do js=1,N_aerosol       
       concentration_gas(js)=gas_conc(js)
       if (aerosol_species_interact(js)>0) then
          concentration_gas_all(aerosol_species_interact(js))=concentration_gas(js)
       endif
    enddo

    do s = 1, N_aerosol_layers
       jesp = List_species(s)
       do j=1,N_size
          init_bin_mass(j,jesp) = init_bin_mass(j,jesp) + concentration_mass(j,s)
          init_mass(jesp)=init_mass(jesp) + aero_conc(j,js)
       enddo
    enddo
    aero_total_mass=sum(aero_conc(:,1:N_aerosol_layers))
    
    do k=1,N_sizebin
       init_bin_number(k)=number(k)
    enddo
   
    call init_parameters()
        
    call init_distributions()   

    do s = 1, N_aerosol_layers
       jesp = List_species(s)
       do j=1,N_size
          concentration_mass(j,s)=aero_conc(j,s)
       enddo
    enddo
    
    do t = 1, nt
       print*,t,sum(concentration_mass(1,93:97))
       
       ! Emissions
       if (tag_emis .ne. 0) call emission(delta_t)
       
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

    end do			! finsh simulation

    do s = 1, N_aerosol
       gas_conc(s)=concentration_gas(s)
    enddo
    do s = 1, N_aerosol_layers
       aero_conc(:,s)=concentration_mass(:,s)
    enddo
    number(:)=concentration_number(:)

    !call delete_empty_file() ! delete empty output files
    !call free_allocated_memory()
    !IF (with_coag.EQ.1) call DeallocateCoefficientRepartition()

  end subroutine launch_ssh_aerosolonly

end module SSHaerosol