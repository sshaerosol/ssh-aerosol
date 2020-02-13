!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------

module SSHaerosol
  implicit none

contains

  subroutine initialisation_repart_coeff(namelist_ssh)

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
    ssh_standalone = .false.

    call dimensions(N_gas, n_reaction, n_photolysis)  

    call read_namelist(namelist_ssh)                  

    tag_emis=0
    i_compute_repart = 1
    i_write_repart = 1
    
    call read_inputs()

    call init_parameters()

    call init_distributions()
    
    call free_allocated_memory()
    IF (with_coag.EQ.1) call DeallocateCoefficientRepartition()
  end subroutine 
  
  
  subroutine launch_ssh_aerosolonly(namelist_ssh,nbins,nspecies,duration,lat,lon,temp,pres,relh,gas_conc,aero_conc,number)

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

    integer :: t, j, k, s, js,jesp,day,nbins,nspecies
    character (len=40) :: namelist_ssh  ! Configuration file
    double precision :: duration,lat,lon,temp,pres,humid,relh
    double precision :: ttmassaero = 0.d0, ttmass = 0.d0, totsulf = 0.d0
    double precision,dimension(nbins) :: number
    double precision,dimension(nspecies) :: gas_conc
    double precision,dimension(nbins,nspecies) :: aero_conc

    double precision :: t_since_update_photolysis

    ! Initialisation: discretization and distribution  
    !call getarg(1, namelist_ssh)

    ssh_standalone = .false.

    ! Read the number of gas-phase species and chemical reactions
    call dimensions(N_gas, n_reaction, n_photolysis)  

    call read_namelist(namelist_ssh)                  

    tag_emis=0
    i_compute_repart = 0
    i_write_repart = 0
    
    call read_inputs_light()

    if (N_aerosol.ne.nspecies) then
       print*,"Number of aerosol species should be equal to : ",N_aerosol
       stop
    endif

    if (N_sizebin.ne.nbins) then
       print*,"Number of bins species should be equal to : ",N_sizebin
       stop
    endif
    
    latitude=lat
    longitude=lon
    Temperature=temp
    Pressure=pres
    Humidity=humid
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
       do k=1, N_sizebin
          init_bin_mass(k,js) = aero_conc(k,js)
          init_mass(js)=init_mass(js) + aero_conc(k,js)
       enddo
    enddo
    aero_total_mass=sum(init_mass(1:N_aerosol)) 

    do k=1,N_sizebin
       init_bin_number(k)=number(k)
    enddo    

    call init_parameters()

    call init_distributions()    
    
    do t = 1, nt
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

    end do			! finsh simulation

    do s = 1, N_aerosol
       gas_conc(s)=concentration_gas(s)
       aero_conc(:,s)=concentration_mass(:,s)
    enddo    
    number(:)=concentration_number(:)

    !call delete_empty_file() ! delete empty output files
    call free_allocated_memory()
    IF (with_coag.EQ.1) call DeallocateCoefficientRepartition()

  end subroutine launch_ssh_aerosolonly

end module SSHaerosol
