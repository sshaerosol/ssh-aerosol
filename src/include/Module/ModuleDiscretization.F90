!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------

MODULE lDiscretization
  use aInitialization
  use bCoefficientRepartition
  use cThermodynamics
  use dPhysicalbalance
  use eRedistribution

  implicit none

contains


  ! ===================================================================
  !
  !          initial all the necessary parameters
  !  ! read species file and initilise pointers/parameters used in the 
  !  ! simulation
  !
  ! =================================================================== 
  subroutine ssh_init_parameters()

    integer i,b,j,f,k,s
    double precision :: totv,tmp

    ! pointer for cond_time
    cond_time_index(1)=ENH4
    cond_time_index(2)=ENO3
    cond_time_index(3)=ECL

    pH = 4.d5
    n_emis=0.d0
    m_emis=0.d0

    !# Liquid water content threshold above which cloud is present (in g/m3).
    lwc_cloud_threshold = 0.05d0

    ! discretization
    if(.not.allocated(N_fracbin)) allocate(N_fracbin(N_sizebin))
    if(N_frac.eq.1)  N_fracbin = 1
    if (.not.allocated(discretization_composition)) call ssh_discretization()
    if (ssh_standalone) write(*,*) "N_size :", N_size
    if (ssh_logger) write(logfile,*) "N_size :", N_size
    if(.not.allocated(concentration_index)) allocate(concentration_index(N_size, 2))
    if(.not.allocated(concentration_index_iv)) allocate(concentration_index_iv(N_sizebin, N_fracmax))
    j = 1
    do k = 1,N_sizebin
       do i = 1, N_fracbin(k)
	  concentration_index(j, 1) = k
	  concentration_index(j, 2) = i
	  concentration_index_iv(k,i) = j
	  j = j + 1
       enddo
    enddo

    ! initialise bin bounds if need.
    if(.not.allocated(diam_bound)) allocate(diam_bound(N_sizebin+1))
    if (tag_dbd == 0) then  ! auto-generate bin bounds (method need to change in order to fit PM10 & PM2.5)
       do i = 1,N_sizebin+1
          diam_bound(i)= diam_input(1) * (diam_input(2)/ diam_input(1))**((i - 1) / dble(N_sizebin))
       enddo  ! set closer bounds to 1, 2.5, 10 if bounds auto-generated ??
       if (ssh_standalone) write(*,*) "sizebin bound is auto-generated."
       if (ssh_logger) write(logfile,*) "sizebin bound is auto-generated."
    else if (tag_dbd == 1) then 
       diam_bound = diam_input
    end if

    !calculate ICUT the corresponding cell index of the cutting diameter if tag_icut = 0
   if (tag_icut .eq. 0) then
     if(Cut_dim.gt.diam_bound(1)) then
       if (cut_dim .gt. diam_bound(N_sizebin + 1)) cut_dim = diam_bound(N_sizebin + 1)
       do k= 1,N_sizebin
          if(diam_bound(k).lt.Cut_dim.and.diam_bound(k+1).ge.Cut_dim) then
             ICUT=concentration_index_iv(k,N_fracbin(k))
          endif
       enddo
     else
       ICUT=0
     endif
     if (ssh_standalone) write(*,*)'Cut_dim :',Cut_dim, 'ICUT :', ICUT
     if (ssh_logger) write(logfile,*)'Cut_dim :',Cut_dim, 'ICUT :', ICUT
   else
     if (ssh_standalone) write(*,*)'tag_icut :',tag_icut,'Cut_dim :',Cut_dim
     if (ssh_logger) write(logfile,*)'tag_icut :',tag_icut,'Cut_dim :',Cut_dim
   endif

    IF(ICUT.GT.0.D0) THEN
       section_pass=concentration_index(ICUT,1)
    ELSE
       section_pass=1
    ENDIF
    
    if (with_cond.EQ.1) then  ! for condensation

       if(.not. allocated(quadratic_speed)) allocate(quadratic_speed(N_aerosol))    
       quadratic_speed=0.D0
       if(.not. allocated(diffusion_coef)) allocate(diffusion_coef(N_aerosol))
       diffusion_coef=0.D0

    end if
    
    ! diameter !
    if(.not.allocated(size_diam_av)) allocate(size_diam_av(N_sizebin))
    size_diam_av = 0.0
    if(.not.allocated(size_mass_av)) allocate(size_mass_av(N_sizebin))
    size_mass_av=0.d0

    do k= 1,N_sizebin
       size_diam_av(k) = (diam_bound(k) * diam_bound(k+1)) ** 0.5
       if(size_diam_av(k) .lt. diam_bound(k)) size_diam_av(k) = diam_bound(k)
    end do

    if(.not. allocated(density_aer_bin)) allocate(density_aer_bin(N_size))
    density_aer_bin = 0.0
    if(.not. allocated(density_aer_size)) allocate(density_aer_size(N_sizebin))
    density_aer_size = 0.0
    if(.not. allocated(rho_wet_cell)) allocate(rho_wet_cell(N_size))
    rho_wet_cell = 0.0

    if (with_fixed_density.eq.1) then
       ! convert from kg/m3 to µg/µm3 or µg/m3
       !rho1 = fixed_density * 1.0d-9    !µg/µm3     
       !rho2 = fixed_density * 1.0d+9    !µg/m3	       
       mass_density(EH2O)=fixed_density
       density_aer_bin=fixed_density
       density_aer_size=fixed_density
       rho_wet_cell = fixed_density
    endif
    do k = 1, N_sizebin
       size_mass_av(k) = (size_diam_av(k)**3)*fixed_density*cst_PI6
    end do
     
    ! relative_humidity
    if (relative_humidity .eq. 0) then
       call ssh_compute_relative_humidity(humidity, Temperature, &
            Pressure, relative_humidity) 
       Relative_Humidity = DMIN1(DMAX1(Relative_Humidity, Threshold_RH_inf), Threshold_RH_sup)
    end if

    ! for gas phase chemistry
    ! This index is modified by adding one later
    ! See ispeclost in Chemistry/common/hetrxn.f
    do k=1,N_gas
      if(species_name(k) == 'HO2') then 
         heterogeneous_reaction_index(1)= k-1
      endif
      if(species_name(k) == 'NO2') then 
         heterogeneous_reaction_index(2)= k-1
      endif
      if(species_name(k) == 'NO3') then 
         heterogeneous_reaction_index(3)= k-1 
      endif
      if(species_name(k) == 'N2O5') then 
         heterogeneous_reaction_index(4)= k-1 
      endif
    enddo

    ns_source = 1
    if(.not.allocated(source_index)) allocate(source_index(ns_source))
    source_index = [1]
    if(.not.allocated(source)) allocate(source(ns_source))
    source = 0.d0
    if(.not.allocated(conversionfactor)) allocate(conversionfactor(n_gas))
    if(.not.allocated(conversionfactorjacobian)) allocate(conversionfactorjacobian(n_gas, n_gas))

    do j = 1, n_gas
       conversionfactor(j) = Navog * 1.d-12 / molecular_weight(j)
       do k = 1, n_gas
          conversionfactorjacobian(j, k) = molecular_weight(j) / &
               molecular_weight(k)
       enddo
    enddo

    ! initialise fraction discretization 

    if(.not.allocated(discretization_mass)) allocate(discretization_mass(N_sizebin+1))
    discretization_mass = 0.d0
    do k = 1,N_sizebin+1
       discretization_mass(k) = fixed_density * pi * diam_bound(k) ** 3 / 6.d0
    end do

    if(.not.allocated(lwc_nsize)) allocate(lwc_nsize(n_size)) ! SOAP !
    lwc_nsize = 0.d0
    if(.not.allocated(ionic_nsize)) allocate(ionic_nsize(n_size)) ! SOAP !
    ionic_nsize = 0.d0
    if(.not.allocated(proton_nsize)) allocate(proton_nsize(n_size))  ! SOAP !
    proton_nsize = 0.d0
    if(.not.allocated(chp_nsize)) allocate(chp_nsize(n_size))  ! SOAP !
    chp_nsize = 0.d0
    if(.not.allocated(liquid_nsize)) allocate(liquid_nsize(12,n_size))   ! SOAP !
    liquid_nsize = 0.d0
    if(.not.allocated(surface_equilibrium_conc_nsize)) allocate(surface_equilibrium_conc_nsize(N_size,N_aerosol)) 
    surface_equilibrium_conc_nsize = 0.d0

    if(.not.allocated(concentration_inti)) allocate(concentration_inti(N_size,N_inside_aer)) ! ModuleCondensation !
    concentration_inti=0.d0

    if(.not.allocated(total_mass)) allocate(total_mass(N_aerosol))
    total_mass=0.d0    

    if(.not.allocated(concentration_number)) allocate(concentration_number(N_size)) 
    concentration_number=0.d0

    if(.not.allocated(concentration_number_tmp)) allocate(concentration_number_tmp(N_size)) ! ModuleAdaptstep
    concentration_number_tmp=0.d0

    if(.not.allocated(concentration_mass)) allocate(concentration_mass(N_size,N_aerosol_layers)) 
    concentration_mass=0.d0

    if(.not.allocated(concentration_mass_tmp)) allocate(concentration_mass_tmp(N_size,N_aerosol_layers)) ! ModuleAdaptstep
    concentration_mass_tmp=0.d0

    if(.not.allocated(emission_rate)) allocate(emission_rate(N_size,N_aerosol_layers)) ! ModuleEmission
    emission_rate=0.d0
    if(.not.allocated(emission_num_rate)) allocate(emission_num_rate(N_size))       ! ModuleEmission
    emission_num_rate=0.d0

    if(.not.allocated(cell_mass)) allocate(cell_mass(N_size))
    cell_mass = 0.d0

    if(.not.allocated(cell_diam_av)) allocate(cell_diam_av(N_size)) ! average diameter of each grid cell
    cell_diam_av=0.d0

    if(.not.allocated(total_aero_mass)) allocate(total_aero_mass(N_aerosol)) ! ModulePhysicalbalance, Emission, Adapstep
    total_aero_mass=0.d0

    if(.not.allocated(mass_total_grid)) allocate(mass_total_grid(N_size))  ! ModulePhysicalbalance
    mass_total_grid=0.d0

    ! wet !
    if(.not.allocated(wet_mass)) allocate(wet_mass(N_size))
    wet_mass=0.d0

    if(.not.allocated(wet_diameter)) allocate(wet_diameter(N_size))
    wet_diameter=0.d0

    if(.not.allocated(wet_volume)) allocate(wet_volume(N_size))
    wet_volume=0.d0

    if(.not.allocated(bin_mass)) allocate(bin_mass(N_sizebin))  ! ModulePhysicalbalance
    bin_mass = 0.d0
    if(.not.allocated(bin_number)) allocate(bin_number(N_sizebin)) ! ModulePhysicalbalance
    bin_number = 0.d0

    ! ModuleBulkequibrium ModuleAdaptstep ModuleCondensation ModuleCongregation ModuleThermodynamics
    if(.not.allocated(ce_kernal_coef)) allocate(ce_kernal_coef(N_size,N_aerosol)) 
    ce_kernal_coef = 0.d0

    if(.not.allocated(frac_grid)) allocate(frac_grid(N_size,N_groups)) ! ModulePhysicalbalance ModuleRedistribution
    frac_grid = 0.d0

    if(.not.allocated(dqdt)) allocate(dqdt(N_size,N_aerosol_layers)) ! ModuleCongregation ModuleCondensation ModuleAdaptstep
    dqdt = 0.d0

    if (ssh_standalone) write(*,*) "=====================finish initialising parameters==================="
    if (ssh_logger) write(logfile,*) "=====================finish initialising parameters==================="
  end subroutine ssh_init_parameters


  subroutine ssh_discretization()
    !------------------------------------------------------------------------
    !
    !     -- DESCRIPTION
    !     This subroutine ssh_automatically computes particle compositions.
    !     Information of particle compositions is saved under "INIT/fractions.txt"
    !
    !------------------------------------------------------------------------
    !
    !     -- INPUT VARIABLES
    !
    !------------------------------------------------------------------------   
    implicit none

    double precision:: sumfrac
    integer,dimension(:), allocatable:: counter
    integer:: Nubvaild
    integer:: i,s,j,k,s1,rankk

    open(unit = 11, file = particles_composition_file)
    Nubvaild=0
    rankk=0
    allocate(counter(N_groups-1))

    !calculate the maximum fraction combinations
    do i = 1, N_frac
       do s = 1, N_groups-1
          counter(s)=1!initial the counter
       enddo
       if(N_groups.gt.2) then
          ! when the index counter of second species reaches its top, move to the N_aerosol fraction bin of first species
          do while(counter(2).le.N_frac)!Traversal all the possible combination
             sumfrac=frac_bound(i)!(i+1)!take the base fraction bounds of current bin of first species
             do s =2, N_groups-1
                rankk=rankk+1
                j=counter(s)!the fraction bin index for species s
                sumfrac=sumfrac+frac_bound(j)!(j+1)!calculate one possible combination
             enddo
             if (sumfrac.lt.1.d0) then
                Nubvaild=Nubvaild+1!get one possible combination
                write(unit=11,FMT=*) frac_bound(i), frac_bound(i+1)!for first species
                do s=2, N_groups-1
                   j=counter(s)!write down possible combinations
                   write(unit=11,FMT=*) frac_bound(j), frac_bound(j+1)
                enddo
                write(unit=11,FMT=*) frac_bound(1), frac_bound(N_frac+1)!for last species
             endif
             !when the second last species hasn't reaches its top,
             if(counter(N_groups-1).le.N_frac) then
                counter(N_groups-1)=counter(N_groups-1)+1!move the index of second last species
             endif
             !!optimized rank method
             do s=3,N_groups-1!check every neighbor counter, form back to forward
                j=N_groups+2-s
                sumfrac=frac_bound(counter(j-1))+frac_bound(counter(j))
                if(sumfrac.ge.1.d0) then
                   do s1=j,N_groups-1
                      counter(s1)=1
                   enddo
                   counter(j-1)=counter(j-1)+1
                endif
             enddo
          enddo
       else!in case of only two/one species
          Nubvaild=N_frac
          do s=1, N_groups-1
             j=counter(s)
             write(unit=11,FMT=*) frac_bound(i), frac_bound(i+1)
          enddo
          write(unit=11,FMT=*) frac_bound(1), frac_bound(N_frac+1)
       endif
    enddo
    CLOSE(11)

    N_fracmax=Nubvaild!get the N_fracmax
    allocate(discretization_composition(N_fracmax, N_groups, 2))

    n_size = 0
    do k = 1,N_sizebin!k is the number of bins
       N_fracbin(k) = N_fracmax
       N_size = N_size + N_fracbin(k)
    enddo

    open(unit = 11, file = particles_composition_file, status = "old")

    do i = 1, N_fracmax
       do s = 1, N_groups
          read(11,*)discretization_composition(i, s, 1), discretization_composition(i, s, 2)
       enddo
    enddo

    CLOSE(11)
    deallocate(counter)
  end subroutine ssh_discretization

  subroutine ssh_init_coag()
    
    integer:: tag_file,i

    !------------------------------------------------------------------------
    !
    !     -- DESCRIPTION
    !     This subroutine initializes coagulation coefficients.
    !
    !------------------------------------------------------------------------
    !
    !     -- INPUT VARIABLES
    !
    !------------------------------------------------------------------------   

    
    if (with_coag.eq.1) then !if coagulation
       if(.not. allocated(kernel_coagulation)) allocate(kernel_coagulation(N_size,N_size))
       kernel_coagulation = 0.d0

       if (i_compute_repart == 0 .or. i_write_repart == 1) then
          do i=1,len(trim(Coefficient_file))!judge the input files
             if(Coefficient_file(i:i)==".")then
                if(Coefficient_file(i+1:i+2)=="nc".or.Coefficient_file(i+1:i+2)=="NC") then
                   tag_file=1
                elseif (Coefficient_file(i+1:i+3)=="bin".or.Coefficient_file(i+1:i+3)=="BIN") then
                   tag_file=0
                elseif (Coefficient_file(i+1:i+3)=="txt".or.Coefficient_file(i+1:i+3)=="TXT") then
                   tag_file=2
                else
                   if (ssh_standalone) write(*,*) "Unsupported input coefficient file type for coagulation."
                   if (ssh_logger) write(logfile,*) "Unsupported input coefficient file type for coagulation."
                   if (i_compute_repart == 0) i_compute_repart = 1
                   if (i_write_repart == 1) i_write_repart = 0
                endif
             endif
          enddo
          if (ssh_standalone) write(*,*) 'Coefficient Repartition Database:',Coefficient_file
          if (ssh_logger) write(logfile,*) 'Coefficient Repartition Database:',Coefficient_file
       endif

       ! Subroutines are defined in ModuleCoefficientRepartition
       if (i_compute_repart == 0) then
          call ssh_ReadCoefficientRepartition(Coefficient_file, tag_file)
       else if (i_compute_repart == 1) then
          if (.not. allocated(repartition_coefficient)) call ssh_ComputeCoefficientRepartition()
       else
          if (ssh_standalone) write(*,*) "Coefficient for coagulation must be read or computed."
          if (ssh_logger) write(logfile,*) "Coefficient for coagulation must be read or computed."
          stop
       endif
       if (i_write_repart == 1) call ssh_WriteCoefficientRepartition(Coefficient_file, tag_file)

       ! Check the quality of coagulation repartition coefficients
       call ssh_check_repart_coeff()

    endif
    
  end subroutine ssh_init_coag
  
  subroutine ssh_Init_distributions()
    !------------------------------------------------------------------------
    !
    !     -- DESCRIPTION
    !     This subroutine ssh_initialize mass and number concentration based
    !     on initialization methods indicate within configuration files.
    !
    !------------------------------------------------------------------------
    !
    !     -- INPUT VARIABLES
    !
    !------------------------------------------------------------------------   
    IMPLICIT NONE
    integer :: j,k,s, i,j1,j2,Czero,f,n,jesp,g,lay
    double precision::mass_frac(N_sizebin),binx_mass(N_sizebin),binx_emis(N_sizebin)
    double precision::speciesfrac(N_aerosol,N_sizebin),tmp

    speciesfrac=1.d0  !adjust mass distribution for different species
    concentration_mass = 0.d0
    concentration_number = 0.d0

    ! calculate the referenced initial mass and number distribution of each bin
    if(tag_init.eq.0) then 
       tmp = 0.d0
       ! Compute the percentage of mass distribution
       do k=1,N_sizebin
          do j = 1, N_aerosol
             tmp = tmp + init_bin_mass(k,j)
	  end do
          mass_frac(k) = tmp / aero_total_mass
	  tmp = 0.d0
       enddo
    endif

    binx_mass=0.d0
    do k=1,N_sizebin
       do s=1, N_aerosol
          binx_mass(k) = binx_mass(k) + init_bin_mass(k,s)
       enddo
    enddo

    binx_emis = 0.d0
    if (tag_emis .ne. 0) then
       do k=1,N_sizebin
          do s=1, N_aerosol
             binx_emis(k)=binx_emis(k) + emis_bin_mass(k,s)
          enddo
       enddo
    end if

    if (N_frac.gt.1) then
       !case of external mixing
       if (ssh_standalone) write(*,*) 'External mixing...'
       if (ssh_logger) write(logfile,*) 'External mixing...'
       do j=1,N_size                     ! j : index of cells
          do s=1, N_aerosol             ! s : index of species
             g=Index_groups(s)           ! g : group number
             k= concentration_index(j, 1)! k : index of size bins
             f= concentration_index(j, 2)! f : index of frac combinations
             if(g.lt.N_groups) then
                if(discretization_composition(f, g ,2).eq.1.d0) then!checking top fraction (particles with single species)

                   if(tag_emis .ne. 0 ) then
                      if(emis_bin_mass(k,s).gt.0.d0) then
                         if(s.LE.N_nonorganics) then
                            emission_rate(j,s)=emis_bin_mass(k,s)
                         else
                            do lay=1,nlayer
                               jesp = index_species(s,lay)
                               emission_rate(j,jesp) = emis_bin_mass(k,s)* Vlayer(lay)
                            enddo
                         endif
                      endif
                   endif

                   !                  if(tag_external.eq.1) then !!in case of external mixed initial condition
                   concentration_number(j) = concentration_number(j)+init_bin_number(k)&
                        *init_bin_mass(k,s)/binx_mass(k)
                   if(s.LE.N_nonorganics) then
                      concentration_mass(j,s) = init_bin_mass(k,s) 
                   else
                      do lay=1,nlayer
                         jesp = index_species(s,lay)
                         concentration_mass(j,jesp) = init_bin_mass(k,s) * Vlayer(lay)
                      enddo
                      if(i_hydrophilic==1) then
                        jesp = index_species(s,nlayer+i_hydrophilic)
                        concentration_mass(j,jesp) = 0.d0 ! Initialise the aqueous phase to zero
                      endif
                   endif
                   !                  endif
                   if(tag_emis .ne. 0 .and. with_emis_num.eq.1) then !!in case of external mixed initial condition
		      emission_num_rate(j) = emis_bin_number(k) &
                           *emis_bin_mass(k,s)/binx_emis(k)   
                   endif
                   !endif
                endif
             else
                Czero=0
                do g=1,N_groups-1
                   if(discretization_composition(f, g,1).gt.0.d0) then
                      Czero=1
                   endif
                enddo
                if(Czero.eq.0.d0) then
                   if(tag_emis .ne. 0) then
                      if(emis_bin_mass(k,s).gt.0.d0) then
                         if(s.LE.N_nonorganics) then
                            emission_rate(j,s)=emis_bin_mass(k,s)
                         else
                            do lay=1,nlayer
                               jesp = index_species(s,lay)
                               emission_rate(j,jesp)=emis_bin_mass(k,s)* Vlayer(lay)
                            enddo
                         endif
                      endif
                   endif
                   !                   if(tag_external.eq.1) then !!in case of external mixed initial condition
                   concentration_number(j) = concentration_number(j)+init_bin_number(k)&
                        *init_bin_mass(k,s)/binx_mass(k)
                   if(s.LE.N_nonorganics) then
                      concentration_mass(j,s) = init_bin_mass(k,s)
                   else
                      do lay=1,nlayer
                         jesp = index_species(s,lay)
                         concentration_mass(j,jesp) = init_bin_mass(k,s)* Vlayer(lay)
                      enddo
                      if(i_hydrophilic==1) then
                        jesp = index_species(s,nlayer+i_hydrophilic)
                        concentration_mass(j,jesp) = 0.d0 ! Initialise the aqueous phase to 0
                      endif
                   endif
                   !                   endif
                   if(tag_emis .ne. 0 .and.with_emis_num.eq.1) then !!in case of external mixed initial condition
                      emission_num_rate(j) = emis_bin_number(k) &
                           *emis_bin_mass(k,s)/binx_emis(k)  
                   endif
		   !endif
                endif
             endif
	  enddo
       enddo
    else ! default & N_frac.eq.1
       !case of internal mixing N_size=N_sizebin
       if (ssh_standalone) write(*,*) "Internal mixing..."
       if (ssh_logger) write(logfile,*) "Internal mixing..."
       do k=1,N_sizebin
          do s=1,N_aerosol
             if(s.LE.N_nonorganics) then
                if(tag_emis .ne. 0) then ! if with emission
                   emission_rate(k,s) = emis_bin_mass(k,s)
                endif
                concentration_mass(k,s) = init_bin_mass(k,s)
             else
                do lay=1,nlayer
                   jesp = index_species(s,lay)
                   if(tag_emis .ne. 0) then ! if with emission
                      emission_rate(k,jesp) = emis_bin_mass(k,s)* Vlayer(lay)
                   endif
                   concentration_mass(k,jesp) = init_bin_mass(k,s)* Vlayer(lay)
                enddo
                if(i_hydrophilic==1) then
                  jesp = index_species(s,nlayer+i_hydrophilic)
                  concentration_mass(k,jesp) = 0.d0 ! Initialise the aqueous phase to 0
                endif
             endif
          enddo
          ! number
          if(tag_emis .ne. 0) emission_num_rate(k) = emis_bin_number(k)
          concentration_number(k) = init_bin_number(k)
       enddo

    endif

    ! ! incase of internal mixed initial condition ! !
    !   if(tag_external.eq.0.and.N_frac.gt.1) then
    !      do k=1,N_sizebin
    !         j=concentration_index_iv(k,1)
    !         concentration_number(j) = init_bin_number(k)
    !         do s=1, N_nonorganics!index of non organic species
    !            concentration_mass(j,s) = init_bin_mass(k,s) 
    !         enddo
    !         do s=N_nonorganics+1, N_aerosol !index of organics and water
    !            do lay=1,nlayer
    !               jesp = index_species(s,lay)
    !               concentration_mass(j,jesp) = init_bin_mass(k,s) * Vlayer(lay)
    !            enddo
    !         enddo
    !      enddo
    !   endif

    if(tag_external.ne.0) then!incase of internal mixed initial condition
       if(N_frac.gt.1) then
          call ssh_redistribution_fraction()!fraction redistribution
       endif
    endif

    ! Initialize the concentrations of aerosol with gas-phase percursors
    ! total_number = 0.0
    total_aero_mass = 0.0
    total_mass = 0.0
    do s = 1, N_aerosol
       if (aerosol_species_interact(s) .gt. 0) then
          concentration_gas(s) = concentration_gas_all(aerosol_species_interact(s)) !µg/m3
       end if
    enddo
    do s = 1, N_aerosol_layers
       jesp=List_species(s)
       do i=1,N_size
          total_aero_mass(jesp)=total_aero_mass(jesp)+concentration_mass(i,s)
       enddo
    enddo
    do s = 1, N_aerosol
       total_mass(s)=total_mass(s) + total_aero_mass(s) + concentration_gas(s)
    end do

    ! done with initialising concentrations 
    ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

    ! initialise density ! after concentration_mass and number concentrations if necessary

!    if (with_fixed_density.NE.1) then
!       call ssh_compute_all_density() ! density_aer_bin(N_size)  density_aer_size(N_sizebin)
!       if (ssh_standalone) write(*,*)"Density is auto-generated."
!       if (ssh_logger) write(logfile,*)"Density is auto-generated."
!    end if

    if (with_init_num.NE.1) then
       if (with_fixed_density.NE.1) then
         call ssh_compute_all_density() ! density_aer_bin(N_size)  density_aer_size(N_sizebin)
         if (ssh_standalone) write(*,*)"Density is auto-generated."
         if (ssh_logger) write(logfile,*)"Density is auto-generated."
       end if
       ! need size_diam_av and conc._mass
       call ssh_compute_number()  ! only for initialisation
       if (ssh_standalone) write(*,*)"Initial PM number concentration is auto-generated."
       if (ssh_logger) write(logfile,*)"Initial PM number concentration is auto-generated."
    end if
    call ssh_compute_average_diameter() ! Compute average_diame for 1st iteration output


    if (ssh_standalone) write(*,*)"=================================finish initial distribution==============================="
    if (ssh_logger) write(logfile,*)"=================================finish initial distribution==============================="

    if (.not. allocated(ratio_water)) allocate(ratio_water(N_size))
    if (.not. allocated(ratio_eqconc)) allocate(ratio_eqconc(4,N_size))
    if (.not. allocated(iter_water)) allocate(iter_water(N_size))
    if (.not. allocated(iter_eqconc)) allocate(iter_eqconc(N_size))  

  end subroutine ssh_Init_distributions

end MODULE lDiscretization
