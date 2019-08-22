!!-----------------------------------------------------------------------
!!     Copyright (C) 2012-2018, ENPC - EDF R&D - INERIS
!!     Author(s): Shupeng Zhu
!!
!!     This file is part of the Size Composition Resolved Aerosol Model (SCRAM), a
!!     component of the SSH-aerosol model.
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
!!    This module read configuration file and initialize all global variables
!!-----------------------------------------------------------------------
module aInitialization
    use iso_c_binding
!    use dPhysicalbalance
    implicit none
    INCLUDE 'CONST.INC'
    INCLUDE 'CONST_A.INC'
    include 'pointer.inc'
    include 'smw.inc'    ! molecular_weight_solid
    include 'solmd.inc'  ! mass_density_solid
    include 'imw.inc'    ! molecular_weight_inside
    include 'CONST_B.inc'

    !!part 1: parameters of system dimension
    Integer :: N_gas        !Complete gas species number
    integer :: N_size       !Total number of size and composition sections
    integer :: N_groups     !Number of groups
    integer :: N_fracmax    !Maximum number of composition sections per size section
    integer :: N_aerosol    !Number of aerosol species
    integer :: N_sizebin    !Number of  size sections
    integer :: N_frac       !Number of fraction of each species
    integer :: N_reaction   !Number of gas-phase reactions
    integer :: N_photolysis !Number of photolyses
    integer :: Nt           !Number of iteration
    integer :: Nmc          
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: N_organics   !Number of organics aerosol species
    integer :: N_inorganic  !Number of inorganic aerosol species
    integer :: N_inert      !number of inert aerosol species
    integer :: N_liquid     !Number of liquid internal species
    integer :: N_solid      !Number of solid internal species
    integer :: N_inside_aer !Number of internal species
    integer :: N_hydrophilic!Number of hydrophilic organics aerosol species
    !parameter (N_aerosol = N_organics + N_inorganic + N_inert + 1)
    parameter (N_organics=27,N_inorganic=5,N_inert=2,N_liquid=12)
    parameter (N_solid=9,N_inside_aer=21)
    parameter (N_hydrophilic=9)

    !!part 2: parameters of system options
  
    integer :: with_fixed_density!IDENS
    integer :: tag_init		! 0 internally mixed; 1 mixing_state resolved
    integer :: with_init_num	! 0 estimated from mass and diameter; 1 number conc. for each bin is read
    integer :: wet_diam_estimation	! 0 = isorropia ?
    integer :: tag_dbd    ! Method for defining particle size bounds (0 auto generated, 1 read)
    integer :: tag_emis	     ! 0 Without emissions 1 with internally-mixed emissions 2 with externally-mixed emissions
    integer :: with_emis_num ! 0 estimated from mass and diameter; 1 number conc. for each bin is read

    integer :: tag_external  ! 0 for internally mixed, 1 for mixing-state resolved
    integer :: kind_composition  ! 1 for auto discretization and 0 for manual discretization

    integer :: tag_chem
    integer :: with_photolysis
    integer :: with_heterogeneous  !Tag of heterogeneous reaction 

    integer :: with_adaptive       !Tag of adaptive time step for chemistry 1 if adaptive time step.
    double precision :: adaptive_time_step_tolerance !Relative tolerance for deciding if the time step is kept
    double precision :: min_adaptive_time_step       !Minimum time step
    double precision :: DTAEROMIN !Minimum time step for aerosol dynamics
    double precision :: epser !  Relative error for time step adjustment
    double precision :: epser_soap !  Relative difference of ros2 in SOAP
    integer :: dynamic_solver = 1 !KDSLV Tag type of solver
    integer :: sulfate_computation = 0 !ISULFCOND tag of sulfate condensation method
    integer :: redistribution_method !tag of redistribution method
    integer :: with_coag   !Tag gCoagulation
    integer :: i_compute_repart ! 0 if repartition coeff are read
    integer :: with_cond   !Tag fCondensation
    integer :: with_nucl   !Tag nucleation
    Integer :: nucl_model  !ITERN !1= Ternary, 0= binary
    integer :: ICUT        !ICUT
    double precision :: Cut_dim  !cuting diameter between equi/dynamic inorganic
    integer :: ISOAPDYN    ! organic equilibrium  = 0 or dynamic = 1
    integer :: with_oligomerization!IOLIGO
    integer :: output_type
    integer :: splitting

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer :: aqueous_module!ICLD
    Integer :: with_incloud_scav!IINCLD
    integer :: with_kelvin_effect!IKELV
    integer :: section_pass
    double precision :: tequilibrium ! time under which equilibrium is assumed

    ! ! part 3: System pointers
    Integer :: E1,E2,G1,G2 !Mark the begin and end of dynamic aerosol (except EH2O)
    ! Number of different species group
    Integer, dimension(:), allocatable :: isorropia_species
    Integer, dimension(:), allocatable :: aec_species
    Integer, dimension(:), allocatable :: pankow_species
    Integer, dimension(:), allocatable :: poa_species
    Integer :: nesp, nesp_isorropia, nesp_aec, nesp_pankow, nesp_pom, nesp_eq_org
    parameter (nesp_isorropia=5,nesp_aec=19,nesp_pankow=1,nesp_pom=6, nesp_eq_org=26)
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Integer :: ENa,ESO4,ENH4,ENO3,ECl,EMD,EBC,EH2O!inorganic pointers
    Integer :: ictmNH3,ictmHNO3,ictmHCl,ictmSO2,ictmH2O2,ictmHCHO,ictmHNO2
    Integer :: ictmO3,ictmOH,ictmHO2,ictmNO3,ictmNO,ictmNO2,ictmPAN,ictmH2SO4
    ! pointers of cloud species.



    !!part 4: System state parameters    
    ! time setting
    double precision :: final_time,dt,time_emis,delta_t, initial_time  

    double precision :: Temperature,Relative_Humidity,Pressure,Humidity, pH, pressure_sat
    double precision :: longitude, latitude
    double precision :: attenuation
    double precision :: fixed_density
    double precision :: lwc_cloud_threshold

    integer :: tag_coag,tag_cond,tag_nucl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    double precision :: viscosity!Dynamic viscosity ([kg/m/s]).
    double precision :: air_free_mean_path
    double precision :: total_water!total mass of water
    double precision :: total_IH!total mass of H+
    double precision :: total_PH!overall PH value
    double precision :: n_emis
    double precision :: m_emis
    double precision :: aero_total_mass 
    Double precision,dimension(:), allocatable :: total_mass!, total_mass_old!total mass of each species
    Double precision,dimension(:), allocatable :: discretization_mass
    Double precision :: p_fact,k_fact
    Double precision :: DQLIMIT

    !!part5: dimension data array    
    integer, dimension(:), allocatable :: Index_groups	!index of which group the species belongs to
    integer, dimension(:), allocatable :: List_species	!read species defined in cfg files
    Integer, dimension(:), allocatable :: aerosol_species_interact
    integer, dimension(:), allocatable :: N_fracbin	!vector of number of composition sections for each section

    Double precision,dimension(:), allocatable :: photolysis
    integer, dimension(:), allocatable :: photolysis_reaction_index

    Double precision,dimension(:), allocatable :: density_aer_bin 	!density of each grid bins
    Double precision,dimension(:), allocatable :: density_aer_size 	!density of each size section
    Double precision , dimension(:), allocatable :: rho_wet_cell

    Double precision,dimension(:), allocatable :: diam_bound	! DBF diameter bounds of each size section
    double precision,dimension(:), allocatable :: diam_input
    double precision,dimension(:), allocatable :: frac_bound
    double precision,dimension(:), allocatable :: frac_input


    Double precision,dimension(:), allocatable :: size_diam_av	!DSF average diameter of each size section
    Double precision,dimension(:), allocatable :: size_mass_av	!MSF average mass of each size section
    !Double precision,dimension(:), allocatable :: size_log_av	!XSF
    Double precision,dimension(:), allocatable :: cell_diam_av	!DSF average diameter of each grid cell
    Double precision,dimension(:), allocatable :: cell_mass_av	!MSF average mass of each grid cell
    Double precision,dimension(:), allocatable :: cell_log_av	!XSF

    DOUBLE PRECISION, dimension(:), allocatable :: concentration_gas_all
    Double precision,dimension(:), allocatable :: concentration_gas	! gas concentration of each species
    integer, dimension(:,:), allocatable :: concentration_index !matrix from grid index to size and composition index
    integer, dimension(:,:), allocatable :: concentration_index_iv !matrix from size and composition to grid index
    Double precision,dimension(:), allocatable :: concentration_number	!number concentration of each grid cell
    double precision , dimension(:,:), allocatable :: concentration_mass

    double precision, dimension(:), allocatable    :: gas_emis !storing Gas consentration (emission) micm^3cm^-3 
    double precision,dimension(:,:), allocatable   :: init_bin_mass
    double precision,dimension(:), allocatable     :: init_mass
    double precision,dimension(:,:), allocatable   :: emis_bin_mass
    double precision,dimension(:), allocatable     :: init_bin_number 
    double precision,dimension(:), allocatable     :: emis_bin_number

    double precision , dimension(:,:), allocatable :: emission_rate
    double precision , dimension(:), allocatable   :: emission_num_rate

    Double precision,dimension(:), allocatable :: wet_diameter	!Aerosol wet diameter (\B5m). of each grid cell
    Double precision,dimension(:), allocatable :: wet_mass	!Aerosol wet mass (\B5g). of each grid cell
    Double precision,dimension(:), allocatable :: wet_volume	!Aerosol wet volume (\B5m^3). of each grid cell
    double precision , dimension(:,:,:), allocatable :: discretization_composition! multi-array storing discretization of composition

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Double precision,dimension(:), allocatable :: mass_bound! MBF
    Double precision,dimension(:), allocatable :: log_bound!XBF
    Double precision,dimension(:), allocatable :: total_bin_mass!total mass of each size section
    Double precision,dimension(:), allocatable :: size_sect!HSF log size of each section


    Double precision,dimension(:), allocatable :: mass_total_grid!total mass of each grid cell
    Double precision,dimension(:), allocatable :: total_aero_mass!total aerosol mass of each species
    Double precision,dimension(:), allocatable :: bin_mass!mass concentration of each size section
    Double precision,dimension(:), allocatable :: bin_number!number concentration of each size section
    Double precision,dimension(:), allocatable :: concentration_number_tmp!first order approximation of number

    Double precision , dimension(:), allocatable :: cell_mass
    double precision, dimension(:), allocatable :: number_init  !for each sizebin
    !double precision, dimension(:), allocatable :: mass_init    !for each sizebin


    double precision, dimension(:), allocatable :: per_mass_init!initial percentage of each species within aerosol



    double precision,dimension(:), allocatable:: gas_mass_init


    double precision , dimension(:,:), allocatable :: kernel_coagulation
    double precision , dimension(:,:), allocatable :: ce_kernal_coef!c/e kernal
    double precision , dimension(:,:), allocatable :: Kelvin_effect_ext!kelvin effect
    double precision , dimension(:,:), allocatable :: frac_grid !excat fraction of each species in each grid

    double precision , dimension(:,:), allocatable :: concentration_mass_tmp!first order apporximation
    double precision , dimension(:,:), allocatable :: concentration_inti!internal inorganic aerosol concentration ([�g.m-3]).
    double precision , dimension(:,:), allocatable :: dqdt



    !! part 7: basic physical and chemical parameters
    !double precision :: SMD(SNaNO3:SLC) = ()!molar weight of internal solids species
    !double precision :: IMW(N_liquid)!molar weight of inorganic species in aqueous_phase
    !double precision :: SMW(SNaNO3:SLC)!molar weight of solids
    double precision ,dimension(:), allocatable :: accomodation_coefficient
    double precision ,dimension(:), allocatable :: surface_tension
    double precision ,dimension(:), allocatable :: molecular_weight_aer! (\B5g/mol)
    double precision ,dimension(:), allocatable :: molecular_diameter
    double precision ,dimension(:), allocatable :: collision_factor_aer
    double precision ,dimension(:), allocatable :: mass_density!(\B5g/m3) liquid mass density
    double precision ,dimension(:), allocatable :: quadratic_speed! (m.s-1)
    double precision ,dimension(:), allocatable :: diffusion_coef! (m2.s-1)
    double precision ,dimension(:), allocatable :: soa_sat_conc! (\B5g.m-3)
    double precision ,dimension(:), allocatable :: soa_part_coef!(m3/microg)
    double precision ,dimension(:), allocatable :: molecular_weight! (\B5g/mol) gas=phase

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, dimension(:), allocatable :: Ncoefficient, index_first, index_second
    double precision, dimension(:), allocatable :: coefficient
    integer :: coef_size
    double precision :: surface_tension_inorg, surface_tension_aq, surface_tension_org

    double precision :: dorg
    integer :: coupled_phases
    integer :: nlayer
    integer :: activity_model
    double precision, dimension(:), allocatable :: lwc_nsize, &
         ionic_nsize, proton_nsize, chp_nsize
    double precision, dimension(:,:), allocatable :: liquid_nsize
    
    !! part 8: divers parameters (species, I/O)
    character (len=80) :: Coefficient_file ! repartition coefficient file
    character (len=80) :: init_aero_conc_mass_file ! File for aeroslos initial mass concentrations
    character (len=80) :: init_aero_conc_num_file ! File for aerosols initial number concentrations
    character (len=80) :: init_gas_conc_file ! File for gas-phase initial conc.
    character (len=80) :: species_list_file ! File for species list.
    character (len=80) :: aerosol_species_list_file ! File for species list.
    character (len=80) :: namelist_species ! Namelist file for species list.
    character (len=80) :: emis_gas_file
    character (len=80) :: emis_aero_mass_file
    character (len=80) :: emis_aero_num_file
    character (len=80) :: mineral_dust(10), black_carbon(10)
    character (len=80) :: isorropia_species_name(10)
    character (len=80) :: aec_species_name(30)
    character (len=80) :: pankow_species_name(6)
    character (len=80) :: poa_species_name(10)
    character (len=80) :: PSO4
    character (len=10) :: precursor
    character (len=10), dimension(:), allocatable :: species_name
    character (len=20), dimension(:), allocatable :: aerosol_species_name
    integer :: spec_name_len
    character (len=10), dimension(:), allocatable :: emis_gas_species_name
    character (len=10), dimension(:), allocatable :: emis_aer_species_name
    character (len=100) :: output_directory, output_dir2


    !!part 6: used in ssh-aerosol.f90 chem()
    integer :: ns_source
    integer, dimension(:), allocatable :: source_index
    double precision, dimension(:), allocatable :: source  
    double precision, dimension(:), allocatable :: conversionfactor
    double precision, dimension(:,:), allocatable :: conversionfactorjacobian
    ! Array of chemical volumic emissions at final time ([\mu.g/m^3/s]).
    integer :: heterogeneous_reaction_index(4)
    integer :: ind_jbiper, ind_kbiper   
    ! To take into account BiPER degradation inside the particle

    !! part 7: used in coupling with external tools (ModuleSaturne)
    !
    ! This flag defines if SSH-aerosol is running standalone or not
    !   true if it is running standalone, default
    !   false if it is running with Code_Saturne
    !
    logical(kind=c_bool), save :: ssh_standalone = .true.
    !
    ! This flag defines if SSH-aerosol is logging information to a file
    !   true if it is writing to a file
    !   false otherwise, default
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!                                                                    !!!!!
    !!!!!                          IMPORTANT                                 !!!!!
    !!!!!                                                                    !!!!!
    !!!!! change the flag to true with the subroutine set_logger             !!!!!
    !!!!!                                                                    !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    logical(kind=c_bool), save :: ssh_logger = .false.
    !
    ! This is the name of the log file
    !
    character(len=80), save :: ssh_logger_file = "ssh-aerosol.log"
    !
    ! This is the unit used to write to the log file
    !
    integer, parameter :: logfile = 99
     
 contains


   subroutine read_namelist(namelist_file)
     

! =============================================================
!
!    read simulation setting from file
!
!    input : namelist.ssh file 
! ============================================================= 

     implicit none
     integer :: i,ierr,tag_file
     character (len=40), intent(in) :: namelist_file


     ! namelists to read namelist.ssh file 

     namelist /setup_meteo/ latitude, longitude, Temperature, Pressure,&
                             Humidity, Relative_Humidity

     namelist /setup_time/ initial_time, final_time, delta_t,time_emis

     namelist /initial_condition/ with_init_num, tag_init, tag_dbd, N_sizebin,&
                                  wet_diam_estimation, init_gas_conc_file,&
                                  init_aero_conc_mass_file, init_aero_conc_num_file
                                  
     namelist /initial_diam_distribution/ diam_input

     namelist /emissions/ tag_emis, with_emis_num, emis_gas_file, &
				emis_aero_mass_file, emis_aero_num_file

     namelist /mixing_state/ tag_external, N_groups, N_frac, kind_composition
			     
     namelist /fraction_distribution/ frac_input

     namelist /gas_phase_species/ species_list_file

     namelist /aerosol_species/ aerosol_species_list_file, mineral_dust,&
 				black_carbon, isorropia_species_name,&
				aec_species_name,pankow_species_name,&
         			poa_species_name

     namelist /physic_gas_chemistry/ tag_chem, attenuation, with_photolysis, &
					with_heterogeneous, with_adaptive, &
					adaptive_time_step_tolerance, min_adaptive_time_step

     namelist /physic_particle_numerical_issues/ DTAEROMIN, redistribution_method,&
		  	     with_fixed_density, fixed_density, splitting
		 
     namelist /physic_coagulation/ with_coag, i_compute_repart, Coefficient_file, Nmc

     namelist /physic_condensation/ with_cond, Cut_dim, ISOAPDYN, nlayer,&
          with_kelvin_effect, tequilibrium,&
          dorg, coupled_phases, activity_model, epser, epser_soap

     namelist /physic_nucleation/ with_nucl, nucl_model

     namelist /physic_organic/ with_oligomerization

     namelist /output/ output_directory, output_type

     if (ssh_standalone) write(*,*) "=========================start read namelist.ssh file======================"
     if (ssh_logger) write(logfile,*) "=========================start read namelist.ssh file======================"
     ! read namelist.ssh file !
     open(unit = 10, file = namelist_file, status = "old")
     
     ! meteorological setup
     read(10, nml = setup_meteo, iostat = ierr)
     if (ierr .ne. 0) then
	write(*,*) "setup_meteo data can not be read."
        stop
     else ! default output meteo data to check
	if (temperature - 273.15e0 == 0) then
		write(*,*) 'temperature should higher than 273.15 K - stop'
		stop
	end if
	pressure_sat = 611.2 * dexp(17.67 * (temperature - 273.15) / (temperature - 29.65))
	Relative_Humidity = DMIN1(DMAX1(Relative_Humidity, Threshold_RH_inf), Threshold_RH_sup)
	humidity =  1/(Pressure/(pressure_sat *0.62197* Relative_Humidity)-1)
	if (ssh_standalone) write(*,*) ''
	if (ssh_logger) write(logfile,*) ''
	if (ssh_standalone) write(*,*) '<<<< Meteorological setup >>>>'
	if (ssh_logger) write(logfile,*) '<<<< Meteorological setup >>>>'
        if (ssh_standalone) write(*,*) 'location', latitude, 'N', longitude,'E'
        if (ssh_logger) write(logfile,*) 'location', latitude, 'N', longitude,'E'
        if (ssh_standalone) write(*,*) 'Temperature', Temperature, 'K'
        if (ssh_logger) write(logfile,*) 'Temperature', Temperature, 'K'
        if (ssh_standalone) write(*,*) 'Pressure', Pressure, 'Pa'
        if (ssh_logger) write(logfile,*) 'Pressure', Pressure, 'Pa'
        if (ssh_standalone) write(*,*) 'Relative Humidity', Relative_Humidity
        if (ssh_logger) write(logfile,*) 'Relative Humidity', Relative_Humidity
        if (ssh_standalone) write(*,*) 'Specific Humidity', Humidity
        if (ssh_logger) write(logfile,*) 'Specific Humidity', Humidity
        if (ssh_standalone) write(*,*) 'Cloud attenuation field', attenuation
        if (ssh_logger) write(logfile,*) 'Cloud attenuation field', attenuation
     end if

     ! time setup
     read(10, nml = setup_time, iostat = ierr)
     if (ierr .ne. 0) then
	write(*,*) "setup_time data can not be read."
        stop
     else
         nt = int((final_time-initial_time) / delta_t) 
	 if (ssh_standalone) write(*,*) ''
	 if (ssh_logger) write(logfile,*) ''
	 if (ssh_standalone) write(*,*) '<<<< Simulation time setup >>>>'
	 if (ssh_logger) write(logfile,*) '<<<< Simulation time setup >>>>'
         if (ssh_standalone) write(*,*) 'Begining time (from Jan. 1st)', initial_time, 's'
         if (ssh_logger) write(logfile,*) 'Begining time (from Jan. 1st)', initial_time, 's'
         if (ssh_standalone) write(*,*) 'Simulation Time', final_time,'s'
         if (ssh_logger) write(logfile,*) 'Simulation Time', final_time,'s'
         if (ssh_standalone) write(*,*) 'Initial Time Step', delta_t,'s'
         if (ssh_logger) write(logfile,*) 'Initial Time Step', delta_t,'s'
         if (ssh_standalone) write(*,*) 'Number of iterations:', nt
         if (ssh_logger) write(logfile,*) 'Number of iterations:', nt
     end if
    
     ! initial_condition
     read(10, nml = initial_condition, iostat = ierr)

     if (ierr .ne. 0) then
	write(*,*) "initial_condition data can not be read."
        stop
     else
	if (ssh_standalone) write(*,*) ''
	if (ssh_logger) write(logfile,*) ''
	if (ssh_standalone) write(*,*) '<<<< Inition condition >>>>'
	if (ssh_logger) write(logfile,*) '<<<< Inition condition >>>>'

	if (tag_init == 1) then 
            write(*,*) ' Mixing state resolved aerosol species are provided for initial condition.' 
            write(*,*) ' -- this option is not yet available.'
            stop
        else
	    tag_init = 0 	    ! default tag_init == 0
	    if (ssh_standalone) write(*,*) 'Internally mixed aerosol species are provided for initial condition.'
	    if (ssh_logger) write(logfile,*) 'Internally mixed aerosol species are provided for initial condition.'
        end if 

	if (ssh_standalone) write(*,*) 'Gas-phase conc. input file :', init_gas_conc_file
	if (ssh_logger) write(logfile,*) 'Gas-phase conc. input file :', init_gas_conc_file
	if (ssh_standalone) write(*,*) 'Particle conc. input file :', init_aero_conc_mass_file
	if (ssh_logger) write(logfile,*) 'Particle conc. input file :', init_aero_conc_mass_file
        if (with_init_num .eq. 1) then
	   if (ssh_standalone) write(*,*) 'Aerosol number conc. is read from file :', init_aero_conc_num_file
	   if (ssh_logger) write(logfile,*) 'Aerosol number conc. is read from file :', init_aero_conc_num_file
	else 
	    with_init_num = 0  ! default with_init_num == 0
	    if (ssh_standalone) write(*,*) ' Aerosol number conc. is estimated from mass and diameter.' 
	    if (ssh_logger) write(logfile,*) ' Aerosol number conc. is estimated from mass and diameter.' 
        end if
     end if

	if (ssh_standalone) write(*,*) 'N_sizebin', N_sizebin
	if (ssh_logger) write(logfile,*) 'N_sizebin', N_sizebin
	   allocate(init_bin_number(N_sizebin))
           init_bin_number = 0.d0

	if (tag_dbd == 1)  then  ! initialisation for sizebin
	   allocate(diam_input(N_sizebin+1))
	   diam_input = 0.d0
        else ! default read the boundary of sizebin (tag_dbd = 0)
	   tag_dbd = 0
	   allocate(diam_input(2))
	   diam_input = 0.d0
        end if 

     ! initial_diam_distribution
     read(10, nml = initial_diam_distribution, iostat = ierr)
     if (ierr .ne. 0) then
	write(*,*) "initial_diam_distribution data can not be read."
        stop
     else
	if (tag_dbd == 1 .and. ssh_standalone) write(*,*) 'Diameter bin bounds are read.'
	if (tag_dbd == 1 .and. ssh_logger) write(logfile,*) 'Diameter bin bounds are read.'
	if (tag_dbd == 0 .and. ssh_standalone) write(*,*) &
		'Lower and higher diameter boundary are read for auto-generated diameter bin bounds.'
	if (tag_dbd == 0 .and. ssh_logger) write(logfile,*) &
		'Lower and higher diameter boundary are read for auto-generated diameter bin bounds.'
	if (ssh_standalone) write(*,*) 'diam_input', diam_input
	if (ssh_logger) write(logfile,*) 'diam_input', diam_input
     end if

     ! emissions
     read(10, nml = emissions, iostat = ierr)
     if (ierr .ne. 0) then
	write(*,*) "emissions data can not be read."
        stop
     else
	if (ssh_standalone) write(*,*) ''
	if (ssh_logger) write(logfile,*) ''
	if (ssh_standalone) write(*,*) '<<<< Emissions >>>>'
	if (ssh_logger) write(logfile,*) '<<<< Emissions >>>>'
	if (tag_emis == 2) then
           write(*,*) 'With externally-mixed emissions. -- this option is not yet available.'
	   stop
        else if (tag_emis == 1) then
           if (ssh_standalone) write(*,*) 'With internally-mixed emissions.'
           if (ssh_logger) write(logfile,*) 'With internally-mixed emissions.'
        else
	   tag_emis = 0
           if (ssh_standalone) write(*,*) 'Without emission.'  ! default tag_emis == 0
           if (ssh_logger) write(logfile,*) 'Without emission.'  ! default tag_emis == 0
        end if

           allocate(emis_bin_number(N_sizebin))
           emis_bin_number = 0.d0
	if (tag_emis .ne. 0) then
	    if (ssh_standalone) write(*,*) 'Gas-phase conc. emission file :', emis_gas_file
	    if (ssh_logger) write(logfile,*) 'Gas-phase conc. emission file :', emis_gas_file
	    if (ssh_standalone) write(*,*) 'Particle conc. emission file :', emis_aero_mass_file
	    if (ssh_logger) write(logfile,*) 'Particle conc. emission file :', emis_aero_mass_file
            if (with_emis_num == 1) then
	       if (ssh_standalone) write(*,*) 'Emitted aerosol number conc. is read from file :', emis_aero_num_file
	       if (ssh_logger) write(logfile,*) 'Emitted aerosol number conc. is read from file :', emis_aero_num_file
            else  ! default with_emis_num == 0 
	       with_emis_num = 0
	       if (ssh_standalone) write(*,*) 'Emitted aerosol number conc. is estimated from mass and diameter.'
	       if (ssh_logger) write(logfile,*) 'Emitted aerosol number conc. is estimated from mass and diameter.'
	    end if
        end if 
     end if


     ! mixing_state
     read(10, nml = mixing_state, iostat = ierr)
     if (ierr .ne. 0) then
	write(*,*) "mixing_state data can not be read - not enough data"
        stop
     else 
	if (ssh_standalone) write(*,*) ''
	if (ssh_logger) write(logfile,*) ''
	if (ssh_standalone) write(*,*) '<<<< Mixing state >>>>'
	if (ssh_logger) write(logfile,*) '<<<< Mixing state >>>>'
	if (tag_external == 1) then
	   if (ssh_standalone) write(*,*) 'simulation is mixing-state resolved.'
	   if (ssh_logger) write(logfile,*) 'simulation is mixing-state resolved.'
           if (ssh_standalone) write(*,*) 'number of composition fractions:', N_frac 
           if (ssh_logger) write(logfile,*) 'number of composition fractions:', N_frac 
	   if (ssh_standalone) write(*,*) 'number of grounps :', N_groups
	   if (ssh_logger) write(logfile,*) 'number of grounps :', N_groups
	else  
	   tag_external = 0   ! default tag_external == 0
	   if (ssh_standalone) write(*,*) 'simulation is internally mixed.' 
	   if (ssh_logger) write(logfile,*) 'simulation is internally mixed.' 
	   N_frac = 1
	   N_groups = 1
	end if
	if (tag_emis == 2 .and. tag_external == 0) then
           write(*,*) " Pb -- internally mixed simulation can not process externally-mixed emissions data."        
	   stop
	end if
	   
        if (kind_composition == 1) then
	   allocate(frac_input(N_frac+1)) 
	else
	   kind_composition = 0 ! default
	   allocate(frac_input(2)) 
	end if
    end if

    ! fraction_distribution
     allocate(frac_bound(N_frac+1))
     read(10, nml = fraction_distribution, iostat = ierr)
     if (ierr .ne. 0) then
	write(*,*) "mixing_state data can not be read."
        stop
     else
	if (kind_composition == 1)  then
	   if (ssh_standalone) write(*,*) 'manual fraction discretization methods.'
	   if (ssh_logger) write(logfile,*) 'manual fraction discretization methods.'
	   frac_bound = frac_input
        else
	   if (frac_input(1) .ge. frac_input(2)) then  ! default
		frac_input(1) = 0.d0
		frac_input(2) = 1.d0
	   end if
           do i = 1, N_frac+1
	      frac_bound(i)= (dble(i-1) / N_frac * (frac_input(2) - frac_input(1)) ) + frac_input(1)
           enddo
           if (ssh_standalone) write(*,*) 'auto fraction discretization methods.'
           if (ssh_logger) write(logfile,*) 'auto fraction discretization methods.'
	end if
	if (ssh_standalone) write(*,*) 'fraction bounds :', frac_bound
	if (ssh_logger) write(logfile,*) 'fraction bounds :', frac_bound
     end if

     ! species
     read(10, nml = gas_phase_species, iostat = ierr)
     if (ierr .ne. 0) then
	write(*,*) "gas_phase_species data can not be read."
        stop
     else
	if (ssh_standalone) write(*,*) ''
	if (ssh_logger) write(logfile,*) ''
	if (ssh_standalone) write(*,*) '<<<< Species lists >>>>'
	if (ssh_logger) write(logfile,*) '<<<< Species lists >>>>'
	if (ssh_standalone) write(*,*) 'gas phase species file :', species_list_file
	if (ssh_logger) write(logfile,*) 'gas phase species file :', species_list_file
     end if

     read(10, nml = aerosol_species, iostat = ierr)
     if (ierr .ne. 0) then
	write(*,*) "aerosol_species data can not be read."
        stop
     else
	if (ssh_standalone) write(*,*) 'particle species file :', aerosol_species_list_file
	if (ssh_logger) write(logfile,*) 'particle species file :', aerosol_species_list_file
     end if

     ! gas chemistry
     read(10, nml = physic_gas_chemistry, iostat = ierr)
     if (ierr .ne. 0) then
	write(*,*) "physic_gas_chemistry data can not be read."
        stop
     else

	if (tag_chem == 0) then
		if (ssh_standalone) write(*,*) ''
		if (ssh_logger) write(logfile,*) ''
		if (ssh_standalone) write(*,*) '<<<< Without Gas-phase chemistry >>>>'
		if (ssh_logger) write(logfile,*) '<<<< Without Gas-phase chemistry >>>>'
	else
		if (ssh_standalone) write(*,*) ''
		if (ssh_logger) write(logfile,*) ''
		if (ssh_standalone) write(*,*) '<<<< Gas-phase chemistry >>>>'
		if (ssh_logger) write(logfile,*) '<<<< Gas-phase chemistry >>>>'
        	if (with_heterogeneous == 1) then
	   		if (ssh_standalone) write(*,*) 'with heterogeneous reaction.'
	   		if (ssh_logger) write(logfile,*) 'with heterogeneous reaction.'
		else  ! default with_heterogeneous == 0
	   		with_heterogeneous = 0
	   		if (ssh_standalone) write(*,*) 'without heterogeneous reaction.' 
	   		if (ssh_logger) write(logfile,*) 'without heterogeneous reaction.' 
		end if
		if (with_adaptive == 1) then
	   		if (ssh_standalone) write(*,*) 'with adaptive step.'
	   		if (ssh_logger) write(logfile,*) 'with adaptive step.'
	   		if (ssh_standalone) write(*,*) 'adaptive time step tolerance', adaptive_time_step_tolerance
	   		if (ssh_logger) write(logfile,*) 'adaptive time step tolerance', adaptive_time_step_tolerance
	   		if (ssh_standalone) write(*,*) 'min adaptive time step', min_adaptive_time_step
	   		if (ssh_logger) write(logfile,*) 'min adaptive time step', min_adaptive_time_step
        	else
	   		with_adaptive = 0
	   		if (ssh_standalone) write(*,*) 'without adaptive step.'
	   		if (ssh_logger) write(logfile,*) 'without adaptive step.'
		end if
	end if
     end if

     ! particle numerical issues
     read(10, nml = physic_particle_numerical_issues, iostat = ierr)
     if (ierr .ne. 0) then
	write(*,*) "physic_particle_numerical_issues data can not be read."
        stop
     else
	if (ssh_standalone) write(*,*) ''
	if (ssh_logger) write(logfile,*) ''
	if (ssh_standalone) write(*,*) '<<<< Particle numerical issues >>>>'
	if (ssh_logger) write(logfile,*) '<<<< Particle numerical issues >>>>'
	if (ssh_standalone) write(*,*) 'minimum time step for aerosol dynamics', DTAEROMIN,'s'
	if (ssh_logger) write(logfile,*) 'minimum time step for aerosol dynamics', DTAEROMIN,'s'

	SELECT CASE (redistribution_method)
	   CASE (3)
	      if (ssh_standalone) write(*,*) 'redistibution method : euler_mass'
	      if (ssh_logger) write(logfile,*) 'redistibution method : euler_mass'
	   CASE (4)
	      if (ssh_standalone) write(*,*) 'redistibution method : euler_number'
	      if (ssh_logger) write(logfile,*) 'redistibution method : euler_number'
	   CASE (5)
	      if (ssh_standalone) write(*,*) 'redistibution method : hemen'
	      if (ssh_logger) write(logfile,*) 'redistibution method : hemen'
	   CASE (6)
	      if (ssh_standalone) write(*,*) 'redistibution method : euler_coupled'
	      if (ssh_logger) write(logfile,*) 'redistibution method : euler_coupled'
	   CASE (10)
	      if (ssh_standalone) write(*,*) 'redistibution method : moving diameter'
	      if (ssh_logger) write(logfile,*) 'redistibution method : moving diameter'
	   CASE (11)
	      if (ssh_standalone) write(*,*) 'redistibution method : siream'
	      if (ssh_logger) write(logfile,*) 'redistibution method : siream'
	   CASE (12)
	      if (ssh_standalone) write(*,*) 'redistibution method : siream - euler coupled'
	      if (ssh_logger) write(logfile,*) 'redistibution method : siream - euler coupled'
	   CASE (13)
	      if (ssh_standalone) write(*,*) 'redistibution method : siream - moving diameter'
	      if (ssh_logger) write(logfile,*) 'redistibution method : siream - moving diameter'
	   CASE DEFAULT ! default redistribution_method = 0
	      redistribution_method = 0
	      if (ssh_standalone) write(*,*) '! ! without redistibution'
	      if (ssh_logger) write(logfile,*) '! ! without redistibution'
	END SELECT

	if (with_fixed_density == 1) then
	   if (ssh_standalone) write(*,*) 'fixed density is used.', fixed_density, 'kg/m^3'
	   if (ssh_logger) write(logfile,*) 'fixed density is used.', fixed_density, 'kg/m^3'
	else
	   if (ssh_standalone) write(*,*) 'real density is computed.'
	   if (ssh_logger) write(logfile,*) 'real density is computed.'
        end if

	if (splitting == 1) then
	   if (ssh_standalone) write(*,*) 'coagulation - cond/evap+nucl are coupled'
	   if (ssh_logger) write(logfile,*) 'coagulation - cond/evap+nucl are coupled'
	else
	   if (ssh_standalone) write(*,*) 'coagulation - cond/evap+nucl are splitted'
	   if (ssh_logger) write(logfile,*) 'coagulation - cond/evap+nucl are splitted'
        end if


	if (wet_diam_estimation == 1) then
	   if (ssh_standalone) write(*,*) 'initial wet diameter is computed from initial water conc. if available'
	   if (ssh_logger) write(logfile,*) 'initial wet diameter is computed from initial water conc. if available'
	else !default wet_diam_estimation == 0
	   wet_diam_estimation = 0
	   if (ssh_standalone) write(*,*) 'initial wet diameter is computed by isorropia.'
	   if (ssh_logger) write(logfile,*) 'initial wet diameter is computed by isorropia.'
	end if
     end if

     ! coagulation
     read(10, nml = physic_coagulation, iostat = ierr)
     if (ierr .ne. 0) then
	write(*,*) "physic_coagulation data can not be read."
        stop
     else 
	if (with_coag == 1) then
	   if (ssh_standalone) write(*,*) '! ! ! with coagulation.'
	   if (ssh_logger) write(logfile,*) '! ! ! with coagulation.'
           if(i_compute_repart == 0) then
             if (ssh_standalone) write(*,*) i_compute_repart,'repartition coefficient are read'
             if (ssh_logger) write(logfile,*) i_compute_repart,'repartition coefficient are read'
	     if (ssh_standalone) write(*,*) 'coefficient file : ', Coefficient_file
	     if (ssh_logger) write(logfile,*) 'coefficient file : ', Coefficient_file
           else
             if (ssh_standalone) write(*,*) i_compute_repart,'repartition coefficient are computed'
             if (ssh_logger) write(logfile,*) i_compute_repart,'repartition coefficient are computed'
	     if (ssh_standalone) write(*,*) 'Nmc = ',Nmc
	     if (ssh_logger) write(logfile,*) 'Nmc = ',Nmc
           endif
	else
	   with_coag = 0
	   if (ssh_standalone) write(*,*) 'without coagulation.'
	   if (ssh_logger) write(logfile,*) 'without coagulation.'
	end if
     end if

     ! condensation/ evaporation
     read(10, nml = physic_condensation, iostat = ierr)
     if (ierr .ne. 0) then
	write(*,*) "physic_condensation data can not be read."
        stop
     else
	if (with_cond == 1)  then ! defalut
	   if (ssh_standalone) write(*,*) '! ! ! with condensation/ evaporation.'
	   if (ssh_logger) write(logfile,*) '! ! ! with condensation/ evaporation.'
		if (Cut_dim .le. diam_input(1)) then
		    if (ssh_standalone) write(*,*) 'inorganic compounds are at dynamic.'
		    if (ssh_logger) write(logfile,*) 'inorganic compounds are at dynamic.'
		else if (Cut_dim .ge. diam_input(size(diam_input))) then
		    if (ssh_standalone) write(*,*) 'inorganic compounds are at equilibrium.'
		    if (ssh_logger) write(logfile,*) 'inorganic compounds are at equilibrium.'
		else 
		    if (ssh_standalone) write(*,*) 'Cutting diameter for equilibrium for inorganic compounds :', Cut_dim
		    if (ssh_logger) write(logfile,*) 'Cutting diameter for equilibrium for inorganic compounds :', Cut_dim

		end if

                if (ISOAPDYN == 1) then
	   	    if (ssh_standalone) write(*,*) 'Dynamic SOA computation are at dynamic.'
	   	    if (ssh_logger) write(logfile,*) 'Dynamic SOA computation are at dynamic.'
		else
	   	    ISOAPDYN = 0   ! defalut
	   	    if (ssh_standalone) write(*,*) 'Dynamic SOA computation are at equilibrium.'
	   	    if (ssh_logger) write(logfile,*) 'Dynamic SOA computation are at equilibrium.'
        	end if
	else
	   with_cond = 0
	   if (ssh_standalone) write(*,*) 'without condensation/ evaporation.'
	   if (ssh_logger) write(logfile,*) 'without condensation/ evaporation.'
	end if
     end if

     ! nucleation
     read(10, nml = physic_nucleation, iostat = ierr)
     if (ierr .ne. 0) then
	write(*,*) "physic_nucleation data can not be read."
        stop
     else
	if (with_nucl == 1) then
	   if (ssh_standalone) write(*,*) '! ! ! with nucleation. -- Need to have lower diameter about 1nm.'
	   if (ssh_logger) write(logfile,*) '! ! ! with nucleation. -- Need to have lower diameter about 1nm.'
	   if (nucl_model == 0) then
	      if (ssh_standalone) write(*,*) 'nucleation model : binary'
	      if (ssh_logger) write(logfile,*) 'nucleation model : binary'
	   else 
              nucl_model = 1
	      if (ssh_standalone) write(*,*) 'nucleation model : ternary'
	      if (ssh_logger) write(logfile,*) 'nucleation model : ternary'
	   end if
	else
	   with_nucl = 0   ! defalut
	   if (ssh_standalone) write(*,*) 'Without nucleation.'
	   if (ssh_logger) write(logfile,*) 'Without nucleation.'
	end if
     end if

     ! organic
     read(10, nml = physic_organic, iostat = ierr)
     if (ierr .ne. 0) then
	write(*,*) "physic_organic data can not be read."
        stop
     else
	if (with_oligomerization == 1) then
	   if (ssh_standalone) write(*,*) 'with oligomerization.'
	   if (ssh_logger) write(logfile,*) 'with oligomerization.'
	else
	   with_oligomerization = 0   ! defalut
	   if (ssh_standalone) write(*,*) 'without oligomerization.'
	   if (ssh_logger) write(logfile,*) 'without oligomerization.'
	end if
     end if

     ! output
     read(10, nml = output, iostat = ierr)
     if (ierr .ne. 0) then
	write(*,*) "output setting can not be read."
        stop
     else
	if (ssh_standalone) write(*,*) ''
	if (ssh_logger) write(logfile,*) ''
	if (ssh_standalone) write(*,*) '<<<< Results output >>>>'
	if (ssh_logger) write(logfile,*) '<<<< Results output >>>>'
	if (output_type == 2) then
	   if (ssh_standalone) write(*,*) 'results are saved in binary files.'
	   if (ssh_logger) write(logfile,*) 'results are saved in binary files.'
	else
	   output_type = 1   ! defalut
	   if (ssh_standalone) write(*,*) 'results are saved in text files.'
	   if (ssh_logger) write(logfile,*) 'results are saved in text files.'
	end if
	if (ssh_standalone) write(*,*) 'output directory :', output_directory
	if (ssh_logger) write(logfile,*) 'output directory :', output_directory
     end if
     
     close(10)
     
     if (ssh_standalone) write(*,*) "=========================finish read namelist.ssh file======================"
     if (ssh_logger) write(logfile,*) "=========================finish read namelist.ssh file======================"


end subroutine read_namelist



! =============================================================
!
!      read initial/emitted species file
!
! 
! ============================================================= 


subroutine read_inputs()
    

	implicit none
	integer :: k,i,j,s,js, ind, count, ierr
	double precision :: tmp
	double precision, dimension(:), allocatable :: tmp_aero
	character (len=40) :: ic_name, sname
     
     ! read gas-phase species namelist ! unit = 11
     allocate(molecular_weight(N_gas))
     allocate(species_name(N_gas))

     open(unit = 11, file = species_list_file, status = "old")
     count = 0
     ierr = 0
     do while(ierr .eq. 0)
        read(11, *, iostat=ierr)
        if (ierr == 0) count = count + 1
     end do

     if (ssh_standalone) write(*,*) 'read gas-phase species list.'
     if (ssh_logger) write(logfile,*) 'read gas-phase species list.'
     if (N_gas == count - 1) then   ! minus the first comment line
	if (ssh_standalone) write(*,*) 'Number of gas-phase species', N_gas
	if (ssh_logger) write(logfile,*) 'Number of gas-phase species', N_gas
     else 
	write(*,*) 'Given gas-phase species list does not fit chem() setting.'
	stop
     end if

     rewind 11
     read(11, *)  ! read the first comment line
     do s = 1, N_gas
        read(11, *) species_name(s), molecular_weight(s)
     enddo
     close(11)

     ! read aerosol species namelist ! unit = 12
     open(unit = 12, file = aerosol_species_list_file, status = "old")
     count = 0
     ierr = 0
     do while(ierr .eq. 0)
        read(12, *, iostat=ierr)
        if (ierr == 0) count = count + 1
     end do
     if (ssh_standalone) write(*,*) 'read aerosol species list.'
     if (ssh_logger) write(logfile,*) 'read aerosol species list.'
     N_aerosol = count - 1  ! minus the first comment line

     allocate(aerosol_species_name(N_aerosol))
     spec_name_len = len(aerosol_species_name(1))
     allocate(Index_groups(N_aerosol))
     ! initialize basic physical and chemical parameters
     allocate(molecular_weight_aer(N_aerosol))
     ! precursor
!!     allocate(saturation_pressure_mass(N_aerosol))        
!!     allocate(saturation_pressure_torr(N_aerosol))
!!     allocate(partition_coefficient(N_aerosol))
!!     allocate(deliquescence_relative_humidity(N_aerosol))
     allocate(collision_factor_aer(N_aerosol))
     allocate(molecular_diameter(N_aerosol)) 
     allocate(surface_tension(N_aerosol))
     allocate(accomodation_coefficient(N_aerosol))
     allocate(mass_density(N_aerosol))
!!     allocate(saturation_pressure(N_aerosol))
!!     allocate(vaporization_enthalpy(N_aerosol))
     ! aerosol species
     allocate(List_species(N_aerosol))
     allocate(isorropia_species(nesp_isorropia))
     allocate(aec_species(nesp_aec))
     allocate(pankow_species(nesp_pankow)) 
     allocate(poa_species(nesp_pom))
     ! relation between Aerosol and GAS
     allocate(aerosol_species_interact(N_aerosol))      
     aerosol_species_interact = 0

     rewind 12
     count = 0
     read(12, *) ! red # comment line
     do s = 1, N_aerosol
        ! read(12, *) aerosol_species_name(s), Index_groups(s), molecular_weight_aer(s), &
        !      precursor, saturation_pressure_mass(s), &
        !      saturation_pressure_torr(s),  partition_coefficient(s), &
        !      deliquescence_relative_humidity(s), &
        !      collision_factor_aer(s), molecular_diameter(s), &
        !      surface_tension(s), accomodation_coefficient(s), &
        !      mass_density(s), saturation_pressure(s), &
        !      vaporization_enthalpy(s)

        !=== Warning (YK/KS) ===
        ! Surface_tension for organic and aqueous phases of organic aerosols
        ! is hardly coded in SOAP/parameters.cxx
        ! And Unit used in SOAP is different to surface_tension (N/m) by 1.e3.
        read(12, *) aerosol_species_name(s), Index_groups(s), molecular_weight_aer(s), &
             precursor, &
             collision_factor_aer(s), molecular_diameter(s), &
             surface_tension(s), accomodation_coefficient(s), &
             mass_density(s)
        
        molecular_weight_aer(s) = molecular_weight_aer(s) * 1.0D06 ! g/mol to \B5g/mol  !!! change later
        list_species(s) = s
        if (aerosol_species_name(s) .eq. "PMD") EMD = s
        if (aerosol_species_name(s) .eq. "PBC") EBC = s
        if (aerosol_species_name(s) .eq. "PNA") ENa = s
        if (aerosol_species_name(s) .eq. "PSO4") ESO4 = s
        if (aerosol_species_name(s) .eq. "PNH4") ENH4 = s
        if (aerosol_species_name(s) .eq. "PNO3") ENO3 = s
        if (aerosol_species_name(s) .eq. "PHCL") ECl = s

        do js = 1, nesp_isorropia
           if (aerosol_species_name(s) .eq. isorropia_species_name(js)) then
              isorropia_species(js) = s
           endif
        enddo
        do js = 1, nesp_aec
           if (aerosol_species_name(s) .eq. aec_species_name(js)) then
              aec_species(js) = s
           endif
        enddo
        do js = 1, nesp_pankow
           if (aerosol_species_name(s) .eq. pankow_species_name(js)) then
              pankow_species(js) = s
           endif
        enddo
        do js = 1, nesp_pom
           if (aerosol_species_name(s) .eq. poa_species_name(js)) then
              poa_species(js) = s
           endif
        enddo

        ind = 0
        do js = 1, N_gas
           if (species_name(js) .eq. trim(precursor)) then
              aerosol_species_interact(s) = js
              count = count + 1
              ind = 1
           endif
           if (ind == 1) exit
        enddo
        !! Check if a precursor name is found in the list of gas-phase species.
        if ((ind .eq. 0) .and. (trim(precursor) .ne. "--")) then
           if (ssh_standalone) write(*,*) "Error: wrong species name is given ", aerosol_species_list_file, trim(precursor)
           if (ssh_logger) write(logfile,*) "Error: wrong species name is given ", aerosol_species_list_file, trim(precursor)
           stop
        endif
        
     enddo
     close(12)

    ! Initialize the concentrations of percursors
     if (ssh_standalone) write(*,*) "Number of precursors: ", count 
     if (ssh_logger) write(logfile,*) "Number of precursors: ", count 


      
     ! read gas-phase initial concentrations unit 21
     ! no comment lines for initial & emitted data
     allocate(concentration_gas_all(N_gas))
     concentration_gas_all = 0.d0 ! set original value to 0

     open(unit = 21, file = init_gas_conc_file, status = "old")
     count = 0 
     ierr = 0
     do while(ierr .eq. 0)
        read(21, *, iostat=ierr)
        if (ierr == 0) count = count + 1
     end do
	count = count - 1 ! minus comment line

     rewind 21

        read(21,*)
     do s= 1, count
        read(21,*) ic_name, tmp
        ind = 0
        do js = 1, N_gas
           if (species_name(js) .eq. ic_name) then
              concentration_gas_all(js) = tmp
              ind = 1
           endif
	   if (ind == 1) exit
        enddo
        if (ind .eq. 0) then
           if (ssh_standalone) write(*,*) "Error: wrong species name is given ", init_gas_conc_file, ic_name
           if (ssh_logger) write(logfile,*) "Error: wrong species name is given ", init_gas_conc_file, ic_name
        endif
     enddo
     close(21)
  if (ssh_standalone) write(*,*) 'gas concentrations have been read'
  if (ssh_logger) write(logfile,*) 'gas concentrations have been read'

  if (tag_init == 0) then  ! change if species list and init are not in the same order
     ! Read aerosol initial mass concentrations unit 22
     ! internally mixed : mass for each sizebin of each species is given
     open(unit = 22, file = init_aero_conc_mass_file, status = "old") 
     count = 0
     ierr = 0
     do while(ierr .eq. 0)
        read(22, *, iostat=ierr)
        if (ierr == 0) count = count + 1
     end do  
	count = count - 1 ! minus comment line
     allocate(init_mass(N_aerosol))   ! aerosol initial mass concentrations inti_mass for each species
     init_mass = 0.d0
     allocate(init_bin_mass(N_sizebin,N_aerosol))
     init_bin_mass = 0.d0
     allocate(tmp_aero(N_sizebin))
     tmp_aero = 0.d0
     aero_total_mass = 0.d0
     rewind 22

        read(22,*)
     do s= 1, count
        read(22,*) ic_name, (tmp_aero(k), k = 1, N_sizebin)
        ind = 0
        do js = 1, N_aerosol
           if (aerosol_species_name(js) .eq. ic_name) then
               do k=1, N_sizebin
	          init_bin_mass(k,js) = tmp_aero(k)
		  init_mass(js)=init_mass(js) + tmp_aero(k)
        	  aero_total_mass = aero_total_mass +  init_mass(js) !total mass
               enddo
               ind = 1
           endif
	   if (ind  == 1) exit
        enddo
        if (ind .eq. 0) then
           if (ssh_standalone) write(*,*) "Error: wrong aerosol species name is given ",init_aero_conc_mass_file, ic_name
           if (ssh_logger) write(logfile,*) "Error: wrong aerosol species name is given ",init_aero_conc_mass_file, ic_name
        endif
     enddo

     close(22)

     if (ssh_standalone) write(*,*) 'initial mass concentrations have been read'
     if (ssh_logger) write(logfile,*) 'initial mass concentrations have been read'
  else if (tag_init == 1) then ! mixing_state resolved
	     ! need to fill in the future
 	write(*,*) "Tag_init = 1, mixing_state resolved - not yet available"
        stop
  end if

     ! Read aerosol initial number concentrations unit 23
     ! auto-generate or read number and size bins distributions

   if (with_init_num == 1) then
      if (tag_init == 0) then
     	 open(unit = 23, file = init_aero_conc_num_file, status = "old") 
	      read(23,*, iostat = ierr) ! read comment line
	      read(23,*, iostat = ierr) ic_name, init_bin_number
         close(23)
	 if (ierr == 0) then
            if (ssh_standalone) write(*,*) "Aerosol number conc. is read."
            if (ssh_logger) write(logfile,*) "Aerosol number conc. is read."
            if (ssh_standalone) write(*,*) "Aerosol initial number conc. distribution :", init_bin_number
            if (ssh_logger) write(logfile,*) "Aerosol initial number conc. distribution :", init_bin_number
	 else 
            write(*,*) "Aerosol number conc. can not be read from file ",init_aero_conc_num_file
	    stop
         end if 
      else if (tag_init == 1) then
               write(*,*) 'with_init_num == 1 .and. tag_init == 1 - not yet build'
	       stop
      end if
   end if 


   ! ! ! ! ! ! 
  if (tag_emis == 1) then  ! with internal emission 

        allocate(emis_bin_mass(N_sizebin,N_aerosol))
	emis_bin_mass = 0.d0
        allocate(gas_emis(N_gas))
	gas_emis = 0.d0
        ! read gas emission concentrations unit 31
        open(unit=31, file = emis_gas_file, status = "old")
        count = 0
        ierr = 0
        do while(ierr .eq. 0)
           read(31, *, iostat=ierr)
        if (ierr == 0) count = count + 1
        end do
	count = count - 1 ! minus comment line
        if (ssh_standalone) write(*,*) "Number of emitted gas-phase species:", count
        if (ssh_logger) write(logfile,*) "Number of emitted gas-phase species:", count
        
        rewind 31
        read(31,*)
        do s = 1, count
           read(31, *) ic_name, tmp
           ind = 0
	   do js = 1, N_gas
		if (species_name(js) .eq. ic_name) then
		   gas_emis(js) = tmp
		   ind = 1
	!if (ssh_standalone) write(*,*) 'gas_emis', species_name(js), gas_emis(js)
	!if (ssh_logger) write(logfile,*) 'gas_emis', species_name(js), gas_emis(js)
                endif
	        if (ind == 1) exit
           enddo
           if (ind .eq. 0) then
              if (ssh_standalone) write(*,*) "Error: wrong species name is given in gas emission", init_gas_conc_file, ic_name
              if (ssh_logger) write(logfile,*) "Error: wrong species name is given in gas emission", init_gas_conc_file, ic_name
           end if
        end do

        close(31)
  
       ! Read aerosol emission concentrations unit 32
	tmp_aero = 0.d0
       open(unit=32, file = emis_aero_mass_file, status = "old")
       count = 0
       ierr = 0
       do while(ierr .eq. 0)
          read(32, *, iostat=ierr)
        if (ierr == 0) count = count + 1
       end do
	count = count - 1 ! minus comment line
       if (ssh_standalone) write(*,*) "Number of emitted aerosols species:", count
       if (ssh_logger) write(logfile,*) "Number of emitted aerosols species:", count

       rewind 32
       read(32,*)
       if (Tag_init .eq. 0) then
 	    do s=1, count
                read(32,*, iostat = ierr) ic_name, (tmp_aero(k),k=1,N_sizebin)
		if (ierr .ne. 0) then
                  write(*,*) "Error when reading ic_name."
                  stop
                endif
		ind = 0
		do js = 1, N_aerosol
			if (aerosol_species_name(js) == ic_name) then
				do k = 1, N_sizebin
					emis_bin_mass(k,js) = tmp_aero(k)
				end do
				ind = 1
			end if
			if (ind ==1) exit
		end do
		if (ind == 0 .and. ssh_standalone) write(*,*) 'Not find the emission species', emis_aero_mass_file, ic_name
		if (ind == 0 .and. ssh_logger) write(logfile,*) 'Not find the emission species', emis_aero_mass_file, ic_name
            enddo
            if (ssh_standalone) write(*,*) "Emission mass conc. has been read."
            if (ssh_logger) write(logfile,*) "Emission mass conc. has been read."
        else if (tag_init .eq. 1) then ! external mixed
		 if (ssh_standalone) write(*,*) "Not yet build -- tag_emis == 1 and tag_init == 1" 
		 if (ssh_logger) write(logfile,*) "Not yet build -- tag_emis == 1 and tag_init == 1" 
                 tag_emis = 0
        endif
       close(32)

      ! Read aerosol number emission concentrations unit 33 if need
      if (with_emis_num == 1) then
         open(unit=33, file = emis_aero_num_file, status = "old")
              read(33,*, iostat=ierr)    
              read(33,*, iostat=ierr) ic_name, emis_bin_number
         close(33)
         if (ierr .eq. 0) then
	     if (ssh_standalone) write(*,*) "Emission number conc. has been read."
	     if (ssh_logger) write(logfile,*) "Emission number conc. has been read."
	     if (ssh_standalone) write(*,*) 'emis_bin_number', emis_bin_number
	     if (ssh_logger) write(logfile,*) 'emis_bin_number', emis_bin_number
         else 
            write(*,*) "can not read aerosol number conc. from ", emis_aero_num_file
            stop
         end if
      end if

   else if  (tag_emis == 2) then
	  write(*,*) "Not yet build -- with externally-mixed emissions -- tag_emis = 2" 
            ! 1 mixing_state resolved !! option 1 not yet available)
          stop
   end if
   if (ssh_standalone) write(*,*) "=========================finish read inputs file======================"
   if (ssh_logger) write(logfile,*) "=========================finish read inputs file======================"
   if (allocated(tmp_aero))  deallocate(tmp_aero)
   end subroutine read_inputs


! ============================================================
!
!  !   ! free allocated memory
!
! ============================================================

  subroutine free_allocated_memory()
    
    integer :: ierr = 0

	! read_namelist()
	if (allocated(init_bin_number))  deallocate(init_bin_number, stat=ierr)
	if (allocated(emis_bin_number))  deallocate(emis_bin_number, stat=ierr)
	if (allocated(diam_input))  deallocate(diam_input, stat=ierr)
	if (allocated(diam_bound))  deallocate(diam_bound, stat=ierr)   ! init_parameter()
	if (allocated(frac_input))  deallocate(frac_input, stat=ierr)
	if (allocated(frac_bound))  deallocate(frac_bound, stat=ierr)
       if (ierr .ne. 0) then
          write(*,*) "Deallocation error"
          stop
       endif
	! read_inputs()
	if (allocated(init_mass))  deallocate(init_mass, stat=ierr)  
	if (allocated(init_bin_mass))  deallocate(init_bin_mass, stat=ierr)
	if (allocated(emis_bin_mass))  deallocate(emis_bin_mass, stat=ierr)
	if (allocated(gas_emis))  deallocate(gas_emis, stat=ierr)

	if (allocated(molecular_weight))  deallocate(molecular_weight, stat=ierr)
	if (allocated(species_name))  deallocate(species_name, stat=ierr)
	if (allocated(aerosol_species_name))  deallocate(aerosol_species_name, stat=ierr)
	if (allocated(Index_groups))  deallocate(Index_groups, stat=ierr)
	if (allocated(molecular_weight_aer))  deallocate(molecular_weight_aer, stat=ierr)
!!	if (allocated(saturation_pressure_mass))  deallocate(saturation_pressure_mass, stat=ierr)
!!	if (allocated(saturation_pressure_torr))  deallocate(saturation_pressure_torr, stat=ierr)
!!	if (allocated(partition_coefficient))  deallocate(partition_coefficient, stat=ierr)
!!	if (allocated(deliquescence_relative_humidity))  deallocate(deliquescence_relative_humidity, stat=ierr)
	if (allocated(collision_factor_aer))  deallocate(collision_factor_aer, stat=ierr)
	if (allocated(molecular_diameter))  deallocate(molecular_diameter, stat=ierr)
	if (allocated(surface_tension))  deallocate(surface_tension, stat=ierr)
	if (allocated(accomodation_coefficient))  deallocate(accomodation_coefficient, stat=ierr)
	if (allocated(mass_density))  deallocate(mass_density, stat=ierr)
!!	if (allocated(saturation_pressure))  deallocate(saturation_pressure, stat=ierr)
!!	if (allocated(vaporization_enthalpy))  deallocate(vaporization_enthalpy, stat=ierr)
	if (allocated(List_species))  deallocate(List_species, stat=ierr)
	if (allocated(isorropia_species))  deallocate(isorropia_species, stat=ierr)
	if (allocated(aec_species))  deallocate(aec_species, stat=ierr)
	if (allocated(pankow_species))  deallocate(pankow_species, stat=ierr)
	if (allocated(poa_species))  deallocate(poa_species, stat=ierr)
	if (allocated(aerosol_species_interact))  deallocate(aerosol_species_interact, stat=ierr)

	if (allocated(concentration_gas_all))  deallocate(concentration_gas_all, stat=ierr)
       if (ierr .ne. 0) then
          write(*,*) "Deallocation error"
          stop
       endif
	! init_parameters()
	if (allocated(N_fracbin))  deallocate(N_fracbin, stat=ierr)
	if (allocated(concentration_index))  deallocate(concentration_index, stat=ierr)
	if (allocated(concentration_index_iv))  deallocate(concentration_index_iv, stat=ierr)
	if (allocated(size_diam_av))  deallocate(size_diam_av, stat=ierr)
	if (allocated(size_mass_av))  deallocate(size_mass_av, stat=ierr)

	if (allocated(per_mass_init))  deallocate(per_mass_init, stat=ierr)

	if (allocated(density_aer_bin))  deallocate(density_aer_bin, stat=ierr)
	if (allocated(density_aer_size))  deallocate(density_aer_size, stat=ierr)
	if (allocated(rho_wet_cell))  deallocate(rho_wet_cell, stat=ierr)

	if (allocated(photolysis))  deallocate(photolysis, stat=ierr)
	if (allocated(photolysis_reaction_index))  deallocate(photolysis_reaction_index, stat=ierr)
	if (allocated(source_index))  deallocate(source_index, stat=ierr)
	if (allocated(source))  deallocate(source, stat=ierr)
	if (allocated(conversionfactor))  deallocate(conversionfactor, stat=ierr)
	if (allocated(conversionfactorjacobian))  deallocate(conversionfactorjacobian, stat=ierr)

	if (allocated(lwc_nsize))  deallocate(lwc_nsize, stat=ierr)
	if (allocated(ionic_nsize))  deallocate(ionic_nsize, stat=ierr)
        if (allocated(proton_nsize))  deallocate(proton_nsize, stat=ierr)
        if (allocated(chp_nsize))  deallocate(chp_nsize, stat=ierr)
	if (allocated(liquid_nsize))  deallocate(liquid_nsize, stat=ierr)
	if (allocated(concentration_inti))  deallocate(concentration_inti, stat=ierr)

	if (allocated(total_mass))  deallocate(total_mass, stat=ierr)
!!	if (allocated(total_mass_old))  deallocate(total_mass_old, stat=ierr)
	if (allocated(concentration_gas))  deallocate(concentration_gas, stat=ierr)
	if (allocated(concentration_number))  deallocate(concentration_number, stat=ierr)
	if (allocated(concentration_number_tmp))  deallocate(concentration_number_tmp, stat=ierr)
	if (allocated(concentration_mass))  deallocate(concentration_mass, stat=ierr)
	if (allocated(concentration_mass_tmp))  deallocate(concentration_mass_tmp, stat=ierr)

	if (allocated(emission_rate))  deallocate(emission_rate, stat=ierr)
	if (allocated(emission_num_rate))  deallocate(emission_num_rate, stat=ierr)

	if (allocated(cell_mass_av))  deallocate(cell_mass_av, stat=ierr)
	if (allocated(cell_mass))  deallocate(cell_mass, stat=ierr)
	if (allocated(cell_diam_av))  deallocate(cell_diam_av, stat=ierr)

	if (allocated(total_aero_mass))  deallocate(total_aero_mass, stat=ierr)
	if (allocated(mass_total_grid))  deallocate(mass_total_grid, stat=ierr)
	if (allocated(wet_mass))  deallocate(wet_mass, stat=ierr)
	if (allocated(wet_diameter))  deallocate(wet_diameter, stat=ierr)
	if (allocated(wet_volume))  deallocate(wet_volume, stat=ierr)
	if (allocated(bin_mass))  deallocate(bin_mass, stat=ierr)
	if (allocated(bin_number))  deallocate(bin_number, stat=ierr)
	if (allocated(ce_kernal_coef))  deallocate(ce_kernal_coef, stat=ierr)
	if (allocated(frac_grid))  deallocate(frac_grid, stat=ierr)
	if (allocated(dqdt))  deallocate(dqdt, stat=ierr)
       if (ierr .ne. 0) then
          write(*,*) "Deallocation error"
          stop
       endif
	! discretization()
	if (allocated(discretization_composition))  deallocate(discretization_composition, stat=ierr)
	! Init_distributions()
	! for coagulation
	if (allocated(kernel_coagulation))  deallocate(kernel_coagulation, stat=ierr)

	! for condensation
	if (allocated(quadratic_speed)) deallocate(quadratic_speed, stat=ierr)
	if (allocated(diffusion_coef))  deallocate(diffusion_coef, stat=ierr)
	! if (allocated(soa_sat_conc)) deallocate(soa_sat_conc, stat=ierr)
	! if (allocated(soa_part_coef)) deallocate(soa_part_coef, stat=ierr)
	if (allocated(discretization_mass)) deallocate(discretization_mass)
       if (ierr .ne. 0) then
          write(*,*) "Deallocation error"
          stop
       endif


  END subroutine free_allocated_memory

! =============================================================
!
! Set the flag to decide if SSH-aerosol is logging informations
!
! Important: This subroutine must be called before read_namelist
!
! input : true if logging to a file, false (default) otherwise
! =============================================================

    subroutine set_logger(flag)

      implicit none

      logical, intent(in) :: flag

      integer :: ierr
      logical :: log_file_exists

      ssh_logger = flag

      ! Create log file if needed
      if (ssh_logger) then
        ! Check if file exists
        inquire(file = trim(ssh_logger_file), exist = log_file_exists, iostat = ierr)
        if (ierr.ne.0) then
          write(*,*) "SSH-aerosol: error when inquiring log file."
          stop
        endif
        ! Open or create the file
        if (log_file_exists) then
          open(unit = logfile, file = trim(ssh_logger_file), access = "append", status = "old", action = "write", iostat = ierr)
        else
          open(unit = logfile, file = trim(ssh_logger_file), status = "new", iostat = ierr)
        endif
        if (ierr.ne.0) then
          write(*,*) "SSH-aerosol: error when creating / opening log file."
          stop
        endif
      endif

    end subroutine set_logger

! =============================================================
!
! Properly close the log file to finalize the logger
!
! =============================================================

    subroutine close_logger()

      implicit none

      integer :: ierr

      ! Close log file if needed
      if (ssh_logger) then
        close(unit = logfile, iostat = ierr)
        if (ierr.ne.0) then
          write(*,*) "SSH-aerosol: error when closing log file."
          stop
        endif
      endif


    end subroutine close_logger

end module aInitialization
