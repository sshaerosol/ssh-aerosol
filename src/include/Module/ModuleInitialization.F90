!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
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
  Integer, save :: N_gas        !Complete gas species number
  integer, save :: N_size       !Total number of size and composition sections
  integer, save :: N_groups     !Number of groups
  integer, save :: N_fracmax    !Maximum number of composition sections per size section
  integer, save :: N_aerosol    !Number of aerosol species
  integer, save :: N_aerosol_layers    !Number of aerosol species in the different layers
  integer, save :: N_sizebin    !Number of  size sections
  integer, save :: N_frac       !Number of fraction of each species
  integer, save :: N_reaction   !Number of gas-phase reactions
  integer, save :: N_photolysis !Number of photolyses
  integer, save :: Nt           !Number of iteration
  integer, save :: Nmc          

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, save :: N_organics   !Number of organics aerosol species
  integer, save :: N_nonorganics   !Number of non-organics aerosol species
  integer, save :: N_inorganic  !Number of inorganic aerosol species
  integer, save :: N_inert      !number of inert aerosol species
  integer :: N_liquid     !Number of liquid internal species
  integer :: N_solid      !Number of solid internal species
  integer :: N_inside_aer !Number of internal species
  !parameter (N_aerosol = N_organics + N_inorganic + N_inert + 1)
!  parameter (N_organics=32,N_inorganic=5,N_inert=2,N_liquid=12) ! YK
  parameter (N_liquid=12)
  parameter (N_solid=9,N_inside_aer=21)

  !!part 2: parameters of system options

  integer, save :: with_fixed_density!IDENS
  integer, save :: tag_init		! 0 internally mixed; 1 mixing_state resolved
  integer, save :: with_init_num	! 0 estimated from mass and diameter; 1 number conc. for each bin is read
  integer, save :: wet_diam_estimation	! 0 = isorropia ?
  integer, save :: tag_dbd    ! Method for defining particle size bounds (0 auto generated, 1 read)
  integer, save :: tag_emis	     ! 0 Without emissions 1 with internally-mixed emissions 2 with externally-mixed emissions
  integer, save :: with_emis_num ! 0 estimated from mass and diameter; 1 number conc. for each bin is read

  integer, save :: tag_external  ! 0 for internally mixed, 1 for mixing-state resolved
  integer, save :: kind_composition  ! 1 for auto discretization and 0 for manual discretization

  integer, save :: tag_chem
  integer, save :: option_photolysis, day0_photolysis ! current day in photolysis files
  double precision, save :: time_update_photolysis
  integer, save :: with_heterogeneous  !Tag of heterogeneous reaction 

  integer, save :: with_adaptive       !Tag of adaptive time step for chemistry 1 if adaptive time step.
  double precision, save :: adaptive_time_step_tolerance !Relative tolerance for deciding if the time step is kept
  double precision, save :: min_adaptive_time_step       !Minimum time step
  double precision, save :: DTAEROMIN !Minimum time step for aerosol dynamics
  double precision, save :: epser !  Relative error for time step adjustment
  double precision, save :: epser_soap !  Relative difference of ros2 in SOAP
  integer, save :: dynamic_solver = 1 !KDSLV Tag type of solver
  integer, save :: redistribution_method !tag of redistribution method
  integer, save :: with_coag   !Tag gCoagulation
  integer, save :: i_compute_repart ! 0 if repartition coeff are read
  integer, save :: i_write_repart ! 1 if repartition coeff are written, 0 otherwise
  integer, save :: with_cond   !Tag fCondensation
  integer, save :: with_nucl   !Tag nucleation
  Integer, save :: nucl_model  !ITERN !1= Ternary, 0= binary
  integer, save :: ICUT        !ICUT
  double precision, save :: Cut_dim  !cuting diameter between equi/dynamic inorganic
  integer, save :: ISOAPDYN    ! organic equilibrium  = 0 or dynamic = 1
  integer, save :: with_oligomerization!IOLIGO
  integer, save :: output_type
  integer, save :: splitting

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Integer, save :: aqueous_module!ICLD
  Integer, save :: with_incloud_scav!IINCLD
  integer, save :: with_kelvin_effect!IKELV
  integer, save :: section_pass
  double precision, save :: tequilibrium ! time under which equilibrium is assumed

  ! ! part 3: System pointers
  Integer, save :: E1,E2,G1,G2 !Mark the begin and end of dynamic aerosol (except EH2O)
  ! Number of different species group
  Integer, dimension(:), allocatable, save :: isorropia_species
  Integer, dimension(:), allocatable, save :: aec_species
  Integer, save :: nesp, nesp_isorropia, nesp_aec, nesp_eq_org

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Integer :: ENa,ESO4,ENH4,ENO3,ECl,EMD,EBC,EH2O!inorganic pointers
  Integer, save :: ictmNH3,ictmHNO3,ictmHCl,ictmSO2,ictmH2O2,ictmHCHO,ictmHNO2
  Integer, save :: ictmO3,ictmOH,ictmHO2,ictmNO3,ictmNO,ictmNO2,ictmPAN,ictmH2SO4
  ! pointers of cloud species.



  !!part 4: System state parameters    
  ! time setting
  double precision, save :: final_time,dt,time_emis,delta_t, initial_time  
  double precision, save :: current_time
  double precision, save :: Temperature,Relative_Humidity,Pressure,Humidity, pH, pressure_sat
  double precision, save :: longitude, latitude
  double precision, save :: attenuation
  double precision, save :: fixed_density
  double precision, save :: lwc_cloud_threshold

  integer, save :: tag_coag,tag_cond,tag_nucl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  double precision, save :: viscosity!Dynamic viscosity ([kg/m/s]).
  double precision, save :: air_free_mean_path
  double precision, save :: total_water!total mass of water
  double precision, save :: total_IH!total mass of H+
  double precision, save :: total_PH!overall PH value
  double precision, save :: n_emis
  double precision, save :: m_emis
  double precision, save :: aero_total_mass 
  Double precision,dimension(:), allocatable, save :: total_mass!, total_mass_old!total mass of each species
  Double precision,dimension(:), allocatable, save :: discretization_mass
  Double precision, save :: p_fact,k_fact
  Double precision, save :: DQLIMIT

  !!part5: dimension data array    
  integer, dimension(:), allocatable, save :: Index_groups	!index of which group the species belongs to
  integer, dimension(:), allocatable, save :: aerosol_type
  integer, dimension(:), allocatable, save :: List_species	!read species defined in cfg files
  integer, dimension(:,:), allocatable, save :: index_species	!index of species if viscosity is taken into account
  Integer, dimension(:), allocatable, save :: aerosol_species_interact
  integer, dimension(:), allocatable, save :: N_fracbin	!vector of number of composition sections for each section

  Double precision,dimension(:), allocatable, save :: photolysis_rate
  integer, dimension(:), allocatable, save :: photolysis_reaction_index

  Double precision,dimension(:), allocatable, save :: density_aer_bin 	!density of each grid bins
  Double precision,dimension(:), allocatable, save :: density_aer_size 	!density of each size section
  Double precision , dimension(:), allocatable, save :: rho_wet_cell

  Double precision,dimension(:), allocatable, save :: diam_bound	! DBF diameter bounds of each size section
  double precision,dimension(:), allocatable, save :: diam_input
  double precision,dimension(:), allocatable, save :: frac_bound
  double precision,dimension(:), allocatable, save :: frac_input


  Double precision,dimension(:), allocatable, save :: size_diam_av	!DSF average diameter of each size section
  Double precision,dimension(:), allocatable, save :: size_mass_av	!MSF average mass of each size section
  !Double precision,dimension(:), allocatable :: size_log_av	!XSF
  Double precision,dimension(:), allocatable, save :: cell_diam_av	!DSF average diameter of each grid cell
  Double precision,dimension(:), allocatable, save :: cell_mass_av	!MSF average mass of each grid cell
  Double precision,dimension(:), allocatable, save :: cell_log_av	!XSF

  DOUBLE PRECISION, dimension(:), allocatable, save :: concentration_gas_all
  Double precision,dimension(:), allocatable, save :: concentration_gas	! gas concentration of each species
  integer, dimension(:,:), allocatable, save :: concentration_index !matrix from grid index to size and composition index
  integer, dimension(:,:), allocatable, save :: concentration_index_iv !matrix from size and composition to grid index
  Double precision,dimension(:), allocatable, save :: concentration_number	!number concentration of each grid cell
  double precision , dimension(:,:), allocatable, save :: concentration_mass

  double precision, dimension(:), allocatable, save :: gas_emis !storing Gas consentration (emission) micm^3cm^-3 
  double precision,dimension(:,:), allocatable, save :: init_bin_mass
  double precision,dimension(:), allocatable, save :: init_mass
  double precision,dimension(:,:), allocatable, save :: emis_bin_mass
  double precision,dimension(:), allocatable, save :: init_bin_number 
  double precision,dimension(:), allocatable, save :: emis_bin_number

  double precision , dimension(:,:), allocatable, save :: emission_rate
  double precision , dimension(:), allocatable, save :: emission_num_rate

  Double precision,dimension(:), allocatable, save :: wet_diameter	!Aerosol wet diameter (\B5m). of each grid cell
  Double precision,dimension(:), allocatable, save :: wet_mass	!Aerosol wet mass (\B5g). of each grid cell
  Double precision,dimension(:), allocatable, save :: wet_volume	!Aerosol wet volume (\B5m^3). of each grid cell
  double precision , dimension(:,:,:), allocatable, save :: discretization_composition! multi-array storing discretization of composition

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Double precision,dimension(:), allocatable, save :: mass_bound! MBF
  Double precision,dimension(:), allocatable, save :: log_bound!XBF
  Double precision,dimension(:), allocatable, save :: total_bin_mass!total mass of each size section
  Double precision,dimension(:), allocatable, save :: size_sect!HSF log size of each section


  Double precision,dimension(:), allocatable, save :: mass_total_grid!total mass of each grid cell
  Double precision,dimension(:), allocatable, save :: total_aero_mass!total aerosol mass of each species
  Double precision,dimension(:), allocatable, save :: bin_mass!mass concentration of each size section
  Double precision,dimension(:), allocatable, save :: bin_number!number concentration of each size section
  Double precision,dimension(:), allocatable, save :: concentration_number_tmp!first order approximation of number

  Double precision , dimension(:), allocatable, save :: cell_mass
  double precision, dimension(:), allocatable, save :: number_init  !for each sizebin
  !double precision, dimension(:), allocatable :: mass_init    !for each sizebin


  double precision, dimension(:), allocatable, save :: per_mass_init!initial percentage of each species within aerosol



  double precision,dimension(:), allocatable:: gas_mass_init


  double precision , dimension(:,:), allocatable, save :: kernel_coagulation
  double precision , dimension(:,:), allocatable, save :: ce_kernal_coef!c/e kernal
  double precision , dimension(:,:), allocatable, save :: Kelvin_effect_ext!kelvin effect
  double precision , dimension(:,:), allocatable, save :: frac_grid !excat fraction of each species in each grid

  double precision , dimension(:,:), allocatable, save :: concentration_mass_tmp!first order apporximation
  double precision , dimension(:,:), allocatable, save :: concentration_inti!internal inorganic aerosol concentration ([�g.m-3]).
  double precision , dimension(:,:), allocatable, save :: dqdt



  !! part 7: basic physical and chemical parameters
  !double precision :: SMD(SNaNO3:SLC) = ()!molar weight of internal solids species
  !double precision :: IMW(N_liquid)!molar weight of inorganic species in aqueous_phase
  !double precision :: SMW(SNaNO3:SLC)!molar weight of solids
  double precision ,dimension(:), allocatable, save :: accomodation_coefficient
  double precision ,dimension(:), allocatable, save :: surface_tension
  double precision ,dimension(:), allocatable, save :: molecular_weight_aer! (\B5g/mol)
  double precision ,dimension(:), allocatable, save :: molecular_diameter
  double precision ,dimension(:), allocatable, save :: collision_factor_aer
  double precision ,dimension(:), allocatable, save :: mass_density!(\B5g/m3) liquid mass density
  double precision ,dimension(:), allocatable, save :: mass_density_layers!(\B5g/m3) liquid mass density
  double precision ,dimension(:), allocatable, save :: quadratic_speed! (m.s-1)
  double precision ,dimension(:), allocatable, save :: diffusion_coef! (m2.s-1)
  double precision ,dimension(:), allocatable, save :: soa_sat_conc! (\B5g.m-3)
  double precision ,dimension(:), allocatable, save :: soa_part_coef!(m3/microg)
  double precision ,dimension(:), allocatable, save :: molecular_weight! (\B5g/mol) gas=phase
  integer ,dimension(:), allocatable, save :: inon_volatile 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, dimension(:), allocatable, save :: Ncoefficient, index_first, index_second
  double precision, dimension(:), allocatable, save :: coefficient
  integer, save :: coef_size
  double precision, save :: surface_tension_inorg, surface_tension_aq, surface_tension_org

  double precision, save :: dorg
  integer, save :: coupled_phases
  integer, save :: nlayer
  integer ,dimension(:), allocatable, save :: layer_number ! Number of the layer
  double precision, dimension(:), allocatable, save :: Vlayer
  integer, save :: activity_model
  double precision, dimension(:), allocatable, save :: lwc_nsize, &
       ionic_nsize, proton_nsize, chp_nsize
  double precision, dimension(:,:), allocatable, save :: liquid_nsize

  !! part 8: divers parameters (species, I/O)
  character (len=80), save :: Coefficient_file ! repartition coefficient file
  character (len=80), save :: init_aero_conc_mass_file ! File for aeroslos initial mass concentrations
  character (len=80), save :: init_aero_conc_num_file ! File for aerosols initial number concentrations
  character (len=80), save :: init_gas_conc_file ! File for gas-phase initial conc.
  character (len=80), save :: species_list_file ! File for species list.
  character (len=80), save :: aerosol_species_list_file ! File for species list.
  character (len=80), save :: namelist_species ! Namelist file for species list.
  character (len=80), save :: particles_composition_file ! File for particles composition
  character (len=80), save :: emis_gas_file
  character (len=80), save :: emis_aero_mass_file
  character (len=80), save :: emis_aero_num_file
  character (len=80), dimension(:), allocatable, save :: isorropia_species_name
  character (len=80), dimension(:), allocatable, save :: aec_species_name
  character (len=10), save :: precursor
  character (len=10), dimension(:), allocatable, save :: species_name
  character (len=20), dimension(:), allocatable, save :: aerosol_species_name
  integer, save :: spec_name_len
  character (len=10), dimension(:), allocatable, save :: emis_gas_species_name
  character (len=10), dimension(:), allocatable, save :: emis_aer_species_name
  character (len=100), save :: output_directory, output_dir2

  ! Photolysis
  character (len=80), save :: photolysis_file ! File for photolysis list.
  character (len=80), save :: photolysis_dir ! Directory for photolysis list.
  character (len=10), dimension(:), allocatable, save :: photolysis_name 
  integer, save :: n_time_angle
  double precision, save :: time_angle_min, delta_time_angle
  integer, save :: n_latitude
  double precision, save :: latitude_min, delta_latitude
  integer, save :: n_altitude
  double precision, save :: altitude_photolysis_input(30)

  !!part 6: used in ssh-aerosol.f90 chem()
  integer, save :: ns_source
  integer, dimension(:), allocatable, save :: source_index
  double precision, dimension(:), allocatable, save :: source  
  double precision, dimension(:), allocatable, save :: conversionfactor
  double precision, dimension(:,:), allocatable, save :: conversionfactorjacobian
  ! Array of chemical volumic emissions at final time ([\mu.g/m^3/s]).
  integer, save :: heterogeneous_reaction_index(4)
  integer, save :: ind_jbiper, ind_kbiper   
  ! To take into account BiPER degradation inside the particle

  !! part 7: used in coupling with external tools (ModuleAPI)
  !
  ! This flag defines if SSH-aerosol is running standalone or not
  !   true if it is running standalone, default
  !   false if it is running with an external tool
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
    integer :: i,ierr,tag_file, nml_out
    character (len=40), intent(in) :: namelist_file
    character (len=40) :: namelist_out

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

    namelist /aerosol_species/ aerosol_species_list_file

    namelist /physic_gas_chemistry/ tag_chem, attenuation, option_photolysis, &
         time_update_photolysis, & 
         with_heterogeneous, with_adaptive, &
         adaptive_time_step_tolerance, min_adaptive_time_step, &
         photolysis_dir, photolysis_file, &
         n_time_angle, time_angle_min, delta_time_angle, &
         n_latitude, latitude_min, delta_latitude, &
         n_altitude, altitude_photolysis_input

    namelist /physic_particle_numerical_issues/ DTAEROMIN, redistribution_method,&
         with_fixed_density, fixed_density, splitting

    namelist /physic_coagulation/ with_coag, i_compute_repart, i_write_repart, Coefficient_file, Nmc

    namelist /physic_condensation/ with_cond, Cut_dim, ISOAPDYN, nlayer,&
         with_kelvin_effect, tequilibrium,&
         dorg, coupled_phases, activity_model, epser, epser_soap

    namelist /physic_nucleation/ with_nucl, nucl_model

    namelist /physic_organic/ with_oligomerization

    namelist /output/ output_directory, output_type, particles_composition_file

    if (ssh_standalone) write(*,*) "=========================start read namelist.ssh file======================"
    if (ssh_logger) write(logfile,*) "=========================start read namelist.ssh file======================"
    ! read namelist.ssh file !
    open(unit = 10, file = namelist_file, status = "old")

    ! Use default values if they are not given in namelist.ssh
    ! And write to namelist.out
    nml_out = 101
    namelist_out = "namelist.out"
    open(nml_out, file = namelist_out)

    ! meteorological setup
    read(10, nml = setup_meteo, iostat = ierr)
    if (ierr .ne. 0) then
       write(*,*) "setup_meteo data can not be read."
       stop
    else ! default output meteo data to check

       if  (Relative_Humidity .gt. 0d0 ) then
          if (      Relative_Humidity.lt.Threshold_RH_inf &
               .or. Relative_Humidity.gt.Threshold_RH_sup) then
             if (ssh_standalone) write(*,*) 'Warning : clipping relative humidity.'
             if (ssh_logger) write(logfile,*) 'Warning : clipping relative humidity.'
             Relative_Humidity = DMIN1(DMAX1(Relative_Humidity, Threshold_RH_inf), Threshold_RH_sup)
          endif
          call compute_psat_sh(Relative_Humidity, temperature, Pressure, pressure_sat, humidity)
       else
          call compute_psat_rh(humidity, temperature, Pressure, pressure_sat, Relative_Humidity)
          if (      Relative_Humidity.lt.Threshold_RH_inf &
               .or. Relative_Humidity.gt.Threshold_RH_sup) then
             if (ssh_standalone) write(*,*) 'Warning : clipping relative humidity.'
             if (ssh_logger) write(logfile,*) 'Warning : clipping relative humidity.'
             Relative_Humidity = DMIN1(DMAX1(Relative_Humidity, Threshold_RH_inf), Threshold_RH_sup)
             call compute_psat_sh(Relative_Humidity, temperature, Pressure, pressure_sat, humidity)
          endif
       end if
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
    wet_diam_estimation = -999
    read(10, nml = initial_condition, iostat = ierr)

    if (ierr .ne. 0) then
       write(*,*) "initial_condition data can not be read."
       stop
    else
       ! wet_diam_estimation = 1 by default
       ! if it is not given in namelist
       if (wet_diam_estimation == -999) then
          wet_diam_estimation = 1
       endif
       
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
       write(*,*) "mixing_state data can not be read.",N_frac
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
    attenuation = -999.d0
    option_photolysis = -999
    time_update_photolysis = -999.d0
    with_adaptive = -999
    adaptive_time_step_tolerance = -999.d0
    min_adaptive_time_step = -999.d0
    photolysis_dir = "---"
    photolysis_file = "---"
    n_time_angle = -999
    time_angle_min = -999.d0
    delta_time_angle = -999.d0
    n_latitude = -999
    latitude_min = -999.d0
    delta_latitude = -999.d0
    n_altitude = -999
    altitude_photolysis_input = -999.d0
    read(10, nml = physic_gas_chemistry, iostat = ierr)
    if (ierr .ne. 0) then
       write(*,*) "physic_gas_chemistry data can not be read."
       stop
    else
       ! attenuation = 1.d0 by default
       ! if it is not given in namelist
       if (attenuation == -999.d0) then
          attenuation = 1.d0
       endif
       ! option_photolysis = 1 by default
       ! if it is not given in namelist
       if (option_photolysis == -999) then
          option_photolysis = 1
       endif
       ! time_update_photolysis = 100000.d0 by default
       ! if it is not given in namelist
       if (time_update_photolysis == -999.d0) then
          time_update_photolysis = 100000.d0
       endif
       ! with_adaptive = 1 by default
       ! if it is not given in namelist
       if (with_adaptive == -999) then
          with_adaptive = 1
       endif
       ! adaptive_time_step_tolerance = 0.001 by default
       ! if it is not given in namelist
       if (adaptive_time_step_tolerance == -999.d0) then
          adaptive_time_step_tolerance = 0.001
       endif
       ! min_adaptive_time_step = 0.001 by default
       ! if it is not given in namelist
       if (min_adaptive_time_step == -999.d0) then
          min_adaptive_time_step = 0.001
       endif
       ! photolysis_dir = "./photolysis/" by default
       ! if it is not given in namelist
       if (trim(photolysis_dir) == "---") then
          photolysis_dir = "./photolysis/"
       endif
       ! photolysis_file = "./photolysis/" by default
       ! if it is not given in namelist
       if (trim(photolysis_file) == "---") then
          photolysis_file = "./photolysis/photolysis-cb05.dat"
       endif       
       ! n_time_angle = 9 by default
       ! if it is not given in namelist
       if (n_time_angle == -999) then
          n_time_angle = 9
       endif
       ! time_angle_min = 0.d0 by default
       ! if it is not given in namelist
       if (time_angle_min == -999.d0) then
          time_angle_min = 0.d0
       endif
       ! delta_time_angle = 1.d0 by default
       ! if it is not given in namelist
       if (delta_time_angle == -999.d0) then
          delta_time_angle = 1.d0
       endif
       ! n_latitude = 10 by default
       ! if it is not given in namelist
       if (n_latitude == -999) then
          n_latitude = 10
       endif
       ! latitude_min = 0.d0 by default
       ! if it is not given in namelist
       if (latitude_min == -999.d0) then
          latitude_min = 0.d0
       endif
       ! delta_latitude = 10.d0 by default
       ! if it is not given in namelist
       if (delta_latitude == -999.d0) then
          delta_latitude = 10.d0
       endif
       ! n_altitude = 9 by default
       ! if it is not given in namelist
       if (n_altitude == -999) then
          n_altitude = 9
       endif
       ! altitude_photolysis_input by default
       ! if it is not given in namelist
       if (altitude_photolysis_input(1) == -999.d0) then
          altitude_photolysis_input(1:9) = [0.0, 1000.0, 2000.0, 3000.0, &
               4000.0, 5000.0, 10000.0, 15000.0, 20000.0]
       endif       
       
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

       if (ssh_standalone) write(*,*) 'Cloud attenuation field', attenuation
       if (ssh_logger) write(logfile,*) 'Cloud attenuation field', attenuation

    end if

    ! particle numerical issues
    dtaeromin = -999.d0
    with_fixed_density = -999
    fixed_density = -999.d0
    splitting = -999
    read(10, nml = physic_particle_numerical_issues, iostat = ierr)
    if (ierr .ne. 0) then
       write(*,*) "physic_particle_numerical_issues data can not be read."
       stop
    else

       ! dtaerominn = 1.d-5 by default
       ! if it is not given in namelist
       if (dtaeromin == -999.d0) then
          dtaeromin = 1.d-5
       endif
       ! with_fixed_density = 0 by default
       ! if it is not given in namelist
       if (with_fixed_density == -999) then
          with_fixed_density = 0
       endif
       ! fixed_density = 1.84d-06 by default
       ! if it is not given in namelist
       if (fixed_density == -999.d0) then
          fixed_density = 1.84d-06
       endif
       ! splitting = 1 by default
       ! if it is not given in namelist
       if (splitting == -999) then
          splitting = 1
       endif

       
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
    nmc = -999
    read(10, nml = physic_coagulation, iostat = ierr)
    if (ierr .ne. 0) then
       write(*,*) "physic_coagulation data can not be read."
       stop
    else 
       ! nmc = 1000000 by default
       ! if it is not given in namelist
       if (nmc == -999) then
          nmc = 1000000
       endif

       if (with_coag == 1) then
          if (ssh_standalone) write(*,*) '! ! ! with coagulation.'
          if (ssh_logger) write(logfile,*) '! ! ! with coagulation.'
          if(i_compute_repart == 0) then
             if (ssh_standalone) write(*,*) i_compute_repart,'repartition coefficient are read'
             if (ssh_logger) write(logfile,*) i_compute_repart,'repartition coefficient are read'
	     if (ssh_standalone) write(*,*) 'coefficient file : ', Coefficient_file
	     if (ssh_logger) write(logfile,*) 'coefficient file : ', Coefficient_file
             i_write_repart = 0 ! Do no write repartition coefficients if they are not computed.
          else
             if (ssh_standalone) write(*,*) i_compute_repart,'repartition coefficient are computed'
             if (ssh_logger) write(logfile,*) i_compute_repart,'repartition coefficient are computed'
	     if (ssh_standalone) write(*,*) 'Nmc = ',Nmc
	     if (ssh_logger) write(logfile,*) 'Nmc = ',Nmc
          endif
          if (i_write_repart == 1) then
             if (ssh_standalone) write(*,*) i_write_repart,'repartition coefficient are written'
             if (ssh_logger) write(logfile,*) i_write_repart,'repartition coefficient are written'
             if (ssh_standalone) write(*,*) 'coefficient file : ', Coefficient_file
             if (ssh_logger) write(logfile,*) 'coefficient file : ', Coefficient_file
          endif
          ! Safety check
          if (i_compute_repart .lt. 0 .or. i_compute_repart .gt. 1) then
             write(*,*) "Wrong value for the parameter i_compute_repart. Abort."
             stop
          endif
          if (i_write_repart .lt. 0 .or. i_write_repart .gt. 1) then
             write(*,*) "Wrong value for the parameter i_write_repart. Abort."
             stop
          endif
       else
          with_coag = 0
          if (ssh_standalone) write(*,*) 'without coagulation.'
          if (ssh_logger) write(logfile,*) 'without coagulation.'
       end if
    end if

    ! condensation/ evaporation
    nlayer = -999
    with_kelvin_effect = -999
    tequilibrium = -999.d0
    dorg = -999.d0
    coupled_phases = -999
    epser = -999.d0
    epser_soap = -999.d0
    read(10, nml = physic_condensation, iostat = ierr)
    if (ierr .ne. 0) then
       write(*,*) "physic_condensation data can not be read."
       stop
    else
       ! nlayer = 1 by default
       ! if it is not given in namelist
       if (nlayer == -999) then
          nlayer = 1
       endif
       ! with_kelvin_effect = 1 by default
       ! if it is not given in namelist
       if (with_kelvin_effect == -999) then
          with_kelvin_effect = 1
       endif
       ! tequilibrium = 0.1d0 by default
       ! if it is not given in namelist
       if (tequilibrium == -999.d0) then
          tequilibrium = 0.1d0
       endif
       ! dorg = 1.d-12 by default
       ! if it is not given in namelist
       if (dorg == -999.d0) then
          dorg = 1.d-12
       endif       
       ! coupled_phases = 0 by default
       ! if it is not given in namelist
       if (coupled_phases == -999) then
          coupled_phases = 0
       endif
       ! epser = 0.01 by default
       ! if it is not given in namelist
       if (epser == -999.d0) then
          epser = 0.01
       endif       
       ! epser_soap = 0.01 by default
       ! if it is not given in namelist
       if (epser_soap == -999.d0) then
          epser_soap = 0.01
       endif

       
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

    if((with_coag == 1).AND.nlayer > 1) then
       print*, 'coagulation can not be switched on if several aerosol layers.'
       stop
    endif

    if((ISOAPDYN == 0).AND.nlayer > 1) then
       print*, 'Organics can not be at equilibrium if several aerosol layers.'
       stop
    endif

    if ((with_coag == 1).AND.(with_cond == 0)) then
          wet_diam_estimation = 0 ! Compute water content during coagulation, as it is not computed during condensation.
    endif
    if((with_cond == 1).AND.(Cut_dim.GE.diam_input(N_sizebin+1))) then
            wet_diam_estimation = 0 ! Compute water content if full equilibrium
    endif
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
          output_type = 1   ! default
          if (ssh_standalone) write(*,*) 'results are saved in text files.'
          if (ssh_logger) write(logfile,*) 'results are saved in text files.'
       end if
       if (ssh_standalone) write(*,*) 'output directory :', output_directory
       if (ssh_logger) write(logfile,*) 'output directory :', output_directory
       if (ssh_standalone) write(*,*) 'Particles composition file : ', particles_composition_file
       if (ssh_logger) write(logfile,*) 'Particles composition file : ', particles_composition_file
    end if

    close(10)

    ! Write to namelist.out
    write(nml_out, setup_meteo)
    write(nml_out, setup_time)
    write(nml_out, initial_condition)
    write(nml_out, initial_diam_distribution)
    write(nml_out, emissions)
    write(nml_out, mixing_state)
    write(nml_out, fraction_distribution)
    write(nml_out, gas_phase_species)
    write(nml_out, aerosol_species)
    write(nml_out, physic_gas_chemistry)
    write(nml_out, physic_particle_numerical_issues)
    write(nml_out, physic_coagulation)
    write(nml_out, physic_condensation)
    write(nml_out, physic_nucleation)
    write(nml_out, physic_organic)
    write(nml_out, output)
    close(nml_out)

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
    integer :: k,i,j,s,js, ind, count, ierr, ilayer, esp_layer, nline
    double precision :: tmp
    double precision, dimension(:), allocatable :: tmp_aero
    character (len=40) :: ic_name, sname, tmp_name

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
    allocate(aerosol_type(N_aerosol))
    allocate(index_species(N_aerosol,nlayer))
    ! initialize basic physical and chemical parameters
    allocate(molecular_weight_aer(N_aerosol))
    allocate(collision_factor_aer(N_aerosol))
    allocate(molecular_diameter(N_aerosol)) 
    allocate(surface_tension(N_aerosol))
    allocate(accomodation_coefficient(N_aerosol))
    allocate(mass_density(N_aerosol))
    allocate(inon_volatile(N_aerosol))
    allocate(Vlayer(nlayer))
    ! relation between Aerosol and GAS
    allocate(aerosol_species_interact(N_aerosol))      
    aerosol_species_interact = 0
    inon_volatile = 0

    ! Read lines from aerosol species file.
    rewind 12
    count = 0
    read(12, *) ! Read a header line (#)
    do s = 1, N_aerosol
       ! Surface_tension for organic and aqueous phases of organic aerosols
       ! is hardly coded in SOAP/parameters.cxx
       ! And Unit used in SOAP is different to surface_tension (N/m) by 1.e3.
       read(12, *) aerosol_species_name(s), aerosol_type(s), &
            Index_groups(s), molecular_weight_aer(s), &
            precursor, &
            collision_factor_aer(s), molecular_diameter(s), &
            surface_tension(s), accomodation_coefficient(s), &
            mass_density(s), inon_volatile(s)
        if((inon_volatile(s).NE.1).AND.(inon_volatile(s).NE.0)) then
            write(*,*) "non_volatile should be 0 or 1", inon_volatile(s),s
            stop
        endif
       ! Find pairs of aerosol species and its precursor.
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
    ! Safety check if index_groups is used
    if (tag_external.eq.1) then
       if (minval(Index_groups(1:N_aerosol_Layers-1)).lt.1) then
          if (ssh_standalone) write(*,*) "Error: Incorrect group index in aerosol_species_list_file."    
          if (ssh_logger) write(logfile,*) "Error: Incorrect group index in aerosol_species_list_file."  
          stop
       endif
       if (maxval(Index_groups(1:N_aerosol_Layers-1)).gt.N_groups) then
          if (ssh_standalone) write(*,*) "Error: Increase N_groups in namelist file."
          if (ssh_logger) write(logfile,*) "Error: Increase N_groups in namelist file."
          stop
       endif
    endif
    if (ssh_standalone) write(*,*) "   --- Number of precursors: ", count 
    if (ssh_logger) write(logfile,*) "   --- Number of precursors: ", count 
    
    ! Count the number of species for each type.
    N_inert = 0
    nesp_isorropia = 0
    nesp_aec = 0
    do s = 1, N_aerosol
       ! Inert aerosols: BC and mineral dust
       if (aerosol_type(s) == 1 .or. aerosol_type(s) == 2) then
          N_inert = N_inert + 1
       ! Inorganic species   
       else if (aerosol_type(s) == 3) then
          nesp_isorropia = nesp_isorropia + 1
       ! Organic species   
       else if (aerosol_type(s) == 4) then
          nesp_aec = nesp_aec + 1
       ! Water   
       else if (aerosol_type(s) == 9) then
          EH2O = s
       end if
    end do
    N_inorganic = nesp_isorropia
    N_organics = nesp_aec
    nesp_eq_org = N_organics
    
    if (ssh_standalone) write(*,*) "   --- Number of inert species:", N_inert
    if (ssh_logger) write(logfile,*) "   --- Number of inert species:", N_inert
    if (ssh_standalone) write(*,*) "   --- Number of inorganic species:", N_inorganic
    if (ssh_logger) write(logfile,*) "   --- Number of inorganic species:", N_inorganic
    if (ssh_standalone) write(*,*) "   --- Number of organic species:", N_organics
    if (ssh_logger) write(logfile,*) "   --- Number of organic species:", N_organics
    if (ssh_standalone) write(*,*) "   --- Index for water:", EH2O
    if (ssh_logger) write(logfile,*) "   --- Index for water:", EH2O

    ! Allocate aerosol arrays
    N_nonorganics = N_aerosol - N_organics -1 ! Remove organics and water
    N_aerosol_layers = N_organics * (nlayer-1) + N_aerosol
    EH2O_layers = N_aerosol_layers
    allocate(mass_density_layers(N_aerosol_layers))
    allocate(List_species(N_aerosol_layers))
    allocate(layer_number(N_aerosol_layers))
    allocate(isorropia_species(nesp_isorropia))
    allocate(isorropia_species_name(nesp_isorropia))
    allocate(aec_species(nesp_aec))
    allocate(aec_species_name(nesp_aec))

    ! Read aerosol species name.
    js = 0
    i = 0
    do s = 1, N_aerosol
       if (aerosol_type(s) == 3) then
          i = i + 1
          isorropia_species_name(i) = aerosol_species_name(s)
          isorropia_species(i) = s
       else if (aerosol_type(s) == 4) then
          js = js + 1
          aec_species_name(js) = aerosol_species_name(s)
          aec_species(js) = s
       endif
    end do
    
    do s = 1, N_aerosol
       ! For non-organic species.
       if (s <= N_nonorganics) then
          mass_density_layers(s) = mass_density(s)
          molecular_weight_aer(s) = molecular_weight_aer(s) * 1.0D06 ! g/mol to \B5g/mol  !!! change later
          List_species(s) = s
          if (aerosol_species_name(s) .eq. "PMD") EMD = s
          if (aerosol_species_name(s) .eq. "PBC") EBC = s
          if (aerosol_species_name(s) .eq. "PNA") ENa = s
          if (aerosol_species_name(s) .eq. "PSO4") ESO4 = s
          if (aerosol_species_name(s) .eq. "PNH4") ENH4 = s
          if (aerosol_species_name(s) .eq. "PNO3") ENO3 = s
          if (aerosol_species_name(s) .eq. "PHCL") ECl = s
          if (aerosol_species_name(s) .eq. "PBiPER") ind_jbiper = s

          do ilayer=1,nlayer
             index_species(s,ilayer) = s
          enddo
          layer_number(s) = 1
       ! For organic species
       else
          molecular_weight_aer(s) = molecular_weight_aer(s) * 1.0D06 ! g/mol to \B5g/mol  !!! change later
          if(s.NE.N_aerosol) then !avoid water
             do ilayer = 0,nlayer-1
                esp_layer = (s-N_nonorganics-1) *(nlayer-1) + s + ilayer
                index_species(s,ilayer+1) = esp_layer
                mass_density_layers(esp_layer) = mass_density(s)
                List_species(esp_layer) = s
                !               molecular_weight_aer(esp_layer) = molecular_weight_aer(s)
                !!aerosol_species_name(esp_layer) = aerosol_species_name(s)
                !Index_groups(esp_layer) = Index_groups(s)
                !               mass_density(esp_layer) = mass_density(s)
                layer_number(esp_layer) = ilayer
             enddo
          else
             List_species(N_aerosol_layers) = s
             do ilayer=1,nlayer
                index_species(N_aerosol,ilayer) = N_aerosol_layers 
             enddo
             layer_number(N_aerosol_layers) = 1
          endif
       endif
    enddo

    if(with_nucl.EQ.1) then
          inon_volatile(ESO4) = 1 ! sulfate needs to be computed dynamically in case of nucleation
    endif
    ! read gas-phase initial concentrations unit 21
    ! no comment lines for initial & emitted data
    allocate(concentration_gas_all(N_gas))
    concentration_gas_all = 0.d0 ! set original value to 0
    allocate(concentration_gas(N_aerosol))
    concentration_gas=0.d0

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

    ! Initialize Vlayer  !! Need to be removed from SOAP
    if(nlayer == 1) then
       Vlayer(1)=1.0
    else 
       if(nlayer == 2) then
          Vlayer(1)=0.99
          Vlayer(2)=0.01
       else 
          if(nlayer == 3) then
             Vlayer(1)=0.6
             Vlayer(2)=0.39
             Vlayer(3)=0.01
          else 
             if(nlayer == 4) then
                Vlayer(1)=0.6
                Vlayer(2)=0.26
                Vlayer(3)=0.13
                Vlayer(4)=0.01
             else 
                if (nlayer == 5) then
                   Vlayer(1)=0.608
                   Vlayer(2)=0.2184165
                   Vlayer(3)=0.12102374
                   Vlayer(4)=0.04255976
                   Vlayer(5)=0.01
                else
                   if (ssh_standalone) write(*,*) "Number of layers not implemented in ssh"
                   if (ssh_logger) write(logfile,*) "Number of layers not implemented in ssh"
                endif
             endif
          endif
       endif
    endif
   
    if(nlayer > 1) then ! Consider non-volatile species in SOAP
        do i = 1, N_aerosol
           if(i.NE.ESO4) inon_volatile(i) = 0
        enddo
    endif

    ! read input file for photolysis rate (unit 34)
    allocate(photolysis_name(n_photolysis))
    allocate(photolysis_reaction_index(n_photolysis))
    allocate(photolysis_rate(n_photolysis))
    photolysis_rate = 0.d0
     
    open(unit = 34, file = photolysis_file, status = "old")
    count = 0 
    ierr = 0
    nline = 0
    do while(ierr .eq. 0)
       read(34, *, iostat=ierr) tmp_name
       if (ierr == 0) then
          if (trim(tmp_name) .eq. "#") then
             nline = nline + 1
          else
             count = count + 1
          end if
       end if
    end do

    if (count .ne. n_photolysis) then
       write(*,*) "Error: number of files for photolysis rate should be ", &
            n_photolysis
       write(*,*) "However, the number of given files is ", count
       stop
    end if
       
    rewind 34
    do s = 1, nline
       read(34, *) ! read comment lines.
    end do
    do s = 1, count
       read(34, *) photolysis_name(s), photolysis_reaction_index(s)
       if (photolysis_name(s) == "BiPER") then
          ind_kbiper = s    ! photolysis index for BiPER
       end if
    end do
    close(34)
    
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
    if (allocated(aerosol_type))  deallocate(aerosol_type, stat=ierr) !! YK
    if (allocated(index_species))  deallocate(index_species, stat=ierr)
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
    if (allocated(mass_density_layers))  deallocate(mass_density_layers, stat=ierr)
    if (allocated(inon_volatile))  deallocate(inon_volatile, stat=ierr)
    if (allocated(layer_number))  deallocate(layer_number, stat=ierr)
    if (allocated(Vlayer))  deallocate(Vlayer, stat=ierr)
    !!	if (allocated(saturation_pressure))  deallocate(saturation_pressure, stat=ierr)
    !!	if (allocated(vaporization_enthalpy))  deallocate(vaporization_enthalpy, stat=ierr)
    if (allocated(List_species))  deallocate(List_species, stat=ierr)
    if (allocated(isorropia_species))  deallocate(isorropia_species, stat=ierr)
    if (allocated(isorropia_species_name))  deallocate(isorropia_species_name, stat=ierr)
    if (allocated(aec_species))  deallocate(aec_species, stat=ierr)
    if (allocated(aec_species_name))  deallocate(aec_species_name, stat=ierr)
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

    if (allocated(photolysis_rate))  deallocate(photolysis_rate, stat=ierr)
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

    ! for photolysis
    if (allocated(photolysis_name)) deallocate(photolysis_name, stat=ierr)
    if (allocated(photolysis_reaction_index)) deallocate(photolysis_reaction_index, stat=ierr)    
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
