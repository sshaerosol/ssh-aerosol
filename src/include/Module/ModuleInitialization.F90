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
  integer, save :: tag_init             ! 0 internally mixed; 1 mixing_state resolved
  integer, save :: with_init_num        ! 0 estimated from mass and diameter; 1 number conc. for each bin is read
  integer, save :: wet_diam_estimation  ! 0 = isorropia ?
  integer, save :: tag_dbd    ! Method for defining particle size bounds (0 auto generated, 1 read)
  integer, save :: tag_emis             ! 0 Without emissions 1 with internally-mixed emissions 2 with externally-mixed emissions
  integer, save :: with_emis_num ! 0 estimated from mass and diameter; 1 number conc. for each bin is read

  integer, save :: tag_external  ! 0 for internally mixed, 1 for mixing-state resolved
  integer, save :: kind_composition  ! 1 for auto discretization and 0 for manual discretization

  integer, save :: tag_chem
  integer, save :: tag_twostep ! 1 for use the two-step time numerical solver to solve chemistry 
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
  integer, save :: nucl_model_binary, nucl_model_ternary, nucl_model_hetero
  integer, save :: nucl_model_org
  double precision, save :: scal_ternary, scal_hetero, scal_org, nexp_org
  double precision, save :: co2_conc_ppm
  integer, save :: nesp_org_h2so4_nucl, nesp_org_nucl
  integer, dimension(2), save :: org_h2so4_nucl_species
  character (len=40), dimension(2) :: name_org_h2so4_nucl_species
  integer, dimension(10), save :: org_nucl_species
  character (len=40), dimension(10) :: name_org_nucl_species
  integer, save :: ISOAPDYN    ! organic equilibrium  = 0 or dynamic = 1
  integer, save :: IMETHOD     ! numerical method for SOAP, 0= explicit, 1= implicit, 2=implicit semi-dynamic
  integer, save :: SOAPlog     ! 0=no text output from SOAP, 1=on screen, 2=written on files
  integer, save :: output_type
  integer, save :: splitting
  integer, save :: soap_inorg,soap_inorg_loc
  integer, save :: niter_water,niter_eqconc
  integer, save :: NACL_IN_THERMODYNAMICS
  integer, dimension(:), allocatable, save :: iter_eqconc,iter_water
  integer, save :: option_cloud
  
  ! cut-off diameter
  integer, save :: ICUT,ICUT_org      ! section index
  Integer, save :: set_icut  ! 0 = need to update ICUT, 1 = keep ICUT unchanged after initialization
  Integer, save :: tag_icut  ! 0 = fixed ICUT, icut is computed in the program using 1 = c/e timescale criteria, 2 = ETR criteria, 3 = QSSA criteria
  integer, save :: cond_time_index(3) ! store the index of ENO3, ECL and ENH4 for icut computation
  double precision, save :: Cut_dim   ! value of the user-chosen parameter. Depending on tag_icut.
  double precision, save :: kwall_gas, kwall_particle, Cwall, surface_volume_ratio, eddy_turbulence, kwp0,radius_chamber

  character (len=200), save :: reaction_soap_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! genoa related parameters
  ! tag for genoa, tag to use keep_gp in twostep
  integer, save :: tag_genoa = 0, keep_gp
  ! tag for RO2 treatment
  integer, save :: tag_RO2, nref_file ! 0 for without RO2, 1 for simulated with generated RO2 only, 2 for with background RO2 only, 3 for with background + generated RO2
  
  ! numbers that needs in the error computation
  integer, save :: nout_err   ! individual species w err compute
  integer, save :: nout_soa   ! a group of species: total concs w err compute
  integer, save :: nout_total ! total num for err computation! default total SOAs
  ! number of output species
  integer, save :: nout_gas, nout_aero
  ! index of initial set - can be updated in ssh-aerosol.f90
  integer, save :: tag_init_set
  ! number of species with constant profiles, RO2 species
  integer, save :: ncst_gas, ncst_aero, nRO2_chem, nRO2_group
  ! RO2 pool index in gas species list
  integer, save :: iRO2_cst
  
  ! index to find corresponding species
  integer, dimension(:), allocatable, save :: output_gas_index, output_aero_index
  integer, dimension(:), allocatable, save :: cst_gas_index, cst_aero_index, RO2out_index
  integer, dimension(:,:), allocatable, save :: RO2index
  integer, dimension(:,:), allocatable, save :: output_err_index ! (gas/aero, index)

  ! concentrations array
  double precision, dimension(:,:), allocatable, save :: cst_gas, total_soa  ! constant concs.
  double precision, dimension(:,:,:), allocatable, save :: ref_soa  ! err compute
  double precision, dimension(:,:,:), allocatable, save :: cst_aero ! constant concs.

  ! file names
  character (len=800), save :: ref_conc_files_in ! reference cases
  character (len=200), save :: RO2_list_file ! File for RO2 species to generate the RO2 pool
  character (len=200), save :: cst_gas_file ! File for species that need to be keep as constants.
  character (len=200), save :: cst_aero_file ! File for constant aerosol species.
  character (len=200), save :: init_species_file ! File for initial sets of SOA precursors
  character (len=200), dimension(:), allocatable, save :: ref_conc_files
  
  ! list of species name
  character (len=40), dimension(:), allocatable, save :: output_err_sps
  character (len=40), dimension(:), allocatable, save :: output_gas_species, output_aero_species


!!!! KINETIC RELATED - GENOA
  double precision, save :: SumMc ! M (in molec/cm3) - constant or change by pressure and temperature
  double precision, save :: YlH2O ! water concentration (molec/cm3) in the gas phase
  integer, dimension(:,:), allocatable, save :: photo_rcn, TB_rcn
  integer, dimension(:), allocatable, save :: fall_rcn, extra_rcn, hetero_rcn
  integer, dimension(:,:), allocatable, save :: index_RCT, index_PDT, index_PDT_ratio
  double precision, dimension(:), allocatable, save :: ratio_PDT
  ! kinetics
  double precision, dimension(:,:), allocatable, save :: Arrhenius,fall_coeff,extra_coeff  ! coefficients for spec reactions
  integer, dimension(:), allocatable, save :: hetero_ind, irdi_ind
  double precision, dimension(:), allocatable, save :: kinetic_rate, chem_prod, chem_loss
  double precision, dimension(:), allocatable, save :: rcn_rate, gas_yield, drv_knt_rate
  ! photolysis
  double precision, dimension(:,:,:), allocatable, save :: photo_ratio
  double precision, dimension(:,:), allocatable, save :: photo_ratio_read
  ! double precision, parameter :: szas(11) = (/0d0,1d1,2d1,3d1,4d1, &
  !                                       5d1,6d1,7d1,7.8d1,8.6d1,9d1/)
  ! integer, parameter :: nsza = 11! sza size 
  integer :: nsza
  double precision, dimension(:), allocatable, save :: szas

  
  ! Wall loss
  integer, dimension(:), allocatable, save :: wall_rcn  
  double precision, dimension(:,:), allocatable, save :: wall_coeff

  ! Irreversible dicarbonyl
  integer, dimension(:), allocatable, save :: irdi_rcn  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Integer, save :: aqueous_module!ICLD
  Integer, save :: with_incloud_scav!IINCLD
  integer, save :: with_kelvin_effect!IKELV
  integer, save :: section_pass
  double precision, save :: tequilibrium ! time under which equilibrium is assumed

  ! Number of different species group
  Integer, dimension(:), allocatable, save :: isorropia_species
  Integer, dimension(:), allocatable, save :: aec_species
  Integer, save :: nesp, nesp_isorropia, nesp_aec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  double precision, save :: cloud_water

  double precision, dimension(:), allocatable, save :: temperature_array
  double precision, dimension(:), allocatable, save :: humidity_array
  double precision, dimension(:), allocatable, save :: pressure_array
  double precision, dimension(:), allocatable, save :: relative_humidity_array

  double precision, dimension(:), allocatable, save :: ratio_water
  double precision, dimension(:,:), allocatable, save :: ratio_eqconc
  
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
  integer, dimension(:), allocatable, save :: Index_groups      !index of which group the species belongs to
  integer, dimension(:), allocatable, save :: aerosol_type
  integer, dimension(:), allocatable, save :: List_species      !read species defined in cfg files
  integer, dimension(:,:), allocatable, save :: index_species   !index of species if viscosity is taken into account
  Integer, dimension(:), allocatable, save :: aerosol_species_interact
  integer, dimension(:), allocatable, save :: N_fracbin !vector of number of composition sections for each section
  integer, dimension(:), allocatable, save :: oligo_index

  Double precision,dimension(:), allocatable, save :: photolysis_rate
  integer, dimension(:), allocatable, save :: photolysis_reaction_index

  Double precision,dimension(:), allocatable, save :: density_aer_bin  !density of each grid bins
  Double precision,dimension(:), allocatable, save :: density_aer_size !density of each size section
  Double precision , dimension(:), allocatable, save :: rho_wet_cell

  Double precision,dimension(:), allocatable, save :: diam_bound       ! DBF diameter bounds of each size section
  double precision,dimension(:), allocatable, save :: diam_input
  double precision,dimension(:), allocatable, save :: frac_bound
  double precision,dimension(:), allocatable, save :: frac_input


  Double precision,dimension(:), allocatable, save :: size_diam_av      !DSF average diameter of each size section
  Double precision,dimension(:), allocatable, save :: size_mass_av      !MSF average mass of each size section
  !Double precision,dimension(:), allocatable :: size_log_av	!XSF
  Double precision,dimension(:), allocatable, save :: cell_diam_av      !DSF average diameter of each grid cell
  Double precision,dimension(:), allocatable, save :: cell_log_av       !XSF

  DOUBLE PRECISION, dimension(:), allocatable, save :: concentration_gas_all
  Double precision,dimension(:), allocatable, save :: concentration_gas	! gas concentration of each species
  Double precision,dimension(:), allocatable, save :: concentration_wall	! wall concentration of each species
  integer, dimension(:,:), allocatable, save :: concentration_index !matrix from grid index to size and composition index
  integer, dimension(:,:), allocatable, save :: concentration_index_iv !matrix from size and composition to grid index
  Double precision,dimension(:), allocatable, save :: concentration_number	!number concentration of each grid cell
  double precision , dimension(:,:), allocatable, save :: concentration_mass
  double precision, dimension(:),allocatable,save :: frac_oligo

  double precision, dimension(:), allocatable, save :: gas_emis !storing Gas consentration (emission) micm^3cm^-3 
  double precision,dimension(:,:), allocatable, save :: init_bin_mass
  double precision,dimension(:), allocatable, save :: init_mass
  double precision,dimension(:,:), allocatable, save :: emis_bin_mass
  double precision,dimension(:), allocatable, save :: init_bin_number 
  double precision,dimension(:), allocatable, save :: emis_bin_number

  double precision , dimension(:,:), allocatable, save :: emission_rate
  double precision , dimension(:), allocatable, save :: emission_num_rate

  Double precision,dimension(:), allocatable, save :: wet_diameter      !Aerosol wet diameter (\B5m). of each grid cell
  Double precision,dimension(:), allocatable, save :: wet_mass          !Aerosol wet mass (\B5g). of each grid cell
  Double precision,dimension(:), allocatable, save :: wet_volume        !Aerosol wet volume (\B5m^3). of each grid cell
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
  double precision ,dimension(:), allocatable, save :: saturation_vapor_pressure! (in torr)
  double precision ,dimension(:), allocatable, save :: enthalpy_vaporization! (in kJ/mol)
  double precision ,dimension(:), allocatable, save :: t_ref ! a reference temperature at which Henry's law constant is given (in K)
  double precision ,dimension(:), allocatable, save :: henry ! Henry's law constant at t_ref (in M/atm)  
  character(len=800),dimension(:), allocatable, save :: smiles ! (\B5g/mol) gas=phase
  integer,dimension(:), allocatable, save :: aerosol_hydrophilic,aerosol_hydrophobic
  
  integer ,dimension(:), allocatable, save :: inon_volatile 

  character(len=4),dimension(:), allocatable, save :: partitioning ! HPHO HPHI BOTH

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, dimension(:), allocatable, save :: Ncoefficient, index_first, index_second
  double precision, dimension(:), allocatable, save :: coefficient
  integer, save :: coef_size
  double precision, save :: surface_tension_inorg, surface_tension_aq, surface_tension_org

  double precision, save :: dorg
  integer, save :: coupled_phases,i_hydrophilic
  integer, save :: nlayer
  integer ,dimension(:), allocatable, save :: layer_number ! Number of the layer
  double precision, dimension(:), allocatable, save :: Vlayer
  integer, save :: activity_model
  double precision, dimension(:), allocatable, save :: lwc_nsize, &
       ionic_nsize, proton_nsize, chp_nsize
  double precision, dimension(:,:), allocatable, save :: liquid_nsize
  double precision, dimension(:,:), allocatable, save :: surface_equilibrium_conc_nsize 

  !! part 8: divers parameters (species, I/O)
  character (len=200), save :: Coefficient_file ! repartition coefficient file
  character (len=200), save :: init_aero_conc_mass_file ! File for aeroslos initial mass concentrations
  character (len=200), save :: init_aero_conc_num_file ! File for aerosols initial number concentrations
  character (len=200), save :: init_gas_conc_file ! File for gas-phase initial conc.
  character (len=200), save :: species_list_file  ! File for species list.
  character (len=200), save :: reaction_list_file ! File for species list.
  character (len=40) , save :: chemID,resID,initID      ! IDs for chemistry/output
  character (len=80), save  :: chemID2 ! IDs to save chemID/chemID
  character (len=200), save :: aerosol_species_list_file ! File for species list.
  character (len=200), save :: aerosol_structure_file ! File for species list.
  character (len=200), save :: namelist_species ! Namelist file for species list.
  character (len=200), save :: particles_composition_file ! File for particles composition
  character (len=200), save :: emis_gas_file
  character (len=200), save :: emis_aero_mass_file
  character (len=200), save :: emis_aero_num_file
  character (len=80), dimension(:), allocatable, save :: isorropia_species_name
  character (len=80), dimension(:), allocatable, save :: aec_species_name
  character (len=40), save :: precursor
  character (len=40), dimension(:), allocatable, save :: species_name
  character (len=40), dimension(:), allocatable, save :: aerosol_species_name
  integer, save :: spec_name_len
  character (len=40), dimension(:), allocatable, save :: emis_gas_species_name
  character (len=40), dimension(:), allocatable, save :: emis_aer_species_name
  character (len=200), save :: output_directory, output_conc_file
  
  ! Photolysis
  character (len=200), save :: photolysis_file ! File for photolysis list.
  character (len=200), save :: photolysis_dir ! Directory for photolysis list.
  character (len=40), dimension(:), allocatable, save :: photolysis_name
  integer, save :: n_time_angle
  double precision, save :: time_angle_min, delta_time_angle
  integer, save :: n_latitude
  double precision, save :: latitude_min, delta_latitude
  integer, save :: n_altitude
  double precision, save :: altitude_photolysis_input(30)

  ! meteo
  character (len=200) , save:: meteo_file
  logical, save :: imeteo 
  
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
  
  ! FINAL PART: OUTPUTS
  Double precision, dimension(:), allocatable, save     :: output_time    !(nt+1)                  
  Double precision, dimension(:,:), allocatable, save   :: output_gas     !(nt+1,N_gas)             
  Double precision, dimension(:,:), allocatable, save   :: output_numb    !(nt+1,N_size+1)         
  Double precision, dimension(:,:), allocatable, save   :: output_diam    !(nt+1,N_size)           
  Double precision, dimension(:,:,:), allocatable, save :: output_aero    !(nt+1,N_aerosol,N_size) 
  Double precision, dimension(:,:,:), allocatable, save :: output_TM      !(nt+1,N_aerosol,3)        
  Double precision, dimension(:,:), allocatable, save   :: output_special !(nt+1,8) 
  Double precision, dimension(:,:), allocatable, save   :: output_pH      !(nt+1,N_size)   
  integer,save                                          :: t_out
  integer,save                                          :: ipm1,ipm25,ipm10,mixing_nb
  
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

  ! Flag for the verbose mode:
  ! When activated, more information is displayed.
  logical :: ssh_verbose = .false.

contains


  subroutine ssh_read_namelist(namelist_file)


    ! =============================================================
    !
    !    read simulation setting from file
    !
    !    input : namelist.ssh file 
    ! ============================================================= 

    implicit none
    integer :: i,ierr, nml_out, s , count, icmt
    character (len=200), intent(in) :: namelist_file
    character (len=200) :: namelist_out

    ! genoa read input lists
    character (len=40)  :: tmpID0, tmpID1, tmp_name
    character (len=200) :: tmp_string
    character (len=400) :: output_aero_list, output_gas_list
    character (len=400) :: err_species_list 

    ! namelists to read namelist.ssh file 

    namelist /setup_meteo/ latitude, longitude, Temperature, Pressure,&
         Humidity, Relative_Humidity, meteo_file, cloud_water

    namelist /setup_time/ initial_time, final_time, delta_t,time_emis

    namelist /initial_condition/ with_init_num, tag_init, tag_dbd, N_sizebin,&
         wet_diam_estimation, init_gas_conc_file,&
         init_aero_conc_mass_file, init_aero_conc_num_file, &
         cst_gas_file, cst_aero_file, init_species_file ! genoa inputs
         
    namelist /initial_diam_distribution/ diam_input

    namelist /emissions/ tag_emis, with_emis_num, emis_gas_file, &
         emis_aero_mass_file, emis_aero_num_file

    namelist /mixing_state/ tag_external, N_groups, N_frac, kind_composition

    namelist /fraction_distribution/ frac_input

    namelist /gas_phase_species/ species_list_file, reaction_list_file

    namelist /aerosol_species/ aerosol_species_list_file, aerosol_structure_file

    namelist /physic_gas_chemistry/ tag_chem, attenuation, option_photolysis, &
         time_update_photolysis, & 
         with_heterogeneous, with_adaptive, &
         adaptive_time_step_tolerance, min_adaptive_time_step, &
         RO2_list_file, tag_RO2, &
         photolysis_dir, photolysis_file, &
         n_time_angle, nsza, time_angle_min, delta_time_angle, &
         n_latitude, latitude_min, delta_latitude, &
         n_altitude, altitude_photolysis_input, & 
         tag_twostep, keep_gp, &
         kwall_gas, kwall_particle, Cwall, eddy_turbulence, surface_volume_ratio, kwp0,radius_chamber, option_cloud

    namelist /physic_particle_numerical_issues/ DTAEROMIN, redistribution_method,&
         with_fixed_density, fixed_density, splitting

    namelist /physic_coagulation/ with_coag, i_compute_repart, i_write_repart, Coefficient_file, Nmc

    namelist /physic_condensation/ with_cond, tag_icut, Cut_dim, ISOAPDYN, IMETHOD, &
         soap_inorg, nlayer,&
         with_kelvin_effect, tequilibrium,&
         dorg, coupled_phases, activity_model, epser, epser_soap, niter_eqconc, niter_water, co2_conc_ppm, NACL_IN_THERMODYNAMICS, &
         SOAPlog, reaction_soap_file

    namelist /physic_nucleation/ with_nucl, nucl_model_binary, nucl_model_ternary, &
         scal_ternary, nucl_model_hetero, scal_hetero, nesp_org_h2so4_nucl, &
         name_org_h2so4_nucl_species, nucl_model_org, scal_org, nexp_org,&
         nesp_org_nucl, name_org_nucl_species 


    namelist /output/ output_directory, output_type, particles_composition_file, &
                      output_aero_list, output_gas_list, & ! genoa inputs
                      err_species_list, nout_soa, & ! genoa inputs
                      ref_conc_files_in !genoa inputs


    ! Update print settings based on tag_genoa
    if (tag_genoa.eq.1) then ! Fast mode
        ! No further print info
        ssh_standalone = .true.
        ssh_logger = .false.
    else
        ! And write to namelist.out
        nml_out = 101
        namelist_out = "namelist.out"
        open(nml_out, file = namelist_out)
    endif


    if (ssh_standalone) write(*,*) "=========================start read namelist.ssh file======================"
    if (ssh_logger) write(logfile,*) "=========================start read namelist.ssh file======================"
    
    ! read namelist.ssh file !
    open(unit = 10, file = namelist_file, status = "old")

    ! Use default values if they are not given in namelist.ssh

    ! meteorological setup

    Temperature = -999.d0
    Pressure = -999.d0
    Relative_humidity = -999.d0
    cloud_water = -999.d0
    humidity = -999.d0
    
    meteo_file = ""
    read(10, nml = setup_meteo, iostat = ierr)

    if (meteo_file == "") then
       if (ssh_standalone) write(*,*) "File for meteorological data is not given."
       if (ssh_logger) write(logfile,*) "File for meteorological data is not given."
       imeteo = .false.
    else
       if (ssh_standalone) write(*,*) "File for meteorological data is read from file",trim(meteo_file)
       if (ssh_logger) write(logfile,*) "File for meteorological data is read from file",trim(meteo_file)
       imeteo = .true.
    endif

    if (ierr .ne. 0) then
       write(*,*) "setup_meteo data can not be read."
       stop
    else ! default output meteo data to check

       ! if  (Relative_Humidity .gt. 0d0 ) then
       !    if (      Relative_Humidity.lt.Threshold_RH_inf &
       !         .or. Relative_Humidity.gt.Threshold_RH_sup) then
       !       if (ssh_standalone) write(*,*) 'Warning : clipping relative humidity.'
       !       if (ssh_logger) write(logfile,*) 'Warning : clipping relative humidity.'
       !       Relative_Humidity = DMIN1(DMAX1(Relative_Humidity, Threshold_RH_inf), Threshold_RH_sup)
       !    endif
       !    call compute_psat_sh(Relative_Humidity, temperature, Pressure, pressure_sat, humidity)
       ! else
       !    call compute_psat_rh(humidity, temperature, Pressure, pressure_sat, Relative_Humidity)
       !    if (      Relative_Humidity.lt.Threshold_RH_inf &
       !         .or. Relative_Humidity.gt.Threshold_RH_sup) then
       !       if (ssh_standalone) write(*,*) 'Warning : clipping relative humidity.'
       !       if (ssh_logger) write(logfile,*) 'Warning : clipping relative humidity.'
       !       Relative_Humidity = DMIN1(DMAX1(Relative_Humidity, Threshold_RH_inf), Threshold_RH_sup)
       !       call compute_psat_sh(Relative_Humidity, temperature, Pressure, pressure_sat, humidity)
       !    endif
       ! end if

       ! cloud_water = 0.d0 by default
       ! if it is not given in namelist
       if (cloud_water == -999.d0) then
          cloud_water = 0.d0 ! in kg/kg
       endif
      
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
       ! if (ssh_standalone) write(*,*) 'Specific Humidity', Humidity
       ! if (ssh_logger) write(logfile,*) 'Specific Humidity', Humidity
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

    ! Allocate meteo data
    if ( .not. allocated(temperature_array)) allocate(temperature_array(nt))
    if ( .not. allocated(pressure_array)) allocate(pressure_array(nt))
    if ( .not. allocated(humidity_array)) allocate(humidity_array(nt))
    if ( .not. allocated(relative_humidity_array)) allocate(relative_humidity_array(nt))

    ! initial_condition
    wet_diam_estimation = -999

    ! genoa init - section initial_condition
    cst_gas_file= "---"
    cst_aero_file = "---"
    init_species_file = "---"
    
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
          tag_init = 0     ! default tag_init == 0
          if (ssh_standalone) write(*,*) 'Internally mixed aerosol species are provided for initial condition.'
          if (ssh_logger) write(logfile,*) 'Internally mixed aerosol species are provided for initial condition.'
       end if

       if (ssh_standalone) write(*,*) 'Gas-phase conc. input file :', trim(init_gas_conc_file)
       if (ssh_logger) write(logfile,*) 'Gas-phase conc. input file :', trim(init_gas_conc_file)
       if (ssh_standalone) write(*,*) 'Particle conc. input file :', trim(init_aero_conc_mass_file)
       if (ssh_logger) write(logfile,*) 'Particle conc. input file :', trim(init_aero_conc_mass_file)
       if (with_init_num .eq. 1) then
          if (ssh_standalone) write(*,*) 'Aerosol number conc. is read from file :', trim(init_aero_conc_num_file)
          if (ssh_logger) write(logfile,*) 'Aerosol number conc. is read from file :', trim(init_aero_conc_num_file)
       else 
          with_init_num = 0  ! default with_init_num == 0
          if (ssh_standalone) write(*,*) ' Aerosol number conc. is estimated from mass and diameter.' 
          if (ssh_logger) write(logfile,*) ' Aerosol number conc. is estimated from mass and diameter.' 
       end if
    end if

    if (ssh_standalone) write(*,*) 'N_sizebin', N_sizebin
    if (ssh_logger) write(logfile,*) 'N_sizebin', N_sizebin    

    if ( .not. allocated(init_bin_number)) allocate(init_bin_number(N_sizebin))
    init_bin_number = 0.d0

    if (tag_dbd == 1)  then  ! initialisation for sizebin
       if ( .not. allocated(diam_input)) allocate(diam_input(N_sizebin+1))
       diam_input = 0.d0
    else ! default read the boundary of sizebin (tag_dbd = 0)
       tag_dbd = 0
       if ( .not. allocated(diam_input)) allocate(diam_input(2))
       diam_input = 0.d0
    end if

    ! GENOA read tag_init_set from initID
    if (trim(initID).ne."-" .and. trim(initID).ne."0") then
        read(initID,*,iostat=ierr) tag_init_set ! read integer
        if (ierr.ne.0) then
            print*, "Error: can not read tag_init_set from initID: ",trim(initID)
            stop
        else
            if (ssh_standalone) write(*,*)  "Read tag_init_set is ",tag_init_set
            if (ssh_logger) write(logfile,*)"Read tag_init_set is ",tag_init_set
        endif
    else
        tag_init_set = 0            
    endif
    
    ! genoa read initial conditions    
    if (init_species_file.ne."---" .and. tag_init_set.ne.0) then
        if (ssh_standalone) write(*,*) 'Read inital conditions for some speices: ',trim(init_species_file)
        if (ssh_logger) write(logfile,*) 'Read inital conditions for some speices: ',trim(init_species_file)
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

       if ( .not. allocated(emis_bin_number)) allocate(emis_bin_number(N_sizebin))
       emis_bin_number = 0.d0
       if (tag_emis .ne. 0) then
          if (ssh_standalone) write(*,*) 'Gas-phase conc. emission file :', trim(emis_gas_file)
          if (ssh_logger) write(logfile,*) 'Gas-phase conc. emission file :', trim(emis_gas_file)
          if (ssh_standalone) write(*,*) 'Particle conc. emission file :', trim(emis_aero_mass_file)
          if (ssh_logger) write(logfile,*) 'Particle conc. emission file :', trim(emis_aero_mass_file)
          if (with_emis_num == 1) then
             if (ssh_standalone) write(*,*) 'Emitted aerosol number conc. is read from file :', &
                 trim(emis_aero_num_file)
             if (ssh_logger) write(logfile,*) 'Emitted aerosol number conc. is read from file :', &
                 trim(emis_aero_num_file)
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
          if ( .not. allocated(frac_input)) allocate(frac_input(N_frac+1)) 
       else
          kind_composition = 0 ! default
          if ( .not. allocated(frac_input)) allocate(frac_input(2)) 
       end if
    end if
    
    ! fraction_distribution
    if ( .not. allocated(frac_bound)) allocate(frac_bound(N_frac+1))
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
    reaction_list_file = "---"
    read(10, nml = gas_phase_species, iostat = ierr)
    if (ierr .ne. 0) then
       write(*,*) "gas_phase_species data can not be read."
       stop
    else
       if (ssh_standalone) write(*,*) ''
       if (ssh_logger) write(logfile,*) ''
       if (ssh_standalone) write(*,*) '<<<< Species lists >>>>'
       if (ssh_logger) write(logfile,*) '<<<< Species lists >>>>'

       ! update species_list_file with chemID
       if (trim(chemID) .ne. "-") species_list_file = trim(adjustl(species_list_file))//"/"//trim(chemID2)//".species"

       if (ssh_standalone) write(*,*) 'gas phase species file :', trim(species_list_file)
       if (ssh_logger) write(logfile,*) 'gas phase species file :', trim(species_list_file)
       
       ! read gas-phase species namelist ! unit = 11
        open(unit = 11, file = species_list_file, status = "old")

        if (ssh_standalone) write(*,*) 'read gas-phase species list.'
        if (ssh_logger) write(logfile,*) 'read gas-phase species list.'
        
        ierr = 0  ! reset ierr
        count = 0 ! count gas-phase species
        icmt = 0  ! count comment lines
        do while(ierr .eq. 0)
           read(11, *, iostat=ierr) tmp_name
           if (ierr /= 0) exit ! no line to read
           tmp_name = adjustl(tmp_name)
           if (trim(tmp_name) == "" .or. tmp_name(1:1) == "%" .or. tmp_name(1:1) == "#") then
             icmt = icmt + 1 ! read comment lines
           else
             count = count + 1 ! read gas-phase species
           endif
        end do

        N_gas = count ! number of gas-phase species
        
        if (ssh_standalone) write(*,*) 'Number of gas-phase species', N_gas
        if (ssh_logger) write(logfile,*) 'Number of gas-phase species', N_gas

        if ( .not. allocated(molecular_weight)) allocate(molecular_weight(N_gas))
        if ( .not. allocated(species_name)) allocate(species_name(N_gas))

        rewind 11
        
        if (icmt .gt. 0) then ! read comment lines if needed
            if (ssh_standalone) write(*,*) 'read comment lines from gas-phase species list',icmt
            if (ssh_logger) write(logfile,*) 'read comment lines from gas-phase species list',icmt
          do s = 1, icmt
            read(11, *)
          enddo
        endif

        do s = 1, N_gas
           read(11, *) species_name(s), molecular_weight(s)
           if (molecular_weight(s).le. 0.d0) then
             print*,'Error: input MWs <= 0',s, species_name(s),&
                     molecular_weight(s)
             stop
           endif
        enddo
        close(11)
    
       ! read & update reaction file
       if (reaction_list_file.ne."---") then

          ! update reaction_list_file with chemID
          if (trim(chemID) .ne. "-") reaction_list_file=trim(adjustl(reaction_list_file))//"/"//trim(chemID2)//".reactions"

          if (ssh_standalone) write(*,*) 'reaction list file:', trim(reaction_list_file)
          if (ssh_logger) write(logfile,*) 'reaction list file', trim(reaction_list_file)
          
          ! get n_reaction, n_photolysis from the first line
          ! call ssh_read_reaction_file()
       else
         reaction_list_file = "./src/include/CHEMISTRY/cb05/CB05_poa_nospack.reactions" 
         print*, "Not reaction file is provided"
         print*, "Default file is used: ", trim(reaction_list_file)
       endif   
    end if

    aerosol_structure_file="---"
    read(10, nml = aerosol_species, iostat = ierr)
    if (ierr .ne. 0) then
       write(*,*) "aerosol_species data can not be read."
       stop
    else
       ! update aerosol_species_list_file with chemID
       if (trim(chemID).ne."-") then
         aerosol_species_list_file = trim(adjustl(aerosol_species_list_file))//"/"//trim(chemID2)//".aer.vec"
       endif
       if (ssh_standalone) write(*,*) 'particle species file :', trim(aerosol_species_list_file)
       if (ssh_logger) write(logfile,*) 'particle species file :', trim(aerosol_species_list_file)
       if (aerosol_structure_file.ne."---") then
          if (ssh_standalone) write(*,*) 'Read aerosol structure from file :', & 
                trim(aerosol_structure_file)
          if (ssh_logger) write(logfile,*) 'Read aerosol structure from file :', & 
                trim(aerosol_structure_file)
       endif
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
    nsza = -999
    time_angle_min = -999.d0
    delta_time_angle = -999.d0
    n_latitude = -999
    latitude_min = -999.d0
    delta_latitude = -999.d0
    n_altitude = -999
    altitude_photolysis_input = -999.d0
    kwall_gas=0.d0
    kwall_particle=0.d0
    Cwall=0.d0
    eddy_turbulence=0.d0
    surface_volume_ratio=0.d0
    kwp0=0.d0
    radius_chamber=0.d0
    option_cloud = -999

    ! default genoa related paramters
    tag_twostep = 1 ! default option
    keep_gp = 0
    
    ! init genoa files: RO2 reaction
    tag_RO2 = 0
    RO2_list_file="---"

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
       ! nsza = 11 by default
       ! if it is not given in namelist
       if (nsza == -999) then
          nsza = 11
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
          ! set default values
          n_reaction = 0
          n_photolysis = 0
       else
          if (ssh_standalone) write(*,*) ''
          if (ssh_logger) write(logfile,*) ''
          if (ssh_standalone) write(*,*) '<<<< Gas-phase chemistry >>>>'
          if (ssh_logger) write(logfile,*) '<<<< Gas-phase chemistry >>>>'

          ! read reaction files
          call ssh_read_reaction_file()
          
          if (with_heterogeneous == 1) then
             if (ssh_standalone) write(*,*) 'with heterogeneous reaction.'
             if (ssh_logger) write(logfile,*) 'with heterogeneous reaction.'
          else  ! default with_heterogeneous == 0
             with_heterogeneous = 0
             if (ssh_standalone) write(*,*) 'without heterogeneous reaction.' 
             if (ssh_logger) write(logfile,*) 'without heterogeneous reaction.' 
          end if

          if (tag_genoa.gt.0.or.keep_gp.eq.1) then
            keep_gp = 1 ! Force keep_gp = 1 if under genoa modes
            if (ssh_standalone) write(*,*) 'keep_gp is activated' 
            if (ssh_logger) write(logfile,*) 'keep_gp is activated' 
            if (tag_twostep.eq.0) then
                if (ssh_standalone) write(*,*) 'keep_gp is only available with two_step. Check namelist' 
                if (ssh_logger) write(logfile,*) 'keep_gp is only available with two_step. Check namelist'
                stop
            endif
          else
            if (ssh_standalone) write(*,*) 'use ROS2 solver.'
            if (ssh_logger) write(logfile,*) 'use ROS2 solver.'
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
       
       !Force to use twostep if wall losses are used
       if (Cwall.eq.0.d0) then
          kwall_gas=0.d0
       else
          tag_twostep=1
       endif
       if (kwall_particle>0.d0) then
          tag_twostep=1
       endif

       ! option_cloud = 0 by default
       ! if it is not given in namelist
       if (option_cloud == -999) then
          option_cloud = 0
       endif
       
    endif

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
          with_fixed_density = 1
       endif
       ! fixed_density = 1.84d-06 by default
       ! if it is not given in namelist
       if (fixed_density == -999.d0) then
          fixed_density = 1.4d-06
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
             if (ssh_standalone) write(*,*) 'coefficient file : ', trim(Coefficient_file)
             if (ssh_logger) write(logfile,*) 'coefficient file : ', trim(Coefficient_file)
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
             if (ssh_standalone) write(*,*) 'coefficient file : ', trim(Coefficient_file)
             if (ssh_logger) write(logfile,*) 'coefficient file : ', trim(Coefficient_file)
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
    ICUT=0 !default in case no input of ICUT in namelist
    tag_icut=0 !default in case no input of tag_icut in namelist
    set_icut = 1 !default fix ICUT in the simulation
    imethod=1 !ROS2 explicit method in SOAP
    SOAPlog=2
    reaction_soap_file="--"
    co2_conc_ppm=410.d0 !Default CO2 concentrations set to 460 ppm
    soap_inorg=0
    niter_eqconc=1
    niter_water=1
    NACL_IN_THERMODYNAMICS=0
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

       
       if (with_cond == 1)  then ! default
          if (ssh_standalone) write(*,*) '! ! ! with condensation/ evaporation.'
          if (ssh_logger) write(logfile,*) '! ! ! with condensation/ evaporation.'

          ! initialise parameters regarding ICUT
          if (tag_icut .ne. 0) then
                set_icut = 0
                if (ssh_standalone) write(*,*) 'inorganics are computed with the hybrid method.'
                if (ssh_logger) write(logfile,*) 'inorganics are computed with the hybrid method.'
                if (ssh_standalone) write(*,*) 'Criteria for the computation of cut-off diameter: '
                if (ssh_logger) write(logfile,*) 'Criteria for the computation of cut-off diameter: '
                if (tag_icut .eq. 1) then ! c/e timescale criteria
                    if (ssh_standalone) write(*,*) '-- c/e characteristic timescales with the cut-off parameter given at ',Cut_dim
                    if (ssh_logger) write(logfile,*) '-- c/e characteristic timescales with the cut-off parameter given at ',Cut_dim
                elseif (tag_icut .eq. 2) then ! solver efficientcy criteria
                    if (ssh_standalone) write(*,*) '-- ETR solver efficiency with the cut-off parameter given at ',Cut_dim
                    if (ssh_logger) write(logfile,*) '-- ETR solver efficiency with the cut-off parameter given at ',Cut_dim
                elseif (tag_icut .eq. 3) then ! modified QSSA criteria
                    if (ssh_standalone) write(*,*) '-- modified QSSA criteria with the cut-off parameter given at ',Cut_dim
                    if (ssh_logger) write(logfile,*) '-- modified QSSA criteria with the cut-off parameter given at ',Cut_dim
                else
                    print*,'tag_icut is not recognize.',tag_icut
                    stop
                endif
          else ! fix ICUT, tag_icut = 0
                set_icut = 1
                if (Cut_dim .lt. diam_input(1)) then !dynamic
                    if (ssh_standalone) write(*,*) 'inorganics are at dynamic.'
                    if (ssh_logger) write(logfile,*) 'inorganics are at dynamic.'
                elseif (Cut_dim .ge. diam_input(size(diam_input))) then !equilibrium
                    if (ssh_standalone) write(*,*) 'inorganic compounds are at equilibrium.'
                    if (ssh_logger) write(logfile,*) 'inorganic compounds are at equilibrium.'
                else ! hybrid method with fixed cut-off diameter
                        if (ssh_standalone) write(*,*) 'inorganics are computed with the hybrid method.'
                        if (ssh_logger) write(logfile,*) 'inorganics are computed with the hybrid method.'
                        if (ssh_standalone) write(*,*) 'The cut-off diameter is fixed at ',Cut_dim
                        if (ssh_logger) write(logfile,*) 'The cut-off diameter is fixed at ',Cut_dim
                endif
          endif

          if (ISOAPDYN == 1) then
             ICUT_org = 0
             if (ssh_standalone) write(*,*) 'Dynamic SOA computation are at dynamic.'
             if (ssh_logger) write(logfile,*) 'Dynamic SOA computation are at dynamic.'
          else
             ISOAPDYN = 0   ! defalut
             ICUT_org = N_sizebin
             if (ssh_standalone) write(*,*) 'Dynamic SOA computation are at equilibrium.'
             if (ssh_logger) write(logfile,*) 'Dynamic SOA computation are at equilibrium.'
          end if
       else
          with_cond = 0
          if (ssh_standalone) write(*,*) 'without condensation/ evaporation.'
          if (ssh_logger) write(logfile,*) 'without condensation/ evaporation.'
       end if
       if(coupled_phases == 1) then
          i_hydrophilic = 1 
          if(ISOAPDYN == 0) then
          !  if (ssh_standalone) write(*,*) 'Phases can not be coupled if ISOAPDYN = 0'
          !  if (ssh_logger) write(logfile,*) 'Phases can not be coupled if ISOADYN = 0'
            i_hydrophilic = 0
          endif
       else 
          i_hydrophilic = 0 
       endif
    end if

    ! nucleation
    nucl_model_binary = 0
    nucl_model_ternary = 0
    scal_ternary = 1.0
    nucl_model_hetero = 0
    scal_hetero = 0.1            
    nesp_org_h2so4_nucl = 0
    name_org_h2so4_nucl_species(1) = 'none'
    name_org_h2so4_nucl_species(2) = 'none'
    nucl_model_org = 0
    scal_org = 0.1            
    nexp_org = 0.1            
    nesp_org_nucl = 0
    !name_org_nucl_species(1) = 'none'
    !name_org_nucl_species(2) = 'none'

    read(10, nml = physic_nucleation, iostat = ierr)
    if (ierr .ne. 0) then
       write(*,*) "physic_nucleation data can not be read."
       stop
    else
       if (with_nucl == 1) then
          if (nucl_model_binary == -999) nucl_model_binary = 0
          if (nucl_model_ternary == -999) nucl_model_ternary = 0
          if (nucl_model_hetero == -999) nucl_model_hetero = 0
          if (nucl_model_org == -999) nucl_model_org = 0
             
          if (ssh_standalone) write(*,*) '! ! ! with nucleation. -- Need to have lower diameter about 1nm.'
          if (ssh_logger) write(logfile,*) '! ! ! with nucleation. -- Need to have lower diameter about 1nm.'
          if (nucl_model_binary == 0) then
             if (ssh_standalone) write(*,*) ' No binary nucleation'
             if (ssh_logger) write(logfile,*) 'No binary nucleation'
          else
             if (nucl_model_binary == 1) then
                if (ssh_standalone) write(*,*) 'nucleation model : binary Veahkamaki'
                if (ssh_logger) write(logfile,*) 'nucleation model : binary Veahkamaki'
             else 
                if (nucl_model_binary == 2) then
                   if (ssh_standalone) write(*,*) 'nucleation model : binary Kuang'
                   if (ssh_logger) write(logfile,*) 'nucleation model : binary Kuang'
                endif
             endif
          endif
          if (nucl_model_ternary == 0) then
             if (ssh_standalone) write(*,*) ' No ternary nucleation'
             if (ssh_logger) write(logfile,*) 'No ternary nucleation'
          else
             if(scal_ternary == -999.d0) scal_ternary = 1
             if (nucl_model_ternary == 1) then
                if (ssh_standalone) write(*,*) 'nucleation model : ternary Napari - scaling:',scal_ternary
                if (ssh_logger) write(logfile,*) 'nucleation model : ternary Napari - scaling:',scal_ternary
             else 
                if (nucl_model_ternary == 2) then
                  if (ssh_standalone) write(*,*) 'nucleation model : ternary Merikanto - scaling:',scal_ternary
                  if (ssh_logger) write(logfile,*) 'nucleation model : ternary Merikanto - scaling:',scal_ternary
               endif
             endif
          endif
          if (nucl_model_hetero == 0) then
             if (ssh_standalone) write(*,*) ' No heteromolecular nucleation'
             if (ssh_logger) write(logfile,*) 'No heteromolecular nucleation'
          else
             if(scal_hetero == -999.d0) scal_hetero = 1
             if (nucl_model_hetero == 1) then
                if(nesp_org_h2so4_nucl.GT.2) then
                  if (ssh_standalone) write(*,*) 'nucleation model : heteromolecular - no more than 2 species',nesp_org_h2so4_nucl
                  if (ssh_logger) write(logfile,*) 'nucleation model : heteromolecular- no more than 2 species',nesp_org_h2so4_nucl
                endif   
                if (ssh_standalone) write(*,*) 'nucleation model : heteromolecular',scal_hetero
                if (ssh_logger) write(logfile,*) 'nucleation model : heteromolecular',scal_hetero
             endif
          endif
          if (nucl_model_org == 0) then
             name_org_nucl_species(1) = 'none'
             if (ssh_standalone) write(*,*) ' No organic nucleation'
             if (ssh_logger) write(logfile,*) 'No organic nucleation'
          else
             if(scal_org == -999.d0) scal_org = 1
             if(nexp_org == -999.d0) nexp_org = 3
             if (nucl_model_org == 1) then
                write(*,*) nesp_org_nucl, 'nesp org nucl'
                !if(nesp_org_nucl.GT.2) then
                !  if (ssh_standalone) write(*,*) 'nucleation model : org - no more than 2 species',nesp_org_nucl
                !  if (ssh_logger) write(logfile,*) 'nucleation model : org - no more than 2 species',nesp_org_nucl
                !endif   
                if (ssh_standalone) write(*,*) 'nucleation model : organic',scal_org,nexp_org
                if (ssh_logger) write(logfile,*) 'nucleation model : organic',scal_org,nexp_org
             endif
          endif
      endif
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
    
    ! genoa ! default values
    nout_soa = 0 ! default no SOA groups
    ref_conc_files_in = "---"
    output_aero_list  = "---"
    output_gas_list   = "---"
    err_species_list  = "---"
    
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
       if (output_type == 3) then
          if (ssh_standalone) write(*,*) 'results are saved in netcdf files.'
          if (ssh_logger) write(logfile,*) 'results are saved in netcdf files.'
       elseif (output_type == 2) then
          if (ssh_standalone) write(*,*) 'results are saved in binary files.'
          if (ssh_logger) write(logfile,*) 'results are saved in binary files.'
       elseif (output_type == 0) then
          if (ssh_standalone) write(*,*) 'No results are output except for the species specified in the namelist.'
          if (ssh_logger) write(logfile,*) 'No results are output except for the species specified in the namelist.'
       else
          output_type = 1   ! default
          if (ssh_standalone) write(*,*) 'results are saved in text files.'
          if (ssh_logger) write(logfile,*) 'results are saved in text files.'
       end if

       output_directory = adjustl(output_directory)
       if (ssh_standalone) write(*,*) 'output directory :', trim(output_directory)
       if (ssh_logger) write(logfile,*) 'output directory :', trim(output_directory)

       output_conc_file = trim(output_directory)//"/concs.txt"
       
       ! update output_conc_file and particles_composition_file with initID,chemID,resID
       ! get suffix
       if(trim(chemID).ne."-") then ! add chemID & resultID
           tmpID0 = "/"//trim(chemID)//"/"//trim(initID)//"."//trim(resID)
           tmpID1 = "/"//trim(initID)//"."//trim(resID)//".concs" ! for pre/ref files
       else if (trim(initID).ne."-") then ! incase only put initID
           tmpID0 = "/"//trim(initID) 
           tmpID1 = "/"//trim(initID)//".concs"
       else
           tmpID0 = "-"
       endif
       if (trim(tmpID0).ne."-") then
           ! detail output by IDs in tmpID0
           output_directory = trim(output_directory)//trim(tmpID0)
           particles_composition_file = trim(output_directory)//'.fac'
           output_conc_file = trim(output_directory)//'.concs'
           ! Create directory if it does not exist.
           call system("mkdir -p "// trim(output_directory))
       endif

       if (ssh_standalone) write(*,*) 'Particles composition file : ', trim(particles_composition_file)
       if (ssh_logger) write(logfile,*) 'Particles composition file : ', trim(particles_composition_file)

       if (ssh_standalone) write(*,*)   'Concs file : ', trim(output_conc_file)
       if (ssh_logger) write(logfile,*) 'Concs file : ', trim(output_conc_file)
       
    end if

    !!!!! genoa
    ! Read output_aero_list, output_gas_list, err_species_list from the user
    ! read a list of output aerosol species
    if (output_aero_list .ne. "---") then
        if (ssh_standalone) write(*,*) 'Read output aerosol speices: ',trim(output_aero_list)
        if (ssh_logger) write(logfile,*) 'Read output aerosol speices: ',trim(output_aero_list)

        ! Count the number of substrings
        call count_delimiter(output_aero_list, ",", nout_aero)
        !print*,"Find number: ", nout_aero
        
        ! Allocate memory
        allocate(output_aero_species(nout_aero))
        allocate(output_aero_index(nout_aero))
        output_aero_index = 0
        
        ! Separate and remove spaces from the input string
        do i = 1, nout_aero
            call get_token(output_aero_list, tmp_name, ",")
            output_aero_species(i) = tmp_name
            !print*,i," now: ",output_aero_list," species taken: ",output_aero_species(i)
        end do
    else
        nout_aero = 0
        allocate(output_aero_species(1))
        allocate(output_aero_index(1))
    end if
    
    ! read a list of output gas-phase species - only used if tag_genoa = 1 or out_type = 0
    if (output_gas_list .ne. "---") then
        if (ssh_standalone) write(*,*) 'Read output gas species: ',trim(output_gas_list)
        if (ssh_logger) write(logfile,*) 'Read output gas species: ',trim(output_gas_list)
        
        ! Count the number of substrings
        call count_delimiter(output_gas_list, ",", nout_gas)
        !print*,"Find number: ", nout_gas
        
        ! Allocate memory
        allocate(output_gas_species(nout_gas))
        allocate(output_gas_index(nout_gas))
        output_gas_index = 0
        
        ! Separate and remove spaces from the input string
        do i = 1, nout_gas
            call get_token(output_gas_list, tmp_name, ",")
            output_gas_species(i) = tmp_name
            !print*,i," now: ",output_gas_list," species taken: ",output_gas(i)
        end do
    else
        nout_gas = 0
        allocate(output_gas_index(1))
        allocate(output_gas_species(1))
    end if
    
    ! read the list of species to compute reduction errors
    if (err_species_list .ne. "---") then
        if (ssh_standalone) write(*,*) 'Read output speices for error computation: ',trim(err_species_list)
        if (ssh_logger) write(logfile,*) 'Read output speices for error computation: ',trim(err_species_list)

        ! Count the number of substrings
        call count_delimiter(err_species_list, ",", nout_err)
        !print*,"Find number: ", nout_err
        
        ! Allocate memory
        allocate(output_err_sps(nout_err))
        allocate(output_err_index(2,nout_err)) ! 1: gas/ 2: aerosol
        output_err_index = 0
        
        ! Separate and remove spaces from the input string
        do i = 1, nout_err
            call get_token(err_species_list, tmp_name, ",")
            output_err_sps(i) = tmp_name
            !print*,i," now: ",err_species_list," species taken: ",output_err(i)
        end do
    else
        nout_err = 0
        allocate(output_err_sps(1))
        allocate(output_err_index(1,1))
    end if

    ! total number of output - 1 for total SOAs
    nout_total = 1 + nout_err + nout_soa

    ! Update ref soa conc file lists
    if (ref_conc_files_in .ne. "---") then
    
        ! Count number of ref files
        call count_delimiter(ref_conc_files_in, ",", nref_file)
        allocate(ref_conc_files(nref_file))
    
        count = 0 ! find output file    
        ! Separate and remove spaces from the input string
        do i = 1, nref_file
            call get_token(ref_conc_files_in, tmp_string, ",")
            
            ! detail paths to references by IDs in tmpID1
            if (trim(tmpID0).ne."-") then
              ref_conc_files(i) = trim(adjustl(tmp_string))//trim(tmpID1)
            else
              ref_conc_files(i) = tmp_string
            endif
            !print*,i," now: ",ref_conc_files_in," species taken: ",ref_conc_files(i)

            ! check if repeat file to output
            if (trim(output_conc_file).eq.trim(ref_conc_files(i))) then
                ref_conc_files(i) = "---" ! reset name
                count = count + 1
            endif
            nref_file = nref_file - count ! update
        end do
        
    else
        nref_file = 0
        allocate(ref_conc_files(1))
    end if
    
    close(10)

    if (tag_genoa.ne.1) then ! Not for genoa fast mode
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
        write(nml_out, output)
        close(nml_out)
    endif ! Not for genoa fast mode

    if (ssh_standalone) write(*,*) "=========================finish read namelist.ssh file======================"
    if (ssh_logger) write(logfile,*) "=========================finish read namelist.ssh file======================"


  end subroutine ssh_read_namelist



  ! =============================================================
  !
  !      read initial/emitted species file
  !
  ! 
  ! =============================================================


  subroutine ssh_read_inputs(nspecies,name_input_species,index_species_ssh)
    
    implicit none
    
    integer:: nspecies   !For chimere interface, number of aerosol species we want to select
    character (len=40),dimension(nspecies),optional :: name_input_species !For Chimere interface, name of the species we want to select
    integer,dimension(nspecies),optional :: index_species_ssh   !For Chimere interface, link between the index of species in CHIMERE and SSH-aerosol
    integer :: k,i,j,s,js, ind, count, ierr, ilayer, esp_layer, nline, N_count,icoun,s2, jesp
    double precision :: tmp
    double precision, dimension(:), allocatable :: tmp_aero
    character (len=40) :: ic_name, sname, tmp_name
    integer :: aerosol_type_tmp,Index_groups_tmp,inon_volatile_tmp,found_spec
    double precision :: molecular_weight_aer_tmp,collision_factor_aer_tmp, molecular_diameter_tmp
    double precision :: surface_tension_tmp, accomodation_coefficient_tmp, mass_density_tmp
    double precision :: saturation_vapor_pressure_tmp,enthalpy_vaporization_tmp
    double precision :: henry_tmp, t_ref_tmp
    character (len=800) :: smiles_tmp
    character (len=40) :: aerosol_species_name_tmp, char1,char2        
    character (len=4) :: partitioning_tmp

    ! currently only used for genoa
    character (len=2) :: ivoc0, ivoc1
    character (len=200) :: tmp_string
    double precision :: tmp_conc
    integer, dimension(:), allocatable :: tmp_index ! (gas/aero, sps index, ref sps index)
    double precision, dimension(:), allocatable :: tmp_read, tmp_fgls
  
    if (nspecies>0) then
       index_species_ssh(:)=-1
    endif
    
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
    N_count = count - 1  ! minus the first comment line
    if (nspecies>0) then
       N_aerosol = nspecies
    else
       N_aerosol = N_count
    endif

    if ( .not. allocated(aerosol_species_name)) allocate(aerosol_species_name(N_aerosol))
    spec_name_len = len(aerosol_species_name(1))
    if ( .not. allocated(Index_groups)) allocate(Index_groups(N_aerosol))
    if ( .not. allocated(aerosol_type)) allocate(aerosol_type(N_aerosol))
    if ( .not. allocated(index_species)) allocate(index_species(N_aerosol,nlayer+i_hydrophilic))
    ! initialize basic physical and chemical parameters
    if ( .not. allocated(molecular_weight_aer)) allocate(molecular_weight_aer(N_aerosol))
    if ( .not. allocated(collision_factor_aer)) allocate(collision_factor_aer(N_aerosol))
    if ( .not. allocated(molecular_diameter)) allocate(molecular_diameter(N_aerosol)) 
    if ( .not. allocated(surface_tension)) allocate(surface_tension(N_aerosol))
    if ( .not. allocated(accomodation_coefficient)) allocate(accomodation_coefficient(N_aerosol))
    if ( .not. allocated(mass_density)) allocate(mass_density(N_aerosol))
    if ( .not. allocated(inon_volatile)) allocate(inon_volatile(N_aerosol))
    if ( .not. allocated(Vlayer)) allocate(Vlayer(nlayer))
    ! relation between Aerosol and GAS
    if ( .not. allocated(aerosol_species_interact)) allocate(aerosol_species_interact(N_aerosol))
    if ( .not. allocated(smiles)) allocate(smiles(N_aerosol))
    if ( .not. allocated(saturation_vapor_pressure)) allocate(saturation_vapor_pressure(N_aerosol))
    if ( .not. allocated(enthalpy_vaporization)) allocate(enthalpy_vaporization(N_aerosol))
    if ( .not. allocated(aerosol_hydrophilic)) allocate(aerosol_hydrophilic(N_aerosol))
    if ( .not. allocated(aerosol_hydrophobic)) allocate(aerosol_hydrophobic(N_aerosol))
    if ( .not. allocated(t_ref)) allocate(t_ref(N_aerosol))
    if ( .not. allocated(henry)) allocate(henry(N_aerosol))
    
    aerosol_species_interact = 0
    inon_volatile = 0
    
    if ( .not. allocated(partitioning)) allocate(partitioning(N_aerosol))

    
    ! Read lines from aerosol species file.
    rewind 12
    count = 0
    read(12, *, iostat = ierr) ! Read a header line (#)

    s=0
    do icoun = 1, N_count
      read(12, *) aerosol_species_name_tmp, aerosol_type_tmp, &
        Index_groups_tmp, molecular_weight_aer_tmp, &
        precursor, &
        collision_factor_aer_tmp, molecular_diameter_tmp, &
        surface_tension_tmp, accomodation_coefficient_tmp, &
        mass_density_tmp, inon_volatile_tmp, partitioning_tmp, smiles_tmp, &
        saturation_vapor_pressure_tmp,enthalpy_vaporization_tmp, &
        henry_tmp, t_ref_tmp

       if (nspecies==0) then
          found_spec=1
       else
          found_spec=0
          do j=1,nspecies          
             if (name_input_species(j)=="H2SO4") name_input_species(j)="SO4"
             if (name_input_species(j)=="NH3") name_input_species(j)="NH4"
             if (name_input_species(j)=="HNO3") name_input_species(j)="NO3"
             if (name_input_species(j)=="WATER") name_input_species(j)="H2O"
             if (trim(aerosol_species_name_tmp)==trim(name_input_species(j)) &
                  .or.trim(aerosol_species_name_tmp)=="P"//trim(name_input_species(j))) then
                found_spec=1
                index_species_ssh(j)=s+1
             endif
          enddo
       endif

       if (found_spec>0) then
          s=s+1
          aerosol_species_name(s)=aerosol_species_name_tmp
          aerosol_type(s)=aerosol_type_tmp
          Index_groups(s)=Index_groups_tmp
          molecular_weight_aer(s)=molecular_weight_aer_tmp
          collision_factor_aer(s)=collision_factor_aer_tmp
          molecular_diameter(s)=molecular_diameter_tmp
          surface_tension(s)=surface_tension_tmp
          accomodation_coefficient(s)=accomodation_coefficient_tmp
          mass_density(s)=mass_density_tmp
          inon_volatile(s)=inon_volatile_tmp
          partitioning(s)=partitioning_tmp
          smiles(s)=smiles_tmp
          saturation_vapor_pressure(s)=saturation_vapor_pressure_tmp
          enthalpy_vaporization(s)=enthalpy_vaporization_tmp
          henry(s)=henry_tmp
          t_ref(s)=t_ref_tmp

          aerosol_hydrophilic(s)=0
          if (trim(partitioning(s))=="HPHI".or.trim(partitioning(s))=="BOTH") then
             aerosol_hydrophilic(s)=1
          endif
          if (aerosol_type(s)==3) then
             aerosol_hydrophilic(s)=1
          endif

          aerosol_hydrophobic(s)=0
          if (trim(partitioning(s))=="HPHO".or.trim(partitioning(s))=="BOTH") then
             aerosol_hydrophobic(s)=1
          endif

          if((inon_volatile(s).NE.1).AND.(inon_volatile(s).NE.0)) then
             write(*,*) "non_volatile should be 0 or 1", inon_volatile(s),s
             stop
          endif
          if((partitioning(s).NE."HPHO").AND.(partitioning(s).NE."HPHI").AND.(partitioning(s).NE."BOTH").AND. &
               (trim(partitioning(s)).NE."--")) then
             write(*,*) "partitioning should be --, HPHO, HPHI or BOTH, aerspec nb ", s," : ", partitioning(s)
             stop
          endif

          if(inon_volatile(s).eq.1.and.partitioning(s).eq."--".and.aerosol_species_name(s).ne."PSO4") then
             write(*,*) trim(aerosol_species_name_tmp)," is non volatile. partitioning should be defined: HPHO, HPHI or BOTH"
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
             if (ssh_standalone) write(*,*) "Error: wrong species name is given ",&
                  trim(aerosol_species_list_file), trim(precursor)
             if (ssh_logger) write(logfile,*) "Error: wrong species name is given ",&
                  trim(aerosol_species_list_file), trim(precursor)
             stop
          endif
       endif
    enddo
    close(12)

    if (nspecies>0) then
       do s=1,nspecies
          if (index_species_ssh(s)<0) print*,trim(name_input_species(s))," not found"
          stop
       enddo
    endif
    
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
    N_aerosol_layers = N_organics * (nlayer-1 + i_hydrophilic) + N_aerosol
    EH2O_layers = N_aerosol_layers
    if ( .not. allocated(mass_density_layers)) allocate(mass_density_layers(N_aerosol_layers))
    if ( .not. allocated(List_species)) allocate(List_species(N_aerosol_layers))
    if ( .not. allocated(layer_number)) allocate(layer_number(N_aerosol_layers))
    if ( .not. allocated(isorropia_species)) allocate(isorropia_species(nesp_isorropia))
    if ( .not. allocated(isorropia_species_name)) allocate(isorropia_species_name(nesp_isorropia))
    if ( .not. allocated(aec_species)) allocate(aec_species(nesp_aec))
    if ( .not. allocated(aec_species_name)) allocate(aec_species_name(nesp_aec))

    
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

    ECO3=-1
    do s = 1, N_aerosol
       ! For non-organic species.
       if (s <= N_nonorganics) then
          mass_density_layers(s) = mass_density(s)
          molecular_weight_aer(s) = molecular_weight_aer(s) * 1.0D06 ! g/mol to \B5g/mol  !!! change later
          List_species(s) = s
          if (aerosol_species_name(s) .eq. "PNA") ENa = s
          if (aerosol_species_name(s) .eq. "PSO4") ESO4 = s
          if (aerosol_species_name(s) .eq. "PNH4") ENH4 = s
          if (aerosol_species_name(s) .eq. "PNO3") ENO3 = s
          if (aerosol_species_name(s) .eq. "PHCL") ECl = s
          if (aerosol_species_name(s) .eq. "PBiPER") ind_jbiper = s
          if (aerosol_species_name(s) .eq. "PCO3") ECO3 = s

          do ilayer=1,(nlayer + i_hydrophilic)
             index_species(s,ilayer) = s
          enddo
          layer_number(s) = 1
       ! For organic species
       else
          molecular_weight_aer(s) = molecular_weight_aer(s) * 1.0D06 ! g/mol to \B5g/mol  !!! change later
          if(s.NE.N_aerosol) then !avoid water
             do ilayer = 0,(nlayer-1 + i_hydrophilic)
                esp_layer = (s-N_nonorganics-1) *(nlayer-1+i_hydrophilic) + s + ilayer
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
             !! Water
             List_species(N_aerosol_layers) = s
             do ilayer=1,(nlayer + i_hydrophilic)
                index_species(N_aerosol,ilayer) = N_aerosol_layers 
             enddo
             layer_number(N_aerosol_layers) = 1
             mass_density_layers(N_aerosol_layers) = mass_density(s)
          endif
       endif
    enddo

    if ( .not. allocated(oligo_index)) allocate(oligo_index(N_aerosol))
    if ( .not. allocated(frac_oligo)) allocate(frac_oligo(N_aerosol))
    oligo_index=0
    do s=1,N_aerosol
       char2=aerosol_species_name(s)       
       if (trim(char2(2:6))=="Oligo") then
          do s2=1,N_aerosol
             char1=aerosol_species_name(s2)
             if (trim(char2(2:))=="Oligo"//trim(char1(2:))) then
                oligo_index(s2)=s                
                exit
             endif
          enddo
       endif
    enddo

    if(with_nucl.EQ.1) then
          inon_volatile(ESO4) = 1 ! sulfate needs to be computed dynamically in case of nucleation
    endif
    
    ! read gas-phase initial concentrations unit 21
    ! no comment lines for initial & emitted data
    if ( .not. allocated(concentration_gas_all)) allocate(concentration_gas_all(N_gas))
    concentration_gas_all = 0.d0 ! set original value to 0
    if ( .not. allocated(concentration_gas)) allocate(concentration_gas(N_aerosol))
    concentration_gas=0.d0
    if ( .not. allocated(concentration_wall)) allocate(concentration_wall(N_gas))
    concentration_wall = 0.d0 ! set original value to 0    

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
             !if (ssh_standalone) write(*,*) 'initialize conc ',ic_name,js,concentration_gas_all(js)
             !if (ssh_logger) write(logfile,*) 'initialize conc ',ic_name,js,concentration_gas_all(js)
             ind = 1
             exit
          endif
       enddo

       if (ind .ne. 1) then ! not found
          if (ssh_standalone) write(*,*) "Error: wrong species name is given in file ",&
              trim(init_gas_conc_file),": " ,trim(ic_name)
          if (ssh_logger) write(logfile,*) "Error: wrong species name is given in file ",&
              trim(init_gas_conc_file),": " ,trim(ic_name)
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
       if ( .not. allocated(init_mass)) allocate(init_mass(N_aerosol))   ! aerosol initial mass concentrations inti_mass for each species
       init_mass = 0.d0
       if ( .not. allocated(init_bin_mass)) allocate(init_bin_mass(N_sizebin,N_aerosol))
       init_bin_mass = 0.d0
       if ( .not. allocated(tmp_aero)) allocate(tmp_aero(N_sizebin))
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
             if (ssh_standalone) write(*,*) "Error: wrong aerosol species name is given ", &
                 trim(init_aero_conc_mass_file), trim(ic_name)
             if (ssh_logger) write(logfile,*) "Error: wrong aerosol species name is given ", &
                 trim(init_aero_conc_mass_file), trim(ic_name)
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
             write(*,*) "Aerosol number conc. can not be read from file ", trim(init_aero_conc_num_file)
             stop
          end if
       else if (tag_init == 1) then
          write(*,*) 'with_init_num == 1 .and. tag_init == 1 - not yet build'
          stop
       end if
    end if


    ! ! ! ! ! ! 
    if (tag_emis == 1) then  ! with internal emission 

       if ( .not. allocated(emis_bin_mass)) allocate(emis_bin_mass(N_sizebin,N_aerosol))
       emis_bin_mass = 0.d0
       if ( .not. allocated(gas_emis)) allocate(gas_emis(N_gas))
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
             if (ssh_standalone) write(*,*) "Error: wrong species name is given in gas emission ", &
               trim(init_gas_conc_file), trim(ic_name)
             if (ssh_logger) write(logfile,*) "Error: wrong species name is given in gas emission ", &
               trim(init_gas_conc_file), trim(ic_name)
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
             if (ind == 0 .and. ssh_standalone) write(*,*) 'Not find the emission species', &
                 trim(emis_aero_mass_file), trim(ic_name)
             if (ind == 0 .and. ssh_logger) write(logfile,*) 'Not find the emission species', &
                 trim(emis_aero_mass_file), trim(ic_name)
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
             write(*,*) "can not read aerosol number conc. from ", trim(emis_aero_num_file)
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
    if ( .not. allocated(photolysis_name)) allocate(photolysis_name(n_photolysis))
    if ( .not. allocated(photolysis_reaction_index)) allocate(photolysis_reaction_index(n_photolysis))
    if ( .not. allocated(photolysis_rate)) allocate(photolysis_rate(n_photolysis))
    photolysis_rate = 0.d0
    
    if (option_photolysis .ne. 1) then !do not read photolysis file if not needed     
      open(unit = 34, file = photolysis_file, status = "old")
      count = 0 
      ierr = 0
      nline = 0
      do while(ierr .eq. 0)
         read(34, *, iostat=ierr) tmp_name
         if (ierr == 0) then
            tmp_name = adjustl(tmp_name)
            if (trim(tmp_name) == "" .or. tmp_name(1:1) == "#") then
               nline = nline + 1 ! Comment line
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
      enddo
      close(34)

    else
      do s = 1, n_photolysis
        photolysis_name(s) = "useless"
        photolysis_reaction_index(s) = s
      enddo
      
    endif
    
    if(nucl_model_hetero == 1) then
       do s=1,nesp_org_h2so4_nucl
          do i = 1,N_aerosol
              jesp= List_species(i)
              if(aerosol_species_name(jesp) == name_org_h2so4_nucl_species(s)) then
                org_h2so4_nucl_species(s) = i
                if (ssh_standalone) write(*,*) "Nucl. species found ",aerosol_species_name(jesp),i
                if (ssh_logger) write(logfile,*) "Nucl. species found ",aerosol_species_name(jesp),i
                exit
             endif
           enddo
       enddo
    endif

    if(nucl_model_org == 1) then
       do s=1,nesp_org_nucl
          do i = 1,N_aerosol
              jesp= List_species(i)
              if(aerosol_species_name(jesp) == name_org_nucl_species(s)) then
                org_nucl_species(s) = i
                if (ssh_standalone) write(*,*) "Org. Nucl. species found ",aerosol_species_name(jesp),i
                if (ssh_logger) write(logfile,*) "Org. Nucl. species found ",aerosol_species_name(jesp),i
                exit
             endif
           enddo
       enddo
    endif

    !!!!!!!!!! genoa related treatments

    if (tag_chem .eq. 1) then
       ! read photolysis files
       call ssh_read_photolysis_file()
    endif
    
    ! genoa assign output aero index
    do s = 1, nout_aero
        do js = 1, N_aerosol
            if (aerosol_species_name(js) .eq. output_aero_species(s)) then
                output_aero_index(s) = js
                if (ssh_standalone) write(*,*) 'Output aerosol: ', trim(output_aero_species(s)), s, js
                if (ssh_logger) write(logfile,*) 'Output aerosol: ', trim(output_aero_species(s)), s, js
            endif
        enddo
        ! check
        if (output_aero_index(s).eq.0) then
            print*,'not found aero in species list: ', s, trim(output_aero_species(s))
            stop
        endif
    enddo

    ! assign output gas index
    do s = 1, nout_gas
        do js = 1, N_gas
            if (species_name(js) .eq. output_gas_species(s)) then
                output_gas_index(s) = js
                if (ssh_standalone) write(*,*) 'Output gas: ', trim(output_gas_species(s)), s, js
                if (ssh_logger) write(logfile,*) 'Output gas: ', trim(output_gas_species(s)), s, js
            endif
        enddo
        ! check
        if (output_gas_index(s).eq.0) then
            print*,'not found output gas name in species list: ', s, trim(output_gas_species(s))
            stop
        endif
    enddo

    ! genoa assign output err species index
    do s = 1, nout_err
        count = 0
        ! check gas speices
        do js = 1, N_gas
            if (species_name(js) .eq. output_err_sps(s)) then
                output_err_index(1,s) = js
                count = count + 1
                if (ssh_standalone) write(*,*) 'Output species (aerosol) for error computation: ', trim(output_err_sps(s)), s, js
                if (ssh_logger) write(logfile,*) 'Output species (aerosol) for error computation: ', trim(output_err_sps(s)), s, js
            endif
        enddo
        ! check aerosol
        do js = 1, N_aerosol
            if (aerosol_species_name(js) .eq. output_err_sps(s)) then
                output_err_index(2,s) = js
                count = count + 1
                if (ssh_standalone) write(*,*) 'Output species (gas) for error computation: ', trim(output_err_sps(s)), s, js
                if (ssh_logger) write(logfile,*) 'Output species (gas) for error computation: ', trim(output_err_sps(s)), s, js
            endif
        enddo
        ! check
        if (count.ne.1 ) then
            print*,'Not found error species in species list: ', s, trim(output_err_sps(s)), count
            stop
        endif
    enddo
    
    ! genoa RO2 treatment settings
    if (tag_RO2.ne.0) then
        if (ssh_standalone) write(*,*) "Treat RO2-RO2 reaction:"
        if (ssh_logger) write(logfile,*) "Treat RO2-RO2 reaction:"

        if (tag_RO2 .eq. 1) then
            if (ssh_standalone) write(*,*) '  Use generated RO2 pool. tag_RO2 = ', tag_RO2
            if (ssh_logger) write(logfile,*) '  Use generated RO2 pool. tag_RO2 = ', tag_RO2
        else if (tag_RO2 .eq. 2) then
            if (ssh_standalone) write(*,*) '  Use background RO2 pool. tag_RO2 = ', tag_RO2
            if (ssh_logger) write(logfile,*) '  Use background RO2 pool. tag_RO2 = ', tag_RO2
        else if (tag_RO2 .eq. 3) then
            if (ssh_standalone) write(*,*) '  Use background + generated RO2 pool. tag_RO2 = ', tag_RO2
            if (ssh_logger) write(logfile,*) '  Use background + generated RO2 pool. tag_RO2 = ', tag_RO2
        else
            if (ssh_standalone) write(*,*) 'Error: unknown input tag_RO2 = ',tag_RO2
            if (ssh_logger) write(logfile,*) 'Error: unknown input tag_RO2 = ',tag_RO2
            stop
        endif
    endif

    ! genoa read initial profiles
    if (init_species_file.ne."---".and.tag_init_set.gt.0) then ! ivoc is used here

      if (ssh_standalone) write(*,*) 'Reading initial sets file ', trim(init_species_file)
      if (ssh_logger) write(logfile,*) 'Reading initial sets file ', trim(init_species_file)

      open(unit = 20, file = init_species_file, status = "old")
        ierr = 0 ! tag read line
        ind = 0  ! tag index - 0/1
        j = 0
        write(ivoc0,'(I0)') tag_init_set   ! start sign
        write(ivoc1,'(I0)') tag_init_set + 1 ! stop sign

        allocate(tmp_index(2)) ! index of sps/ref sps
        tmp_index = 0

        do while(ierr .eq. 0)
            count = 0 ! default tag
            read(20, '(A)', iostat=ierr) tmp_string ! read line
            if (ierr /= 0) exit
            tmp_string = adjustl(tmp_string)
            s = len_trim(tmp_string) ! len
            if (s == 0 .or. tmp_string(1:1) == '#') cycle ! Comment lines
            if (s > 2 ) then
                if (ind /= 1) cycle ! Read line in other set
                ! Read line
                if (ssh_standalone) write(*,*) '  Read init set info: ', trim(tmp_string)
                if (ssh_logger) write(logfile,*)'  Read init set info: ', trim(tmp_string)
                read(tmp_string,*, iostat=j) ic_name, tmp_name, sname
                if (j.eq.0) then ! Read line
                    count = 3
                    ! get factor
                    read(tmp_name,*, iostat=j) tmp_conc
                    !print*,"read 3: ",j,trim(ic_name),tmp_conc, trim(sname)
                else
                    count = 2
                    read(tmp_string,*, iostat=j) ic_name, tmp_conc
                    !print*, "read 2:",j, trim(ic_name), tmp_conc                 
                end if
                if (j.ne.0) then
                    count = 0
                    print*, " Warning! Can not read line from initial sets file: ", trim(tmp_string)
                endif 
            else if (tmp_string(1:2) .eq. ivoc1) then ! find the stop sign, exit
                exit
            else if (tmp_string(1:2) .eq. ivoc0) then ! find the start sign, reading next line
                ind = 1
                if (ssh_standalone) write(*,*) '  Find sign:', trim(tmp_string)
                if (ssh_logger) write(logfile,*)'  Find sign:', trim(tmp_string)
            end if

            if (count.ne.0) then ! read line info
                ! check species name and type
                do js = 1, N_gas
                    if (ic_name.eq.species_name(js)) then
                        tmp_index(1) = js
                        if (ssh_standalone) write(*,*) '  Find gas (initial concs may change): ', trim(ic_name), ind, js
                        if (ssh_logger) write(logfile,*) '  Find gas (initial concs may change): ', trim(ic_name), ind, js
                    endif
                    if (count.eq.3.and.sname.eq.species_name(js)) then
                        tmp_index(2) = js
                        if (ssh_standalone) write(*,*) '  Find gas (used to compute initial concs): ', trim(sname), ind, js
                        if (ssh_logger) write(logfile,*) '  Find gas (used to compute initial concs): ', trim(sname), ind, js
                    endif
                enddo

                ! check
                if (tmp_index(1).eq.0) then
                    print*, "  index of input species is not found", trim(ic_name)
                    stop
                else if (count.eq.3.and.tmp_index(2).eq.0) then
                    print*, "  index of reference species is not found", trim(sname)
                    stop
                endif
                
                ! change concentration
                js = tmp_index(1) ! index in concentration_gas_all
                tmp = concentration_gas_all(js) ! old concs
                if (count.eq.2) then
                    concentration_gas_all(js) = tmp_conc
                else if (count.eq.3) then
                    concentration_gas_all(js) = tmp_conc * concentration_gas_all(tmp_index(2))
                endif
                if (ssh_standalone) write(*,*) "  Old/New concentration: ",trim(species_name(js)),tmp,concentration_gas_all(js)
                if (ssh_logger) write(logfile,*)"  Old/New concentration: ",trim(species_name(js)),tmp,concentration_gas_all(js)
                
            endif
        end do
    endif
    
    ! genoa read constant aerosol profiles
    ncst_aero = 0
    if (cst_aero_file .ne. "---") then ! with input cst aero file

        if (ssh_standalone) write(*,*) 'Particle concentration constant file : ', trim(cst_aero_file)
        if (ssh_logger) write(logfile,*) 'Particle concentration constant file : ', trim(cst_aero_file)

        open(unit = 24, file = cst_aero_file, status = "old")
        count = 0
        ierr = 0
        nline = 0
        do while(ierr .eq. 0)
           read(24, *, iostat=ierr) tmp_name
           if (ierr == 0) then
              tmp_name = adjustl(tmp_name)
              if (trim(tmp_name) == "" .or. tmp_name(1:1) .eq. "#") then
                 nline = nline + 1 ! Comment line
              else
                 count = count + 1
              end if
           end if
        end do

        ! number of constant aerosol species
        ncst_aero = count
        if (.not.allocated(cst_aero_index)) allocate(cst_aero_index(ncst_aero))
        if (.not.allocated(cst_aero)) allocate(cst_aero(ncst_aero,N_sizebin,nt))
        if (.not.allocated(tmp_read)) allocate(tmp_read(nt))
        tmp_read = 0.0

        rewind 24
        do s = 1, nline
           read(24, *) ! read comment lines.
        end do
        do s= 1, count
           ! aerosol name, number of sizebin, concentrations per simulation time step
           read(24,*) ic_name, i, (tmp_read(k), k = 1, nt)
           ind = 0
           ! find index in aerosol list
           do js = 1, N_aerosol
               if (aerosol_species_name(js) .eq. ic_name) then
                   cst_aero_index(s) = js
                   do k = 1, nt
                        ! for now assume N_size == N_sizebin.
                        ! Need to be modified if use for external mixing
                        if (i.gt.N_sizebin) then
                            if (ssh_standalone) write(*,*) 'Error: Input number of bin >= N_sizebin', &
                                  i,N_sizebin,ic_name,' ',trim(cst_aero_file)
                            if (ssh_logger) write(logfile,*) 'Error: Input number of bin >= N_sizebin', &
                                  i,N_sizebin,ic_name,' ',trim(cst_aero_file)
                            stop
                        endif
                        cst_aero(s,i,k) = tmp_read(k)
                   enddo
                   ind = 1
              endif

              if (ind  == 1) then
                if (ssh_standalone) write(*,*) 'aerosol ',ic_name,' with index',js , &
                    ' Concentrations are given at each time step.'
              if (ssh_logger) write(logfile,*) 'aerosol ',ic_name,' with index',js , &
                    ' Concentrations are given at each time step.'
                exit
              endif
           enddo
           if (ind .eq. 0) then
              if (ssh_standalone) write(*,*) "Error: wrong aerosol species name is given ",cst_aero_file, ic_name
              if (ssh_logger) write(logfile,*) "Error: wrong aerosol species name is given ",cst_aero_file, ic_name
           endif
        enddo
        close(24)
    endif

    ! genoa read aerosol structure file
    if (aerosol_structure_file.ne."---") then
      open(unit = 25, file = aerosol_structure_file, status = "old")
        ! count comment lines and the number of input aerosol structures
        count = 0 
        ierr = 0
        nline = 0
        do while(ierr .eq. 0)
           read(25, *, iostat=ierr) tmp_name
           if (ierr == 0) then
              tmp_name = adjustl(tmp_name)
              if (trim(tmp_name) == "" .or. tmp_name(1:1) .eq. "#") then
                 nline = nline + 1 ! Comment line
              else
                 count = count + 1
              end if
           end if
        end do
        ! read
        rewind 25
        do s = 1, nline
           read(25, *) ! read comment lines.
        end do
        if (.not.allocated(tmp_fgls)) allocate(tmp_fgls(60)) !read unifac founctional groups (60)
        do s= 1, count
           ! name, unifac groups
           read(25,*) ic_name, (tmp_fgls(k), k = 1, 60)
           do js=1, N_aerosol
              if (aerosol_species_name(js).eq.ic_name) then
                 ! update smiles in the format &1.23E+02&1.00E+01&...
                 smiles(js)='' !init
                 do k=1, 60 !update
                    write(char1,'(A1,ES8.2)') "&",tmp_fgls(k)
                    smiles(js)=trim(smiles(js))//trim(char1)
                 enddo
                 !print*,aerosol_species_name(js),smiles(js)
              endif
           enddo
        enddo
      close(25)
    endif

    ! init
    iRO2_cst = 0  ! index of RO2 pool
    ncst_gas = 0  ! number of gas species that are constants
    nRO2_chem = 0 ! number of RO2 species
    nRO2_group = 1! number of RO2 group
    
    ! genoa read contant gasoes profile
    if (cst_gas_file.ne."---") then
        open(unit = 34, file = cst_gas_file, status = "old")

        if (ssh_standalone) write(*,*) 'Gas-phase concentration constant file : ', trim(cst_gas_file)
        if (ssh_logger) write(logfile,*) 'Gas-phase concentration constant file : ', trim(cst_gas_file)

        count = 0 
        ierr = 0
        nline = 0
        do while(ierr .eq. 0)
           read(34, *, iostat=ierr) tmp_name
           if (ierr == 0) then
              tmp_name = adjustl(tmp_name)
              if (trim(tmp_name) == "" .or. tmp_name(1:1) .eq. "#") then
                 nline = nline + 1 ! Comment line
              else
                 count = count + 1
              end if
           end if
        end do

        ! number of species that keep as constants
        ncst_gas =count
        if (.not.allocated(cst_gas_index)) allocate(cst_gas_index(ncst_gas)) ! index of unchanged species
        if (.not.allocated(cst_gas)) allocate(cst_gas(ncst_gas,nt))
        if (.not.allocated(tmp_read)) allocate(tmp_read(nt))
        tmp_read = 0.0

        rewind 34
        do s = 1, nline
           read(34, *) ! read comment lines.
        end do
        do s= 1, count
           read(34,*) ic_name, (tmp_read(k), k = 1, nt)
           ind = 0
           do js = 1, N_gas

              if (species_name(js) .eq. ic_name) then
                 cst_gas_index(s)=js

                 if (ic_name .eq. 'RO2pool') then
                    iRO2_cst = s
                    if (ssh_standalone) write(*,*) 'Index for cst background RO2 pool in cst_gas list', iRO2_cst
                    if (ssh_logger) write(logfile,*) 'Index for cst background RO2 pool in cst_gas list', iRO2_cst
                 endif

                 do k =1 ,nt
                    !cst_gas(s)=concentration_gas_all(js)
                    cst_gas(s,k) = tmp_read(k)
                 enddo
                 ind = 1
                 !print*, ic_name,tmp_read
              endif

              if (ind == 1) then
                if (ssh_standalone) write(*,*) ic_name,'with index',cst_gas_index(s), & 
                    ' Concentrations are given at each time step.'
                if (ssh_logger) write(logfile,*) ic_name,'with index',cst_gas_index(s), &
                    ' Concentrations are given at each time step.'

                exit
              endif
           enddo

           if (ind .eq. 0) then
              if (ssh_standalone) write(*,*) "Error: wrong species name is given in ",&
                trim(cst_gas_file),', species ', trim(ic_name)
              if (ssh_logger) write(logfile,*) "Error: wrong species name is given in ",&
                trim(cst_gas_file),', species ', trim(ic_name)
              stop
           endif
        enddo
        close(34)
        
        ! check iRO2_cst
        if (iRO2_cst.eq.0) then ! no input
            if (tag_RO2.eq.2 .or. tag_RO2.eq.3) then
                print*, "Modify tag_RO2 cuz not iRO2_cst", iRO2_cst, tag_RO2
                tag_RO2 = tag_RO2 - 2
                if (tag_RO2.eq.0) then
                    print*, "Update: No RO2-RO2 reactions", tag_RO2
                else ! tag_RO2.eq.1
                    print*, "Update: Only generated RO2-RO2 reactions", tag_RO2
                endif
            endif
        endif
    else
        !print*,'No input cst gas file.'
        if ( .not. allocated(cst_gas_index)) allocate(cst_gas_index(ncst_gas)) ! index of unchanged species
        if ( .not. allocated(cst_gas)) allocate(cst_gas(ncst_gas,nt))
    endif

    ! genoa RO2 species list - build RO2pool
    ! update species_list_file with chemID
    if (RO2_list_file.ne."---".and.chemID.ne."-") RO2_list_file = trim(adjustl(RO2_list_file))//"/"//trim(chemID2)//".RO2"
       
    if (tag_RO2 .eq. 1 .or. tag_RO2 .eq. 3) then ! need to have the list
      if (RO2_list_file.ne."---") then
      
        if (ssh_standalone) write(*,*) 'Read RO2 list file : ', trim(RO2_list_file)," tag_RO2 = ",tag_RO2
        if (ssh_logger) write(logfile,*) 'Read RO2 list file : ', trim(RO2_list_file)," tag_RO2 = ",tag_RO2

        open(unit = 34, file = RO2_list_file, status = "old")
        count = 0 
        ierr = 0
        nline = 0
        do while(ierr .eq. 0)
           read(34, *, iostat=ierr) tmp_name
           if (ierr == 0) then
              tmp_name = adjustl(tmp_name)
              if (trim(tmp_name) == "" .or. tmp_name(1:1) .eq. "#") then
                 nline = nline + 1 ! Comment line
              else
                 count = count + 1
              end if
           end if
        end do

        nRO2_chem = count
        if ( .not. allocated(RO2index)) allocate(RO2index(nRO2_chem,2)) ! index of RO2 and group id
        
        rewind 34
        do s = 1, nline
           read(34, *) ! read comment lines.
        end do
        do s= 1, count
           read(34,*) ic_name, k
           ind = 0
           do js = 1, N_gas
              if (species_name(js) .eq. ic_name) then
                 RO2index(s,1)=js ! isps
                 RO2index(s,2)=k  ! group id
                 ! update no. group if need
                 if (k.gt.nRO2_group) nRO2_group = k
                 ind = 1
              endif
              if (ind == 1) exit
           enddo
           if (ind .eq. 0) then
              if (ssh_standalone) write(*,*) "Error: wrong species name is given in ",&
                trim(RO2_list_file)," ",trim(ic_name)
              if (ssh_logger) write(logfile,*) "Error: wrong species name is given in ",&
                trim(RO2_list_file)," ", trim(ic_name)
              stop
           endif
        enddo
        close(34)
        
        ! check nRO2_group with reaction list
        if (nRO2_chem.gt.0) then
          ind = minval(TB_rcn(:,2)) !number read from reaction list
          if (ind.gt.0 .and. abs(ind).gt.nRO2_group) then
            if (ssh_standalone) write(*,*)"Error: RO2 group no. < RO2 index read from reaction file", &
                                 abs(ind), nRO2_group
            if (ssh_logger) write(logfile,*) "Error: RO2 group no. < RO2 index read from reaction file", &
                                 abs(ind), nRO2_group
            stop
          endif
        endif

        ! check nRO2_group with iRO2_cst
        if (tag_RO2.gt.1.and.nRO2_group.gt.1.and.iRO2_cst.gt.0) then
          if (ssh_standalone) write(*,*)   "Error: RO2 group no. > 1 & iRO2_cst exists", nRO2_group
          if (ssh_logger) write(logfile,*) "Error: RO2 group no. > 1 & iRO2_cst exists", nRO2_group
          stop
        endif
        
        ! for output
        if ( .not. allocated(RO2out_index)) allocate(RO2out_index(nRO2_group))
        
      else ! no RO2 file can be read
          if (ssh_standalone) write(*,*) "Error: RO2 list info is needed but not given."
          if (ssh_standalone) write(*,*) "Please check RO2_list_file in the namelist. tag_RO2 = ",&
              tag_RO2, " ",trim(RO2_list_file)
          if (ssh_logger) write(logfile,*) "Error: RO2 list info is needed but not given."
          if (ssh_logger) write(logfile,*) "Please check RO2_list_file in the namelist. tag_RO2 = ",&
              tag_RO2, " ",trim(RO2_list_file)
          stop
      endif
      
    else ! no need to read RO2 file  
      if (ssh_standalone) write(*,*) "RO2 list is given but not read because tag_RO2 = ",&
        tag_RO2, " ", trim(RO2_list_file)
      if (ssh_logger) write(logfile,*) "RO2 list is given but not read because tag_RO2 = ",&
        tag_RO2, " ", trim(RO2_list_file)
    endif
    
    ! get RO2 output index
    if (tag_RO2.ne.0) then
        ! YK
        if ( .not. allocated(RO2out_index)) allocate(RO2out_index(nRO2_group))
        RO2out_index = 0 ! init
        ! find RO2 index
        do k = 1, nRO2_group
          ! convert the integer to a string => RO2pool1 to RO2pooln, RO2pool => total RO2 pool
          write(ic_name, '("RO2pool", I0)') k
          
          ! check index in species list
          do js = 1, N_gas
            if (INDEX(species_name(js),trim(ic_name)) /= 0) then
              RO2out_index(k) = js
              if (ssh_standalone) write(*,*) 'RO2 pool id found for group ',k, RO2out_index(k)
              if (ssh_logger) write(logfile,*) 'RO2 pool id found for group ',k, RO2out_index(k)
              exit
            endif
          enddo
          ! check RO2 index - not found
          if ( RO2out_index(k) .eq. 0) then
            if (ssh_standalone) write(*,*)   "Not found RO2 pool id in species list: ", k, trim(ic_name), trim(species_list_file)
            if (ssh_logger) write(logfile,*) "Not found RO2 pool id in species list: ", k, trim(ic_name), trim(species_list_file)
          end if
        enddo
    endif

    ! for twostep solver input
    if ( .not. allocated(RO2index)) allocate(RO2index(nRO2_chem,2)) ! index of RO2
    if ( .not. allocated(RO2out_index)) allocate(RO2out_index(nRO2_group))
    
    ! genoa read ref and pre concentrations
    if (nref_file.gt.0) then

      allocate(ref_soa(nref_file,nout_total,nt)) ! save ref soa conc
      ref_soa = 0.d0
        
      do ind = 1, nref_file
        if (trim(ref_conc_files(ind)).ne."---") then
            if (ssh_standalone) write(*,*) ind,'Read ref SOA file: ',trim(ref_conc_files(ind))
            if (ssh_logger) write(logfile,*) ind,'Read ref SOA file: ',trim(ref_conc_files(ind))
            
            open(unit = 35, file = trim(ref_conc_files(ind)), status = "old")

            ! check number of values in the list
            count = 1 
            ierr = 0
            
            read(35,*, iostat=ierr) ! read the init conc
            do while(ierr.eq.0.and.count.le.nt) 
               read(35,*, iostat=ierr) (ref_soa(ind,k,count), k=1, nout_total)
               if (ierr == 0) count = count + 1
            enddo
            
            if (ierr.ne.0) then
                write(*,*) "Error: can not read the ref conc. ", count," ",trim(ref_conc_files(ind))
                stop
            endif

            ! check if count number == nt
            if (nt.ne.count-1) then
                write(*,*) "Error: not the same number in ref conc. ", nt, count, " ",trim(ref_conc_files(ind))
                stop
            endif
            
            ! print info
            if (ssh_standalone) write(*,*) "Read ref soa concentration.", nt, count," ",trim(ref_conc_files(ind))
            if (ssh_logger) write(logfile,*) "Read ref soa concentration.", nt, count," ",trim(ref_conc_files(ind))
            
            close(35)
        end if
      end do
    else
        if (ssh_standalone) write(*,*) "Not read ref soa concentration."
        if (ssh_logger) write(logfile,*) "Not read ref soa concentration."
    endif
    
    ! save for error computation
    if (tag_genoa.ne.0.or.nout_total.gt.0) then                  
        allocate(total_soa(nt+1,nout_total))
        total_soa = 0d0
    endif
    
! genoa


    if (ssh_standalone) write(*,*) "=========================finish read inputs file======================"
    if (ssh_logger) write(logfile,*) "=========================finish read inputs file======================"

    if (allocated(tmp_aero))  deallocate(tmp_aero)
    ! genoa
    if (allocated(tmp_read))  deallocate(tmp_read)
    if (allocated(tmp_fgls))  deallocate(tmp_fgls)
    if (allocated(tmp_index)) deallocate(tmp_index)
    
  end subroutine ssh_read_inputs


  subroutine ssh_read_reaction_file()
  
  implicit none
      
  integer :: i, j, k, js, ierr, iarrow, nrct, npdt
  integer :: ircn, iknc, irct, ipdt, ipho, ipho_t, itb
  integer :: iext, ifall, iro2, ipdt_r, ipdt_f, ihet, iwall, iirdi !counts
  integer :: start, finish
  integer :: finish_plus, finish_minus
  integer :: index_start
  double precision, dimension(3) :: ratios ! use to read ratio function
  
  character (len=800) :: line, subline, combinedline, previousline, nextline
  character (len=80) :: tmp_line
  character (len=80) :: sname
  
  integer, parameter :: ntmp = 12 ! max number of coefficients in reaction file
  double precision :: a_tmp(ntmp), tmp_in ! read coefficients
  character(3), parameter :: TBs(5) = ["O2 ", "H2O", "M  ", "N2 ", "H2 "] ! do not change the order ! used in kinetic
  double precision :: xlw
  
  ! set up parameters for gas-phase chem - GECKO
  ! Use perfect gas law :sumc=(pres*6.022E+23)/(8.32*temp)
  ! SumMc = 2.5E19 ! default
  SumMc = Pressure * 7.243d16 / Temperature
  ! Water concentration (molec/cm3)
  ! YlH2O = 0.d0
  xlw = humidity
  YlH2O = 29.d0*SumMc*xlw/(18.d0+11.d0*xlw)
  
  ! Open reactions file
  open(unit=11, file=reaction_list_file, status='old')
    ! get size
    ircn = 0 ! count the number of reactions
    irct = 0 ! count the number of reactants
    ipdt = 0 ! count the number of products
    ipdt_r= 0 ! count the number of products with ratio
    ipdt_f= 0 ! count the number of products with ratio as a function
    iknc = 0  ! count the number of kinetics
    ipho  = 0 ! count the number of photolysis read from file
    ipho_t= 0 ! count the number of photolysis read from line
    itb   = 0 ! count the number of TB reactions
    iro2  = 0 ! count the number of RO2-RO2 reactions
    ifall = 0 ! count the number of falloff reactions 
    iext  = 0 ! count the number of other types of reactions
    ihet  = 0 ! count the number of heterogeneous reactions
    iwall  = 0 ! count the number of wall loss
    iirdi = 0 ! count the number of irreversible dicarbonyl
    
    do
        read(11, '(A)', iostat=ierr) line
        if (ierr /= 0) exit ! no line
        line = adjustl(line)
        i = len_trim(line)
        if(i < 3 .or. line(1:1) == '%' .or. line(1:1) == '#' .or. line(1:3) == 'END') cycle ! read comments
        iarrow = index(line, '->') ! arrow position
        if(iarrow > 0) then ! Find a reaction
            ircn = ircn + 1     ! no.reaction
            call count_delimiter(line(1:iarrow - 1), ' + ',i)
            irct = i + irct     ! no.reactant
            call count_products(line(iarrow + 2:), i, j, k)
            ipdt = i + ipdt     ! no.product
            ipdt_r = ipdt_r + j ! no.product w ratio
            ipdt_f = ipdt_f + k ! no.ratio as a function
        ! find all kinetic keywords
        else if (line(1:7) == "KINETIC") then
            iknc = iknc + 1 !no.kinetic
            read(line(8:), '(A20)') sname ! keywords
            ! Check 1st keyword (two keywords in a line)
            if (index(sname,'TB') /= 0) then
                itb = itb + 1 ! third body
            else if (index(sname,'RO2') /= 0) then ! RO2-RO2
                iro2 = iro2 + 1
            endif
            ! Check 2nd keyword or the only keyword
            if(index(sname,'PHOT') /= 0) then
                ipho = ipho + 1
                ! read in-line tabulation for photolysis
                if(index(sname,'PHOTOLYSIS') /= 0) ipho_t = ipho_t + 1
            else if(index(sname,'FALLOFF') /= 0) then ! falloff
                ifall = ifall + 1
            else if(index(sname,'EXTRA') /= 0) then ! extra and specified kinetics
                iext = iext + 1
            else if(index(sname,'HETERO') /= 0) then ! heterogenous reactions
                ihet = ihet + 1
            else if(index(sname,'WALL') /= 0) then ! wall loss
               iwall = iwall + 1
            else if(index(sname,'IRDICARB') /= 0) then ! irreversible dicarbonyl
                iirdi = iirdi + 1               
                
            endif 
        end if
    end do
    
    n_reaction = ircn
    n_photolysis = ipho

    ! if (ssh_standalone) then
        print*, "Finish 1st time reading reactions. Read No.reactions: ",n_reaction
        if (ircn /= iknc) print*, "No.kinetics /= No.reactions: ", iknc, n_reaction
    ! endif
    
    ! allocate arrays
    allocate(index_RCT(irct,2))   ! index of species in species list and reaction list
    allocate(index_PDT(ipdt,3))   ! index of species in species list, reaction list, ratio list
    allocate(index_PDT_ratio(ipdt_f,3))
    allocate(ratio_PDT(ipdt_r))   ! branching ratios
    
    ! kinetics
    allocate(photo_rcn(ipho,2)) ! code of photolysis

    allocate(photo_ratio_read(ipho_t,nsza)) ! read ratio
    allocate(szas(nsza))
    
    allocate(TB_rcn(itb+iro2,2))  ! TBs & RO2-RO2 2. index(TB+,RO2-); 1. reaction index
    allocate(Arrhenius(iknc,3))   ! kinetic rate contants
    allocate(fall_rcn(ifall))     ! 1. spec index; 2: reaction index
    allocate(fall_coeff(ifall,8)) ! coefficients for troe reactions
    allocate(extra_rcn(iext))     ! code and reaction index for EXTRA type of reactions
    allocate(extra_coeff(iext,8)) ! coefficients for extra reactions
    allocate(hetero_rcn(ihet))    ! heterogenous reactions
    allocate(hetero_ind(ihet))    ! coefficients for heterogeneous reactions
    allocate(wall_rcn(iwall))     ! wall loss
    allocate(wall_coeff(iwall,5)) ! coefficients for wall loss
    allocate(irdi_rcn(iirdi))     ! irreversible dicarbonyl
    allocate(irdi_ind(iirdi))      ! index for irreversible dicarbonyl
    
    ! for two-step: y,w,rk,prod,loss,dw
    allocate(gas_yield(n_gas))
    allocate(rcn_rate(iknc))
    allocate(kinetic_rate(iknc))
    allocate(chem_prod(n_gas))
    allocate(chem_loss(n_gas))
    allocate(drv_knt_rate(irct))
    
    if (ssh_standalone) print*, "Finish allocating arrays for reactions."

    ! init
    index_RCT = 0
    index_PDT = 0
    index_PDT_ratio = 0
    ratio_PDT = 0.d0
    
    photo_rcn = 0
    TB_rcn = 0
    fall_rcn = 0
    extra_rcn = 0
    hetero_rcn = 0
    
    Arrhenius = 0.d0
    fall_coeff = 0.d0
    extra_coeff = 0.d0
    hetero_ind = 0
    irdi_ind = 0    

    !!!! (To do) Move to namelist
    if (nsza == 2) then
       szas = [0.0d0, 9.0d1]
    else if (nsza == 11) then
       szas = [0d0,1d1,2d1,3d1,4d1, &
            5.0d1,6d1,7d1,7.8d1,8.6d1,9d1]
    endif

    rewind(11)  ! Reset the file position

    ! recount size
    ircn = 0 ! count the number of reactions
    iknc = 0  ! count the number of kinetics

    irct = 0 ! count the number of reactants
    ipdt = 0 ! count the number of products
    ipdt_r =0 ! count the number of products with prodcut ratio
    ipdt_f =0 ! count the number of products with prodcut ratio as a function

    ipho  = 0 ! count the number of photolysis read from a file
    ipho_t= 0 ! count the number of photolysis read from a tabulation
    itb   = 0 ! count the number of TB reactions
    iro2  = 0 ! count the number of RO2-RO2 reactions
    ifall = 0 ! count the number of falloff reactions 
    iext  = 0 ! count the number of other types of reactions
    ihet = 0  ! count the number of heterogenous reactions
    iwall = 0  ! count the number of wall loss
    iirdi = 0  ! count the number of irreversible dicarbonyl
    
    ! Read reactions and their rates
    do
        read(11, '(A)', iostat=ierr) line
        if (ierr /= 0) exit ! no line

        !! YK TEST
        ! ! Combine the lines having the line break symbol "//"
        ! previousline = trim(line)
        ! do while (index(line, "//") .ne. 0)

        !    index_start = index(line, "//")
           
        !    ! Read a next line
        !    read(11, '(A)', iostat=ierr) line

        !    combinedline = trim(previousline(1:index_start - 1)) // trim(line)
        !    previousline = combinedline
        !    line = combinedline
        ! enddo
        
        line = adjustl(line)
        i = len_trim(line)
        if(i < 3 .or. line(1:1) == '%' .or. line(1:1) == '#' .or. line(1:3) == 'END') cycle ! read comments
        iarrow = index(line, '->') ! arrow position
        if (iarrow > 0) then ! read reactions
            ircn = ircn + 1 ! reaction index
            nrct = 0 ! Reset count reactants
            npdt = 0 ! Reset count products
            !if (ssh_standalone) print*, "Reading reaction: ", ircn, " / ", n_reaction

            ! Extract reactants
            subline = line(1:iarrow-1)
            start = 1
            do
                ! get sname
                finish = index(subline(start:), ' + ')
                if (finish == 0) then ! last species
                    sname = adjustl(subline(start:))
                else
                    sname = adjustl(subline(start:start+finish-2))
                    start = start + finish + 1
                endif
                
                ! get sname index
                irct = irct + 1 ! index in index_RCT
                nrct = nrct + 1 ! count
                ! get index in species list
                do js = 1, N_gas
                  if (species_name(js) .eq. sname) then
                    index_RCT(irct,1)=ircn ! no.reaction
                    index_RCT(irct,2)=js   ! no.species
                    exit
                  endif
                enddo
                ! check index
                if (index_RCT(irct,1) .eq. 0) then
                    print*, 'reactants not found in species list: ',trim(sname)
                    stop
                endif    
                ! check exit
                if (finish == 0) exit
            end do

            ! Extract products & branching ratios
            subline = line(iarrow+2:)
            ! If there is no product - continue
            if (len_trim(subline) == 0 ) cycle
            
            ! Start to process products
            start = 1
            do
                ! get sname in the format "xxx" or "5E-2 xxx" 
                !  finish = index(subline(start:), ' + ')
                !! read both '+' and '-' signs (YK)
                finish_plus = index(subline(start:), ' + ')
                finish_minus = index(subline(start:), ' - ')
                if (finish_plus == 0 .and. finish_minus == 0) then
                   finish = 0
                else if (finish_plus == 0 .and. finish_minus .ne. 0) then
                   finish = finish_minus
                else if (finish_plus .ne. 0 .and. finish_minus == 0) then
                   finish = finish_plus
                else
                   finish = min(finish_plus, finish_minus)
                endif !! YK
                
                if (finish == 0) then ! last species
                    tmp_line = adjustl(subline(start:))
                else
                    tmp_line = adjustl(subline(start:start+finish-2))
                    start = start + finish + 1
                endif

                ipdt = ipdt + 1 ! index in index_RCT
                npdt = npdt + 1 ! count
                ! get species num and index
                if (len_trim(tmp_line) == 0) then ! update for no product
                    index_PDT(ipdt,1) = ircn
                else ! get index in species list
                    call get_prodcut_and_ratios(tmp_line, ratios, sname, i)
                    if (i /= 0) then ! Find ratio
                        ipdt_r = ipdt_r + 1 ! count
                        ! update index
                        if (i == 3) then ! Find ratio, k1, k2 - use a function
                            ipdt_f = ipdt_f + 1 ! count function
                            index_PDT(ipdt,3) = -ipdt_f  ! Record index
                            index_PDT_ratio(ipdt_f,1) = ipdt_r ! ratio index
                            index_PDT_ratio(ipdt_f,2) = INT(ratios(2))+ircn-1 ! k1
                            index_PDT_ratio(ipdt_f,3) = INT(ratios(3))+ircn-1 ! k2
                        else
                            index_PDT(ipdt,3) = ipdt_r  ! ratio as a number
                        endif
                        ! update ratio
                        ratio_PDT(ipdt_r) = ratios(1)
                    endif
                    do js = 1, N_gas
                      if (species_name(js) .eq. sname) then
                        index_PDT(ipdt,1)=ircn
                        index_PDT(ipdt,2)=js
                        exit
                      endif
                    enddo
                    ! check index
                    if (index_PDT(ipdt,1) .eq. 0 .and. &
                         (trim(sname) .ne. "O2" .and. &
                         trim(sname) .ne. "H2") ) then
                       ! if the species name is not found in species_list
                       ! an error message raises.
                       ! 'O2' and 'H2' are allowed even though they are not
                       ! listed.
                       print*, trim(line)
                       print*, 'product not found in species list: ', &
                            trim(sname), ' ', trim(tmp_line)
                       stop
                    endif
                endif
                ! check exit
                if (finish == 0) exit
            end do
                        
        elseif (line(1:7) == "KINETIC") then ! read kinetics

            iknc = iknc + 1 ! Count kinetics
            ! Update ircn if need
            if (ircn /= iknc) then
                !print*, 'Update index: ',ircn,iknc
                ! Update ircn index in index_RCT & index_PDT
                do js=1, nrct
                    index_RCT(irct-nrct+js,1) = iknc
                enddo
                do js=1, npdt
                    index_PDT(ipdt-npdt+js,1) = iknc
                enddo
                ircn = iknc
            endif

            subline = adjustl(line(8:))
            do js=1, len_trim(subline)
                if (subline(js:js) .eq. " ") exit ! found space to locate 1st keyword
            enddo 
            
            read(subline(1:js), '(A20)', iostat=ierr) sname ! keyword
            js = js + 1 ! cursor position
            
            ! read 1st keyword TBs & RO2s
            if (trim(sname) == "TB".or.trim(sname) == "RO2") then
            
                ! read info related to TBs & RO2s
                itb = itb + 1 ! both TB & RO2
                TB_rcn(itb,1) = iknc
                
                ! read second input: for TB is a string, for RO2 is a number
                subline = adjustl(subline(js:))
                do js=1, len_trim(subline)
                    if (subline(js:js) .eq. " ") exit ! found " " loc
                enddo
                
                if (trim(sname) == "TB") then ! read TB
                    read(subline(1:js), '(A20)', iostat=ierr) sname
                    if (ierr /= 0) then
                        print*, "Can not read TB name from line", trim(subline), iknc
                        stop
                    endif
                    ! get index in TBs
                    j = 0
                    do i = 1, 5 ! Need to update here the size of TBs
                        if (trim(TBs(i)) == trim(sname)) then
                            j = i
                            exit
                        end if
                    end do
                    ! check & record index
                    if (j.eq.0) then
                        print*, 'Not find TB in TBs: ', trim(sname)
                        stop
                    else
                        TB_rcn(itb, 2) = j ! index of TB
                    endif
                else ! read RO2 index
                    iro2 = iro2 + 1
                    read(subline(1:js), *, iostat=ierr) j
                    if (ierr /= 0) then
                        print*, "Can not read RO2 integer index from line", trim(subline(1:js)), &
                              js
                        stop
                    endif
                    TB_rcn(itb, 2) = -j ! index of RO2 - negative!
                endif
                
                ! read 2nd keyword
                js = js + 1
                subline = adjustl(subline(js:))
                do js=1, len_trim(subline)
                    if (subline(js:js) .eq. " ") exit ! found " "
                enddo
                read(subline(1:js), '(A20)', iostat=ierr) sname
            end if
            
            ! Read all coefficients - up to No. ntmp
            a_tmp = 0d0
            do i = 1, ntmp
                ! adjust the rest of line
                subline = adjustl(subline(js:))
                ! get 1st space loc
                js = index(subline, ' ')
                if (js /= 0) then ! find space
                    ! check length of str
                    if (len_trim(adjustl(subline(1:js))).eq.0) exit
                    ! read number
                    read(subline(1:js), *, iostat=ierr) tmp_in
                    if (ierr .eq. 0) a_tmp(i) = tmp_in ! read
                else
                    exit ! no space
                endif
            enddo
            
            ! check no. read coefficients
            if (i-1.lt.1) then
                print*, 'read no kinetic coefficient', iknc
                stop
            else
                finish = i - 1
            endif
            
            ! assign values for different types of kinetics
            if (trim(sname) == "ARR") then
                ! read arrhenius contants - can be 3 or 4 (with a ratio from GENOA reduction)
                if (finish.eq.3 .or. finish.eq.4) then
                    do i = 1, 3
                        Arrhenius(iknc,i) = a_tmp(i)
                    enddo
                    if (finish.eq.4) then ! multiply read ratio
                        Arrhenius(iknc,1) = Arrhenius(iknc,1) * a_tmp(4)
                    endif
                else
                    print*, "Error:  ARR read no. coeff should be 3 or 4", &
                            iknc, finish, a_tmp(1:finish)
                    stop
                endif
             elseif (trim(sname) == "PHOT") then
                ! read photolysis
                ipho = ipho + 1
                photo_rcn(ipho,1) = iknc

                if (finish .eq. 2) then ! read only index & ratio
                    photo_rcn(ipho,2) = int(a_tmp(1))
                    Arrhenius(iknc,1) = a_tmp(2)
                    
                else if (finish.eq.nsza .or. finish.eq.nsza+1) then ! in-list tabulation
                    ! id: negative in photo_ratio_read
                    ipho_t = ipho_t + 1
                    photo_rcn(ipho,2) = -ipho_t
                    
                    ! read tabulation
                    do i = 1, nsza
                        photo_ratio_read(ipho_t,i) = a_tmp(i)
                    enddo
                    
                    ! ratio
                    if (finish.eq.nsza+1) then ! read
                        Arrhenius(iknc,1) = a_tmp(nsza+1)
                    else ! default
                        Arrhenius(iknc,1) = 1d0
                    endif
                else
                    print*, "Error: PHOT read no. coeff not 2." &
                            , iknc, finish, a_tmp(1:finish)
                    stop
                endif

             elseif (trim(sname) == "PHOTOLYSIS") then
                ! read photolysis
                ipho = ipho + 1
                photo_rcn(ipho,1) = iknc

                ! id: negative in photo_ratio_read
                ipho_t = ipho_t + 1
                photo_rcn(ipho,2) = -ipho_t

                if (finish.eq.nsza .or. finish.eq.nsza+1) then ! in-list tabulation
  
                    ! read tabulation
                    do i = 1, nsza
                        photo_ratio_read(ipho_t,i) = a_tmp(i)
                    enddo
                    
                    ! ratio
                    if (finish.eq.nsza+1) then ! read
                        Arrhenius(iknc,1) = a_tmp(nsza+1)
                    else ! default
                        Arrhenius(iknc,1) = 1d0
                    endif
                else
                    print*, "Warning: PHOTOLYSIS read no. coeff not nsza/nsza+1." &
                         , iknc, finish, a_tmp(1:finish), &
                         "nsza in namelist): ", nsza
                    if (finish .lt. nsza) then
                       ! put the 1st value from the tabulation to all angles 
                       do i = 1, nsza
                          photo_ratio_read(ipho_t,i) = a_tmp(1)
                       enddo
                       Arrhenius(iknc,1) = 1d0
                    else
                       write(*,*) "Error: something wrong in nsza."
                       stop
                    endif
                endif

                
            elseif (trim(sname) == "FALLOFF") then
                ! read falloff 
                ifall = ifall + 1
                fall_rcn(ifall) = iknc

                ! arrhenius(3) + FALLOFF(7) + ratio = 11
                if (finish.eq.10 .or. finish.eq.11) then 
                    do i = 1, 3 ! arrhenius coefficents
                        Arrhenius(iknc,i) = a_tmp(i)
                    enddo
                    do i = 1, 7 ! falloff coefficents
                        fall_coeff(ifall,i) = a_tmp(i+3)
                    enddo
                    ! ratio
                    if (finish .eq. 11) then
                        fall_coeff(ifall,8) = a_tmp(11)
                    else ! default
                        fall_coeff(ifall,8) = 1d0
                    endif
                else
                    print*, "Error: FALLOFF read no.coeff should be 10 or 11", &
                            iknc, finish, a_tmp(1:finish)
                    stop
                endif

            elseif (trim(sname) == "EXTRA") then
              ! read extra or specified kinetics
              iext = iext + 1
              extra_rcn(iext) = iknc
                
              js = int(a_tmp(1)) ! get label
              extra_coeff(iext,1) =  js

              SELECT CASE (js)

                CASE(91) !MCM: 91 C1 C2 C3 C4(optional)
                  if (finish.eq.4 .or. finish.eq.5) then
                    do i = 2, 4
                        extra_coeff(iext,i) = a_tmp(i)
                    enddo
                    ! ratio
                    if (finish.eq.5) then
                        Arrhenius(iknc,1) = a_tmp(5)
                    else ! default
                        Arrhenius(iknc,1) = 1d0
                    endif
                  else
                    print*, "Error: EXTRA 91 read no. not 4/5 for MCM photolysis", &
                             iknc, finish,a_tmp(1:finish)
                    stop
                  endif
                  
                CASE(92, 93) !MCM: 92 C1 C2(optional)
                  if (finish.eq.2 .or. finish.eq.3) then
                    extra_coeff(iext,2) = a_tmp(2)
                    ! ratio
                    if (finish.eq.3) then
                        Arrhenius(iknc,1) = a_tmp(3)
                    else ! default
                        Arrhenius(iknc,1) = 1d0
                    endif
                  else
                    print*, "Error: EXTRA 92,93 read no. not 2/3 for MCM generic, complex rate coefficents", &
                              iknc,finish,a_tmp(1:finish)
                    stop
                  endif

                CASE(100, 200, 500, 501, 510, 502, 550) ! GECKO labels
                    start = finish ! get length - no. coff remain unread
                    if (finish.ge.4) then
                        do i = 2, 4
                            Arrhenius(iknc,i-1) = a_tmp(i)
                        enddo
                        start = start - 4
                    else
                        print*, "Error: EXTRA GECKO read no. < 4", &
                                iknc,finish,a_tmp(1:finish)
                        stop
                    endif
                    
                    j = 0 ! no. coeff need to be read
                    if (js.eq.200.and.start.ge.5) then ! 5 coef
                        j = 5
                    else if (js.eq.501.and.start.ge.3) then ! 3 coef
                        j = 3
                    else if (js.eq.550.and.start.ge.6) then ! 6 coef
                        j = 6
                    endif
                    ! read coef
                    do i = 1, j
                        extra_coeff(iext,i+1) = a_tmp(i+4) 
                    enddo
                    start = start - j
                    
                    ! read one ratio
                    if (start.eq.1) then
                        extra_coeff(iext,8) = a_tmp(j+5) ! last number
                        start = start - 1
                    else ! default
                        extra_coeff(iext,8) = 1d0
                    endif
                    
                    ! check if read all
                    if (start.ne.0) then
                        print*, "Error: EXTRA GECKO read no. coefs not match required: ", js, j, &
                                 iknc,finish,start,a_tmp(1:finish)
                        stop
                    endif
                    
                CASE(99, 10, 20) ! 99/10 C1 C2 (optional) 
                  ! Additional calculations that can be updated by the user.
                  ! Please ensure the label has not been used before if adding new labels. 
                  if (finish.eq.2 .or. finish.eq.3) then
                    extra_coeff(iext,2) = a_tmp(2)
                    ! ratio
                    if (finish.eq.3) then
                        Arrhenius(iknc,1) = a_tmp(3)
                    else ! default
                        Arrhenius(iknc,1) = 1d0
                    endif
                  else
                    print*, "Error: EXTRA read no. coeffs not 2/3", js, &
                             iknc,finish,a_tmp(1:finish)
                    stop
                  endif

                CASE DEFAULT
                  print*,"EXTRA label unknown: ",js
                  STOP
              END SELECT

            ! Find heterogenous reactions.  
            elseif (trim(sname) == "HETERO") then

              ihet = ihet + 1
              hetero_rcn(ihet) = iknc
              js = int(a_tmp(1)) ! get label
              hetero_ind(ihet) =  js

            ! Find wall loss
            elseif (trim(sname) == "WALL") then

              iwall = iwall + 1
              wall_rcn(iwall) = iknc
              do i = 1, 5 ! wall loss coefficents
                 wall_coeff(iwall,i) = a_tmp(i)
              enddo

            ! Find irreversible dicarbonyl
            elseif (trim(sname) == "IRDICARB") then

              iirdi = iirdi + 1
              irdi_rcn(iirdi) = iknc
              js = int(a_tmp(1)) ! get label
              irdi_ind(iirdi) = js
              
            else
              print*,"Error: Unknown kinetic keyword detected: ",trim(sname)
              stop
            END if
        else ! neither reaction nor kinetics, not read
            ! print to check
            print*, 'Error: Can not read line from reaction file: ',trim(line)
            stop
        endif
    end do

    close(11)
    
    print*,'Finished reading reaction file. Find No.rcn: ',ircn,' No.rct: ', &
          irct,' No.pdt: ',ipdt, ' No.pdt w ratio: ',ipdt_r, &
          ' No.pdt w ratio as a function: ', ipdt_f, &
          ' No.pho: ',ipho,' No.pho_tab: ',ipho_t ,' No.tb: ', &
          itb-iro2, ' No.fall: ',ifall, ' No.ext: ',iext, ' No.RO2-RO2: ', iro2
    write (*,*) "No. hetero: ", ihet
    write (*,*) "No. irdi: ", iirdi   
    
  end subroutine ssh_read_reaction_file
  
  subroutine ssh_read_photolysis_file()
  
    USE mod_cubicspline
  
    implicit none

    integer, dimension(:), allocatable :: ratio_id
    integer :: i,j,k,s,n,ierr, isize
    integer :: npho_r,npho_t, npho_tot
    integer :: tag, ind
    character (len=200) :: line
    character (len=40) :: ic_name, sname
    double precision :: a_tmp(nsza), tmp, tmp_inter(4,nsza-1)
    

    ! size    
    npho_t = size(photo_ratio_read,1) ! read from table
    npho_tot = size(photo_rcn,1) ! all photolysis
    npho_r = npho_tot - npho_t ! read from file

    ! allocate
    allocate(photo_ratio(npho_tot,nsza-1,4))  
    allocate(ratio_id(npho_r))
    ratio_id = 0 ! init
    photo_ratio = 0.d0
    
    ! orginize input id to get ratio_id
    isize = 0 ! size of current ratio_id
    n = 0 ! number in photo_ratio_read
    do i=1, npho_tot ! 
        !check if contains - tag = 1, if not tag = 0
        tag = 0
        j = photo_rcn(i,2) ! photolysis index
        if (j .gt. 0) then ! read ratio_id
          tag = 0
          do k = 1, isize
            if (j .eq. ratio_id(k)) then
                tag = 1 ! find id
                exit
            endif
          enddo
          if (tag.eq.0) then ! record
            isize = isize + 1
            ratio_id(isize) = j
          endif
        else !write photo_ratio
            n = n + 1
            ! interpolation
            call ssh_SPL3(nsza,szas,photo_ratio_read(-j,:),tmp_inter)
            ! save in photo_ratio
            do s=1, nsza-1
              do k=1,4
                photo_ratio(n,s,k) =  tmp_inter(k,s)
              enddo
            enddo
            ! update index in photo_rcn !!! negative
            photo_rcn(i,2) = -n
        endif
    enddo
    ! check no. tabulation
    if (npho_t .ne. n) then
        print*,'Error in ssh_read_photolysis_file(): tabulation number not match', &
                npho_t, n, npho_r, npho_tot
        stop
    endif
    ! print to check
    if (ssh_standalone) write(*,*) "Read in-line tabulation: ", npho_t
    if (ssh_logger) write(logfile,*) "Read in-line tabulation: ", npho_t

    ! read photolysis from file

    if (isize .gt. 0) then
      if (ssh_standalone) write(*,*) "Reading No.phot ",npho_r," from file: ", trim(photolysis_file)
      if (ssh_logger) write(logfile,*) "Reading No.phot ",npho_r," from file: ", trim(photolysis_file)
      
      n = 0 ! init for total read number
      
      open(unit=11, file=photolysis_file, status='old')
        ! Read file line by line
        do
          read(11, '(A)', iostat=ierr) line
          if (ierr /= 0) exit ! end file
            
          if (index(line, 'PHOT') > 0) then ! Check if the line contains 'PHOT'
            
            read(line, *, iostat=ierr) ic_name, sname, ind, i

            if (ierr == 0 .and. trim(ic_name) == 'PHOT') then
                !print*, 'read ', ic_name, sname, ind, i
                ! check npoints
                if (i .ne. nsza) then
                    print*, 'photolysis read, number of points not match: ',i,' ', nsza
                    stop
                endif
                
                ! check ind_tmp if in ratio_id
                do s=1, isize
                    if (ind .eq. ratio_id(s)) then ! found. id in ratio_id is s
                        n = n + 1 ! record read
                        do j=1, nsza ! read data points for this species 
                            read(11, *) tmp, a_tmp(j) ! sza,j
                            !print*, tmp,' ', a_tmp(j)
                            ! check sza
                            if (tmp .ne. szas(j)) then
                                print*," sza not match: ",tmp, szas(j),j
                                stop
                            endif
                            ! set for search later
                            ratio_id(s) = -1
                        enddo
                        
                        ! interpolation
                        !call ssh_SPL3(nsza,szas,a_tmp,tmp_inter)
                        call ssh_gck_cspline(nsza,szas,a_tmp,tmp_inter)
                        
                        ! save in photo_ratio
                        do j=1, nsza-1
                          do k=1,4
                            photo_ratio(s+npho_t,j,k) =  tmp_inter(k,j)
                          enddo
                        enddo

                        ! update photo index in photo_rcn
                        tag = 0 ! no.reaction
                        do j=1, npho_tot
                            if (photo_rcn(j,2) .eq. ind) then
                                photo_rcn(j,2) = s + npho_t ! set index in photo_ratio
                                tag = tag + 1
                            endif
                        enddo
                        if (tag.eq.0) then ! check read if complete
                            print*,"Not found index in photo_rcn"
                            stop
                        endif
                        exit
                    endif
                enddo
            end if
          end if  
        end do    
      close(11)
    
      ! after reading, check numbers
      if (n.ne.isize) then
        print*,"Photolysis: not found all points from file. find: ", &
                n," need: ",isize
        do s=1, isize
          if (ratio_id(s) /= -1) print*,"PHOT id: ",ratio_id(s)
        enddo
        stop
      endif
    endif
      
    if (allocated(ratio_id)) deallocate(ratio_id)
    
  end subroutine ssh_read_photolysis_file

!------------------------------------------------------------------------
       
  ! ============================================================
  !
  !  !   ! free allocated memory
  !
  ! ============================================================

  subroutine ssh_free_allocated_memory()

    integer :: ierr = 0

    ! genoa
    if (allocated(RO2index)) deallocate(RO2index)
    if (allocated(RO2out_index)) deallocate(RO2out_index)
    if (allocated(cst_gas_index)) deallocate(cst_gas_index)
    if (allocated(cst_aero_index)) deallocate(cst_aero_index)
    if (allocated(output_gas_index)) deallocate(output_gas_index)
    if (allocated(output_aero_index)) deallocate(output_aero_index)
    if (allocated(output_err_index)) deallocate(output_err_index)
        
    if (allocated(cst_gas)) deallocate(cst_gas)
    if (allocated(cst_aero)) deallocate(cst_aero)
    if (allocated(ref_soa)) deallocate(ref_soa)
    if (allocated(total_soa)) deallocate(total_soa)
    
    if (allocated(output_aero_species)) deallocate(output_aero_species)
    if (allocated(output_gas_species)) deallocate(output_gas_species)
    if (allocated(output_err_sps)) deallocate(output_err_sps)
    if (allocated(ref_conc_files)) deallocate(ref_conc_files)
    
    ! genoa - kinetic
    if (allocated(photo_rcn)) deallocate(photo_rcn)
    if (allocated(TB_rcn)) deallocate(TB_rcn)
    if (allocated(fall_rcn)) deallocate(fall_rcn)
    if (allocated(hetero_rcn)) deallocate(hetero_rcn)
    if (allocated(irdi_rcn)) deallocate(irdi_rcn)    
    if (allocated(extra_rcn)) deallocate(extra_rcn)
    if (allocated(index_RCT)) deallocate(index_RCT)
    if (allocated(index_PDT)) deallocate(index_PDT)
    if (allocated(index_PDT_ratio)) deallocate(index_PDT_ratio)
    if (allocated(Arrhenius)) deallocate(Arrhenius)
    if (allocated(fall_coeff)) deallocate(fall_coeff)
    if (allocated(extra_coeff)) deallocate(extra_coeff)
    if (allocated(hetero_ind)) deallocate(hetero_ind)
    if (allocated(irdi_ind)) deallocate(irdi_ind)
    if (allocated(ratio_PDT)) deallocate(ratio_PDT)
    if (allocated(photo_ratio)) deallocate(photo_ratio)
    if (allocated(photo_ratio_read)) deallocate(photo_ratio_read)
    if (allocated(gas_yield)) deallocate(gas_yield)
    if (allocated(rcn_rate)) deallocate(rcn_rate)
    if (allocated(kinetic_rate)) deallocate(kinetic_rate)
    if (allocated(chem_prod)) deallocate(chem_prod)
    if (allocated(chem_loss)) deallocate(chem_loss)
    if (allocated(drv_knt_rate)) deallocate(drv_knt_rate)
    
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
    if (allocated(partitioning)) deallocate(partitioning, stat=ierr)
    if (allocated(smiles)) deallocate(smiles, stat=ierr)
    if (allocated(aerosol_hydrophilic)) deallocate(aerosol_hydrophilic, stat=ierr)
    if (allocated(aerosol_hydrophobic)) deallocate(aerosol_hydrophobic, stat=ierr) 
    if (allocated(saturation_vapor_pressure))  deallocate(saturation_vapor_pressure, stat=ierr)
    if (allocated(enthalpy_vaporization))  deallocate(enthalpy_vaporization, stat=ierr)    
    if (allocated(mass_density))  deallocate(mass_density, stat=ierr)
    if (allocated(mass_density_layers))  deallocate(mass_density_layers, stat=ierr)
    if (allocated(inon_volatile))  deallocate(inon_volatile, stat=ierr)
    if (allocated(layer_number))  deallocate(layer_number, stat=ierr)
    if (allocated(Vlayer))  deallocate(Vlayer, stat=ierr)
    if (allocated(t_ref))  deallocate(t_ref, stat=ierr)
    if (allocated(henry))  deallocate(henry, stat=ierr)
    !!	if (allocated(saturation_pressure))  deallocate(saturation_pressure, stat=ierr)
    !!	if (allocated(vaporization_enthalpy))  deallocate(vaporization_enthalpy, stat=ierr)
    if (allocated(List_species))  deallocate(List_species, stat=ierr)
    if (allocated(isorropia_species))  deallocate(isorropia_species, stat=ierr)
    if (allocated(isorropia_species_name))  deallocate(isorropia_species_name, stat=ierr)
    if (allocated(aec_species))  deallocate(aec_species, stat=ierr)
    if (allocated(aec_species_name))  deallocate(aec_species_name, stat=ierr)
    if (allocated(aerosol_species_interact))  deallocate(aerosol_species_interact, stat=ierr)
    if (allocated(oligo_index))  deallocate(oligo_index, stat=ierr)
    if (allocated(frac_oligo))  deallocate(frac_oligo, stat=ierr)
    
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
    if (allocated(surface_equilibrium_conc_nsize))  deallocate(surface_equilibrium_conc_nsize, stat=ierr)
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

    ! meteo
    if (allocated(temperature_array)) deallocate(temperature_array)
    if (allocated(pressure_array)) deallocate(pressure_array)
    if (allocated(humidity_array)) deallocate(humidity_array)
    if (allocated(relative_humidity_array)) deallocate(relative_humidity_array)
    
    ! FINAL PART: OUTPUTS
    if (allocated(output_time)) deallocate(output_time)       !(nt+1)                  
    if (allocated(output_gas)) deallocate(output_gas)         !(nt+1,N_gas)             
    if (allocated(output_numb)) deallocate(output_numb)       !(nt+1,N_size+1)         
    if (allocated(output_diam)) deallocate(output_diam)       !(nt+1,N_size)           
    if (allocated(output_aero)) deallocate(output_aero)       !(nt+1,N_aerosol,N_size) 
    if (allocated(output_TM)) deallocate(output_TM)           !(nt+1,N_aerosol,3)        
    if (allocated(output_special)) deallocate(output_special) !(nt+1,8)
    if (allocated(output_pH)) deallocate(output_pH)           !(nt+1,N_size)
    
  END subroutine ssh_free_allocated_memory

  ! =============================================================
  !
  ! Set the flag to decide if SSH-aerosol is logging informations
  !
  ! Important: This subroutine must be called before read_namelist
  !
  ! input : true if logging to a file, false (default) otherwise
  ! =============================================================

  subroutine ssh_set_logger(flag)

    implicit none

    logical, intent(in) :: flag

    integer :: ierr

    ssh_logger = flag

    ! Create log file if needed
    if (ssh_logger) then
       ! Create or overwrite the logflie.
       open(unit = logfile, file = trim(ssh_logger_file), status = "replace", iostat = ierr)
       
       if (ierr.ne.0) then
          write(*,*) "SSH-aerosol: error when creating / opening log file."
          stop
       endif
    endif

  end subroutine ssh_set_logger

  ! =============================================================
  !
  ! Properly close the log file to finalize the logger
  !
  ! =============================================================

  subroutine ssh_close_logger()

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


  end subroutine ssh_close_logger


! ========================================================
! genoa used for read strings from namelist

subroutine get_token(input_string, token, delimiter)

    character(*), intent(inout) :: input_string
    character(*), intent(out) :: token
    character(len=1), intent(in) :: delimiter
    integer :: i ! comma position

    input_string = adjustl(input_string)
    i = index(input_string, delimiter) ! Find delimiter
    if (i == 0) then ! No delimiter
        token = trim(input_string)
        input_string = ""
        return
    end if

    token = trim(input_string(1:i-1)) ! Get token before delimiter
    input_string = input_string(i+1:) ! Update input_string

end subroutine get_token

subroutine get_prodcut_and_ratios(s, numbers, str_val, nrt)

    character(len=*), intent(in) :: s
    double precision, dimension(3) :: numbers
    character(len=len(s)) :: str_val
    integer :: i, ierr, nrt

    numbers = 0.d0
    nrt = 0 ! number of ratios <= 3

    str_val = adjustl(s) ! Adjust the string
    
    i = 1  ! current position of a string
    do while (i <= len_trim(str_val))
        if (str_val(i:i) == ' '.and. i > 1) then ! Find a space (i) and a string (1:i-1)
            nrt = nrt + 1
            if (nrt <= 3) then ! Max. 3 ratios
                read(str_val(1:i-1), *, iostat=ierr) numbers(nrt)
                if (ierr == 0) then ! Read a number
                    str_val = adjustl(str_val(i+1:)) ! Update the string
                    i = 1 ! Reset the loop counter
                else
                    print*, nrt,'number can not read in line: ',trim(str_val),' from ',trim(s)
                    stop
                endif
            else
                print*, 'Too many numbers in line: ', trim(s)
                stop
            endif
        endif
        i = i + 1
    end do

    !print*, 'Finish reading ',trim(s)," find: ", numbers, nrt, trim(str_val)

end subroutine get_prodcut_and_ratios

subroutine count_delimiter(str, delimiter, res)
    character(*), intent(in) :: str
    character(len=*), intent(in) :: delimiter
    integer, intent(out) :: res
    integer :: i, n, m

    res = 1
    n = len_trim(str)
    m = len_trim(delimiter)

    if (n >= m .and. m > 0) then
        do i = 1, n - m + 1
            if (str(i:i+m-1) == delimiter) then
                res = res + 1
            end if
        end do
    end if

end subroutine count_delimiter

subroutine count_products(s, count_pdt, count_ratio, count_set)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: str_val
    integer :: count, count_ratio, count_pdt, count_set, i, start, len_s
    character(len=1) :: ch

    str_val = adjustl(s)
    len_s = len_trim(str_val)
    count = 0       ! number of elements
    count_pdt = 0   ! number of products (default one, incremented with '+')
    count_ratio = 0 ! number of ratios
    count_set = 0   ! number of ratios as a set of three numbers 
    start = 1       ! start of the word

    ! Process each character in the string
    do i = 1, len_s
        ch = str_val(i:i)
        if (ch == ' '.or. i == len_s) then ! Find a space
            if (i == len_s .and. ch /= ' ') count = count + 1 ! End of the list
            if (i > start .or. i == len_s) then ! Find a valid string !! YK
               ! if (i == len_s .or. (start == i-1 .and. str_val(start:start) == '+')) then
               if (i == len_s .or. (start == i-1 .and. &
                    (str_val(start:start) == '+' .or. str_val(start:start) == '-'))) then   
                  
                  ! write(*,*) " Found ---> ", i, start, str_val(start:start) !! YK
                  count_pdt = count_pdt + 1 ! Add one product
                    if (count > 1) then ! count == 1: sname
                        count_ratio = count_ratio + 1
                        if (count == 4) then ! count == 3: ratio, k1, k2, sname
                            count_set = count_set + 1
                        else if (count /= 2) then ! count == 2: ratio, sname
                            print*, count, 'count number should be 1, 2, or 4. Check: ', trim(s)
                            stop
                        endif
                    endif
                    count = 0 ! Reset
                else  ! Find an element: can be a ratio or species name
                    count = count + 1
                endif
            endif
            start = i + 1
        endif
    enddo
end subroutine count_products

! ===========================================
end module aInitialization
