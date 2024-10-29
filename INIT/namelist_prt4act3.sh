&setup_meteo
latitude = 48.2,                ! Latitude
longitude = 2.22,               ! Longitude
Temperature = 298.0,    ! TEMP
Pressure = 1.01325D05,         ! PRES
Humidity = 5.33039774746E-3,      ! Specific humidity
Relative_Humidity = 0.8	! if 0, compute RH from specific humidity. if not, use it
/


&setup_time
initial_time = 19040400,             
final_time = 19058400.0,             
delta_t = 60,                 
/

&initial_condition
with_init_num = 0,                ! 0 estimated from mass and diameter; 1 number conc. for each bin is read 
tag_init = 0,                     ! initial method for aerosol species (0 internally mixed; 1 mixing_state resolved !! option 1 not yet available)
wet_diam_estimation = 1       	  ! 0 = isorropia
tag_dbd = 1,	                  ! Method for defining particle size bounds (0 if they are auto generated, 1 if bounds are read)
N_sizebin = 1,                    ! Number of size bin
init_gas_conc_file = "./inputs/inputs-prt4act3/init_gas.dat",                     ! Initial data for gas
init_aero_conc_mass_file = "./inputs/inputs-prt4act3/init_aero.dat",        !  Initial data for aero mass  
/

&initial_diam_distribution       
diam_input = 0.1585 0.4
/

&emissions
tag_emis = 0,                ! 0 Without emissions, 1 with internally-mixed emissions, 2 with externally-mixed emissions  !!KS Check when reading the mixing state that if tag_emis = 2 then tag_external = 1. The program needs to stop if tag_emis=2 and tag_external = 0  !!option 2 not yet available
with_emis_num = 0,    ! 0 estimated from mass and diameter; 1 number conc. for each bin is read
emis_gas_file = "./inputs/emis_gas.dat",
emis_aero_mass_file = "./inputs/emis_aero.dat", 
emis_aero_num_file = "./inputs/emis_aero_num.dat" 
/

&mixing_state
tag_external = 0,             !Mixing state(0 for internally mixed, 1 for mixing-state resolved) !!! Initial concentrations are assumed to be internally mixed.
N_groups = 1,                  !Nb of Specie group (internal within group, external between groups)
N_frac = 1, 	                  ! Mass fraction sections for each species
kind_composition = 1,          !Fraction discretization methods 1 for auto discretization and 0 for manual discretization
/

&fraction_distribution
frac_input= 0.0 1.0,          !Manual discretization, set fraction bounds manully
/

&gas_phase_species
species_list_file = "./species-list/species-list-cb05.dat",                    !Species list
reaction_list_file = "./src/include/CHEMISTRY/cb05/CB05.reactions",
/

&aerosol_species
aerosol_species_list_file = "./species-list/species-list-aer.dat",        !Aerosol species list
/

&physic_gas_chemistry
tag_chem = 1,				!Tag of gas-phase chemistry
attenuation = 1.d0,               ! Cloud attenuation field ([fraction]) - 1 = no attenuation
option_photolysis = 1,			!1 if photolysis rates are estimated in the program, 2 read from binary files.
time_update_photolysis = 100000.        ! if photolysis are read, time in seconds between two reads
with_heterogeneous = 0,                 !Tag of heterogeneous reaction 
with_adaptive = 1,                      !Tag of adaptive time step for chemistry 1 if adaptive time step.
adaptive_time_step_tolerance = 0.001,    !Relative tolerance for deciding if the time step is kept
min_adaptive_time_step = 0.001,          !Minimum time step
photolysis_dir = "./photolysis/",  ! Directory where binary files are located.
photolysis_file = "./photolysis/photolysis-cb05.dat",  ! File for photolysis list
n_time_angle = 9,
time_angle_min = 0.d0,
delta_time_angle = 1.d0,
n_latitude = 10,
latitude_min = 0.d0,
delta_latitude = 10.d0,
n_altitude = 9,
altitude_photolysis_input = 0.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 10000.0, 15000.0, 20000.0,
/


&physic_particle_numerical_issues
DTAEROMIN = 1E-5,              	!Minimum time step
redistribution_method = 0,   	!tag of Redistibution methods /3 euler_mass 4 euler_number 5 hemen 6 euler_coupled, 10 Moving Diameter, 11 SIREAM, 12 euler_couple_siream,  0 without redistribution
with_fixed_density = 1,       	!1 if density is fixed - 0 else
fixed_density = 1.4D3,        
splitting = 0,             ! 0 if coagulation and (condensation/evaporation+nucleation) are splitted - 1 if they are not
/

&physic_coagulation
with_coag = 0,                 		!Tag of coagulation
i_compute_repart = 1,                   ! 0 if repartition coeff are not computed but read from file, 1 if they are 
i_write_repart = 0,                     ! 1 to write repartition coeff file, 0 otherwise
Coefficient_file = "coef_s1_f1_b6.nc",  !Coagulation coefficient file coef_s5_f3_b7.nc for external mixing, coef_s1_f1_b7.nc for internal mixing
Nmc = 10000                             ! Number of Monte Carlo points if repartition coeff are computed
/

&physic_condensation
with_cond = 1,               !Tag of condensation/evaporation
Cut_dim = 0.0,             ! KCUT	ICUT = 0
ISOAPDYN = 1,                !0 = equilibrium, 1 = dynamic
nlayer = 1,
with_kelvin_effect = 1,                 ! 1 if kelvin effect is taken into account.
tequilibrium = 1.e-15,                     ! time under which equilibrium is assumed.
dorg = 1.d-12,                          ! diffusion coefficient in the organic phase.
coupled_phases = 1,                       ! 1 if aqueous and organic phases are coupled
activity_model = 1,               ! 1. "ideal", 2. "unifac", 3. "aiomfac"
epser = 0.01,                     ! relative error precision for time step adjustment
epser_soap = 0.01,                 ! relative difference of ros2 in SOAP  
SOAPlog = 1,                  !0=no text output (Henry + smiles decomposition) from SOAP, 1=on screen, 2=in written files
reaction_soap_file="species-list/REACTIONS.particle.soap",
/

&physic_nucleation
with_nucl = 0,                !Tag of nucleation - Need to have lower diameter about 1nm.
/

&output
output_directory = "results/aerosol-prt4act3/ref/", ! modify result saving paths : basic/ cond/ coag/ icut0/ ISOAPDYN/
output_type = 1               ! 1: text, 2: binary
particles_composition_file = "fractions.txt"
/
