o---------------------o
| SSH-AEROSOL HISTORY |
o---------------------o

# SSH-aerosol version 2.0 release (2024-10-XX)

## Model

- Added a wall loss module
- Added a parameterization of Kelvin effect
- Computed diffusion coefficient based on the viscosity
- Added option to take into account multiple RO2 pool for GENOA
- Added MCM and GECKO kinetic rates
- New reactions between organic and inorganic ions
- Modified ELVOC formation from monoterpenes
- Updated SMILES to check the decomposition of a species into functional groups
- Added new biogenic organic aerosol species


## Test case

- Added test cases for the oxidation of toluene, monoterpene, naphthalene, sesqui-terpene
- Added test cases for AIOMFAC-VISC
- Added GENOA and MCM test cases
- Added a test case for IEPOX with particle-phase reaction
- More size sections for test cases
- Added a test case for HNO3 condensation onto CaCO3
- Added a test case for soap_inorg = 1

## Fixed bugs

- Corrected redistribution option 12
- Corrected mass of nucleated particles for heteromolecular nucleation

## Interface

- SPACK is not required to generate chemistry files.
- Added netCDF-style output

## Architecture

- Installation with Ducker
- O2 optimization used instead of Ofast

## Papers



# SSH-aerosol version 1.3 release (2022-06-27)

## Model

- Added the possibility to consider organic compounds as both hydrophilic and hydrophobic.
- Added the possibility to use MCM chemical mechanism, with beta-carophylene as an example.
- Added Melchior2 chemical mechanism.
- Added RACM2-2020 chemical mechanism (Lannuque et al. 2021).
- Added heteromolecular nucleation scheme for biogenic organics and sulfuric acid and the possibility to take into account different types of nucleation simultaneously.
- Improved convergence in SOAP.
- Added an alternative solver (used when the number of iterations exceeds 50) in SOAP in the coupled mode used ideal conditions.
- Corrected smiles for PAN as C(=O)OON(=O)(=O)
- Added parameters niter_eqconc and niter_water to call ISOROPIA every iteration and gain computation time.

## Papers illustrating the new functionalities:
- Lannuque V., D'Anna B., Couvidat F., Valorso R., Sartelet K. (2021), Improvement in modeling of OH and HO2 radical concentrations during toluene and xylene oxidation with RACM2 using MCM/GECKO-A. Atmosphere, 12, 732, doi:10.3390/atmos12060732.
- Sartelet, K., Kim, Y., Couvidat, F., Merkel, M., Petäjä, T., Sciare, J., and Wiedensohler, A.: Influence of emission size distribution and nucleation on number concentrations over Greater Paris, Atmos. Chem. Phys. Discuss. [preprint], https://doi.org/10.5194/acp-2022-22, accepted, 2022.
- Wang, Z., Couvidat, F., and Sartelet, K.: GENerator of reduced Organic Aerosol mechanism (GENOA v1.0): An automatic generation tool of semi-explicit mechanisms, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2022-245, 2022.


## Interface 

- Added partitioning variable in species-list-aer.
- Added files for CHIMERE species
- Added the possibility to have as many inert species as required by the user.
- Added option with_fixed_density=2. Density is then nonlinear and depends on the concentrations of inorganics.
- Added option to read aerosol species structures in UNIFAC functional group format.
- Added option to run simulations when applying fixed concentration profiles per time step to certain species in gas and aerosol phases.
- Added script to refine the discretization of the size distribution, while conserving mass and number. 

## Architecture

- Updated compiling by scons for Intel compiler.
- Updated Makefile to include external libraries


## Fixed bugs

- Correct a bug when organics can be both hydrophilic and hydrophobic and at equilibrium.
- Force the gain term of coagulation to zero in case the number of particles is zero.
- Unusual increase in the mass concentration when the particles less than the first size bin are taken into account.
- Wrong initialization for inorganic aerosols using the dynamic approach in SOAP.
- Modify the redistribution in the lowest bin for the redistribution 12.


## Test case

- New test case for mass and number conservative upsampling scheme.
- New test case to use MCM chemical mechanism.



# SSH-aerosol version 1.2.2 release (2021-07-02)

## Model

- Improved the reproductability of the generation of random numbers
  when calculating the coagulation partition coefficient.
- For particles with several organic layers, condense low-volatility compounds
  only on the layer which is the most outside.

## Interface

- Added a detailed example for using the shared library
  for coupling SSH-aerosol to other models in FORTRAN.
- Names of gas-phase species are now available in the API.
- The option to use the two-step solver for gaseous chemistry is
  now available in the API.

## Architecture

- scons compile script is now compatible with python 3.x

## Fixed bugs

- Added initialization of wet diameter which was missing in some cases.
- Improved checking on division by tiny number or zero.
- Improved checking on mass conservation.

## Test case

- In the viscosity test case, special treatment is not applied for
  low-volatility compounds.




# SSH-aerosol version 1.2.1 release (2021-03-01)


## Model

- Add the two-step time numerical solver for gas-phase chemistry
- Correct the reallocation of the vectors in unifac.cxx

## Interface 

- Add a python example for the shared library
- Added new functions in API including get_tag_external

## Fixed bugs

- Correct bug on computation of coagulation
- Correct bug on initialization when aqueous and organic phases are coupled
- Fixed missing initialization of some variables
- Correct bug for water content of sections lower than ICUT


# SSH-aerosol version 1.2 release (2020-11-30)

## Model

- Added generic organic aerosol species. A generic species can be used when a particular species, which needs to be modeled, is not found in SSH-aerosol. Using some provided information about the species (smiles code, saturation vapor pressure, enthalpy of vaporization, molar mass, etc...), SSH-aerosol will automatically create a new generic species during the computation. The generic species can be used to easily add a new organic species in SSH-aerosol. The smiles code will be used to automatically estimate the decomposition into UNIFAC functionnal groups. 
- ICUT (to distinguish bulk equilibrium and dynamic approach) differentiated form inorganics and organics (see the variable ICUT_org).
- Allow coagulation if aerosols are made of different layers.
- Added the peroxyacetyl acid, PAN, peroxide groups for the computation of activity coefficients in SOAP
- For organic compounds, added the possibility to store the aqueous and the organic concentrations separately (see  the variable i_hydrophilic).
- Multiple changes regarding the computation of ICUT using the c/e time scale, the ETR solver efficiency, and the weighted QSSA methods.
- Oligomerization of organics species can be taken into account (with a second-order reversible reaction unfavored by humid conditions) based on [Couvidat et al. (2018)](https://acp.copernicus.org/articles/18/15743/2018/).

## Test case

- New test case for oligomerization.
- Improved CPU time for cases where dynamic condensation and coagulation are considered.

## Interface 

- Modified ModuleAPI for the implementation to Polair3d.
- Modified and corrected fortran Interface for the implementation to any fortran programs.
- Improved CPU time for coagulation.
- Added the possibility to gain CPU time by computing surface equilibrium concentrations only once per time (for second order solver).
- Added input file for meteorological data.

## Fixed bugs

- Fixed condensation rate for hydrophobic species.
- Corrected bug for computation of wet diameter for mixing-state resolved particles with a cutoff diameter between dynamic and equilibrium sections.


* Version 1.1 release (2020-05-03)

  + Model
    - Improved redistribution methods.
    - Improved the computation of wet diameter, especially for cases of
      coagulation only.
    - Added the possibility to condense low-volatility compounds dynamically,
      even if full equilibrium is assumed for semi-volatile compounds.
      That improves the computation of the growth of ultrafine particles.
    - Improved mass conservation for evaporation.

  + Test case
    - New test case with cutoff diameter between dynamic and equilibrium size
      sections.
    - Nucleation with condensation of a low-volatility organic.
    - Exhaust test case with ammonia emission, with or without isoprene forming SOA.
      
  + Interface
    - Added implementation examples to use shared library and Fortran 90 module.
    - Makefile added to compile with make.
    - Reduced namelist.ssh to hide some options and input data which should
      not be modified.

* Version 1.0 release (2019-10-09)

