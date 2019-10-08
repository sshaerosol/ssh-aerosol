!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------
!!
!!     -- DESCRIPTION
!!    This module contains methods related to emissions.
!!-----------------------------------------------------------------------
MODULE mEmissions
  use aInitialization
  use dPhysicalbalance

  implicit none

contains

  subroutine emission(emis_dt)
    !------------------------------------------------------------------------
    !
    !     -- DESCRIPTION
    !     This subroutine adds the emission into aerosol and gas phase
    !
    !------------------------------------------------------------------------
    !
    !     -- INPUT VARIABLES
    !
    !     current_time: current time when emission start(s)
    !     time_step: emission time step(s)
    !
    !------------------------------------------------------------------------        
    implicit none

    integer::s,i,j,jesp,k
    ! record emission time
    double precision :: emis_dt 
    double precision :: qemis
    double precision :: subrho,rho_emis 
    !double precision, dimension() :: density_aer_emissions

    ! For gas-phase emissions.
    do s = 1, N_aerosol
       if (aerosol_species_interact(s) .gt. 0) then
          if (gas_emis(aerosol_species_interact(s)) .gt. 0.d0) then 
             qemis = gas_emis(aerosol_species_interact(s)) * emis_dt
             total_mass(s) = total_mass(s) + qemis
             concentration_gas(s) = concentration_gas(s) + qemis
          end if
       end if
    end do
    do i = 1, N_gas
       if (gas_emis(i) .gt. 0.d0) then 
          qemis = gas_emis(i) * emis_dt
          concentration_gas_all(i) = concentration_gas_all(i) + qemis
       end if
    end do


    ! For aerosols emissions.

    do j = 1,N_size
       m_emis = 0.d0
       subrho = 0.d0
       k = concentration_index(j, 1) ! number for sizebin
       do s = 1, N_aerosol_layers
          jesp = List_species(s)   
          qemis = emission_rate(j,s)
          ! add conc. mass
          concentration_mass(k,s)=concentration_mass(k,s)+qemis*emis_dt
          total_mass(jesp) = total_mass(jesp) + qemis*emis_dt

          m_emis = m_emis + qemis ! total emission mass
          subrho = subrho +  emission_rate(j,s) / mass_density(jesp)   ! emission density

       end do

       if(m_emis.gt.0.d0) then
          if (with_emis_num == 0 ) then

             if(subrho.gt.0.d0) then
                rho_emis = m_emis / subrho
             else  
                rho_emis = fixed_density
             endif

             if(size_diam_av(k).gt.0.d0) then
                emission_num_rate(j)=(6.d0*m_emis)/((size_diam_av(k)**3.d0)*pi*rho_emis)
             else
                print*,"Emission number calculation : wrong size_diam_av",k,size_diam_av(k)
             endif
          end if
       end if
       if (emission_num_rate(j) .gt. 0.d0) then
          concentration_number(j) = concentration_number(j) + &
               emission_num_rate(j)* emis_dt 
       endif

    end do

    !	! re-calculate total_mass(N_aerosol) because mass change due to emission   
    total_aero_mass = 0.d0
    do s = 1, N_aerosol_layers
       jesp = List_species(s)   
       do j=1,N_size
          total_aero_mass(jesp) = total_aero_mass(jesp) + concentration_mass(j,s)
       enddo
    end do

end subroutine emission


end MODULE mEmissions
