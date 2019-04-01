!!-----------------------------------------------------------------------
!!     Copyright (C) 2012-2018, ENPC - EDF R&D - INERIS
!!     Author(s): Shupeng Zhu, Karine Sartelet
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
!!    This module contains methods related to emissions.
!!-----------------------------------------------------------------------
MODULE mEmissions
  use aInitialization

  implicit none

contains
 
  subroutine emission(current_time,time_step)
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
    double precision::qemis,current_time,time_step
    double precision::time_emis,emis_dt
    double precision :: elapsed_time

    ! if(nucl_model.eq.5) then
    ! time_emis=60.d0*60.d0*12.d0
    ! else
    time_emis=2.64376d3
    ! endif
    !print*,'emis'
    ! if((current_time+time_step).lt.time_emis) then
    !   emis_dt=time_step
    ! else
    !   emis_dt=time_emis-current_time
    !   if(emis_dt.lt.0.d0) emis_dt=0.d0
    ! endif

    elapsed_time = current_time - initial_time + time_step

    if(elapsed_time.lt.time_emis) then
      emis_dt=time_step
    else
      emis_dt=time_emis-(elapsed_time - time_step)
    endif

    ! For gas-phase emissions.
    do i = 1, ngas_emis
       do j = 1, n_gas
          if (species_name(j) == emis_species_name(i)) then
             qemis = gas_emis(i) * emis_dt
             concentration_gas_all(j) = concentration_gas_all(j) + qemis
          end if
       end do
    end do
       

      !for non C/E species
    do j=1,N_size
      k=concentration_index(j, 1)
      do jesp=EMD,EBC
	if(emission_rate(j,jesp).gt.0.d0) then
	  if(concentration_mass(j,jesp).ne.concentration_mass(j,jesp)) then
	    print*,j,jesp,concentration_mass(j,jesp),'emission 1'
	  endif
	    qemis=emission_rate(j,jesp)*emis_dt!some precision will be lost due to the different maganitude between concentration_mass and qemis
	    !print*,'current_sub_time',current_time,time_step
	  concentration_mass(j,jesp)=concentration_mass(j,jesp)+qemis
	  concentration_number(j)=concentration_number(j)+qemis/size_mass_av(k)!cell_mass_av(j)
	  n_emis=n_emis+qemis/size_mass_av(k)
	  m_emis=m_emis+qemis
	endif
      enddo
    enddo

    do jesp=EMD,EBC
      total_aero_mass(jesp)=0.d0
    enddo

    do j=1,N_size
      do jesp=EMD,EBC
	if(concentration_mass(j,jesp).ne.concentration_mass(j,jesp)) then
	    print*,j,jesp,concentration_mass(j,jesp),'emission 2'
	    stop
	  endif
	total_aero_mass(jesp)=total_aero_mass(jesp)+concentration_mass(j,jesp)
	total_mass(jesp)=total_aero_mass(jesp)
      enddo
    enddo

    !     !for C/E species
    ! do jesp=ENA,N_aerosol
    !   concentration_gas(jesp)=concentration_gas(jesp)+gas_emission_rate(jesp)*emis_dt
    !   total_mass(jesp)=total_mass(jesp)+gas_emission_rate(jesp)*emis_dt
    ! enddo
  end subroutine emission


end MODULE mEmissions
