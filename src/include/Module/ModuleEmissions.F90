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
  use dPhysicalbalance

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
    ! record emission time
    double precision :: current_time,time_step, emis_dt, elapsed_time 
    double precision :: qemis
    double precision :: subrho,rho_emis 
    !double precision, dimension() :: density_aer_emissions

    elapsed_time = current_time - initial_time + time_step

    if(elapsed_time.lt.time_emis) then
      emis_dt=time_step
    else
      emis_dt=time_emis-(elapsed_time - time_step)
    end if

      if (emis_dt .le. 0.d0) return 

	print*, 'called emission()	emis_dt', emis_dt
    ! For gas-phase emissions.
     do i = 1, N_gas
       if (gas_emis(i) .gt. 0.d0) then 
             qemis = gas_emis(i) * emis_dt
             concentration_gas_all(i) = concentration_gas_all(i) + qemis
        end if
     end do 
       

    ! For aerosols emissions.

    ! some precision will be lost 
    ! due to the different maganitude between concentration_mass and qemis


    !!!! emis_bin_number() -> emission_num_rate(N_size, N_aerosol)
    !!!! emis_bin_massn() -> emission_rate(N_size, N_aerosol)

    do j = 1,N_size
        m_emis = 0.d0
        subrho = 0.d0
        k = concentration_index(j, 1) ! number for sizebin
         do s = 1, N_aerosol  ! 
	         qemis = emission_rate(j,s)
                 ! add conc. mass
	         concentration_mass(j,s)=concentration_mass(j,s)+qemis*emis_dt
		 ! total_aero_mass(s) = total_aero_mass(s) + qemis*emis_dt

                 m_emis = m_emis + qemis ! total emission mass
                 subrho = subrho +  emission_rate(j,s) / mass_density(s)   ! emission density

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



	! re-calculate total_mass(N_aerosol) because mass change due to emission   
    total_aero_mass = 0.d0
    total_mass = 0.d0
    do s = 1, N_aerosol
       do j=1,N_size
         total_aero_mass(s) = total_aero_mass(s) + concentration_mass(j,s)
       enddo
	! update mass conc. of aerosol precursors
	! concentration_gas_all(precursor_index) -> concentration_gas(n_aerosol)
       if (aerosol_species_interact(s) .gt. 0) then
          concentration_gas(s) = concentration_gas_all(aerosol_species_interact(s))
       end if
          total_mass(s) = total_mass(s) + concentration_gas(s) + total_aero_mass(s)
    end do

  
   !do i = 1, naero_emis
   !   jesp = 0
	!do j =1, N_aerosol
	 !  if (jesp == 1) exit 
	  ! if (aerosol_species_name(j) == emis_aer_species_name(i)) then
	!	jesp = 1
	!	do k = 1, n_sizebin !!! NEED TO CHECK WHEN MIXING STATE RESOLVED
         !          qemis = init_bin_emission(k,i) * emis_dt 	    
	!	   concentration_mass(k,j)=concentration_mass(k,j)+qemis
	!	end do 
         !  end if
        !end do 
    !end do

    ! init_bin_num_emis(k)



  !  if (with_emis_num == 0) then

  !    do s= 1, N_size
   !      volum_cell=0.d0
    !     k=concentration_index(s, 1)!size bins
    !     do j= 1, (N_aerosol-1)
    !       do i = 1, naero_emis
  	!     if (aerosol_species_name(j) == emis_aer_species_name(i)) then
	!       jesp = emis_aer_species_name(i)
	!       volum_cell = volum_cell +  init_bin_emission(s,jesp)
     !          if(mass_density(j).gt.0.d0) then
     !             subrho = subrho +  init_bin_emission(s,jesp)/mass_density(j)
      !         endif
     !          if (volum_cell.eq.0.d0 .OR. subrho.eq.0.d0) then
     !             density_aer_emissions = fixed_density
      !         else
     !             if(subrho.gt.0.d0) density_aer_emissions = volume_cell/subrho
     !          endif
     !       endif
     !      enddo
    !     enddo
    !    if(density_aer_emissions.gt.0.d0) then
   !        volum_cell=volum_cell/density_aer_emissions
    !    else
   !        write(*,*) 'Pb in computing density at emissions',density_aer_emissions
  !         STOP
  !      endif

 !       if(size_diam_av(k).gt.0.d0) then
!	  init_bin_num_emis(s)=(6.d0*volum_cell)/((size_diam_av(k)**3.d0)*pi)
  !      else
!	  print*,"Wrong size_diam_av",k,size_diam_av(k)
    !    endif
   !   enddo



  end subroutine emission


end MODULE mEmissions
