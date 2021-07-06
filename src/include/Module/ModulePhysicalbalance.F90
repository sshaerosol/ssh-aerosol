 !!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------
!!
!!     -- DESCRIPTION
!!    This module contains methods related to particles physical characteristics
!!-----------------------------------------------------------------------
MODULE dPhysicalbalance
  use aInitialization

  implicit none

contains
  subroutine ssh_compute_average_diameter()
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine compute particle dry diameter of each
!     bins
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!------------------------------------------------------------------------   
    implicit none
    integer::k,s,j,jesp
    double precision::mass_total
    double precision::volum_cell,av_volum_cl

    !call compute_all_density()!renew the density of each bin

    do j= 1, N_size
       k=concentration_index(j, 1)!size bins
       cell_diam_av(j)=size_diam_av(k)
       volum_cell=0.d0
       mass_total =0.d0
       do s=1,N_aerosol_layers
	jesp=List_species(s)
        !if(aerosol_species_name(jesp).NE.'PH2O') then
          mass_total = mass_total + concentration_mass(j,s)
          mass_total_grid (j) =mass_total
	  if(mass_density(jesp).gt.0.d0) then
             volum_cell=volum_cell+concentration_mass(j,s)/mass_density(jesp)
	  endif
	 !endif
	enddo
	if (concentration_number(j).gt.0.d0) then
          av_volum_cl=dble(volum_cell)/dble(concentration_number(j))
	  cell_diam_av(j)  = (6.d0*av_volum_cl/pi)**(1.D0/3.D0) ! Âµm
	else
	  concentration_number(j)=0.d0
	  cell_diam_av(j)=size_diam_av(k)
	endif
	if(cell_diam_av(j).eq.0.d0) then
	  cell_diam_av(j)=size_diam_av(k)
	else if(cell_diam_av(j).lt.diam_bound(k)) then
	  cell_diam_av(j)=diam_bound(k)
	endif
    enddo

  end  subroutine ssh_compute_average_diameter

  subroutine ssh_compute_number()
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine compute particle based on particle mass
!     and fixed diameter.
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!------------------------------------------------------------------------  
    implicit none
    integer::k,s,j,jesp
    double precision::volum_cell

    do j= 1, N_size
      volum_cell=0.d0
      k=concentration_index(j, 1)!size bins
      do s= 1, N_aerosol_layers
         if(s.NE.EH2O_layers) then !aerosol_species_name(List_species(s)).NE.'PH2O') then
  	    volum_cell=volum_cell+concentration_mass(j,s)
         endif
      enddo
      
      if(density_aer_size(k).gt.0.d0) then
       volum_cell=volum_cell/density_aer_size(k)
      endif
      if(size_diam_av(k).gt.0.d0) then
         concentration_number(j)= volum_cell /((size_diam_av(k)**3.d0)*cst_PI6)
      else
	  print*,"Wrong size_diam_av",k,size_diam_av(k)
      endif
    enddo
  end  subroutine ssh_compute_number

  subroutine ssh_compute_all_density()
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine update overall density of each bin and size bin
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!------------------------------------------------------------------------  
    implicit none

    integer::k,s,j,f,jesp
    double precision:: subrho(N_size),masstot(N_size)
    double precision:: subrho2(N_sizebin),masstot2(N_sizebin)

    subrho = 0.d0
    masstot = 0.d0
    subrho2 = 0.d0
    masstot2 = 0.d0
    !calculate rho of each grid point
    do k=1, N_size
       do s = 1, N_aerosol_layers
          jesp=List_species(s)
          if(s.NE.EH2O_layers) then
             if(mass_density(jesp).gt.0.d0) then
                subrho(k) = subrho(k) + concentration_mass(k,s)/mass_density(jesp)
             endif
             masstot(k) = masstot(k) + concentration_mass(k,s)
          endif
          if (masstot(k).eq.0.d0 .OR. subrho(k).eq.0.d0) then
             density_aer_bin(k) = fixed_density
          else
             if(subrho(k).gt.0.d0) density_aer_bin(k) = masstot(k)/subrho(k)
          endif
       enddo
    enddo

    !calculate everage rho of each size bins
    do j =1, N_sizebin
      do f = 1, N_fracmax
        k=concentration_index_iv(j,f)
        subrho2(j)=subrho2(j)+ subrho(k)
        masstot2(j) = masstot2(j) +  masstot(k)
      enddo
      if (masstot2(j).eq.0d0 .OR. subrho2(j).eq.0d0) then
        density_aer_size(j) = fixed_density
      else
        if(subrho2(j).gt.0.d0) density_aer_size(j) = masstot2(j)/subrho2(j)
      endif
    enddo

  end  subroutine ssh_compute_all_density

  subroutine ssh_mass_conservation(c_mass,c_number,c_gas, t_mass)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine computes the gas-phase concentration by mass
!     conservation.
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     c_mass: aerosol mass concentration(µg/m^3)
!     c_number: aerosol number concentration(#/m^3)
!
!     -- INPUT/OUTPUT VARIABLES

!     c_gas: aerosol gas phase concentration(µg/m^3)
!     t_mass: total mass of each species both gas phase and aerosol(µg/m^3)
!
!------------------------------------------------------------------------     
    implicit none

    integer::j,jesp,k,s,lay
    double precision:: c_number(N_size)
    double precision:: c_mass(N_size,N_aerosol_layers)
    double precision:: c_gas(N_aerosol)
    double precision:: t_mass(N_aerosol)
    double precision:: tmp_cell,total_ms

    !double precision:: total_number,total_mass_t

    do jesp=1,N_aerosol
       total_aero_mass(jesp)=0.d0
       if (frac_oligo(jesp)>0.) then
          t_mass(jesp)=t_mass(jesp)+t_mass(oligo_index(jesp))
          t_mass(oligo_index(jesp))=0.d0
       endif
    enddo
    !bin_mass=0.d0
    !cell_mass=0.d0
    !bin_number=0.d0
    !total_number=0.d0
    !total_mass_t=0.d0
    !check all not negtive value is allowed   
    
    do j=1,N_size
      tmp_cell=0.d0
      !k=concentration_index(j, 1)
      !total_number=total_number+c_number(j)
      !bin_number(k)=bin_number(k)+c_number(j)
      do s=1,N_aerosol_layers-1
         jesp=List_species(s)
         !if(s.NE.EH2O_layers) then
         !total_mass_t=total_mass_t+c_mass(j,s)
         tmp_cell=tmp_cell+c_mass(j,s)
         !bin_mass(k)=bin_mass(k)+c_mass(j,s)
         !cell_mass(j)=cell_mass(j)+c_mass(j,s)
         !endif
      enddo
      if(tmp_cell.eq.0.d0) then!
	c_number(j)=0.d0
	do jesp=1,N_aerosol_layers
	  c_mass(j,jesp)=0.d0
	enddo
      endif
    enddo

    !if(total_number.eq.0.d0.or.total_mass_t.eq.0.d0) then
    !  do j=1,N_size
	!c_number(j)=0.d0
	!do jesp=1,N_aerosol_layers
	!  c_mass(j,jesp)=0.d0
	!enddo
      !enddo
    !endif

    do j=1,N_size
      do s=1,N_aerosol_layers-1
         jesp = List_species(s)
         !if(s.NE.EH2O_layers) then
         total_aero_mass(jesp)=total_aero_mass(jesp)+c_mass(j,s)
         !endif
      enddo
    enddo

    do jesp=1,(N_aerosol-1)
      if (aerosol_species_interact(jesp).gt.0) then      
	!renew the gas_concentration (reduce into aerosol)
         c_gas(jesp)=t_mass(jesp)-total_aero_mass(jesp)
         if ((c_gas(jesp) .GT. 0.d0) .and. &
              (c_gas(jesp) .LT. TINYM)) c_gas(jesp) = 0.d0
      else
	  c_gas(jesp)=0.d0
      endif
    enddo

    do jesp=1,N_nonorganics
       if(c_gas(jesp).lt.0.d0) then
	  if(total_aero_mass(jesp).gt.0.d0) then
             do j=1,N_size
                   c_mass(j,jesp)=c_mass(j,jesp)+c_gas(jesp)*(c_mass(j,jesp)/total_aero_mass(jesp))
             enddo
          endif
          total_aero_mass(jesp)=t_mass(jesp)
          c_gas(jesp)=0.d0
       endif
       if(total_aero_mass(jesp).lt.0.d0) then
          total_aero_mass(jesp) =0.d0
          c_gas(jesp)=t_mass(jesp)
          do j=1,N_size
             c_mass(j,jesp)=0.d0
          enddo
       endif
    enddo
    do s=N_nonorganics+1, N_aerosol_layers-1
       jesp = List_species(s)
       if(c_gas(jesp).lt.0.d0) then
	  if(total_aero_mass(jesp).gt.0.d0) then
             do j=1,N_size
                   c_mass(j,s)=c_mass(j,s)+c_gas(jesp)*(c_mass(j,s)/total_aero_mass(jesp))
             enddo
          endif
          total_aero_mass(jesp)=t_mass(jesp)
          c_gas(jesp)=0.d0
       endif
       if(total_aero_mass(jesp).lt.0.d0) then
          total_aero_mass(jesp) =0.d0
          c_gas(jesp)=t_mass(jesp)
          do j=1,N_size
             c_mass(j,s)=0.d0
          enddo
       endif
    enddo
    
    do s=N_nonorganics+1, N_aerosol_layers-1
       jesp = List_species(s)
       if (frac_oligo(jesp)>0.) then          
          do j=1,N_size             
             c_mass(j,oligo_index(jesp))=frac_oligo(jesp)*c_mass(j,s)
             c_mass(j,s)=(1.d0-frac_oligo(jesp))*c_mass(j,s)
             t_mass(oligo_index(jesp))=t_mass(oligo_index(jesp))+c_mass(j,oligo_index(jesp))
             t_mass(jesp)=t_mass(jesp)-c_mass(j,oligo_index(jesp))           
          enddo
       endif                     
    enddo

    
    

  end subroutine ssh_mass_conservation

end MODULE dPhysicalbalance
