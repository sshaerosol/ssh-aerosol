!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------
!!
!!     -- DESCRIPTION
!!    This module contains methods for for mass and number redistribution
!!    for fixed euler scheme in 3D application.
!!-----------------------------------------------------------------------
Module eRedistribution
  use aInitialization
  use dPhysicalbalance
  use cThermodynamics
  implicit none

contains
  subroutine ssh_redistribution_size(scheme)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine provides entries for different size redistribution methods
!     the chosen numerical solver.
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     scheme: type of chosen redistribution method
!
!------------------------------------------------------------------------
    implicit none
    integer::k,j,f,jesp,s,b,lay
    integer:: scheme !redistribution scheme 3 = euler_mass 4 = euler_number 5 = hemen 6 = euler_coupled
    integer:: section_pass!bin include 100nm
    double precision:: Qesp(N_sizebin, N_aerosol_layers)!Temperature mass concentration on each fraction
    double precision:: N(N_sizebin)!Temperature mass concentration on each fraction
    double precision:: totQ(N_sizebin)
    double precision:: d(N_sizebin),d_after_cond(N_sizebin)
    double precision:: tmp_n,mass_redis_layer(N_aerosol_layers)
    Qesp=0.d0
    tmp_n=0.d0

    if(scheme.EQ.5) then
      do k=1,N_sizebin
        if(diam_bound(k).lt.1.d-2.and.diam_bound(k+1).ge.1.d-2) then
         section_pass=k
        endif
      enddo
    endif
    !see the distribution in one fixed fraction as the same case of internal mixing
    do f=1,N_fracmax
      !extrac the mass and number distribution of each fraction out
      do k=1,N_sizebin
	do s=1,N_aerosol_layers
	  Qesp(k,s)=0.d0
        enddo 
      enddo
      do k=1,N_sizebin
	j=concentration_index_iv(k,f)
	N(k)=concentration_number(j)
	!totQ(k)=cell_mass(j)
	do s=1,N_aerosol_layers
	  Qesp(k,s)=Qesp(k,s)+concentration_mass(j,s)
	enddo
	d(k)=size_diam_av(k)
	d_after_cond(k)=cell_diam_av(j)
      enddo

      call ssh_redistribution(N_sizebin,N_aerosol_layers,EH2O_layers,diam_bound, d, scheme, &
      section_pass, mass_density_layers, DQLIMIT, Qesp, N, totQ,&
      with_fixed_density, fixed_density,d_after_cond)


      !update distribution of each bin in current fraction section size_diam_av
      do k=1,N_sizebin
	j=concentration_index_iv(k,f)
	concentration_number(j)=N(k)
	cell_mass(j)=totQ(k)
	cell_diam_av(j)=d_after_cond(k)
	do s=1,N_aerosol_layers
	  concentration_mass(j,s)=Qesp(k,s)
	enddo
      enddo
    enddo

	  
  end subroutine ssh_redistribution_size

  subroutine ssh_redistribution_fraction()
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine redistribute mass and number
!     based on the new fraction composition of aerosol
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!------------------------------------------------------------------------
    implicit none
    integer::k,i,i1,g,dj,j,jesp,s
    double precision:: mass_groups(N_size,N_groups)!mass by groups
    double precision:: mass_var(N_size,(N_aerosol_layers-1))!mass variation map
    double precision:: numb_var(N_size)!number variation map
    double precision::mass_total,f1,f2
    !calculate redistribution map
!! Initialize to zero the variation in mass and number concentrations due to 
!! redistribution

    do j= 1, N_size  
       numb_var(j)=0.d0
       do jesp=1,(N_aerosol_layers-1)
       	  mass_var(j,jesp)=0.d0
       enddo
       do g= 1, N_groups
	  frac_grid(j,g)=0.d0!Fraction of group g in grid j
	  mass_groups(j,g)=0.d0
       enddo
    enddo

    do j= 1, N_size !!for each grid the number of size sections x fraction sections
      mass_total_grid (j)=0.d0
!      Compute fraction of each species in the grid point before redistribution      
      do s=1,N_aerosol_layers-1
         jesp= List_species(s)
         if(s.NE.EH2O_layers) then
            g=index_groups(jesp)
            mass_groups(j,g)=mass_groups(j,g)+concentration_mass(j,s)
            mass_total_grid (j)=mass_total_grid (j) + concentration_mass(j,s)
         endif
      enddo

      do g= 1, N_groups
	if(mass_total_grid (j).gt.0d0) then
	  frac_grid(j,g) =mass_groups(j,g)/mass_total_grid (j)
	endif
      enddo
      
      i1 =concentration_index(j, 1)!size bin index

      ! Loop on fraction combinations in the size bin i1 
      ! to find out where the 
      ! original fraction section (before condensation/evaporation) has moved to.
      do i=1,N_fracmax
         dj=0
         ! Check for each species whether 
         ! this fraction bin is the correct one
         do g=1,N_groups-1
            f1=discretization_composition( i, g, 1)
            f2=discretization_composition( i, g, 2)
            if(f1.eq.0.d0) then
               if(frac_grid(j,g).ge.f1.and.&
                    frac_grid(j,g).le.f2) then
                  dj=dj+1
               endif
            else
               if(frac_grid(j,g).gt.f1.and.&
                    frac_grid(j,g).le.f2) then
                  dj=dj+1!dj is the number of matched groups
               endif
            endif
         enddo
	
         !  The fraction of each species is identified. 
         ! i is  then the correct fraction combination (after redistribution). 
         if(dj.eq.(N_groups-1)) then
            k=concentration_index_iv(i1,i)
            numb_var(k)=numb_var(k)+concentration_number(j)
            concentration_number(j)=0.d0!minimize numerical disarrange due to difference scale
            do s=1,N_aerosol_layers-1
               jesp=List_species(s)
               if(aerosol_species_name(jesp).NE.'PH2O') then
                  mass_var(k,s)=mass_var(k,s)+ &
                       concentration_mass(j,s)
                  concentration_mass(j,s)=0.d0
               endif
            enddo
         endif
      enddo

    enddo
    
    !Update mass and number concentrations after redistribution
    do j= 1, N_size !!for each grid point
       concentration_number(j)=concentration_number(j)+numb_var(j)
       mass_total_grid (j)=0.d0
       do g= 1, N_groups
          frac_grid(j,g)=0.d0!Fraction of group g in grid j
          mass_groups(j,g)=0.d0
       enddo
       !renew Temperature frac of each grid
       do s= 1, (N_aerosol_layers-1)
          jesp=List_species(s)
          concentration_mass(j,s)=concentration_mass(j,s)+mass_var(j,s)
          if(jesp.ne.EH2O) then
             g=index_groups(jesp)
             mass_total_grid (j) = mass_total_grid (j) + concentration_mass(j,s)
             mass_groups(j,g)=mass_groups(j,g)+concentration_mass(j,s)
          endif
       enddo

       do g= 1, N_groups
          if(mass_total_grid (j).gt.0d0) then
             frac_grid(j,g) =mass_groups(j,g)/mass_total_grid (j)
          endif
       enddo
    enddo
    
    
  end subroutine ssh_redistribution_fraction

  subroutine ssh_redistribution_lwc(lwc,ionic,proton,liquid,iredist,end_bin)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine redistribute liquid water content (LWC)
!     based on the fraction of inorganic aerosols.
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!    lwc: liquid water content (ug/m3)
!
!------------------------------------------------------------------------
    implicit none

    integer :: iredist
    double precision :: inorg_total, inorg_bin(N_size)
    integer :: jesp, js,lay,end_bin
    double precision :: lwc,ionic, proton,liquid(12)
!    double precision :: lwc_Nsize(N_size),proton_Nsize(N_size)
!    double precision :: proton_Nsize(N_size)

!    double precision :: ionic_Nsize(N_size),liquid_Nsize(12,N_size)

    inorg_total = 0.D0
    inorg_bin = 0.D0
 
    do jesp=ENa,ECl
      do js=1,N_size
        if(concentration_index(js, 1) <= end_bin) then
          inorg_total = inorg_total + concentration_mass(js, jesp)
          inorg_bin(js) = inorg_bin(js) + concentration_mass(js, jesp)
        endif
      enddo
    enddo

    do js=1,N_size
     if(concentration_index(js, 1) <= end_bin) then
       if (inorg_total .gt. 0.D0) then
          if(iredist.EQ.0) concentration_mass(js, EH2O_layers) = lwc * inorg_bin(js)/inorg_total
          lwc_Nsize(js) = lwc * inorg_bin(js) / inorg_total
          proton_Nsize(js) = proton * inorg_bin(js) / inorg_total
          ionic_Nsize(js) = ionic 
          do jesp=1,12
            liquid_Nsize(jesp,js) = liquid(jesp)
          enddo
       else
          if(iredist.EQ.0) concentration_mass(js, EH2O_layers) = 0.d0
       end if
     endif
    enddo
  end subroutine ssh_redistribution_lwc
  
  subroutine ssh_redistribution_lwcorg(lwcorg,lwcorg_Nsize)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine redistribute liquid water content (LWCORG)
!     based on the fraction of organic aerosols.
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!    lwcorg: liquid water content (ug/m3)
!
!------------------------------------------------------------------------
    implicit none

    double precision :: org_total, org_bin(N_size)
    integer :: jesp, js, jesp2
    double precision :: lwcorg,lwcorg_Nsize(N_size)

    org_total = 0.D0
    org_bin = 0.D0

    do jesp=1,nesp_aec
       jesp2 = aec_species(jesp)
       do js=1,N_size
          org_total = org_total + concentration_mass(js, jesp2)
          org_bin(js) = org_bin(js) + concentration_mass(js, jesp2)
       enddo
    enddo

    do js=1,N_size
       if (org_total .gt. 0.D0) then
          concentration_mass(js, EH2O_layers) = concentration_mass(js, EH2O_layers) + &
                                      lwcorg * org_bin(js) / org_total
          lwcorg_Nsize(js) = lwcorg * org_bin(js) / org_total
       end if
    enddo

  end subroutine ssh_redistribution_lwcorg
  
end Module eRedistribution
