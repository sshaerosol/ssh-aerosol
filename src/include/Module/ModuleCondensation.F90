!!-----------------------------------------------------------------------
!!     Copyright (C) 2012-2018, ENPC - EDF R&D - INERIS
!!     Author(s): Karine Sartelet,  Shupeng Zhu
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
!!    This module contains methods to compute particle condensation rate
!!-----------------------------------------------------------------------
 Module fCondensation
  use aInitialization
  use cThermodynamics
  implicit none
 
contains
  subroutine SULFDYN(Q1,Q,N1,N,c_gas,dqdt,dtx,time_step_sulf)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine solve sulferic acid condensation by explicit solution
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     c_gas: aerosol gas phase concentration(µg/m^3)
!     dtx: condensation time step (s)
!
!     -- INPUT/OUTPUT VARIABLES
!
!     Q: aerosol mass concentration of second order evaluation(µg/m^3)
!     N: aerosol number concentration of second order evaluation(#/m^3)
!
!     -- OUTPUT VARIABLES
!
!     Q1: aerosol mass concentration of first order evaluation(µg/m^3)
!     N1: aerosol number concentration of first order evaluation(#/m^3)
!     dqdt: mass derivation(µg/m^3/s)
!
!------------------------------------------------------------------------
    implicit none
    integer::j,jesp
    double precision :: dqdt(N_size,N_aerosol)
    double precision ::c_gas(N_aerosol)!micg/m^-3
    double precision :: ce_kernal_coef_tot ! c/e kernel coef (m3.s-1)
    double precision :: Q1(N_size,(N_aerosol)) ! Mass concentration
    double precision :: Q(N_size,(N_aerosol)) ! Mass concentration  !!BUG
    double precision :: N1(N_size) ! Number concentration
    double precision :: N(N_size) ! Number concentration
    double precision :: dtx,tmp,cond_so4!Time steps
    double precision :: dexploc,n2err,tmp_n2err,time_step_sulf

    jesp=ESO4!Pointer
    ce_kernal_coef_tot = 0.0d0
    cond_so4 = 0.0d0
    n2err = 0.d0

    do j = 1,N_size! Reassigned distribution by mass of each species
      call compute_condensation_transfer_rate(diffusion_coef(jesp), &
      quadratic_speed(jesp), accomodation_coefficient(jesp), &
      wet_diameter(j), dqdt(j,jesp))
      ce_kernal_coef_tot=ce_kernal_coef_tot+N(j)*dqdt(j,jesp)
    enddo

    do j = 1,N_size! Reassigned distribution by mass of each species
      if(ce_kernal_coef_tot.ne.0.d0) then
        dexploc = DEXP(-ce_kernal_coef_tot*dtx)
	tmp=(dqdt(j,jesp)*N(j)/ce_kernal_coef_tot)*&
	      (1.d0-dexploc) * c_gas(jesp)
        if (Q(j,jesp).GT.TINYM) then
          tmp_n2err = (dqdt(j,jesp)*N(j)*dtx*c_gas(jesp))/Q(j,jesp)
          n2err = n2err + tmp_n2err*tmp_n2err
        endif
        Q(j,jesp) = Q(j,jesp)+tmp!renew mass
        Q1(j,jesp) = Q(j,jesp)
        cond_so4 = cond_so4+tmp
        N1(j)=N(j)
        if(dtx.gt.0.d0) dqdt(j,jesp)=tmp/dtx!for redistribution
      endif
    enddo
    c_gas(jesp) = DMAX1(c_gas(jesp)*dexploc,0.d0)

    n2err=DSQRT(n2err)
    if (n2err > 0.0) then
      tmp = 1.D6
      n2err=DMIN1(n2err,EPSER*tmp)
      n2err=DMAX1(EPSER/tmp,n2err)
      ! New time step new time step
      time_step_sulf = dtx*DSQRT(EPSER/n2err)
    else
      time_step_sulf = 0.0
    endif
  end subroutine SULFDYN

  subroutine KERCOND(qn,q,c_gas,Wet_diam,Temp,ce_kernel,ce_kernal_coef_i,jj,&
                     lwc,ionic,proton,liquid)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine computes particle condensation kernels
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     c_gas: aerosol gas phase concentration(µg/m^3)
!     Wet_diam: wet diameters (s)
!     q: aerosol mass concentration (µg/m^3)
!     qn: aerosol number concentration (#/m^3)
!     Temperature: temperature
!     ce_kernal_coef_i: c/e kernel coefficient          ([m3.s-1]).
!     jj: current bin index 
!
!     -- OUTPUT VARIABLES
!
!     ce_kernel: particle condensation kernels (µg/m^3/s)
!
!------------------------------------------------------------------------
    implicit none
    integer:: jesp,init,jj,s
    double precision:: qn,qext(N_aerosol),init_bulk_gas(N_aerosol)
    double precision:: qinti(N_inside_aer),ce_kernal_coef_i(N_aerosol)
    double precision:: surface_equilibrium_conc(N_aerosol),ce_kernel(N_aerosol)
    double precision:: Kelvin_effect(N_aerosol),Wet_vol
    double precision:: Wet_diam,Wet_diam_used,rhop
    double precision:: c_gas(N_aerosol)!micg/m^-3
    double precision:: q(N_aerosol)!mass concentration in current grid point
    double precision:: qih,emw_tmp,rhop_tmp,Temp
    double precision:: lwc,ionic,proton, liquid(12)
  

!!     ******Initialization to zero
    do s=1,N_aerosol
      jesp=List_species(s)
      Kelvin_effect(jesp)=1.D0
      ce_kernel(jesp)=0.D0
      init_bulk_gas(jesp)=0.D0
      surface_equilibrium_conc(jesp)=0.D0 !surface equilibrium concentration
    end do

    do s=1,N_aerosol
      jesp=List_species(s)
      init_bulk_gas(jesp)=c_gas(jesp)!initial bulk gas conc (µg.m-3)
      qext(jesp)=q(jesp)
    enddo

!     Aerosol wet density in µ g.µ m -3
    rhop = 0.d0
    do s=1,N_aerosol
      jesp=List_species(s)
      rhop=rhop+qext(jesp)
    enddo

    if(rhop.gt.0.d0) then

!calculate the equilibrium between aerosols and gas-phase
      call surface_eq_conc(qext,qinti,surface_equilibrium_conc,lwc,ionic,proton,liquid)

      ! we prevent evaporation when conc
      ! are too near from zero
      do s=1,(N_aerosol-1)
	jesp=List_species(s)
	if (qext(jesp).LE.TINYM) then
	  surface_equilibrium_conc(jesp)=0.D0
	endif
	if (surface_equilibrium_conc(jesp).lt.0.D0) then
	  surface_equilibrium_conc(jesp)=0.D0
	endif
      end do
  !     ******c/e kernel coefficient
      do s=1,nesp_isorropia
	jesp=isorropia_species(s)
	if (aerosol_species_interact(jesp).gt.0) then
	  call COMPUTE_CONDENSATION_TRANSFER_RATE(&
		diffusion_coef(jesp), &! diffusion coef (m2.s-1)
		quadratic_speed(jesp),& ! quadratic mean speed (m.s-1)
		accomodation_coefficient(jesp),& ! accomadation coef (adim)
		Wet_diam,   & ! wet aero diameter (Âµm)
		ce_kernal_coef_i(jesp) ) ! c/e kernel coef (m3.s-1)
	endif
      enddo

      Wet_vol=Wet_diam**3.d0*cst_pi6
      rhop_tmp = rho_wet_cell(jj) * 1.d9 !1400. ! kg/m3
      Wet_diam_used =DMAX1(Wet_diam,diam_bound(1)) 
      if (wet_diam_used .lt. 1.d-3) then
         write(*,*) "kercond: too small wet_diameter",wet_diam_used
         stop
      endif !! YK
      do s=1,(N_aerosol-1)
	jesp=List_species(s)
	if (aerosol_species_interact(jesp).gt.0) then!& .and.rhop.gt.1.d3
	  emw_tmp = molecular_weight_aer(jesp) * 1.D-6 ! g/mol
	  call COMPUTE_KELVIN_COEFFICIENT(&
		  Temp,&          ! temperature (Kelvin)
		  emw_tmp,&       ! ext mol weight (g.mol-1)
		  surface_tension(jesp),&   ! surface tension (N.m-1) from INC
		  Wet_diam_used,&         ! wet aero diameter (µm)
		  rhop_tmp,&      ! aerosol density (kg.m-3)
		  Kelvin_effect(jesp) )   ! kelvin effect coef (adim)
	endif
      enddo

  !!     ******Not limited c/e kernels

      do s=1,nesp_isorropia
	jesp=isorropia_species(s)
	if (aerosol_species_interact(jesp).gt.0) then
	  ce_kernel(jesp)= ce_kernal_coef_i(jesp)& ! kernel coef (m3.s-1)
	      *( init_bulk_gas(jesp)&    ! bulk gas conc (µg.m-3)
	      -surface_equilibrium_conc(jesp)&      ! equi gas conc (µg.m-3) unknow
	      *Kelvin_effect(jesp) )     ! kelvin coef (adim)
	endif
      enddo

      init=0
      do s=1,nesp_isorropia
	  jesp=isorropia_species(s)
	  if (qext(jesp).gt.TINYM) then
	    init=1
	  endif
      enddo

      if (qext(EH2O).eq.0.D0.AND.init.eq.1) then ! solid
	call DRYIN( Temp,&   ! local temperature (Kelvin)
	    qinti,&         ! int sld inorg conc (µg.m-3)
	    N_aerosol,&         ! size of vectors following below
	    init_bulk_gas,&         ! bulk gas conc (µg.m-3)
	    ce_kernal_coef_i,&           ! kernel coef (m3.s-1)
	    Kelvin_effect,&          ! kelvin coef (adim)
	    surface_equilibrium_conc,&         ! equi gas conc (µg.m-3)
	    ce_kernel )       ! modified c/e kernel
				! liq or mix : H+ limitation flux
      else
	if (qn.gt.0.d0) qih=qinti(IH)/qn    ! µg ,qn is number (qih is H+ per particl)
	!print*,'water go',qext(EH2O)
	call HPLFLIM( ALFHP,& ! percentage of H+ allowed to c/e(0.1)
	    qih,&            ! int H+ conc (µg)
	    N_aerosol,&          ! size of vectors following below
	    init_bulk_gas,&          ! bulk gas conc (µg.m-3)
	    ce_kernal_coef_i,&            ! kernel coef (m3.s-1)
	    Kelvin_effect,&           ! kelvin coef (adim)
	    surface_equilibrium_conc,&          ! equi gas conc (µg.m-3)
	    ce_kernel )        ! modified c/e kernel
      endif
      concentration_mass(jj,EH2O)=qext(EH2O)!water updated here
      do jesp=1,N_inside_aer
	concentration_inti(JJ,jesp)=qinti(jesp)
      enddo
   else
     concentration_number(jj)=0.d0
     concentration_mass(jj,EH2O)=0.d0
      do jesp=1,N_inside_aer
	concentration_inti(JJ,jesp)=0.d0
      enddo     
   endif
   
  end subroutine KERCOND
   
  subroutine surface_eq_conc(qext,qinti,surface_equilibrium_conc,lwc,ionic,proton,liquid)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine computes the local (surface) aerosol equilibrium within each bin.
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     qext: aerosol mass concentration(µg/m^3)
!     qinti: aerosol internal species mass concentration(µg/m^3)
!
!     -- OUTPUT VARIABLES
!
!     surface_equilibrium_conc: surface equilibrium concentration of aerosol species
!
!------------------------------------------------------------------------
    implicit none
    integer:: jesp,s
    double precision:: surface_equilibrium_conc(N_aerosol)!surface_equilibrium_conc : equilibrium gas concentration    ([\mu.g.m^-3]).
    double precision:: qext(N_aerosol)!QEXT : external aerosol concentration ([\mu.g.m^-3]).
    double precision:: qinti(N_inside_aer)!QINTI : internal inorganic concentration ([\mu.g.m^-3]).not used
    double precision:: qtinorg
    double precision:: lwc,ionic,proton,liquid(12)

!!     ******zero init
    do jesp=1,N_aerosol
      surface_equilibrium_conc(jesp)=0.D0
    end do

    !!     ******organics et inorganics thermodynamics
    do jesp=1,N_inside_aer!N_inside_aer=21
      qinti(jesp)=0.d0
    enddo
    !qext(EH2O)=0.d0
			      ! sum of inorganic mass
    qtinorg=0.D0
    do s=1,nesp_isorropia
      jesp=isorropia_species(s)
      qtinorg=qtinorg+qext(jesp)
    end do

    call EQINORG( N_aerosol,qext,&         ! ext inorg aero conc (µg.m-3)
		qinti,&         ! int inorg aero conc (µg.m-3)
		surface_equilibrium_conc,lwc,ionic,proton,liquid)        ! inorg eq gas conc (µg.m-3)

    do s=1,nesp_isorropia
      jesp=isorropia_species(s)
      if(surface_equilibrium_conc(jesp).lt.0.0) surface_equilibrium_conc(jesp)=0.d0
    enddo    
  end subroutine surface_eq_conc
   
End module fCondensation
