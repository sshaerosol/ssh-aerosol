!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
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
  subroutine ssh_SULFDYN(Q1,Q,N1,N,c_gas,dtx,time_step_sulf)
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
    double precision ::c_gas(N_aerosol)!micg/m^-3
    double precision :: ce_kernal_coef_tot ! c/e kernel coef (m3.s-1)
    double precision :: Q1(N_size,N_aerosol_layers) ! Mass concentration
    double precision :: Q(N_size,N_aerosol_layers) ! Mass concentration  !!BUG
    double precision :: N1(N_size) ! Number concentration
    double precision :: N(N_size) ! Number concentration
    double precision :: dtx,tmp,cond_so4!Time steps
    double precision :: dexploc,n2err,tmp_n2err,time_step_sulf
    double precision :: dqdt(N_size)

    jesp=ESO4!Pointer
    ce_kernal_coef_tot = 0.0d0
    cond_so4 = 0.0d0
    n2err = 0.d0

    do j = 1,N_size! Reassigned distribution by mass of each species
      call ssh_compute_condensation_transfer_rate(diffusion_coef(jesp), &
      quadratic_speed(jesp), accomodation_coefficient(jesp), &
      wet_diameter(j), dqdt(j))
      ce_kernal_coef_tot=ce_kernal_coef_tot+N(j)*dqdt(j)
    enddo

    do j = 1,N_size! Reassigned distribution by mass of each species
      if(ce_kernal_coef_tot.ne.0.d0) then
        dexploc = DEXP(-ce_kernal_coef_tot*dtx)
	tmp=(dqdt(j)*N(j)/ce_kernal_coef_tot)*&
	      (1.d0-dexploc) * c_gas(jesp)
        if (Q(j,jesp).GT.TINYM) then
          tmp_n2err = (dqdt(j)*N(j)*dtx*c_gas(jesp))/Q(j,jesp)
          n2err = n2err + tmp_n2err*tmp_n2err
        endif
        Q(j,jesp) = Q(j,jesp)+tmp!renew mass
        Q1(j,jesp) = Q(j,jesp)
        cond_so4 = cond_so4+tmp
        N1(j)=N(j)
        if(dtx.gt.0.d0) dqdt(j)=tmp/dtx!for redistribution
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
  end subroutine ssh_SULFDYN

  subroutine ssh_KERCOND(c_mass,c_number,qn,q,c_gas,Wet_diam,wet_mass,Temp,ce_kernel,ce_kernal_coef_i,jj,&
                     lwc,ionic,proton,liquid,qtot,cond_time,iker)
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
!     cond_time: for the computation of cut-off diameter
!		 if tag_icut=1 return the c/e timescales (s)
!		 if tag_icut=3 return the QSSA value for ENH4, EHNO3, and EHCL 
!
!------------------------------------------------------------------------
    implicit none
    integer:: jesp,init,jj,s,iker
    double precision:: c_mass(N_size,N_aerosol_layers)
    double precision:: c_number(N_size)
    double precision:: qn,qext(N_aerosol),init_bulk_gas(N_aerosol)
    double precision:: qinti(N_inside_aer),ce_kernal_coef_i(N_aerosol)
    double precision:: surface_equilibrium_conc(N_aerosol),ce_kernel(N_aerosol)
    double precision:: Kelvin_effect(N_aerosol),Wet_vol
    double precision:: Wet_diam,Wet_diam_used,rhop,wet_mass
    double precision:: c_gas(N_aerosol)!micg/m^-3
    double precision:: q(N_aerosol)!mass concentration in current grid point
    double precision:: qih,emw_tmp,rhop_tmp,Temp,dry_diam
    double precision:: lwc,ionic,proton, liquid(12),qtot
    double precision:: vad,rhoaer,cond_time(3)
    double precision:: pplusl,pminusl
    double precision:: dry_density,wet_density,dry_mass,wmass,dry_to_wet

!!     ******Initialization to zero
    do jesp=1,N_aerosol
      Kelvin_effect(jesp)=1.D0
      ce_kernel(jesp)=0.D0
      init_bulk_gas(jesp)=0.D0
      surface_equilibrium_conc(jesp)=0.D0 !surface equilibrium concentration
    end do

    do jesp=1,N_aerosol
      init_bulk_gas(jesp)=c_gas(jesp)!initial bulk gas conc (µg.m-3)
      qext(jesp)=q(jesp)
    enddo

!     Aerosol wet density in µ g.µ m -3
    rhop = 0.d0
    do jesp=1,N_aerosol
      rhop=rhop+qext(jesp)
    enddo

    if(rhop.gt.0.d0) then

!calculate the equilibrium between aerosols and gas-phase
      if(iker.EQ.0) then
        call ssh_surface_eq_conc(qext,qinti,surface_equilibrium_conc,lwc,ionic,proton,liquid,jj)
        concentration_mass(jj,EH2O_layers)=qext(EH2O)!water updated here
        Do jesp=1,N_aerosol
          surface_equilibrium_conc_nsize(jj,jesp) = surface_equilibrium_conc(jesp)
        Enddo 
      else
        qext(EH2O) = concentration_mass(jj,EH2O_layers)
        Do jesp=1,N_aerosol
          surface_equilibrium_conc(jesp) = surface_equilibrium_conc_nsize(jj,jesp) 
        Enddo 
      endif
      !     Compute wet_diam and wet_mass
      if (with_fixed_density == 2) then
         call ssh_get_nonlinear_density(qext,dry_density,wet_density,dry_mass,wmass,dry_to_wet)
         rhoaer = wet_density
         vad=dry_mass/dry_density !qti/rhoaer!qti total dry mass            
         rho_wet_cell(jj)=wet_density
             
         wet_vol=wmass/wet_density !: wet volume aerosol concentration (µm3/m3).         
         dry_diam=(vad/c_number(jj)/cst_pi6)**cst_FRAC3 ! dry aerosol dimaeter µm             
         wet_diam=dry_to_wet*dry_diam ! wet aerosol diameter µm
         wet_mass=wmass/c_number(jj) !single wet mass (µg)
      else
         if (with_fixed_density == 0) then
            call ssh_compute_density(N_size,N_aerosol_layers,EH2O_layers,TINYM,c_mass,&
                 mass_density_layers,jj,rhoaer)             
         else
            rhoaer = fixed_density 
         endif
         if(rhoaer.gt.0.d0) then
            vad=qtot/rhoaer!qtot total dry mass
            wet_vol=vad+qext(EH2O)/rhoaer!: wet volume aerosol concentration (µm3/m3).
            dry_diam=(vad/c_number(jj)/cst_pi6)**cst_FRAC3 ! dry aerosol dimaeter µm
            wet_diam=((wet_vol)/c_number(jj)/cst_pi6)**cst_FRAC3 ! wet aerosol diameter µm
            wet_diam=DMAX1(wet_diam,dry_diam)!wet diameter is always larger than dry diameter
            wet_mass=(qtot+qext(EH2O))/c_number(jj) ! single wet mass (µg)
         endif
      endif
      ! we prevent evaporation when conc
      ! are too near from zero
      do s=1,nesp_isorropia
	jesp=isorropia_species(s)
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
	  call ssh_COMPUTE_CONDENSATION_TRANSFER_RATE(&
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
      endif
      do s=1,nesp_isorropia
	jesp=isorropia_species(s)
	if (aerosol_species_interact(jesp).gt.0) then!& .and.rhop.gt.1.d3
	  emw_tmp = molecular_weight_aer(jesp) * 1.D-6 ! g/mol
	  call ssh_COMPUTE_KELVIN_COEFFICIENT(&
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

      if(iker.EQ.0) then
        init=0
        do s=1,nesp_isorropia
	  jesp=isorropia_species(s)
	  if (qext(jesp).gt.TINYM) then
	    init=1
	  endif
        enddo

        if (qext(EH2O).eq.0.D0.AND.init.eq.1) then ! solid
	  call ssh_DRYIN( Temp,&   ! local temperature (Kelvin)
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
	   call ssh_HPLFLIM( ALFHP,& ! percentage of H+ allowed to c/e(0.1)
	    qih,&            ! int H+ conc (µg)
	    N_aerosol,&          ! size of vectors following below
	    init_bulk_gas,&          ! bulk gas conc (µg.m-3)
	    ce_kernal_coef_i,&            ! kernel coef (m3.s-1)
	    Kelvin_effect,&           ! kelvin coef (adim)
	    surface_equilibrium_conc,&          ! equi gas conc (µg.m-3)
	    ce_kernel )        ! modified c/e kernel
        endif
        do jesp=1,N_inside_aer
  	  concentration_inti(JJ,jesp)=qinti(jesp)
        enddo
        Do jesp=1,N_aerosol
          surface_equilibrium_conc_nsize(jj,jesp) = surface_equilibrium_conc(jesp)
        Enddo 
      endif
    else
!!       concentration_number(jj)=0.d0
        concentration_mass(jj,EH2O_layers)=0.d0
        do jesp=1,N_inside_aer
  	  concentration_inti(JJ,jesp)=0.d0
        enddo     
    endif

   ! compute c/e timescales (tag_icut=1) or QSSA factors (tag_icut=3) for cut-off diameter
   if (tag_icut.eq.1 .or. tag_icut.eq.3) then
      do s=1,3
	   ! get index of ENH4, ENO3, ECL
	   jesp=cond_time_index(s)
	   ! c/e timescales ! unit: seconds
            if (tag_icut.eq.1) then
		pminusl=concentration_mass(jj,jesp)
		pplusl=kelvin_effect(jesp)*surface_equilibrium_conc(jesp)*&
			       ce_kernal_coef_i(jesp)*concentration_number(jj)
	   ! modified QSSA criteria
            elseif (tag_icut.eq.3) then
	        pminusl = dabs(init_bulk_gas(jesp)-kelvin_effect(jesp)*&
				surface_equilibrium_conc(jesp))
	        pplusl = init_bulk_gas(jesp)+kelvin_effect(jesp)*&
				surface_equilibrium_conc(jesp)
	    endif
	   ! compute cond_time
           if ((pminusl.gt.TINYM).and.(pplusl.gt.TINYM)) then
	    	cond_time(s) = pminusl/pplusl
           else
	    	cond_time(s) = 0.0
           endif
      enddo
   endif

  end subroutine ssh_KERCOND
   
  subroutine ssh_surface_eq_conc(qext,qinti,surface_equilibrium_conc,lwc,ionic,proton,liquid,jbin)
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

    integer :: jbin

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

    call ssh_EQINORG( N_aerosol,qext,&         ! ext inorg aero conc (µg.m-3)
		qinti,&         ! int inorg aero conc (µg.m-3)
		surface_equilibrium_conc,lwc,ionic,proton,liquid,jbin)        ! inorg eq gas conc (µg.m-3)

    do s=1,nesp_isorropia
      jesp=isorropia_species(s)
      if(surface_equilibrium_conc(jesp).lt.0.0) surface_equilibrium_conc(jesp)=0.d0
    enddo    
  end subroutine ssh_surface_eq_conc
   
End module fCondensation
