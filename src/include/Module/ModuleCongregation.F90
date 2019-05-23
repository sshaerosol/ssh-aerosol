!!-----------------------------------------------------------------------
!!     Copyright (C) 2012-2018, ENPC - EDF R&D - INERIS
!!     Author(s): Karine Sartelet, Edouard Debry, Shupeng Zhu
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
!!    This module contains methods for for mass and number redistribution
!!    for fixed euler scheme in 3D application.
!!-----------------------------------------------------------------------
Module hCongregation
  use dPhysicalbalance
  use bCoefficientRepartition
  use aInitialization
  use fCondensation
  use gCoagulation

  implicit none

contains
  ! subroutine fgde(c_mass,c_number,c_gas,dqdt,dndt,ce_kernal_coef, &
  !                 lwc_Nsize,ionic_Nsize,proton_Nsize,liquid_Nsize)
  ! subroutine fgde(c_mass,c_number,c_gas,dqdt,dndt,ce_kernal_coef, &
  !                 ionic_Nsize,proton_Nsize,liquid_Nsize)
  subroutine fgde(c_mass,c_number,c_gas,dqdt,dndt,ce_kernal_coef)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine provides entries for different aerosol dynamic process
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     c_gas: aerosol gas phase concentration(µg/m^3)
!     c_mass: aerosol mass concentration (µg/m^3)
!     c_number: aerosol number concentration (#/m^3)
!     ce_kernal_coef: c/e kernel coefficient          ([m3.s-1]).
!
!     -- OUTPUT VARIABLES
!
!     dqdt: particle mass derivation(µg/m^3/s)
!     dndt: particle number derivation(µg/m^3/s)
!
!------------------------------------------------------------------------
    implicit none
    double precision ::c_number(N_size)
    double precision ::c_mass(N_size,N_aerosol)
    double precision ::dqdt(N_size,N_aerosol)
    double precision ::dndt(N_size)
    double precision ::c_gas(N_aerosol)
    double precision ::ce_kernal_coef(N_size,N_aerosol)
!    double precision:: lwc_Nsize(N_size),ionic_Nsize(N_size)
!    double precision:: ionic_Nsize(N_size)
!    double precision:: proton_Nsize(N_size),liquid_Nsize(12,N_size)

    if (tag_nucl.eq.1) then
      call fgde_nucl(c_mass,c_number,c_gas,dqdt,dndt)
    endif

    if (tag_coag.eq.1) then
      call fgde_coag (c_mass,c_number,dqdt,dndt)
    endif

    if (tag_cond.eq.1) then
     ! call fgde_cond(c_mass,c_number,c_gas,dqdt,dndt,ce_kernal_coef,&
     !               ionic_Nsize,proton_Nsize,liquid_Nsize)
     call fgde_cond(c_mass,c_number,c_gas,dqdt,dndt,ce_kernal_coef)
    endif

  end subroutine fgde

  ! subroutine fgde_cond(c_mass,c_number,c_gas,dqdt,dndt,ce_kernal_coef,&
  !                     lwc_Nsize,ionic_Nsize,proton_Nsize,liquid_Nsize)
  ! subroutine fgde_cond(c_mass,c_number,c_gas,dqdt,dndt,ce_kernal_coef,&
  !                     ionic_Nsize,proton_Nsize,liquid_Nsize)
  subroutine fgde_cond(c_mass,c_number,c_gas,dqdt,dndt,ce_kernal_coef)

!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine computes first order derivatives related to C/E process
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     c_gas: aerosol gas phase concentration(µg/m^3)
!     c_mass: aerosol mass concentration (µg/m^3)
!     c_number: aerosol number concentration (#/m^3)
!     ce_kernal_coef: c/e kernel coefficient          ([m3.s-1]).
!
!     -- OUTPUT VARIABLES
!
!     dqdt: particle mass derivation(µg/m^3/s)
!     dndt: particle number derivation(µg/m^3/s)
!
!------------------------------------------------------------------------
    implicit none
    integer::j,jesp,s
    double precision:: dqdt(N_size,N_aerosol)
    double precision:: dndt(N_size)
    double precision :: ce_kernel(N_aerosol)
    double precision :: ce_kernal_coef(N_size,N_aerosol)
    double precision :: ce_kernal_coef_i(N_aerosol)
    double precision :: qn		!number concentration in current grid point
    double precision :: q(N_aerosol)	!mass concentration in current grid point
    double precision :: c_gas(N_aerosol)
    double precision:: c_mass(N_size,N_aerosol)
    double precision:: c_number(N_size),wet_mass(N_size)
    double precision:: wet_diam(N_size),wet_vol(N_size),cell_diam(N_size)
!    double precision:: lwc_Nsize(N_size),ionic_Nsize(N_size)
!    double precision:: ionic_Nsize(N_size)

!    double precision:: proton_Nsize(N_size),liquid(12)
    double precision:: liquid(12)
!    double precision:: liquid_Nsize(12,N_size)

     call update_wet_diameter_liquid(1,N_size,c_mass,c_number, &
          wet_mass,wet_diam,wet_vol,cell_diam)

    dqdt=0.d0
    dndt=0.d0
    call mass_conservation(c_mass,c_number,c_gas,total_mass)

    do j =(ICUT+1), N_size
      dndt(j)=0.d0
      qn=c_number(j)!initial number and mass
      !if(qn.gt.TINYN) then!skip empty cell
	do s=1,N_aerosol
	    jesp=List_species(s)
	    q(jesp)=c_mass(j,jesp)
	    ce_kernal_coef_i(jesp)=ce_kernal_coef(j,jesp)
	enddo
	call KERCOND(qn,q,c_gas,wet_diam(j),temperature,ce_kernel,ce_kernal_coef_i,j, &
                     lwc_Nsize(j),ionic_Nsize(j),proton_Nsize(j),liquid)

	do s=1,12
           liquid_Nsize(s,j) = liquid(s)
        enddo 
	!calculate the C/E kernal
	do s=1, nesp_isorropia
	  jesp = isorropia_species(s)
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
          IF (jesp.NE.ECl) THEN
#endif
	    ce_kernal_coef(j,jesp)=ce_kernal_coef_i(jesp)
	    dqdt(j,jesp)=c_number(j)*ce_kernel(jesp)


#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
          ENDIF
#endif
	  if(sulfate_computation.eq.1) then!do not take sulfate into account
	    dqdt(j,ESO4)=0.d0
	  endif
	enddo
     ! endif
    enddo

    do j=1, N_size
       do s=1,N_aerosol
          if (isnan(dqdt(j,s))) then
             write(*,*) j,s,dqdt(j,s)
             stop
          endif
       enddo
    enddo

  end subroutine fgde_cond

  subroutine fgde_coag (c_mass,c_number,rate_mass,rate_number)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine computes the coagulation rate
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     c_gas: aerosol gas phase concentration(µg/m^3)
!     c_mass: aerosol mass concentration (µg/m^3)
!     c_number: aerosol number concentration (#/m^3)
!
!     -- OUTPUT VARIABLES
!
!     rate_mass: particle mass derivation(µg/m^3/s)
!     rate_number: particle number derivation(µg/m^3/s)
!
!------------------------------------------------------------------------
    implicit none
    integer :: j
    double precision ::c_number(N_size)
    double precision ::c_mass(N_size,N_aerosol)
    double precision ::rate_number(N_size)
    double precision ::rate_mass(N_size,N_aerosol)
    
    rate_number=0.d0
    rate_mass=0.d0

    do j = 1,N_size
       total_number=total_number+c_number(j)
    enddo !! Added by YK

    if(total_mass_t*total_number.ne.0.d0) then
      
!       call compute_average_diameter()
! 
!       do j1 = 1, N_size
! 	do j2 = 1, N_size
! 	  call compute_bidisperse_coagulation_kernel(Temperature,air_free_mean_path,&
! 		    cell_diam_av(j1),cell_diam_av(j2),&
! 		    cell_mass_av(j1),cell_mass_av(j2), kernel_coagulation(j1,j2))
! 	enddo
!       enddo

      call  Rate(rate_number,rate_mass,c_number,c_mass)

    endif

    !check rate diameter
    !call check_diam_fraction(rate_mass,rate_number)

  end subroutine fgde_coag
  
  subroutine fgde_nucl(c_mass,c_number,c_gas,dqdt,dndt)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine computes source terms for the system of
!     Ordinary Differential Equations defined by nucleation.
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     c_gas: aerosol gas phase concentration(µg/m^3)
!     c_mass: aerosol mass concentration (µg/m^3)
!     c_number: aerosol number concentration (#/m^3)
!
!     -- OUTPUT VARIABLES
!
!     dqdt: particle mass derivation(µg/m^3/s)
!     dndt: particle number derivation(µg/m^3/s)
!
!------------------------------------------------------------------------
    implicit none
    double precision:: dqdt(N_size,N_aerosol)
    double precision:: dndt(N_size)
    double precision :: c_gas(N_aerosol)
    double precision :: c_mass(N_size,N_aerosol)
    double precision :: c_number(N_size)
    double precision	:: mr,na
    double precision	:: jnucl,ntot,ntotnh3
    double precision	:: dpnucl,nanucl,qanucl
    double precision	:: xstar,mSO4
!     Compute gas mass conservation

      call mass_conservation(c_mass,c_number,c_gas, total_mass)
      
      if(nucl_model.eq.1) then            ! sulfuric-acid-ammonia-water nucl'n
!     mr should be in ppt
         mr = 10.d12 * c_gas(ENH4)/(Pressure*MMair/RGAS*Temperature)
         
         na= c_gas(ESO4)*1.D-06&       ! convert to µg.cm-3
	        /molecular_weight_aer(ESO4)&        ! to mol.m-3
                *Navog            ! to #molec.m-3      
                
         call COMPUTE_TERNARY_NUCLEATION(Relative_Humidity,&     ! relative humidity 
		Temperature,&             ! temperature (Kelvin)
		na,&               ! gas h2so4 conc (#molec.cm-3)
		mr,&		      !Mixing ratio of NH3 (ppt).
		jnucl,&           ! nucleation rate (#part.cm-3.s-1)
		ntot,&             ! number of molec of h2so4 in nucleus 
		ntotnh3,&          ! number of molec of nh3 in nucleus 
		dpnucl )          ! nucleation diameter (nm)    
		
                                ! nucleation rate (#part.m-3.s-1)
         jnucl=jnucl*1.D06
         if(Navog.ne.0.d0.and.(.not.IsNaN(jnucl*0.d0))) then
	  dndt(1) =dndt(1) +jnucl ! #part.m-3.s-1
	  dqdt(1,ESO4)=dqdt(1,ESO4)+jnucl*ntot/Navog*molecular_weight_aer(ESO4) ! µg.m-3.s-1
	  dqdt(1,ENH4)=dqdt(1,ENH4)+jnucl*ntotnh3/Navog*molecular_weight_aer(ENH4)
         endif

         if(IsNaN(dndt(1)*0.d0)) then
	   dndt(1)=0.d0
	   dqdt(1,ESO4)=0.d0
	   dqdt(1,ENH4)=0.d0
	 endif
         
      else                      !sulfuric-acid-water nucl'n    !we should use this  
!     Compute H2SO4 threshold concentration

         call NA_THRESHOLD_VEAHKAMAKI(Relative_Humidity,Temperature,nanucl) !#molec.cm-3

         qanucl= nanucl*1.D06&   ! convert to #molec.m-3
		/Navog&            ! to mol.m-3
		*molecular_weight_aer(ESO4)        ! to µg.mol-1      
     
!     Compute nucleation kernel if qSO4 exceed qanucl  

         if (c_gas(ESO4).GE.qanucl) then
            
            na= c_gas(ESO4)*1.D-06&    ! convert to µg.cm-3
		  /molecular_weight_aer(ESO4)&     ! to mol.m-3
		  *Navog         ! to #molec.m-3

            call COMPUTE_BINARY_NUCLEATION_KERNEL( Relative_Humidity,& ! relative humidity 
		    Temperature,&          ! temperature (Kelvin)
		    na,&            ! gas h2so4 conc (#molec.cm-3)
		    jnucl,&         ! nucleation rate (#part.cm-3.s-1)
		    ntot,&          ! num of molec in nucleus 
		    xstar,&         ! mol fraction of h2so4
		    dpnucl )       ! nucleation diameter (nm)
            
                                ! nucleation rate (#part.m-3.s-1)
            jnucl=jnucl*1.D06
	  if(Navog.ne.0.d0.and.(.not.IsNaN(jnucl*0.d0))) then
				  ! h2so4 mass in nucleus (µg)
	      mSO4= ntot&          ! #molec
		    /Navog&         ! Avogadro number (adim)
		    *xstar&         ! mol fraction of h2so4
		    *molecular_weight_aer(ESO4)     ! mol weight µg.mol-1
		    
	    dndt(1) =dndt(1) +jnucl ! #part.m-3.s-1
	    dqdt(1,ESO4)=dqdt(1,ESO4)+jnucl*mSO4! µg.m-3.s-1
	  endif

         if(IsNaN(dndt(1)*0.d0)) then
	   dndt(1)=0.d0
	   dqdt(1,ESO4)=0.d0
	 endif
	  
         endif

      endif     
   end subroutine fgde_nucl
  
end module hCongregation
