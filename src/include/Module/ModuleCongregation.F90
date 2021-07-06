!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
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
  subroutine ssh_fgde(c_mass,c_number,c_gas,dqdt,dndt,ce_kernal_coef,qH2O,cond_time,iker)
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
    !     qH2O: mass of water absorbed by inorganics
    !     cond_time: the condensation/evaporation characteristic 
    !			timescales or the QSSA values for inorganic species
    !
    !------------------------------------------------------------------------
    implicit none
    double precision ::c_number(N_size),qH2O(N_size)
    double precision ::c_mass(N_size,N_aerosol_layers)
    double precision ::dqdt(N_size,N_aerosol_layers)
    double precision ::dndt(N_size)
    double precision ::c_gas(N_aerosol),cond_time(N_size,3)
    double precision ::ce_kernal_coef(N_size,N_aerosol)
    double precision:: wet_diam(N_size),wet_mass(N_size)
    double precision:: wet_vol(N_size),cell_diam(N_size)
    integer :: j,iker

    dqdt = 0.d0
    dndt = 0.d0
    ! Compute wet_mass_diameter using isorropia if wet_diam_estimation == 0 but condensation is not called, 
    ! or if coagulation is splitted (that means that condensation and coagulation are not solved simultaneously
    ! For dynamic sections, the wet_mass_diameter is computed in fgde_cond

    Do j=1,N_size
       if(concentration_index(j, 1) <= ICUT) then
          c_mass(j,N_aerosol_layers) = 0.d0
       endif
    Enddo
    call ssh_update_wet_diameter_liquid(ICUT,c_mass,c_number, &
                                      wet_mass,wet_diam,wet_vol,cell_diam)
    if(((wet_diam_estimation.eq.0).AND.(tag_cond.eq.0)).OR. &
            ((splitting.eq.0).AND.(tag_coag.eq.1))) then !with isorropia
        call ssh_compute_wet_mass_diameter(1,N_size,c_mass,c_number,&
                                   concentration_inti,wet_mass,wet_diam,wet_vol)
    endif
    !call ssh_mass_conservation(c_mass,c_number,c_gas,total_mass)

    if (tag_nucl.eq.1) then
       call ssh_fgde_nucl(c_mass,c_number,c_gas,dqdt,dndt)
    endif

    if (tag_cond.eq.1) then
       call ssh_fgde_cond(c_mass,c_number,c_gas,dqdt,dndt,ce_kernal_coef,qH2O,wet_diam,wet_mass,cond_time,iker)
    endif

    if (tag_coag.eq.1) then
       call ssh_fgde_coag (c_mass,c_number,dqdt,dndt,wet_diam,wet_mass)
    endif

  end subroutine ssh_fgde

  subroutine ssh_fgde_cond(c_mass,c_number,c_gas,dqdt,dndt,ce_kernal_coef,qH2O,wet_diam,wet_mass,cond_time,iker)

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
    !     qH2O: mass of water absorbed by inorganics
    !     cond_time: the condensation/evaporation characteristic 
    !			timescales or the QSSA values for inorganic species
    !
    !------------------------------------------------------------------------
    implicit none
    integer:: iker,j,jesp,s,lay
    double precision:: dqdt(N_size,N_aerosol_layers),qH2O(N_size)
    double precision:: dndt(N_size)
    double precision :: ce_kernel(N_aerosol)
    double precision :: ce_kernal_coef(N_size,N_aerosol)
    double precision :: ce_kernal_coef_i(N_aerosol)
    double precision :: qn		!number concentration in current grid point
    double precision :: q(N_aerosol)	!mass concentration in current grid point
    double precision :: c_gas(N_aerosol)
    double precision:: c_mass(N_size,N_aerosol_layers)
    double precision:: c_number(N_size),wet_mass(N_size)
    double precision:: wet_diam(N_size),wet_vol(N_size),cell_diam(N_size)
    double precision:: cond_time_char(3),cond_time(N_size,3)
    double precision:: liquid(12),qtot

    cond_time_char = 0.0
    
    !init cond_time for c/e timescale and QSSA
    if (tag_icut.eq.1.or.tag_icut.eq.3) cond_time=0.d0

    if (soap_inorg.eq.0) then
       do j =1, N_size
          if(concentration_index(j, 1) > ICUT) then! k : index of size bins
             qn=c_number(j)!initial number and mass
             if(qn.GT.TINYN) then
                do s=1,N_aerosol
                   q(s) = 0.d0
                enddo
                qtot = 0.0
                do s=1,N_aerosol_layers
                   jesp=List_species(s)
                   q(jesp)=q(jesp)+ c_mass(j,s)
                   qtot =qtot + c_mass(j,s)
                   ce_kernal_coef_i(jesp) = 0.d0 !ce_kernal_coef(j,jesp)
                enddo
                ! do s=N_nonorganics+1,N_aerosol
                !   q(s) = 0.0
                !   ce_kernal_coef_i(s) = 0.0
                ! enddo
                call ssh_KERCOND(c_mass,c_number,qn,q,c_gas,wet_diam(j),wet_mass(j),temperature, &
                     ce_kernel,ce_kernal_coef_i,j, &
                     lwc_Nsize(j),ionic_Nsize(j),proton_Nsize(j),liquid,qtot,cond_time_char,iker)

                ! store condensation/evaporation characteristic timescales or modified QSSA
                if (tag_icut.eq.1 .or. tag_icut.eq.3) then
                   do s=1,3
                      cond_time(j,s) = cond_time_char(s)
                   enddo
                endif

                qH2O(j) = q(EH2O)
                if(iker.EQ.0) then
                   do s=1,12
                      liquid_Nsize(s,j) = liquid(s)
                   enddo
                endif
                !calculate the C/E kernal
                do s=1, nesp_isorropia
                   jesp = isorropia_species(s)
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
                   IF (jesp.NE.ECl) THEN
#endif
                      ce_kernal_coef(j,jesp)=ce_kernal_coef_i(jesp)
                      if (jesp.NE.ESO4) then ! SO4 is computed either with non volatile species or in bulk equilibrium
                         dqdt(j,jesp)=dqdt(j,jesp)+c_number(j)*ce_kernel(jesp)
                      endif

#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
                   ENDIF
#endif
                enddo
             endif
          endif
       enddo
    endif

   ! Compute dynamically low-volatility organics
    do j=1,N_size
       qn=c_number(j)!initial number and mass
       do s=1,N_aerosol
          if(s.LE.N_nonorganics) then
             lay = 1 
          else
             lay = nlayer
          endif
          jesp = index_species(s,lay)
          if((inon_volatile(s).EQ.1).OR.  &
              ((inon_volatile(s).EQ.0).AND.(concentration_index(j,1)>ICUT).AND.(s.EQ.ESO4))) then
             call ssh_COMPUTE_CONDENSATION_TRANSFER_RATE(&
                diffusion_coef(s), &! diffusion coef (m2.s-1)
                quadratic_speed(s),& ! quadratic mean speed (m.s-1)
                accomodation_coefficient(s),& ! accomadation coef (adim)
                wet_diam(j),   & ! wet aero diameter (µm)
                ce_kernal_coef_i(s) ) ! c/e kernel coef (m3.s-1)
             ce_kernal_coef(j,s)=ce_kernal_coef_i(s)    ! bulk gas conc (ug.m-3)
             ce_kernel(s)=ce_kernal_coef_i(s) * c_gas(s)    ! bulk gas conc (ug.m-3)
             dqdt(j,jesp)=dqdt(j,jesp)+c_number(j)*ce_kernel(s)
          endif
       enddo
    enddo

    do j=1, N_size
       do s=1,N_aerosol_layers
          if (isnan(dqdt(j,s))) then
             write(*,*) "ModuleCongregation", j,s,dqdt(j,s)
             stop
          endif
       enddo
    enddo

  end subroutine ssh_fgde_cond

  subroutine ssh_fgde_coag (c_mass,c_number,rate_mass,rate_number,wet_diam,wet_mass)
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
    integer :: j,j1,j2
    double precision ::c_number(N_size)
    double precision ::c_mass(N_size,N_aerosol_layers)
    double precision ::rate_number(N_size)
    double precision ::rate_mass(N_size,N_aerosol)
    double precision:: wet_mass(N_size)
    double precision:: wet_diam(N_size),wet_vol(N_size),cell_diam(N_size)

    double precision :: total_number

    total_number = 0.d0

    do j = 1,N_size
       total_number=total_number+c_number(j)
    enddo

    if(total_number > 0.d0) then

       do j1 = 1, N_size
          do j2 = 1, N_size
             call ssh_compute_bidisperse_coagulation_kernel(Temperature,air_free_mean_path,&
                  wet_diam(j1),wet_diam(j2),&
                  wet_mass(j1),wet_mass(j2), kernel_coagulation(j1,j2))
          enddo
       enddo

       call ssh_Rate(rate_number,rate_mass,c_number,c_mass)

    endif

    !check rate diameter
    !call ssh_check_diam_fraction(rate_mass,rate_number)

  end subroutine ssh_fgde_coag

  subroutine ssh_fgde_nucl(c_mass,c_number,c_gas,dqdt,dndt)
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
    double precision:: dqdt(N_size,N_aerosol_layers)
    double precision:: dndt(N_size)
    double precision :: c_gas(N_aerosol)
    double precision :: c_mass(N_size,N_aerosol_layers)
    double precision :: c_number(N_size)
    double precision	:: mr,na
    double precision	:: jnucl,ntot,ntotnh3
    double precision	:: dpnucl,nanucl,qanucl
    double precision	:: xstar,mSO4,nanh3,dmdt
    double precision :: perc_so4, perc_nh4,tot_so4nh4
    integer :: iterp,s,isection
    !     Compute gas mass conservation

    if(nucl_model.eq.1.OR.nucl_model.eq.2) then            ! sulfuric-acid-ammonia-water nucl'n
       !     mr should be in ppt
       mr = c_gas(ENH4)/17. * 22.41d0 *273.d0/Temperature * Pressure/101300.d0

       na= c_gas(ESO4)*1.D-06&       ! convert to µg.cm-3
            /98.d6 &        ! to mol.cm-3
            *Navog            ! to #molec.cm-3      

       if(nucl_model.eq.1) then            ! Napari
          call ssh_COMPUTE_TERNARY_NUCLEATION(Relative_Humidity,&     ! relative humidity 
               Temperature,&             ! temperature (Kelvin)
               na,&               ! gas h2so4 conc (#molec.cm-3)
               mr,&		      !Mixing ratio of NH3 (ppt).
               jnucl,&           ! nucleation rate (#part.cm-3.s-1)
               ntot,&             ! number of molec of h2so4 in nucleus 
               ntotnh3,&          ! number of molec of nh3 in nucleus 
               dpnucl )          ! nucleation diameter (nm)    
       else ! Merikanto
          call ssh_COMPUTE_TERNARY_NUCLEATION_MERIKANTO(Relative_Humidity,&     ! relative humidity 
               Temperature,&             ! temperature (Kelvin)
               na,&               ! gas h2so4 conc (#molec.cm-3)
               mr,&		      !Mixing ratio of NH3 (ppt).
               jnucl,&           ! nucleation rate (#part.cm-3.s-1)
               ntot,&             ! number of molec of h2so4 in nucleus 
               ntotnh3,&          ! number of molec of nh3 in nucleus 
               dpnucl )          ! nucleation diameter (nm)    
       endif
       ! nucleation rate (#part.m-3.s-1)
       jnucl=jnucl*1.D06
       if(jnucl.GT.0.0) then
          dpnucl = dpnucl * 0.001
          isection = 1 !0
          !    if(dpnucl.LT.diam_bound(1)) then
          !       isection=1  ! Consider diameters even if below diam_bound(1)
          !       !jnucl = 0
          !    else
          !       s = 2
          !       Do while ((isection == 0).AND.(s < N_sizebin+1))
          !          if(dpnucl.LT.diam_bound(s)) then
          !             isection = s-1
          !          endif
          !          s = s+1
          !       Enddo
          !    endif

          if(jnucl.gt.0.d0.and.(.not.IsNaN(jnucl*0.d0))) then
            if(ntot.GT.0.d0.OR.ntotnh3.GT.0.0) then
	     dndt(isection) =dndt(isection) +jnucl ! #part.m-3.s-1
             dpnucl = size_diam_av(1) 
             dmdt = jnucl * PI/6.0 * dpnucl**3
             perc_so4 = ntot * molecular_weight_aer(ESO4)
             perc_nh4 = ntotnh3 * molecular_weight_aer(ENH4)
             tot_so4nh4 = perc_so4 + perc_nh4
             perc_so4 = perc_so4/tot_so4nh4
             perc_nh4 = perc_nh4/tot_so4nh4
             !dqdt(isection,ESO4)=dqdt(isection,ESO4)+jnucl*ntot/Navog*molecular_weight_aer(ESO4) ! µg.m-3.s-1
	     !dqdt(isection,ENH4)=dqdt(isection,ENH4)+jnucl*ntotnh3/Navog*molecular_weight_aer(ENH4)
	     dqdt(isection,ESO4)=dqdt(isection,ESO4)+ dmdt * perc_so4 * fixed_density ! µg.m-3.s-1
	     dqdt(isection,ENH4)=dqdt(isection,ENH4)+ dmdt * perc_nh4 * fixed_density

             if(IsNaN(dndt(isection)*0.d0)) then
                dndt(isection)=0.d0
                dqdt(isection,ESO4)=0.d0
                dqdt(isection,ENH4)=0.d0
             endif
            endif
          endif

       endif

    else                     
       if(nucl_model.eq.0) then   !sulfuric-acid-water nucl'n   
          !     Compute H2SO4 threshold concentration

          call ssh_NA_THRESHOLD_VEAHKAMAKI(Relative_Humidity,Temperature,nanucl) !#molec.cm-3

          qanucl= nanucl*1.D06&   ! convert to #molec.m-3
               /Navog&            ! to mol.m-3
               *molecular_weight_aer(ESO4)        ! to µg.mol-1      

          !     Compute nucleation kernel if qSO4 exceed qanucl  

          if (c_gas(ESO4).GE.qanucl) then

             na= c_gas(ESO4)*1.D-06&    ! convert to µg.cm-3
                  /molecular_weight_aer(ESO4)&     ! to mol.m-3
                  *Navog         ! to #molec.m-3

             call ssh_COMPUTE_BINARY_NUCLEATION_KERNEL( Relative_Humidity,& ! relative humidity 
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
       else
          write(*,*) 'nucleation scheme not implemented',nucl_model
          stop
       endif

    endif
  end subroutine ssh_fgde_nucl

end module hCongregation
