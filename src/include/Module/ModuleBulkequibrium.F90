!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------
!!
!!     -- DESCRIPTION
!!    This module contains methods to solve particle condensation/evaporation
!!    based on the gas/aerosol bulk equilibrium method.
!!-----------------------------------------------------------------------
Module iBulkequibrium
  use aInitialization
  use cThermodynamics
  use kSOAP
  implicit none
contains

  subroutine ssh_bulkequi_org(nesp_eq, &
       lwc, watorg, ionic, proton, liquid, delta_t)
    !------------------------------------------------------------------------
    !
    !     -- DESCRIPTION
    !     This subroutine solves gas/aerosol bulk equilibrium for organic species
    !     based on H2O method.
    !     equilibrium will be established between all size bin and bulk gas
    !     mass flux will be redistributed into each cell based on rates
    !
    !------------------------------------------------------------------------
    !
    !     -- INPUT VARIABLES
    !------------------------------------------------------------------------   
    implicit none

    integer::j,jesp,iter,s,k
    integer ::nesp_eq!number of species at equilibirum    
    double precision::dq(N_aerosol),qext(N_aerosol),qextold(N_aerosol)
    double precision::qaero(N_aerosol),qgas(N_aerosol),qgasold(N_aerosol)
    double precision:: Kelvin_effect(N_size,N_aerosol)
    double precision::ce_kernal_coef_tot(N_aerosol)
    double precision::ce_kernal_coef_totho(N_aerosol),ce_kernal_coef_ho(N_size,N_aerosol)
    double precision::ce_kernal_coef_tothi(N_aerosol),ce_kernal_coef_hi(N_size,N_aerosol)
    double precision organion, watorg, proton, lwc, d_ms
    double precision rhop_tmp,emw_tmp,wet_diam,aatoteq,qgasa,qgasi
    integer :: eq_species(nesp_eq)! species compute with equilibrium
    integer :: eq_species2(nesp_eq+nesp_isorropia)
    double precision:: total_ms(N_aerosol)
    double precision:: delta_t, produced, tot_mass,amm_sulfate,amm_to_be_redist,dqamm

    !**** SOAP ****
    double precision :: liquid(12), ionic, other(6)
    !!    DOUBLE PRECISION :: q(N_size*(1+N_aerosol)+N_aerosol)
    double precision qaerona,qaerocl,qgascl

    ! qaero contains the concentrations of inert, inorganic and organic species
    ! if i_hydrophilic ==1, the organic concentrations are hydrophobic
    ! the hydrophilic ones are then stored in qaq
    double precision qaq(N_aerosol),qextaq(N_aerosol),dqaq(N_aerosol)
    double precision qextaqold(N_aerosol)
    double precision orgmass_ho(N_size), orgmass_hi(N_size) 
    double precision coef_ho,coef_hi

    !!     ******zero init
    do jesp=1,(N_aerosol-1)
       dq(jesp)=0.D0
       qgas(jesp)=0.d0
       qgasold(jesp)=0.d0 
       qext(jesp)=0.D0
       qextold(jesp)=0.D0
       qaero(jesp)=0.d0
       qextaq(jesp)=0.D0
       qextaqold(jesp)=0.D0
       qaq(jesp)=0.d0
       ce_kernal_coef_totho(jesp)=0.D0
       ce_kernal_coef_tothi(jesp)=0.D0
       total_ms(jesp)=0.d0
    end do
    Kelvin_effect=1.0

    !     ****** allocate equilibrium species list      

    do s=1,nesp_eq
       eq_species(s)=aec_species(s)
    enddo

    if (soap_inorg==1) then
       do s=1,nesp_isorropia
          eq_species2(s)=isorropia_species(s)
       enddo
    endif

    do s=1,(N_aerosol-1)
       jesp=index_species(s,1) !This index is used to differentiate hydrophobic from hydrophilic part
       do j=1,N_size
          if(concentration_index(j, 1) <= ICUT_org) then
             qaero(s)=qaero(s)+concentration_mass(j,jesp)
          endif
       enddo
       if (inon_volatile(s).EQ.0 .and. aerosol_species_interact(s).GT.0) then
          qgas(s)=concentration_gas(s)
       else
          qgas(s) = 0.d0
       endif
       total_ms(s)=qaero(s)+qgas(s)
       qextold(s)=qaero(s)
    end do

    if(i_hydrophilic == 1) then ! Initialise the aqueous part of organics
       do s=1,nesp_eq
          jesp=index_species(eq_species(s),2) !This index is used to differentiate hydrophobic from hydrophilic part
          do j=1,N_size
             qaq(eq_species(s))=qaq(eq_species(s))+concentration_mass(j,jesp)
          enddo
          total_ms(eq_species(s))=total_ms(eq_species(s))+qaq(eq_species(s))
          qextaqold(eq_species(s))=qaq(eq_species(s))
       end do
    endif

    do jesp=1,(N_aerosol-1)
       if (aerosol_species_interact(jesp).GT.0) then
          emw_tmp = molecular_weight_aer(jesp) * 1.D-6 ! g/mol
          do j = 1,N_size
             wet_diam=wet_diameter(j)
             rhop_tmp = rho_wet_cell (j) * 1.d9 !1400 ! kg/m3
             if (with_kelvin_effect == 1) then
                call ssh_COMPUTE_KELVIN_COEFFICIENT(&
                     Temperature,&          ! temperature (Kelvin)
                     emw_tmp,&       ! ext mol weight (g.mol-1)
                     surface_tension(jesp),&   ! surface tension (N.m-1) from INC
                     wet_diam,&         ! wet aero diameter (µm)
                     rhop_tmp,&      ! aerosol density (kg.m-3)
                     Kelvin_effect(j,jesp) )   ! kelvin effect coef (adim)
             endif
             call ssh_COMPUTE_CONDENSATION_TRANSFER_RATE(&
                  diffusion_coef(jesp), &! diffusion coef (m2.s-1)
                  quadratic_speed(jesp),& ! quadratic mean speed (m.s-1)
                  accomodation_coefficient(jesp),& ! accomadation coef (adim)
                  wet_diameter(j),   & ! wet aero diameter (Âµm)
                  ce_kernal_coef(j,jesp) ) ! c/e kernel coef (m3.s-1)
          enddo
       endif
    enddo

    if (soap_inorg==1) then
       jesp=isorropia_species(2)
       qaero(jesp)=0.d0
       do j=1,N_size
          if(concentration_index(j, 1) <= ICUT_org) then
             qaero(jesp)=qaero(jesp)+concentration_mass(j,jesp)
          endif
       enddo
       qgasi = concentration_gas(jesp)
       if (inon_volatile(jesp).EQ.0) then
          !compute apparent gas concentration of sulfate
          !     i.e. the SO4 gas conc actually
          !     seen by equilibrium aerosols
          aatoteq = 0.d0
          ce_kernal_coef_tot(jesp) = 0.d0
          do j =1, N_size
             if(concentration_index(j, 1) <= ICUT_org) then
                ce_kernal_coef_tot(jesp)= ce_kernal_coef_tot(jesp)&
                     +ce_kernal_coef(j,jesp)*concentration_number(j)
                aatoteq=aatoteq+ce_kernal_coef(j,jesp)*concentration_number(j)
             endif
          enddo
          if(ce_kernal_coef_tot(jesp).gt.0.d0) then ! gas concentration to be condensed
             qgas(jesp)=concentration_gas(jesp)*aatoteq/ce_kernal_coef_tot(jesp)
          endif
          qgasa=qgas(jesp)
       else
          qgas(jesp)=0.d0
          qgasa = 0.d0
       endif
    endif

    jesp=isorropia_species(2)  ! case of sulfate
    qgas(jesp)=0.d0
    qgas(EH2O)=0.0
    if (soap_inorg==1) then
       qaero(EH2O)=0.d0
       do j=1,N_size
          if(concentration_index(j, 1) <= ICUT_org) then
             qaero(EH2O)=qaero(EH2O)+concentration_mass(j,EH2O)
          endif
       enddo
    endif
    organion = 0.D0
    watorg = 0.D0
    proton = 0.D0
    ionic = 0.D0
    lwc = 0.D0

    if (NACL_IN_THERMODYNAMICS==0) then
       qaerona = qaero(ENa)
       qaerocl = qaero(ECl)
       qgascl = qgas(ECl)
       qaero(ENa) = 0.D0
       qaero(ECl) = 0.D0
       qgas(ECl) = 0.D0
    endif

    if (soap_inorg==0) then
       call ssh_isoropia_drv(N_aerosol,&
            qaero,qgas,organion, watorg, ionic, proton, lwc, Relative_Humidity, Temperature, &
            liquid,other)
       do s=1,nesp_isorropia
          jesp=isorropia_species(s)
          qextold(jesp)=qaero(jesp)
          qgasold(jesp)=qgas(jesp)
       enddo
    endif

    if (ISOAPDYN.eq.0.or.soap_inorg==1) then
       call ssh_soap_eq(watorg, lwc, Relative_Humidity, ionic, proton, &
            Temperature, qaero, qgas, liquid, delta_t, qaq)
    endif

    ! compute total organic mass and partitioning coefficients 
    !if (with_kelvin_effect == 1) then
    do j=1,N_size
       orgmass_ho(j) = 0.d0
       orgmass_hi(j) = 0.d0
       do s=1,nesp_eq
          jesp=index_species(eq_species(s),1) !This index is used to differentiate hydrophobic from hydrophilic part
          orgmass_ho(j) = orgmass_ho(j) + concentration_mass(j,jesp)
       enddo
       if(i_hydrophilic == 1) then ! Initialise the aqueous part of organics
          do s=1,nesp_eq
             jesp=index_species(eq_species(s),2) !This index is used to differentiate hydrophobic from hydrophilic part
             orgmass_hi(j) = orgmass_hi(j) + concentration_mass(j,jesp)
          enddo
       endif
       if (orgmass_ho(j)==0.d0) then
          orgmass_ho(j)=1.d-10
       endif
       if (orgmass_hi(j)==0.d0) then
          orgmass_hi(j)=1.d-10
       endif
    enddo
    !else
    !   if (aerosol_hydrophobic(s)=1)
    !endif

    !compute only c/e coefficients
    do jesp=1,(N_aerosol-1)
       if (aerosol_species_interact(jesp).GT.0) then
          ce_kernal_coef_tot(jesp)=0.d0
          !if (with_kelvin_effect == 1) then
          do j = 1,N_size
             coef_ho =  (Kelvin_effect(j,jesp) - 0.999999999d0) * orgmass_ho(j)
             if ((coef_ho.GT.0D0).AND.(inon_volatile(jesp).EQ.0)) then
                ce_kernal_coef_ho(j,jesp) = ce_kernal_coef(j,jesp) * coef_ho
             else if (inon_volatile(jesp).EQ.1) then
                ce_kernal_coef_ho(j,jesp) = ce_kernal_coef(j,jesp)
             else
                ce_kernal_coef_ho(j,jesp) = 0.D0
             endif
             ce_kernal_coef_totho(jesp)= ce_kernal_coef_totho(jesp)&! compute total ce_kernal_coef coef
                  +ce_kernal_coef_ho(j,jesp)*concentration_number(j)
             ce_kernal_coef_tot(jesp)= ce_kernal_coef_tot(jesp)&! compute total ce_kernal_coef coef
                  +ce_kernal_coef(j,jesp)*concentration_number(j) !! for soap_inorg == 1
             if(i_hydrophilic == 1) then
                coef_hi =  (Kelvin_effect(j,jesp) - 0.999999999d0) * orgmass_hi(j)
                if ((coef_hi.GT.0D0).AND.(inon_volatile(jesp).EQ.0)) then
                   ce_kernal_coef_hi(j,jesp) = ce_kernal_coef(j,jesp) * coef_hi
                else if (inon_volatile(jesp).EQ.1) then
                   ce_kernal_coef_hi(j,jesp) = ce_kernal_coef(j,jesp)
                else
                   ce_kernal_coef_hi(j,jesp) = 0.D0
                endif
                ce_kernal_coef_tothi(jesp)= ce_kernal_coef_tothi(jesp)&! compute total ce_kernal_coef coef
                     +ce_kernal_coef_hi(j,jesp)*concentration_number(j) 
             endif
          enddo
          !else
          !    do j = 1,N_size
          !        ce_kernal_coef_tot(jesp)= ce_kernal_coef_tot(jesp)&! compute total ce_kernal_coef coef
          !             +ce_kernal_coef(j,jesp)*concentration_number(j)!&
          !     enddo
          !endif
       endif
    enddo

!!!!!! Reaction intra-particle  !!!!!!!!!!!!!!!!!!
    do s=1,nesp_isorropia
       jesp=isorropia_species(s)
       if (soap_inorg==0.or.inon_volatile(jesp).EQ.1) then
          produced=qaero(jesp)-qextold(jesp)+qgas(jesp)-qgasold(jesp)
          tot_mass=0.
          if(aerosol_species_interact(jesp).GT.0) then
             if (inon_volatile(jesp).EQ.0) then
                tot_mass=concentration_gas(jesp) !+qgas(jesp)-qgasold(jesp),0.)
             endif
          endif
          do j=1,N_size
             tot_mass=tot_mass+concentration_mass(j,jesp)
          enddo
          produced=max(produced,-tot_mass)

          if (tot_mass.ne.0) then
             if(aerosol_species_interact(jesp).GT.0) then
                if (inon_volatile(jesp).EQ.0) then
                   concentration_gas(jesp)=concentration_gas(jesp)*(1.d0+produced/tot_mass) !qgas(jesp)-qgasold(jesp),0.)
                endif
             endif
             do j=1,N_size
                concentration_mass(j,jesp)=concentration_mass(j,jesp)*(1.d0+produced/tot_mass) !*concentration_mass(j,jesp)
             enddo
          endif
       endif
    enddo
!!!!!! End reaction intra-particle  !!!!!!!!!!!!!!!!!!

    if (NACL_IN_THERMODYNAMICS==0) then
       qaero(ENa) = qaerona
       qaero(ECl) = qaerocl
       qgas(ECl) = qgascl
    endif

    qgas(EH2O)=0.0

    !for organic
    do s=1,nesp_eq
       jesp=eq_species(s)
       
       if(aerosol_species_interact(jesp).GT.0) then
          if (inon_volatile(jesp).EQ.0) &
               concentration_gas(jesp)=qgas(jesp)!new qgas is used in N_aerosol bin
          qext(jesp)=qaero(jesp)
          if(qext(jesp).gt.0.d0) then
             dq(jesp)=qext(jesp)-qextold(jesp)! compute delta aero conc
          else
             dq(jesp)=-qextold(jesp)
          endif
       endif
    enddo

    if (soap_inorg==1) then
       if (ECO3>0) then
          qext(ECO3)=qaero(ECO3)
          concentration_gas(ECO3)=0.
          if(qext(ECO3).gt.0.d0) then
             dq(ECO3)=qext(ECO3)-qextold(ECO3)! compute delta aero conc
          else
             dq(ECO3)=-qextold(ECO3)! compute delta aero conc
          endif
       endif

       !     ******redistribute on each cell according to Rates
       ! Ammonium that neutralises sulfate initially
       amm_sulfate = min((2.0*qextold(eq_species2(2))/96. - qextold(eq_species2(1))/22.989769)*17., qextold(eq_species2(3)))
       ! Additional ammonium necessary to neutralise sulfate
       amm_to_be_redist = max(17.* (2.*qaero(eq_species2(2))/96.0 - qaero(eq_species2(1))/22.989769) - amm_sulfate,0.0)
       amm_to_be_redist = min(concentration_gas(eq_species2(3))-qgas(eq_species2(3)),amm_to_be_redist)
       do s=1, nesp_isorropia
          jesp=isorropia_species(s)
          if(aerosol_species_interact(jesp).GT.0) then      
             qext(jesp)=qaero(jesp)
             if(jesp.ne.ECl.or.NACL_IN_THERMODYNAMICS==1) then
                concentration_gas(jesp)=qgas(jesp)!new qgas is used in N_aerosol bin
                if(qext(jesp).gt.0.d0) then
                   dq(jesp)=qext(jesp)-qextold(jesp)! compute delta aero conc
                else
                   dq(jesp)=-qextold(jesp)! compute delta aero conc
                endif
             endif
          endif
       enddo
       !print*,amm_to_be_redist,dq(eq_species2(3))
       !dqamm = min(amm_to_be_redist,dq(eq_species2(3)))  !Ammonium to be redistributed according to mass
       !dq(eq_species2(3)) = dq(eq_species2(3)) - dqamm   ! Ammonium to be redistributed according to cond. dynamics
       !print*,dqamm,dq(eq_species2(3))
   
       ! give back initial SO4 gas conc
       ! minus that consumed by equi bins
       jesp=isorropia_species(2)
       concentration_gas(jesp)=qgasi-qgasa
    endif

    if(i_hydrophilic == 1) then
       do s=1,nesp_eq
          jesp = eq_species(s)
          qextaq(jesp)=qaq(jesp)
          if(qextaq(jesp).gt.0.d0) then
             dqaq(jesp)=qextaq(jesp)-qextaqold(jesp)! compute delta aero conc
          else
             dqaq(jesp)=-qextaqold(jesp)
          endif
       enddo
    endif

    !     ******redistribute on each cell according to Rates
    if (soap_inorg==0) then
       call ssh_bulkequi_redistribution(concentration_number,concentration_mass,&
            nesp_eq,eq_species,N_size,dq,ce_kernal_coef_ho,ce_kernal_coef_hi,ce_kernal_coef_totho,&
            ce_kernal_coef_tothi,i_hydrophilic,dqaq,Kelvin_effect)
    else       
       call ssh_bulkequi_redistribution(concentration_number,concentration_mass,&
            nesp_eq,eq_species,ICUT_org,dq,ce_kernal_coef_ho,ce_kernal_coef_hi,ce_kernal_coef_totho,&
            ce_kernal_coef_tothi,i_hydrophilic,dqaq,Kelvin_effect)
       call ssh_bulkequi_redistribution(concentration_number,concentration_mass,&
            nesp_isorropia,eq_species2,ICUT_org,dq,ce_kernal_coef,ce_kernal_coef,ce_kernal_coef_tot,&
            ce_kernal_coef_tot,0,dqaq,Kelvin_effect)
       !call ssh_redistribution_amm(dqamm,ICUT_org)

       if (ECO3>0.and.NACL_IN_THERMODYNAMICS==1) then
          call ssh_redistribution_co3(qaero(ECO3),ICUT_org)
       endif

    endif

  end subroutine ssh_bulkequi_org

  subroutine ssh_bulkequi_inorg(nesp_eq, lwc, ionic, proton, liquid)
    !------------------------------------------------------------------------
    !
    !     -- DESCRIPTION
    !     This subroutine solves gas/aerosol bulk equilibrium for inorganic species
    !     based on ISORROPIA method.
    !     equilibrium will be established between all size bin < ICUT and bulk gas
    !     mass flux will be redistributed into each cell based on rates
    !
    !------------------------------------------------------------------------
    !
    !     -- INPUT VARIABLES
    !------------------------------------------------------------------------   
    implicit none
    integer::j,jesp,s
    integer ::nesp_eq!number of species at equilibirum    
    double precision::dq(N_aerosol),qext(N_aerosol),qextold(N_aerosol)
    double precision::qaero(N_aerosol),qgas(N_aerosol)
    double precision::ce_kernal_coef_tot(N_aerosol)
    double precision::aatoteq,qgasa,qgasi
    double precision:: Kelvin_effect(N_size,N_aerosol)
    double precision organion, watorg, proton, lwc, d_ms
    double precision rhop_tmp,emw_tmp,wet_diam    
    integer :: eq_species(nesp_eq)! species compute with equilibrium
    double precision:: total_ms(N_aerosol)
    double precision :: liquid(12), ionic, other(6)
    double precision :: daq(N_aerosol),amm_to_be_redist,dqamm
    double precision :: inorg_mass(N_size), amm_sulfate

    !!     ******zero init
    do jesp=1, N_aerosol !nesp_isorropia
       dq(jesp)=0.D0
       qgas(jesp)=0.d0
       qext(jesp)=0.D0
       qextold(jesp)=0.D0
       qaero(jesp)=0.d0
       ce_kernal_coef_tot(jesp)=0.D0
       total_ms(jesp)=0.d0
    end do
    Kelvin_effect=1.0
    !     ****** if sulfate computed dynamically avoid it
    !     ****** allocate equilibrium species list
    do jesp=1,nesp_isorropia
       eq_species(jesp)=isorropia_species(jesp)
    enddo

    do j=1,N_size
       inorg_mass(j) = 0.d0
       do s=1, nesp_isorropia
          jesp=isorropia_species(s)
          inorg_mass(j) = inorg_mass(j) + concentration_mass(j,jesp)
       enddo
    enddo
    !compute only c/e coefficients

    ce_kernal_coef=0.d0

    do s=1, nesp_isorropia
       jesp=isorropia_species(s)
       ce_kernal_coef_tot(jesp) = 0.d0

       if (aerosol_species_interact(jesp).GT.0) then
          ! compute total ce_kernal_coef coef (s-1)	
          IF (jesp.NE.ECl.or.NACL_IN_THERMODYNAMICS==1) THEN
             do j =1, N_size
                if(concentration_index(j, 1) <= ICUT) then
                   rhop_tmp = rho_wet_cell (j) * 1.d9 !1400 ! kg/m3
                   emw_tmp = molecular_weight_aer(jesp) * 1.D-6 ! g/mol

                   wet_diam=wet_diameter(j)

                   if (with_kelvin_effect == 1)  then
                      call ssh_COMPUTE_KELVIN_COEFFICIENT(&
                           Temperature,&          ! temperature (Kelvin)
                           emw_tmp,&       ! ext mol weight (g.mol-1)
                           surface_tension(jesp),&   ! surface tension (N.m-1) from INC
                           wet_diam,&         ! wet aero diameter (µm)
                           rhop_tmp,&      ! aerosol density (kg.m-3)
                           Kelvin_effect(j,jesp) )   ! kelvin effect coef (adim)
                   else
                      Kelvin_effect(j,jesp) = 1.00d0
                   endif
                   call ssh_COMPUTE_CONDENSATION_TRANSFER_RATE(&
                        diffusion_coef(jesp), &! diffusion coef (m2.s-1)
                        quadratic_speed(jesp),& ! quadratic mean speed (m.s-1)
                        accomodation_coefficient(jesp),& ! accomadation coef (adim)
                        wet_diameter(j),   & ! wet aero diameter (Âµm)
                        ce_kernal_coef(j,jesp) ) ! c/e kernel coef (m3.s-1)
                   if (with_kelvin_effect == 1) then
                      ce_kernal_coef(j,jesp) = ce_kernal_coef(j,jesp) * (Kelvin_effect(j,jesp)-0.999999999d0)*inorg_mass(j)
                      ce_kernal_coef_tot(jesp)= ce_kernal_coef_tot(jesp)&! compute total ce_kernal_coef coef
                           +ce_kernal_coef(j,jesp)*concentration_number(j)                
                   else
                      ce_kernal_coef_tot(jesp)= ce_kernal_coef_tot(jesp)&! compute total ce_kernal_coef coef
                           +ce_kernal_coef(j,jesp)*concentration_number(j)!&
                   endif
                endif
             enddo
          ENDIF
          !	endif
       endif
    enddo

    !compute total mass of each species
    do s=1,nesp_isorropia!inorganic
       jesp=isorropia_species(s)
       do j =1, N_size
          if(concentration_index(j, 1) <= ICUT) then
             qaero(jesp)=qaero(jesp)+concentration_mass(j,jesp)
          endif
       enddo
       qgas(jesp)=concentration_gas(jesp)
       total_ms(jesp)=qaero(jesp)+qgas(jesp)
       qextold(jesp)=qaero(jesp)
    end do

    jesp = isorropia_species(2)
    qgasi = concentration_gas(jesp)
    if (inon_volatile(jesp).EQ.0) then
       !compute apparent gas concentration of sulfate
       !     i.e. the SO4 gas conc actually
       !     seen by equilibrium aerosols
       aatoteq = 0.d0
       ce_kernal_coef_tot(jesp) = 0.d0
       do j =1, N_size
          if(concentration_index(j, 1) <= ICUT) then
             ce_kernal_coef_tot(jesp)= ce_kernal_coef_tot(jesp)&
                  +ce_kernal_coef(j,jesp)*concentration_number(j)
             aatoteq=aatoteq+ce_kernal_coef(j,jesp)*concentration_number(j)
          endif
       enddo
       if(ce_kernal_coef_tot(jesp).gt.0.d0) then ! gas concentration to be condensed
          qgas(jesp)=concentration_gas(jesp)*aatoteq/ce_kernal_coef_tot(jesp)
       endif
       qgasa=qgas(jesp)
    else
       qgas(jesp)=0.d0
       qgasa = 0.d0
    endif

    if (NACL_IN_THERMODYNAMICS==0) then
       qaero(ENa) = 0.D0
       qaero(ECl) = 0.D0
       qgas(ECl) = 0.D0
    endif

    organion= 0.D0
    watorg = 0.D0
    proton= 0.D0
    ionic = 0.D0
    lwc= 0.0

    !eqilibirum for inorganic
    call ssh_isoropia_drv(N_aerosol,&
         qaero,qgas,organion, watorg, ionic, proton, lwc, Relative_Humidity, Temperature, &
         liquid,other)

    qaero(N_aerosol)=lwc
    qgas(EH2O)=0.0

    !     ******redistribute on each cell according to Rates
    ! Ammonium that neutralises sulfate initially
    amm_sulfate = min((2.0*qextold(eq_species(2))/96. - qextold(eq_species(1))/22.989769)*17., qextold(eq_species(3)))
    ! Additional ammonium necessary to neutralise sulfate
    amm_to_be_redist = max(17.* (2.*qaero(eq_species(2))/96.0 - qaero(eq_species(1))/22.989769) - amm_sulfate,0.0)
    amm_to_be_redist = min(concentration_gas(eq_species(3))-qgas(eq_species(3)),amm_to_be_redist)
    do s=1, nesp_isorropia
       jesp=eq_species(s)
       if(aerosol_species_interact(jesp).GT.0) then      
          qext(jesp)=qaero(jesp)
          if(jesp.ne.ECl.or.NACL_IN_THERMODYNAMICS==1) then
             concentration_gas(jesp)=qgas(jesp)!new qgas is used in N_aerosol bin
             if(qext(jesp).gt.0.d0) then
                dq(jesp)=qext(jesp)-qextold(jesp)! compute delta aero conc
             else
                dq(jesp)=-qextold(jesp)! compute delta aero conc
             endif
          endif
       endif
    enddo
    dqamm = min(amm_to_be_redist,dq(eq_species(3)))  !Ammonium to be redistributed according to mass
    dq(eq_species(3)) = dq(eq_species(3)) - dqamm   ! Ammonium to be redistributed according to cond. dynamics

    ! give back initial SO4 gas conc
    ! minus that consumed by equi bins
    jesp=isorropia_species(2)
    concentration_gas(jesp)=qgasi-qgasa

    call ssh_bulkequi_redistribution(concentration_number,concentration_mass,&
         nesp_isorropia,eq_species,ICUT,dq,ce_kernal_coef,ce_kernal_coef,ce_kernal_coef_tot, &
         ce_kernal_coef_tot,0,daq,Kelvin_effect)
    call ssh_redistribution_amm(dqamm,ICUT)

  end subroutine ssh_bulkequi_inorg



  subroutine ssh_bulkequi_redistribution(c_number,c_mass,nesp_eq,eq_species,end_bin,dq, &
       AAiho,AAihi,ce_kernal_coef_totho,ce_kernal_coef_tothi,i_hydrophilic_tmp,dqaq,ck_ef)
    !------------------------------------------------------------------------
    !
    !     -- DESCRIPTION
    !     This subroutine redistribute bulk mass variations into each bins
    !
    !------------------------------------------------------------------------
    !
    !     -- INPUT VARIABLES
    !
    !     c_number: aerosol number concentration(#/m^3)
    !     c_mass: aerosol mass concentration(µg/m^3)
    !     nesp_eq: number of species at equilibirum
    !     eq_species: the list of species pointer at equilibirum
    !     end_bin: marks the number of bins concerned during the redistribution
    !     dq: bulk mass variations
    !     AAi: c/e kernel coefficient          ([m3.s-1]).
    !     ce_kernal_coef_tot: sum of c/e kernel coefficient          ([m3.s-1]).
    !------------------------------------------------------------------------
    implicit none
    integer::j,s,jesp,end_bin,iclip,i_hydrophilic_tmp
    integer::nesp_eq
    integer::eq_species(nesp_eq)
    double precision::totaer,temp_mass,totaa
    double precision::ce_kernal_coef_totho(N_aerosol)
    double precision::ce_kernal_coef_tothi(N_aerosol)
    double precision::ck_ef(N_size,N_aerosol)
    double precision::dq(N_aerosol)
    double precision::AAiho(N_size,N_aerosol)
    double precision::AAihi(N_size,N_aerosol)
    double precision::frac(N_size,N_aerosol)
    double precision::c_number(N_size)
    double precision::c_mass(N_size,N_aerosol_layers)
    double precision :: temp_mass_aq,dqaq(N_aerosol)
    integer::k,jespmass,jespmass2,iclipaq

    do s=1, nesp_eq
       jesp=eq_species(s)!index of the species in the N_aerosol list
       jespmass = index_species(jesp,1) !index of the species in the N_aerosol_layer list
       if(i_hydrophilic_tmp == 1) then
          jespmass2 = index_species(jesp,2) !index of the species in the N_aerosol_layer list
       endif
       if (aerosol_species_interact(jesp).GT.0) then
          !if (inon_volatile(jesp).EQ.0) then ! Do not redistribute non-volatile species   !Due to chemistry, concentrations of nonvolatile compounds can change
             IF ((jesp.NE.ECl .and. jesp.ne.ENa).or.NACL_IN_THERMODYNAMICS==1) THEN
                iclip=0
                iclipaq=0

                do j =1, N_size
                   if(concentration_index(j, 1) <= end_bin) then
                      if(ce_kernal_coef_totho(jesp).gt.0.d0) then
                         frac(j,jesp)= AAiho(j,jesp)*c_number(j) /& 
                              ce_kernal_coef_totho(jesp) 
                         temp_mass=c_mass(j,jespmass)+dq(jesp)*frac(j,jesp)
                         if(i_hydrophilic_tmp == 1) then
                            frac(j,jesp)= AAihi(j,jesp)*c_number(j) /& 
                                 ce_kernal_coef_tothi(jesp) 
                            temp_mass_aq=c_mass(j,jespmass2)+dqaq(jesp)*frac(j,jesp)
                            if(temp_mass_aq.lt.0.d0) iclipaq=1!case of over evaporation
                         endif
                         if(temp_mass.lt.0.d0) iclip=1!case of over evaporation
                      endif
                   endif
                enddo

                if(iclip.eq.1) then !over evaporate
                   totaer=0.d0
                   do j =1, N_size
                      if(concentration_index(j, 1) <= end_bin) then
                         totaer=totaer+c_mass(j,jespmass)
                      endif
                   enddo
                   if(totaer.gt.0.d0) then
                      do j =1, N_size
                         if(concentration_index(j, 1) <= end_bin) then
                            frac(j,jesp)=c_mass(j,jespmass)/totaer
                            c_mass(j,jespmass)=c_mass(j,jespmass)+dq(jesp)*frac(j,jesp)
                         endif
                      enddo
                   endif
                else !normal case

                   do j =1, N_size
                      if(concentration_index(j, 1) <= end_bin) then
                         if(ce_kernal_coef_totho(jesp).gt.0.d0) then  
                            frac(j,jesp)= AAiho(j,jesp)*c_number(j)/ &
                                 ce_kernal_coef_totho(jesp) 
                            c_mass(j,jespmass)=c_mass(j,jespmass)+dq(jesp)*frac(j,jesp)
                         endif
                      endif
                   enddo
                endif

                if(i_hydrophilic_tmp==1) then
                   if(iclipaq.eq.1) then ! overevaporate the hydrophilic concentrations if exist
                      totaer=0.d0
                      do j =1, N_size
                         if(concentration_index(j, 1) <= end_bin) then
                            totaer=totaer + c_mass(j,jespmass2)
                         endif
                      enddo
                      if(totaer.gt.0.d0) then
                         do j =1, N_size
                            if(concentration_index(j, 1) <= end_bin) then
                               frac(j,jesp)=c_mass(j,jespmass2)/totaer
                               c_mass(j,jespmass2)=c_mass(j,jespmass2)+dqaq(jesp)*frac(j,jesp)
                            endif
                         enddo
                      endif
                   else !normal case
                      do j =1, N_size
                         if(concentration_index(j, 1) <= end_bin) then
                            if(ce_kernal_coef_tothi(jesp).gt.0.d0) then
                               frac(j,jesp)= AAihi(j,jesp)*c_number(j)/ce_kernal_coef_tothi(jesp)
                               c_mass(j,jespmass2)=c_mass(j,jespmass2)+dqaq(jesp)*frac(j,jesp)
                            endif
                         endif
                      enddo
                   endif
                endif
             ENDIF
          !endif
       endif
    enddo

  end subroutine ssh_bulkequi_redistribution

  subroutine ssh_redistribution_co3(co3,end_bin)
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
    double precision :: co3

    inorg_total = 0.D0
    inorg_bin = 0.D0

    do jesp=1,N_aerosol
       if (aerosol_hydrophilic(jesp)==1.and.jesp.ne.ECO3.and.jesp.ne.EH2O) then
          do js=1,N_size
             if(concentration_index(js, 1) <= end_bin) then
                inorg_total = inorg_total + concentration_mass(js, jesp)
                inorg_bin(js) = inorg_bin(js) + concentration_mass(js, jesp)
             endif
          enddo
       endif
    enddo

    do js=1,N_size
       if(concentration_index(js, 1) <= end_bin) then
          if (inorg_total .gt. 0.D0) then
             concentration_mass(js, index_species(ECO3,1)) = co3 * inorg_bin(js)/inorg_total
          else
             concentration_mass(js, index_species(ECO3,1)) = 0.d0
          end if
       endif
    enddo
  end subroutine ssh_redistribution_co3

  subroutine ssh_redistribution_amm(dqamm,end_bin)
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

    integer :: iredist,iso4,ina
    double precision :: inorg_total, inorg_bin(N_size)
    integer :: jesp, js,lay,end_bin
    double precision :: dqamm

    inorg_total = 0.D0
    inorg_bin = 0.D0

    !do jesp=1,N_aerosol
    do js=1,N_size
       if(concentration_index(js, 1) <= end_bin) then
          inorg_total = inorg_total + concentration_mass(js, ESO4)
          inorg_bin(js) = inorg_bin(js) + concentration_mass(js, ESO4)
          if (NACL_IN_THERMODYNAMICS==1) then
             inorg_total = inorg_total + concentration_mass(js, ENa)
             inorg_bin(js) = inorg_bin(js) + concentration_mass(js, ENa)
          endif
       endif
    enddo
    !enddo

    do js=1,N_size
       if(concentration_index(js, 1) <= end_bin) then
          if (inorg_total .gt. 0.D0) then
             concentration_mass(js, ENH4) = dqamm * inorg_bin(js)/inorg_total
          else
             concentration_mass(js, ENH4) = 0.d0
          end if
       endif
    enddo
  end subroutine ssh_redistribution_amm

end Module iBulkequibrium


