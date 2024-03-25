module mod_wall

  use aInitialization

  implicit none

contains

  subroutine ssh_init_wall_gas_loss(klossg, krevg)

!-----------------------------------------------------------------------
!                                                                       
!     -- DESCRIPTION                                                    
!                                                                       
!     This routine computes ----- .                                        
!                                                                       
!-----------------------------------------------------------------------

    implicit none


    ! Wall loss
    integer :: jsp, jsp2
    integer :: ns, ns_aer
    DOUBLE PRECISION, allocatable, intent(out) :: klossg(:), krevg(:)
    DOUBLE PRECISION, allocatable :: psat(:), tref(:), dhvap(:), dif_air(:)
    double precision, allocatable :: velocity(:)
    double precision :: psat_loc, dltemp
    double precision :: kploc, gamloc, accom
    double precision, allocatable :: Wmol(:)
    
    !=================
    !=== Wall loss ===
    !=================
     
    ns = n_gas
    ns_aer = n_aerosol

    ! allocate arrays
    allocate (klossg(ns_aer))
    allocate (psat(ns_aer))
    allocate (tref(ns_aer))
    allocate (dhvap(ns_aer))    
    allocate (krevg(ns_aer))
    allocate (dif_air(ns_aer))
    allocate (velocity(ns_aer))   
    allocate (wmol(ns))

    ! global variables to local
    psat = saturation_vapor_pressure
    tref = t_ref
    dhvap = enthalpy_vaporization
    dltemp = temperature
    dif_air = diffusion_coef
    wmol = molecular_weight
    velocity = quadratic_speed
    
    ! 
    IF (kwall_gas>0.d0) THEN
       klossg(:)=0.d0
       DO Jsp=1,ns_aer
          IF (aerosol_type(Jsp)==4 .and. psat(Jsp)>0) THEN
             Jsp2=aerosol_type(Jsp)
             if (Tref(Jsp)>0) then
                psat_loc=psat(Jsp) &
                     *exp(-1000.0*dhvap(Jsp)/8.314d0* &
                     (1.0/DLtemp-1.0/Tref(Jsp)))
             else 
                psat_loc=psat(Jsp) & !*
                     *exp(-1000.0*dhvap(Jsp)/8.314d0* &
                     (1.0/DLtemp-1.0/298.d0))
             endif
             kploc=760.d0*8.202d-5*DLtemp/(1.d6*Wmol(Jsp2)*psat_loc)
             gamloc=10.d0**(3.299d0)*kploc**0.6407d0
             krevg(Jsp)=gamloc/kploc/Cwall*200.d0/Wmol(Jsp2)

             klossg(Jsp)=kwall_gas
             if (eddy_turbulence>0.d0.and.Cwall>0.d0.and. &
                  surface_volume_ratio>0.d0) then
                accom=min(10.d0**(-2.744d0)*kploc**1.407d0,1.0d0)
                klossg(Jsp)=surface_volume_ratio/(pi/2*1.d0/ &
                     (eddy_turbulence*dif_air(Jsp))**0.5d0+ &
                     4.d0/(accom*velocity(Jsp)))
             endif
             krevg(Jsp)=krevg(Jsp)*klossg(Jsp)
             kwall_gas=maxval(klossg)
          ENDIF
       ENDDO
    ENDIF
    
    ! deallocate
    if (allocated(psat)) deallocate(psat)
    if (allocated(tref)) deallocate(tref)
    if (allocated(dhvap)) deallocate(dhvap)      
    if (allocated(dif_air)) deallocate(dif_air)
    if (allocated(velocity)) deallocate(velocity)      
    if (allocated(wmol)) deallocate(wmol)            
    
  end subroutine ssh_init_wall_gas_loss


  

  subroutine ssh_compute_zc_wall(klossg, ratio_gas, &
       krevg, zcwall, zctot, tstep, ind, dt_ratio, concz, &
       conczz, gam)

!-----------------------------------------------------------------------
!                                                                       
!     -- DESCRIPTION                                                    
!                                                                       
!     This routine computes the wall loss for gas-phase species.                                        
!                                                                       
!-----------------------------------------------------------------------

    implicit none
   
    integer :: ind, jsp, jsp2
    integer :: ns, ns_aer
    double precision :: tstep, dt_ratio, gam
    double precision :: dun, dzero
    DOUBLE PRECISION, allocatable :: psat(:)
    DOUBLE PRECISION, allocatable :: chpr0(:),chlo0(:)
    DOUBLE PRECISION, allocatable, intent(in) :: klossg(:), krevg(:)
    DOUBLE PRECISION, intent(in) :: ratio_gas(:), zctot(:)
    DOUBLE PRECISION, intent(inout) :: zcwall(:)    
    double precision, allocatable :: concz(:), conczz(:)
    
    PARAMETER (dun=1.d0) 
    PARAMETER (dzero=0.d0)
    
    ns = n_gas
    ns_aer = n_aerosol

    
    ! allocate arrays
    allocate (psat(ns_aer))
    allocate (chpr0(ns))
    allocate (chlo0(ns))

    ! global variables to local
    psat = saturation_vapor_pressure
    chlo0 = chem_loss
    chpr0 = chem_prod

    IF (kwall_gas>0.d0) THEN 
       do Jsp=1,ns_aer
          if (aerosol_species_interact(Jsp)>0) then
             if (aerosol_type(Jsp) == 4 .and. psat(Jsp)>0) then
                Jsp2=aerosol_species_interact(Jsp)
                if (ind .eq. 1) then
                   chlo0(Jsp2)=chlo0(Jsp2)+klossg(Jsp)*ratio_gas(Jsp2)
                   chpr0(Jsp2)=chpr0(Jsp2)+krevg(Jsp)*ZCwall(Jsp2)
                   ZCwall(Jsp2)=(concz(Jsp2)+tstep*klossg(Jsp) &
                        *ZCtot(Jsp2)*ratio_gas(Jsp2)) &
                        /(dun+tstep*krevg(Jsp))
                   
                else if (ind .eq. 2 .or. ind .eq. 3) then
                   chlo0(Jsp2)=chlo0(Jsp2)+klossg(Jsp)
                   chpr0(Jsp2)=chpr0(Jsp2)+krevg(Jsp)*ZCwall(Jsp2)
                   ZCwall(Jsp2)=(((dt_ratio+dun)*(dt_ratio+dun)*concz(Jsp2)- &
                        conczz(Jsp2))/(dt_ratio*dt_ratio+2.d0*dt_ratio) + &
                        gam*tstep*klossg(Jsp)*ZCtot(Jsp2) &
                        * ratio_gas(Jsp2))/(dun+ &
                        gam*tstep*krevg(Jsp))
                endif
                if(ZCwall(Jsp2)<dzero) ZCwall(Jsp2)=dzero
                !     print*,ZCwall(Jsp2),klossg(Jsp),concz(Jsp2),ZCtot(Jsp2),
                !     &                 kpart(Jsp),Cwall,klossg(Jsp)/(kpart(Jsp)*Cwall)
             endif
          endif
       enddo
    ENDIF


    chem_loss = chlo0
    chem_prod = chpr0
   
    ! deallocate
    if (allocated(psat)) deallocate(psat)
    if (allocated(chpr0)) deallocate(chpr0)
    if (allocated(chlo0)) deallocate(chlo0)
    
  end subroutine ssh_compute_zc_wall



  subroutine ssh_compute_wall_particle_loss(fixed_density_aer, IDENS, &
       nbin_aer, ns_aer, dlconc_aer, mass_density_aer, delta_t, &
       DLnumconc_aer, ncomp_aer)

!-----------------------------------------------------------------------
!                                                                       
!     -- DESCRIPTION                                                    
!                                                                       
!     This routine computes the wall loss for particles .                                        
!                                                                       
!-----------------------------------------------------------------------

    implicit none

    INCLUDE 'paraero.inc'
    
    integer :: jb, idens, nbin_aer, ns_aer, jsp, ncomp_aer
    INTEGER :: idx_bs(nbin_aer)
    double precision :: dltemp, dlpress
    double precision :: dp, rho_tmp, rhoa, fixed_density_aer
    double precision :: rho_dry(nbin_aer), vset(nbin_aer), DLnumconc_aer(nbin_aer)
    DOUBLE PRECISION :: DLconc_aer(nbin_aer,ns_aer)
    double precision :: mass_density_aer(ns_aer)
    double precision :: CC, dif_part, debye, kwall_particle
    double precision :: wloss, delta_t
    
    ! global variables to local
    dltemp = temperature
    dlpress = pressure

    !C     Aerosol density converted in microg / microm^3.
    RHOA = fixed_density_aer * 1.D-09
    
    !C     Compute aerosol density
    rho_dry = RHOA  
    IF (IDENS.EQ.1) THEN   ! for varying density
       DO Jb=1,nbin_aer
          CALL SSH_COMPUTE_DENSITY(nbin_aer,ns_aer, ns_aer, TINYM, &
               DLconc_aer, &
               mass_density_aer,Jb,rho_dry(Jb)) 
       ENDDO
    ENDIF

    
    if (eddy_turbulence>0.d0.and.kwp0>0.d0 &
         .and.radius_chamber>0.d0) then
       call ssh_COMPUTE_AIR_FREE_MEAN_PATH(DLtemp,DLpress, &
            air_free_mean_path,viscosity)
       air_free_mean_path=air_free_mean_path*1.0e-6


       !C     relations between bin idx and size idx
       DO Jb=1,nbin_aer
          idx_bs(Jb)=(Jb-1)/ncomp_aer+1
       ENDDO
       
       do Jb=1,nbin_aer
          if (wet_diameter(Jb)>0.d0) then
             DP=wet_diameter(Jb)*1.0e-6
          else
             DP=diam_bound(idx_bs(Jb))*diam_bound(idx_bs(Jb)+1)**0.5*1.0e-6
          endif
          rho_tmp=rho_dry(Jb)*1.d18
          call ssh_compute_CC(air_free_mean_path,DP,CC)         
          call ssh_compute_VSTOKES(DP,rho_tmp,CC,viscosity,vset(Jb))
          dif_part=1.39d-23*DLtemp/(3.d0*pi*viscosity*DP)
          call SSH_DEBYE_N(pi*vset(Jb) &
               /(2*(eddy_turbulence*dif_part)**0.5) &
               ,1,debye)
          kwall_particle=kwp0+vset(Jb)*3.d0/4.d0/8.314d0+debye* &
               6*(eddy_turbulence*dif_part)**0.5/pi/radius_chamber
          !print*,kwall_particle
          wloss=exp(-kwall_particle*delta_t)
          DLnumconc_aer(Jb)=DLnumconc_aer(Jb)*wloss
          DO Jsp = 1, Ns_aer
             DLconc_aer(Jb,Jsp)=DLconc_aer(Jb,Jsp)*wloss
          ENDDO
            
       enddo
    ELSE IF (kwall_particle>0.d0) THEN
       wloss=exp(-kwall_particle*delta_t)
       DO Jb = 1, nbin_aer
          DLnumconc_aer(Jb)=DLnumconc_aer(Jb)*wloss
          DO Jsp = 1, Ns_aer
             DLconc_aer(Jb,Jsp)=DLconc_aer(Jb,Jsp)*wloss
          ENDDO
       ENDDO
    ENDIF

  end subroutine ssh_compute_wall_particle_loss

end module mod_wall
