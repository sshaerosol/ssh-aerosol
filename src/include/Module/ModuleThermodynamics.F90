!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------

Module cThermodynamics
  use aInitialization
  implicit none
contains
  subroutine ssh_compute_nonlinear_density(watconc,so4,nh4,no3,na,cl,massol,soldensity, &
       wet_density,dry_density,massdry,masswet,dry_to_wet)

    double precision :: watconc,so4,nh4,no3,na,cl,wet_density,dry_density,dry_to_wet
    double precision :: availso4, availna, availnh4, availno3, availcl
    double precision :: summass, waq_hcl, waq_hno3, waq_h2so4, waq_nacl
    double precision :: waq_nh4cl, waq_nano3, waq_nh4no3, waq_nh4nh4so4, waq_na2so4
    double precision :: Tr, temp_value, temp_value2, temp_rho, w, massol, soldensity
    double precision :: Tr2, waq_nh3, waq_na, masswet, massdry, volwet, voldry

    waq_na2so4=0.d0
    availso4=so4/98.0d0
    availnh4=nh4/17.0d0
    availno3=no3/63.0d0
    availna=na/23.0d0
    availcl=cl/36.5d0

    if (availso4<0.5d0*availna) then
       waq_na2so4=142d0*availso4
    else
       waq_na2so4=142d0*0.5d0*availna
    endif
    availso4=availso4-waq_na2so4/142d0
    availna=availna-2d0*waq_na2so4/142d0

    waq_nh4nh4so4=0.d0
    if (availso4<0.5d0*availnh4) then
       waq_nh4nh4so4=132d0*availso4
    else
       waq_nh4nh4so4=132d0*0.5d0*availnh4
    endif
    availso4=availso4-waq_nh4nh4so4/132d0
    availnh4=availnh4-2d0*waq_nh4nh4so4/132d0

    waq_nano3=0.d0
    if (availno3<availna) then
       waq_nano3=availno3*85d0
    else
       waq_nano3=availna*85d0
    endif
    availna=availna-waq_nano3/85d0
    availno3=availno3-waq_nano3/85d0

    waq_nh4no3=0.d0
    if (availno3<availnh4) then
       waq_nh4no3=availno3*80d0
    else
       waq_nh4no3=availnh4*80d0
    endif
    availnh4=availnh4-waq_nh4no3/80d0
    availno3=availno3-waq_nh4no3/80d0

    waq_nacl=0.d0
    if (availcl<availna) then
       waq_nacl=availno3*58.5d0
    else
       waq_nacl=availna*58.5d0
    endif
    availna=availna-waq_nacl/58.5d0
    availcl=availcl-waq_nacl/58.5d0

    waq_nh4cl=0.d0
    if (availcl<availnh4) then
       waq_nh4cl=availcl*53.5d0
    else
       waq_nh4cl=availnh4*53.5d0
    endif
    availnh4=availnh4-waq_nh4cl/53.5d0
    availcl=availcl-waq_nh4cl/53.5d0

    waq_hcl=availcl*36.5d0
    waq_hno3=availno3*63d0
    waq_h2so4=availso4*98d0
    waq_nh3=availnh4*17d0
    waq_na=availna*23d0

    summass=watconc+waq_hcl+waq_na+waq_nh3+waq_hno3+waq_h2so4 &
         &      +waq_nacl+waq_nh4cl+waq_nano3+waq_nh4no3+waq_nh4nh4so4+waq_na2so4

    temp_value=0.d0
    temp_value2=0.d0
    Tr=(Temperature-273.15d0)/273.15d0
    Tr2=(293.d0-273.15d0)/273.d0
    if (summass-watconc>0.d0) then
       if (summass>1.0d-10) then

          !     H2SO4
          w=waq_h2So4/summass    ! wet
          temp_rho=98d0/(36.2249d0+5.46891d0*w+8.73468d0*w*w+34.4556d0*Tr-65.0688d0*Tr*w &
               &        +41.7174d0*Tr*w*w+13.6245d0*Tr*Tr-39.1920d0*Tr*Tr*w+25.2329d0*Tr*Tr*w*w)*1000d0
          temp_value=temp_value+w/temp_rho

          w=waq_h2So4/(summass-watconc) !dry
          temp_rho=98d0/(36.2249d0+5.46891d0*w+8.73468d0*w*w+34.4556d0*Tr2-65.0688d0*Tr2*w &
               &        +41.7174d0*Tr2*w*w+13.6245d0*Tr2*Tr2-39.1920d0*Tr2*Tr2*w+25.2329d0*Tr2*Tr2*w*w)*1000d0
          temp_value2=temp_value2+w/temp_rho

          !     (NH4)2SO4
          w=waq_nh4nh4so4/summass ! wet
          temp_rho=132d0/(55.8085d0+40.62d0*w-11.0797d0*w*w+13.5837d0*Tr-15.8075d0*Tr*w)*1000d0
          temp_value=temp_value+w/temp_rho

          w=waq_nh4nh4so4/(summass-watconc) ! dry
          temp_rho=132d0/(55.8085d0+40.62d0*w-11.0797d0*w*w+13.5837d0*Tr2-15.8075d0*Tr2*w)*1000d0
          temp_value2=temp_value2+w/temp_rho

          !     HNO3
          w=waq_hno3/summass     ! wet
          temp_rho=63d0/(28.6889d0-1.354d0*w+15.9166d0*w*w+16.9281d0*Tr-1.36496d0*Tr*w)*1000d0
          temp_value=temp_value+w/temp_rho

          w=waq_hno3/(summass-watconc) ! dry
          temp_rho=63d0/(28.6889d0-1.354d0*w+15.9166d0*w*w+16.9281d0*Tr2-1.36496d0*Tr2*w)*1000d0
          temp_value2=temp_value2+w/temp_rho

          !     NH4NO3
          w=waq_nh4no3/summass
          temp_rho=80d0/(50.3779d0-6.91069d0*w+11.3024d0*w*w-28.7432d0*Tr+155.419d0*Tr*w-157.049d0*Tr*w*w &
               &       +161.094d0*Tr*Tr-585.88d0*Tr*Tr*w+541.843d0*Tr*Tr*w*w)*1000d0
          temp_value=temp_value+w/temp_rho

          w=waq_nh4no3/(summass-watconc)
          temp_rho=80d0/(50.3779d0-6.91069d0*w+11.3024d0*w*w-28.7432d0*Tr2+155.419d0*Tr2*w-157.049d0*Tr2*w*w &
               &       +161.094d0*Tr2*Tr2-585.88d0*Tr2*Tr2*w+541.843d0*Tr2*Tr2*w*w)*1000d0
          temp_value2=temp_value2+w/temp_rho

          !     Na2SO4
          temp_value=temp_value+waq_na2so4/2700d0/summass
          temp_value2=temp_value2+waq_na2so4/2700d0/(summass-watconc)

          !     NaNO3
          temp_value=temp_value+waq_nano3/2260d0/summass
          temp_value2=temp_value2+waq_nano3/2260d0/(summass-watconc)

          !     NaCl
          temp_value=temp_value+waq_nacl/2165d0/summass
          temp_value2=temp_value2+waq_nacl/2165d0/(summass-watconc)

          !     NH4Cl
          temp_value=temp_value+waq_nh4cl/1530d0/summass
          temp_value2=temp_value2+waq_nh4cl/1530d0/(summass-watconc)

          !     HCl
          temp_value=temp_value+waq_hcl/1150d0/summass
          temp_value2=temp_value2+waq_hcl/1150d0/(summass-watconc)

          !     Na
          temp_value=temp_value+waq_na/970d0/summass
          temp_value2=temp_value2+waq_na/970d0/(summass-watconc)

          !     NH3
          temp_value=temp_value+waq_nh3/910d0/summass
          temp_value2=temp_value2+waq_nh3/910d0/(summass-watconc)

          !     H20
          temp_value=temp_value+watconc/1000d0/summass

          masswet=summass+massol
          massdry=masswet-watconc
          !print*,"la",massol,summass+massol,summass

          wet_density=1.d0/(massol/masswet*1.d0/soldensity+1.d9*temp_value*(1.d0-massol/masswet))
          dry_density=1.d0/(massol/massdry*1.d0/soldensity+1.d9*temp_value2*(1.d0-massol/massdry))
          !print*,"loc: ", wet_density,dry_density,soldensity,fixed_density,massol,masswet,massdry,temp_value,temp_value2
          !stop

          if (wet_density.le.900d-9) then
             wet_density=900d-9
          endif
       else
          masswet=summass+massol
          massdry=masswet-watconc

          wet_density=soldensity
          dry_density=soldensity
       endif
    else
       masswet=summass+massol
       massdry=masswet-watconc

       wet_density=soldensity
       dry_density=soldensity
    endif

    volwet=masswet/wet_density
    voldry=massdry/dry_density

    if (voldry.gt.0.) then
       dry_to_wet=(volwet/voldry)**(1d0/3d0)
    else
       dry_to_wet=1.d0
    endif

    if (dry_to_wet.le.0.8d0) then
       dry_to_wet=0.8d0
    endif
    !rint*,"ok stop",wet_density,dry_density
    !stop

  end subroutine ssh_compute_nonlinear_density

  subroutine ssh_get_nonlinear_density(qext,dry_density,wet_density,dry_mass,wmass,dry_to_wet)

    double precision :: qext(N_aerosol)
    double precision :: massol,soldensity,dry_density,wet_density,dry_mass,wmass,dry_to_wet
    integer :: is_inorg,jesp,jesp2,i

    massol=0.d0
    soldensity=0.d0
    do jesp=1,N_aerosol
       if (jesp.ne.EH2O.and.jesp.ne.ESO4.and.jesp.ne.ENH4.and.jesp.ne.ENO3.and.jesp.ne.ENa.and.jesp.ne.ECl) then
          massol=massol+qext(jesp)
          soldensity=soldensity+qext(jesp)/mass_density(jesp)
       endif
    enddo

    if (massol>0.d0) then       
       soldensity=1.d0/(soldensity/massol)       
    else
       soldensity=fixed_density
    endif

    !print*,"ici",massol,massol+qext(EH2O)+qext(ESO4)+qext(ENH4)+qext(ENO3)+qext(ENa)+qext(ECl)
    call ssh_compute_nonlinear_density(qext(EH2O),qext(ESO4),qext(ENH4),qext(ENO3),qext(ENa),qext(ECl), &
         massol,soldensity,wet_density,dry_density,dry_mass,wmass,dry_to_wet)
    !print*,wet_density

  end subroutine ssh_get_nonlinear_density

  subroutine ssh_compute_wet_mass_diameter(start_bin,end_bin,c_mass,c_number,c_inti, &
       wet_m,wet_d,wet_v)
    !------------------------------------------------------------------------
    !
    !     -- DESCRIPTION
    !     This subroutine computes particle wet diameter for each
    !     bin between start_bin and end_bin as well as their water content
    !
    !------------------------------------------------------------------------
    !
    !     -- INPUT VARIABLES
    !
    !     start_bin: index of the first bin need the computation
    !     end_bin: index of the last bin need the computation
    !     c_number: aerosol number concentration(#/m^3)
    !     c_mass: aerosol mass concentration(µg/m^3)
    !
    !     -- OUTPU VARIABLES
    !
    !     c_inti: aerosol internal species mass concentration(µg/m^3)
    !     wet_m: particle wet mass
    !     wet_d: particle wet diameter
    !     wet_v: particle wet volume
    !     dry_d: particle dry diameter
    !------------------------------------------------------------------------
    implicit none
    integer :: jesp,s,j,k,i,start_bin,end_bin
    double precision :: c_number(N_size)
    double precision :: c_mass(N_size,N_aerosol_layers)
    double precision :: c_inti(N_size,N_inside_aer)
    double precision :: wet_m(N_size)!wet mass
    double precision :: wet_d(N_size)!wet diameter
    double precision :: wet_v(N_size)!wet volume
    double precision :: dry_d(N_size)!dry diameter
    double precision :: qext(N_aerosol)
    double precision :: qinti(N_inside_aer)
    double precision :: lwc,qti,vad
    double precision :: aero(5)
    double precision :: rhoaer!aerosol density based on the internal composition
    double precision :: rhoaer_dry
    double precision :: dry_density,wet_density,dry_mass,wmass,dry_to_wet

    rhoaer = fixed_density

    do j=start_bin,end_bin
       !initialization
       do jesp=1,N_aerosol
          qext(jesp) = 0.d0
       enddo
       do jesp=1,N_aerosol_layers
          s = List_species(jesp)
          qext(s) = qext(s) + c_mass(j,jesp)
       enddo
       do jesp=1,N_inside_aer
          qinti(jesp)=0.d0
       enddo
       !     ******total dry mass
       qti=0.D0
       do s=1,(N_aerosol-1)
          qti=qti+qext(s)     !total dry mass µg.m-3
       end do
       ! The threshold 0.0 for the minimum aerosol concentration 
       ! should be avoided because it leads to a too small particle diameter.
       if ((qti.gt.TINYM).AND.(c_number(j).gt.TINYN)) then ! No water  initially - compute it
          do i=1,nesp_isorropia
             jesp=isorropia_species(i)
             aero(i)=qext(jesp)
          enddo
          call ssh_calculatewater(aero,qinti,lwc,j)             
          qext(EH2O)=lwc
          do jesp=1,N_inside_aer
             c_inti(j,jesp)=qinti(jesp)
          enddo

          if (with_fixed_density == 2) then                
             call ssh_get_nonlinear_density(qext,dry_density,wet_density,dry_mass,wmass,dry_to_wet)                

             vad=dry_mass/dry_density !qti/rhoaer!qti total dry mass            
             rho_wet_cell(j)=wet_density
             rhoaer=rho_wet_cell(j)
             wet_v(j)=wmass/wet_density !vad+qext(EH2O)/rhoaer!: wet volume aerosol concentration (µm3/m3).                            
             dry_d(j)=(vad/c_number(j)/cst_pi6)**cst_FRAC3 ! dry aerosol dimaeter µm                
             wet_d(j)=dry_to_wet*dry_d(j) ! (wet_v(j)/c_number(j)/cst_pi6)**cst_FRAC3 ! wet aerosol diameter µm
             wet_m(j)=wmass/c_number(j) !(qti+qext(EH2O))/c_number(j)
             c_mass(j,EH2O_layers)=qext(EH2O)
          else
             if(rhoaer.le.0.d0) then
                rhoaer=density_aer_bin(j)
             endif

             if(rhoaer.gt.0.d0) then
                vad=qti/rhoaer!qti total dry mass
                wet_v(j)=vad+qext(EH2O)/rhoaer!: wet volume aerosol concentration (µm3/m3).
                ! aerosol diameter
                ! qn cannot be zero, checked in eqpart routine                   
                dry_d(j)=(vad/c_number(j)/cst_pi6)**cst_FRAC3 ! dry aerosol dimaeter µm
                k=concentration_index(j, 1)
                wet_d(j)=(wet_v(j)/c_number(j)/cst_pi6)**cst_FRAC3 ! wet aerosol diameter µm
                wet_d(j)=DMAX1(wet_d(j),dry_d(j))!wet diameter is always larger than dry diameter
                c_mass(j,EH2O_layers)=qext(EH2O)
                wet_m(j)=(qti+qext(EH2O))/c_number(j) ! single wet mass (µg)                   
             else
                k=concentration_index(j, 1)
                wet_d(j)=size_diam_av(k)
                wet_m(j)=size_mass_av(k)
             endif
          endif
          !if (j==50) then
          !   print*,rhoaer,wet_v(j),wet_m(j),dry_d(j),wet_d(j)
          !   stop
          !endif
       else
          c_number(j) = 0.d0
          do jesp=1,N_aerosol_layers
             c_mass(j,s) = 0.d0
          enddo
          k=concentration_index(j, 1)
          wet_d(j)=size_diam_av(k)
          wet_m(j)=size_mass_av(k)
       endif
    enddo

  end subroutine ssh_compute_wet_mass_diameter

  subroutine ssh_calculatewater(aero,qinti,lwc,jbin)
    !------------------------------------------------------------------------
    !
    !     -- DESCRIPTION
    !     This subroutine computes particle water content
    !
    !------------------------------------------------------------------------
    !
    !     -- INPUT VARIABLES
    !
    !     aero: aerosol mass concentration(µg/m^3)
    !
    !     -- OUTPU VARIABLES
    !
    !     qinti: aerosol internal species mass concentration(µg/m^3)
    !     lwc: particle water content
    !------------------------------------------------------------------------
    implicit none
    integer jesp,j
    double precision aero(5)
    double precision organion, watorg, proton
    double precision lwc,ionic

    double precision wi(5),w(5),gas(3),cntrl(2), other(6)
    double precision liquid(N_liquid),solid(N_solid)
    double precision qinti(N_inside_aer)

    integer jbin

    organion = 0.D0
    watorg = 0.D0
    proton = 0.D0

    cntrl(1) = 1.D0!reverse mode
    cntrl(2) = 1.D0

    gas(1)=0.d0
    gas(2)=0.d0
    gas(3)=0.d0

    !     conversion unit for isorropia needed in mol.m-3
    do j=1,nesp_isorropia
       jesp=isorropia_species(j)
       if(molecular_weight_aer(jesp).gt.0.d0) then
	  wi(j)= aero(j)& ! µg.m-3
               /molecular_weight_aer(jesp)!&  ! µg.mol-1
       endif
    end do

    !     call isorropia fortran routine
    if (iter_water(jbin)==0) then
       call SSH_ISOROPIA(wi, Relative_Humidity, Temperature, cntrl, w, gas,&
            liquid, solid, other, organion, watorg)
       !     clipping to tinym

       !     Aqueous phase total liquid water content and pH (proton) concentration
       ionic = other(5)
       proton = liquid(IH) * imw(IH) !* gammaH  ! microg.m-3 but equivalent to micromol.m-3

       if (gas(1).lt.0.d0) gas(1)=tinym
       if (gas(2).lt.0.d0) gas(2)=tinym
       if (gas(3).lt.0.d0) gas(3)=tinym

       !     Aqueous phase total liquid water content and pH (proton) concentration
       do jesp=IH,IOH
          qinti(jesp)= DMAX1(liquid(jesp),0.D0)*imw(jesp)   ! moles to µg MOLAR WEIGHT
       end do
       ! solid inorg aerosol
       do jesp=SNaNO3,SLC
          qinti(jesp)= DMAX1(solid(jesp-12),0.D0)&
               *smw(jesp)        ! moles to µg !molecular_weight_solid(jesp)
       end do
       ! liquid water content
       lwc= qinti(IH2O)+qinti(IOH)*1.05882352941D0 ! mwh2o/mwioh
       if(lwc < 1.1d-12) lwc = 0.d0 !Minimum lwc is arbitrary fixed in ISORROPIA. Remove it.
       if (sum(wi)>0.d0) then
          ratio_water(jbin)=lwc/sum(wi)
          iter_water(jbin)=iter_water(jbin)+1
       endif
    else
       iter_water(jbin)=iter_water(jbin)+1
       lwc=ratio_water(jbin)*sum(wi)
       do jesp=SNaNO3,SLC
          qinti(jesp)=0.D0      !FCo: Warning if delisquescent is used one day      
       end do
    endif
    if (iter_water(jbin)==niter_water) iter_water(jbin)=0

  end subroutine ssh_calculatewater

  subroutine ssh_update_wet_diameter(end_bin,c_mass,c_inti,c_number,wet_m,&
       wet_d,wet_v,dry_d)
    !------------------------------------------------------------------------
    !
    !     -- DESCRIPTION
    !     This subroutine computes particle wet diameter for each
    !     bin between start_bin and end_bin based on known water content
    !     It also update the total water and pH of the system
    !
    !------------------------------------------------------------------------
    !
    !     -- INPUT VARIABLES
    !
    !     end_bin: index of the last bin need the computation
    !     c_number: aerosol number concentration(#/m^3)
    !     c_mass: aerosol mass concentration(µg/m^3)
    !     c_inti: aerosol internal species mass concentration(µg/m^3)
    !
    !     -- OUTPUT VARIABLES
    !
    !     wet_m: particle wet mass
    !     wet_d: particle wet diameter
    !     wet_v: particle wet volume
    !     dry_d: particle dry diameter
    !------------------------------------------------------------------------
    implicit none
    integer :: jesp,s,j,k,start_bin,end_bin
    double precision ::c_mass(N_size,N_aerosol_layers)
    double precision ::c_number(N_size)
    double precision :: c_inti(N_size,N_inside_aer)
    double precision :: wet_m(N_size)!wet mass
    double precision :: wet_d(N_size)!wet diameter
    double precision :: wet_v(N_size)!wet volume
    double precision :: dry_d(N_size)!dry diameter
    double precision :: qext(N_aerosol)
    double precision :: qinti(N_inside_aer)
    double precision :: qti,vad,rhoaer
    double precision :: rhoaer_dry ! aerosol density without water (YK)

    total_water=0.d0
    total_IH=0.d0
    Do j=1,N_size
       if(concentration_index(j, 1) <= end_bin) then
          rho_wet_cell(j)=fixed_density
          total_water=total_water+c_mass(j,EH2O)
          total_IH=total_IH+c_inti(j,IH)
          qti=0.D0
          do s=1,N_aerosol_layers-1
             qti=qti+c_mass(j,s)     ! µg.m-3
          end do
          do s=1,N_aerosol
             qext(s) = 0.d0
          enddo
          if (c_number(j).gt.TINYN.and.qti.gt.TINYM ) then
             do jesp=1,N_aerosol_layers
                s = List_species(jesp)
                qext(s) = qext(s) + c_mass(j,jesp)
             enddo
             do jesp=1,N_inside_aer
                qinti(jesp)=c_inti(j,jesp)
             enddo
             call ssh_VOLAERO(N_aerosol,qext,qinti,rhoaer,rhoaer_dry)

             if(rhoaer.gt.0.d0) then
                rho_wet_cell(j)=rhoaer
             else
                rhoaer=density_aer_bin(j)
             endif
             if(rhoaer.gt.0.d0) then
                ! vad=qti/rhoaer!qti total dry mass
                vad=qti/rhoaer_dry!qti total dry mass YK: aerosol density without water.
                wet_v(j)=vad+qext(EH2O)/rhoaer!: wet volume aerosol concentration (µm3/m3).
                ! aerosol diameter
                ! qn cannot be zero, checked in eqpart routine
                dry_d(j)=(vad/c_number(j)/cst_pi6)**cst_FRAC3 ! dry aerosol dimaeter µm
                wet_d(j)=(wet_v(j)/c_number(j)/cst_pi6)**cst_FRAC3 ! wet aerosol diameter µm
                wet_d(j)=DMAX1(wet_d(j),dry_d(j))!wet diameter is always larger than dry diameter
                wet_m(j)=(qti+qext(EH2O))/c_number(j) ! single wet mass (µg)
             endif
          else
             ! if too few aerosols or too few mass
             ! we set variables of given bins as
             ! its initial fixed ones,
             ! thus avoiding zero values
             k=concentration_index(j, 1)
             wet_d(j)=size_diam_av(k)
             wet_m(j)=size_mass_av(k)
          endif
       endif
    enddo

  end subroutine ssh_update_wet_diameter

  subroutine ssh_update_wet_diameter_liquid(end_bin,c_mass,c_number,wet_m,&
       wet_d,wet_v,dry_d)

    !------------------------------------------------------------------------
    !
    !     -- DESCRIPTION
    !     This subroutine computes particle wet diameter for each
    !     bin between start_bin and end_bin based on known water content
    !     It also update the total water and pH of the system
    !
    !------------------------------------------------------------------------
    !
    !     -- INPUT VARIABLES
    !
    !     end_bin: index of the last bin need the computation
    !     c_number: aerosol number concentration(#/m^3)
    !     c_mass: aerosol mass concentration(µg/m^3)
    !     c_inti: aerosol internal species mass concentration(µg/m^3) from isoroppia
    !     lwcorg = organic liquid water content
    !
    !     -- OUTPUT VARIABLES
    !
    !     wet_m: particle wet mass
    !     wet_d: particle wet diameter
    !     wet_v: particle wet volume
    !     dry_d: particle dry diameter
    !------------------------------------------------------------------------
    implicit none
    integer :: jesp,s,j,k,end_bin

    double precision ::c_mass(N_size,N_aerosol_layers)
    double precision ::c_number(N_size)
    double precision :: c_inti(N_size,N_inside_aer)
    double precision :: wet_m(N_size)!wet mass
    double precision :: wet_d(N_size)!wet diameter
    double precision :: wet_v(N_size)!wet volume
    double precision :: dry_d(N_size)!dry diameter
    double precision :: qext(N_aerosol)
    double precision :: qinti(N_inside_aer)
    double precision :: qti,vad,rhoaer
    double precision :: dry_density,wet_density,dry_mass,wmass,dry_to_wet

    total_water=0.d0
    !total_IH=0.d0
    Do j=1,N_size
       if(concentration_index(j, 1) <= end_bin) then
          rho_wet_cell(j)=fixed_density
          qti=0.D0
          do s=1,N_aerosol_layers-1
             qti=qti+c_mass(j,s)     ! µg.m-3
          end do
          total_water=total_water+c_mass(j,N_aerosol_layers)

          if (c_number(j).gt.TINYN.AND.qti.gt.TINYM) then
             do jesp=1,N_aerosol
                qext(jesp)=0.d0
             enddo
             do jesp=1,N_aerosol_layers
                s = List_species(jesp)
                qext(s) = qext(s) + c_mass(j,jesp)
             enddo
             if (with_fixed_density == 2 .and. c_number(j)>0.d0) then
                call ssh_get_nonlinear_density(qext,dry_density,wet_density,dry_mass,wmass,dry_to_wet)             

                vad=dry_mass/dry_density !qti/rhoaer!qti total dry mass            
                rho_wet_cell(j)=wet_density

                wet_v(j)=wmass/wet_density !vad+qext(EH2O)/rhoaer!: wet volume aerosol concentration (µm3/m3).
                ! aerosol diameter
                ! qn cannot be zero, checked in eqpart routine
                dry_d(j)=(vad/c_number(j)/cst_pi6)**cst_FRAC3 ! dry aerosol dimaeter µm
                !k=concentration_index(j, 1)
                wet_d(j)=dry_to_wet*dry_d(j) ! (wet_v(j)/c_number(j)/cst_pi6)**cst_FRAC3 ! wet aerosol diameter µm
                ! wet_d(j)= !DMAX1(wet_d(j),dry_d(j))!wet diameter is always larger than dry diameter
                !c_mass(j,EH2O)=qext(EH2O)
                wet_m(j)=wmass/c_number(j) !(qti+qext(EH2O))/c_number(j) 
             else
                if (with_fixed_density == 0) then
                   call ssh_compute_density(N_size,N_aerosol_layers,EH2O_layers,TINYM,c_mass,&
                        mass_density_layers,j,rhoaer)
                else
                   rhoaer = fixed_density 
                endif
                if(rhoaer.gt.0.d0) then
                   rho_wet_cell(j)=rhoaer
                else
                   rhoaer=density_aer_bin(j)
                   rho_wet_cell(j)=rhoaer
                endif
                if(rhoaer.gt.0.d0) then
                   vad=qti/rhoaer!qti total dry mass
                   wet_v(j)=vad+qext(EH2O)/rhoaer!: wet volume aerosol concentration (µm3/m3).
                   ! aerosol diameter
                   ! qn cannot be zero, checked in eqpart routine
                   dry_d(j)=(vad/c_number(j)/cst_pi6)**cst_FRAC3 ! dry aerosol dimaeter µm
                   wet_d(j)=(wet_v(j)/c_number(j)/cst_pi6)**cst_FRAC3 ! wet aerosol diameter µm
                   wet_d(j)=DMAX1(wet_d(j),dry_d(j))!wet diameter is always larger than dry diameter
                   wet_m(j)=(qti+qext(EH2O))/c_number(j) ! single wet mass (µg)
                endif
             endif
          else
             ! if too few aerosols or too few mass
             ! we set variables of given bins as
             ! its initial fixed ones,
             ! thus avoiding zero values
             k=concentration_index(j, 1)
             wet_d(j)=size_diam_av(k)
             wet_m(j)=size_mass_av(k)
             dry_d(j)=size_diam_av(k)
             wet_v(j)=cst_pi6*(wet_d(j)**3.d0)
          endif
       else
          k=concentration_index(j, 1)
          wet_d(j)=size_diam_av(k)
          wet_m(j)=size_mass_av(k)
          dry_d(j)=size_diam_av(k)
          wet_v(j)=cst_pi6*(wet_d(j)**3.d0)
       endif
    enddo

  end subroutine ssh_update_wet_diameter_liquid

  subroutine ssh_VOLAERO(nesp_aer,qext,qinti,rhoaer, rhoaer_dry)
    !------------------------------------------------------------------------
    !
    !     -- DESCRIPTION
    !     This subroutine computes the dry and wet aerosol volumes and the
    !     aerosol density according to the internal composition.
    !
    !------------------------------------------------------------------------
    !
    !     -- INPUT VARIABLES
    !
    !     nesp_aer: number of aerosol species
    !     qext: aerosol mass concentration(µg/m^3)
    !     qinti: aerosol internal species mass concentration(µg/m^3)
    !
    !     -- OUTPU VARIABLES
    !
    !     rhoaer: aerosol density according to the internal composition.
    !------------------------------------------------------------------------
    implicit none
    integer nesp_aer
    double precision qext(nesp_aer),qinti(N_inside_aer),rhoaer

    integer jesp,s
    double precision vid,viw,vis,vil,vod,vad,vaw,sumint
    double precision rhoaer_dry, sum_dry !! YK

!!!     ******inorganic volume
    vid=0.D0
    vis=0.D0
    vil=0.D0
    viw=0.D0
    vod=0.D0
    sumint=0.d0
    sum_dry = 0.d0
    !print*,mass_density
    !!     Mineral Dust and BC
    do jesp=1,N_inert
      vis = vis + qext(jesp)/mass_density(jesp)
      sumint = sumint + qext(jesp)
    enddo

    !!     Inorganic
    ! compute solid aerosol volume
    do jesp=SNaNO3,SLC
       vis=vis+qinti(jesp)/SMD(jesp)
    end do
    ! compute liquid aerosol volume
    ! sodium volume
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
    vil=vil + qext(ENa)/mass_density(ENa)
#endif
    ! 
#ifndef WITHOUT_NACL_IN_THERMODYNAMICS
    vil=vil + qinti(INa)/mass_density(ENa)
#endif
    ! 
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
    vil = vil + qext(ECl)/mass_density(ECl) ! HCl volume
#endif
    ! ammonium volume
    vil=vil+( qinti(INH4)&
         +qinti(INH4)&
         *0.944444444444D0&    ! mwnh3/mwinh4
         )/mass_density(ENH4)          ! µg.µm-3

    ! nitric acid volume
    vil=vil+( qinti(IHNO3)&
         +qinti(INO3)&
         *1.01612903226D0&     ! mwhno3/mwino3
         )/mass_density(ENO3)          ! µg.µm-3

    ! chlorhydric acid volume
#ifndef WITHOUT_NACL_IN_THERMODYNAMICS
    vil=vil+( qinti(IHCl)&
         +qinti(ICl)&
         *1.02816901408D0&     ! mwhcl/mwicl
         )/mass_density(ECl)           ! µg.µm-3
#endif
    ! sulfuric acid volume
    vil=vil+( qinti(IHSO4)&
         *1.01030927835D0&     ! mwh2so4/mwihso4
         +qinti(ISO4)&
         *1.02083333333D0&     ! mwh2so4/mwiso4
         )/mass_density(ESO4)          ! µg.µm-3

    vid=vil+vis               ! dry inorg vol µm3.m-3

    ! water volume
    !print*,'qext(EH2O)',qext(EH2O),mass_density(EH2O)
    viw=qext(EH2O)/mass_density(EH2O)

    ! total inorg internal mass
    do jesp=1,N_inside_aer
       sumint=sumint+qinti(jesp)
    enddo
    ! correction for water
    sumint=sumint-qinti(IH2O)-qinti(IOH)+qext(EH2O)

#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
    ! 			      ! correction when no NaCl in internal composition
    sumint = sumint - qinti(INa) + qext(ENa)
    sumint = sumint - qinti(IHCl) - qinti(ICl)*1.02816901408D0&
         + qext(ECl)
#endif

    !!     dry organic volume and total internal mass
    do s=1,nesp_aec
       jesp=aec_species(s)
       vod=vod+qext(jesp)/mass_density(jesp)
       sumint=sumint+qext(jesp)
    enddo

    !!     ******tot dry and wet vol, µm3.m-3
    vad=vod+vid
    vaw=vad+viw

    !!     Notice that the density is based on the internal composition
    !!     and as such bigger than the minimal density and less that the maximal
    !!     one
    if ((vaw.gt.0.d0).AND.(sumint.gt.0.D0)) then
       rhoaer=sumint/vaw
       rhoaer_dry = (sumint - qext(EH2O)) / vad
    else
       rhoaer=fixed_density
       rhoaer_dry = fixed_density
    endif
  end subroutine ssh_VOLAERO

  subroutine ssh_EQINORG(nesp_aer,qext,qinti,surface_equilibrium_conc,lwc,ionic,proton,aerliq,jbin)
    !------------------------------------------------------------------------
    !
    !     -- DESCRIPTION
    !     This subroutine computes the equilibrium between inorganic aerosols
    !     and gas-phase (reverse mode).
    !     It calls ISORROPIA by Nenes et al.
    !
    !------------------------------------------------------------------------
    !
    !     -- INPUT VARIABLES
    !
    !     nesp_aer: number of aerosol species
    !     qext: aerosol mass concentration(µg/m^3)
    !     qinti: aerosol internal species mass concentration(µg/m^3)
    !
    !     -- OUTPU VARIABLES
    !
    !     surface_equilibrium_conc: surface equilibrium concentration of aerosol species
    !------------------------------------------------------------------------
    implicit none
    double precision GREAT
    parameter (GREAT=100.D0)

    integer nesp_aer
    double precision qext(nesp_aer),qinti(N_inside_aer)
    double precision surface_equilibrium_conc(nesp_aer)
    double precision lwc,ionic,proton

    integer jesp,j
    double precision wi(nesp_isorropia),w(nesp_isorropia)
    double precision aerliq(N_liquid),aersld(N_solid)
    double precision gas(3),cntrl(2)
    double precision other(6),coef
    !      CHARACTER*15 scase
    double precision organion,watorg

    integer jbin

    organion = 0.D0
    watorg = 0.D0
    !     Inputs  for Isorropia
    gas(1)=0.d0
    gas(2)=0.d0
    gas(3)=0.d0
    cntrl(1)=1.D0             ! reverse mode
    cntrl(2)=1.D0            ! metastable option
    ! convert µg to moles
    !     conversion unit for isorropia needed in mol.m-3
    do j=1,nesp_isorropia
       jesp=isorropia_species(j)
       wi(j)= qext(jesp)& ! µg.m-3
            /molecular_weight_aer(jesp)!&  ! µg.mol-1
    end do

#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
    wi(1) = 0.D0              !Do not consider sea salt in isoropia
    wi(5) = 0.D0
#endif

    if (iter_eqconc(jbin)==0) then
       call ssh_ISOROPIA(wi,Relative_Humidity,Temperature,cntrl,w,gas,&
            aerliq,aersld,other,organion,watorg)
       !     clipping to tinym
       if (gas(1).lt.0.d0) gas(1)=tinym
       if (gas(2).lt.0.d0) gas(2)=tinym
       if (gas(3).lt.0.d0) gas(3)=tinym
       !     Aqueous phase total liquid water content and pH (proton) concentration       

       ionic = other(5)
       proton = aerliq(IH) * imw(IH) !* gammaH  ! microg.m-3 but equivalent to micromol.m-3      

       ! liquid inorg aerosol
       do jesp=IH,IOH
          qinti(jesp)= DMAX1(aerliq(jesp),0.D0)*imw(jesp)  ! moles to µg

       end do
       ! liquid water content
       qext(EH2O)= qinti(IH2O)&
            +qinti(IOH)*1.05882352941D0 ! mwh2o/mwioh
       if(qext(EH2O) < 1.1d-12) qext(EH2O) = 0.d0 !Minimum lwc is arbitrary fixed in ISORROPIA. Remove it.
       ! solid inorg aerosol
       do jesp=SNaNO3,SLC
          qinti(jesp)= DMAX1(aersld(jesp-12),0.D0)&
               *smw(jesp)        ! moles to µg  
          !molecular_weight_solid(jesp)
       end do

       if (wi(3)>0.d0) then
          ratio_eqconc(1,jbin)=gas(1)/wi(3)
       else
          ratio_eqconc(1,jbin)=0.d0
       endif
       if (wi(4)>0.d0) then
          ratio_eqconc(2,jbin)=gas(2)/wi(4)
       else
          ratio_eqconc(2,jbin)=0.d0
       endif
       if (wi(5)>0.d0) then
          ratio_eqconc(3,jbin)=gas(3)/wi(5)
       else
          ratio_eqconc(3,jbin)=0.d0
       endif
       if (sum(wi)>0.d0) then
          ratio_eqconc(4,jbin)=aerliq(IH2O)/sum(wi)
          iter_eqconc(jbin)=iter_eqconc(jbin)+1
       else
          ratio_eqconc(4,jbin)=0.d0
       endif

    else
       gas(1)=wi(3)*ratio_eqconc(1,jbin)
       gas(2)=wi(4)*ratio_eqconc(2,jbin)
       gas(3)=wi(5)*ratio_eqconc(3,jbin)
       aerliq(:)=0.d0
       aersld(:)=0.d0
       other(:)=0.d0
       aerliq(IH2O)=ratio_eqconc(4,jbin)*sum(wi)
       iter_eqconc(jbin)=iter_eqconc(jbin)+1
       qext(EH2O)=aerliq(IH2O) * imw(IH2O) 
       if(qext(EH2O) < 1.1d-12) qext(EH2O) = 0.d0 !Minimum lwc is arbitrary fixed in ISORROPIA. Remove it.
    endif

    !     Outputs isorropia
    ! sulfate surf conc always 0. µg.m-3
    surface_equilibrium_conc(ESO4)=0.D0

    ! convert moles.m-3 to µg.m-3
    surface_equilibrium_conc(ENH4)=gas(1)*molecular_weight_aer(ENH4)
    surface_equilibrium_conc(ENO3)=gas(2)*molecular_weight_aer(ENO3)
    surface_equilibrium_conc(ECl) =gas(3)*molecular_weight_aer(ECl)
#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
    surface_equilibrium_conc(ECl) = 0.d0 
#endif
    lwc = aerliq(IH2O) * imw(IH2O) ! microg.m-3 

    if (iter_eqconc(jbin)==niter_eqconc) iter_eqconc(jbin)=0

  end subroutine ssh_EQINORG

  subroutine ssh_KLIMIT(q,c_gas,k,ce_kernal_coef)
    !------------------------------------------------------------------------
    !
    !     -- DESCRIPTION
    !     This subroutine limits the condensation/evaporation rates for
    !     aerosol and gases in order to avoid clippings.
    !     Two kinds of limitation are performed:
    !     - the 1st is aerosol clipping : as it may reduce evaporation then
    !     enlarge condensation;, it is done before condensation limitation.
    !     - the 2nd is gas clipping : in practice it may reduce, per species,
    !     aerosol condensation only in bins that leads to gas clipping.
    !
    !------------------------------------------------------------------------
    !
    !     -- INPUT VARIABLES
    !
    !     q: aerosol mass concentration(µg/m^3)
    !     c_gas: aerosol gas concentration(µg/m^3)
    !     ce_kernal_coef: c/e kernel coefficient          ([m3.s-1]).
    !
    !     -- OUTPU VARIABLES
    !
    !     k: particle mass derivation(µg/m^3/s)
    !------------------------------------------------------------------------
    implicit none
    double precision:: q(N_size,N_aerosol_layers)!1th order mass concentration
    double precision:: k(N_size,N_aerosol_layers)
    double precision:: c_gas(N_aerosol)
    double precision:: ce_kernal_coef(N_size,N_aerosol)

    integer:: jesp,j,s
    double precision:: ksum,klim,ktlim
    double precision:: qsum,ce_kernal_coef_tot
    double precision:: frac,qnew

    !     ****** prevent aerosol clipping

    do s=1,nesp_isorropia
       jesp=isorropia_species(s)
       do j =1, N_size
          if(concentration_index(j, 1) > ICUT) then! k : index of size bins
             ! only when evaporation
             if (k(j,jesp).lt.0.D0) then
                ! if q(j,s) is <=TINYM or =0
                ! then k should be >=0, but
                ! due to bad matrix inversion
                ! this case may occur
                if (q(j,jesp).lt.TINYM) k(j,jesp)=0.D0
                ! test clipping in other cases
                qnew=q(j,jesp)+k(j,jesp)
                if (qnew.lt.0.D0) then
                   ! we are sure that q>=tinym
                   ! otherwise k would be = 0
                   ! from previous case
                   k(j,jesp)=(TINYM-q(j,jesp)) !/timestep_splitting ! <=0 µg.m-3.s-1
                   ! we force q to be a
                   ! 'little' more than TINYM
                   k(j,jesp)=0.99D0*k(j,jesp)
                endif
             endif
          endif
       enddo
    enddo
    !     ****** prevent gas clipping
    do s=1,nesp_isorropia
       jesp=isorropia_species(s)
       ! compute total mass rate per species
       ksum=0.D0
       do j =1, N_size
          if(concentration_index(j, 1) > ICUT) then! k : index of size bins
             ksum=ksum+k(j,jesp)     !µg.m-3.s-1
          endif
       enddo
       ! we perform limiting in
       ! case of condensation only
       if(ksum.gt.0.D0) then
          ! this is the total lumped mass
          ! to perserve from clipping
          qsum=c_gas(jesp)
          do j =1, N_size
             if(concentration_index(j, 1) <= ICUT) then! k : index of size bins
                qsum=qsum+q(j,jesp)  ! µg.m-3
             endif
          enddo
          ! test if clipping occurs
          ! then perform the limitation
          if (ksum.gt.qsum) then
             ! sum of ce_kernal_coef(*) c/e coefficient
             ce_kernal_coef_tot=0.D0
             do j= 1, N_size
                ce_kernal_coef_tot=ce_kernal_coef_tot+ce_kernal_coef(j,jesp)
             enddo
             ! tot max rate, µg.m-3.s-1
             ktlim=qsum       !/timestep_splitting
             do j =1, N_size
                if(concentration_index(j, 1) > ICUT) then! k : index of size bins
                   ! fraction, adim
                   if (ce_kernal_coef_tot .ne. 0.d0) then !! YK
                      frac=ce_kernal_coef(j,jesp)/ce_kernal_coef_tot ! ce_kernal_coef_tot != 0
                   else
                      frac = 0.d0
                   endif !! YK

                   ! we allow only a given fraction
                   ! of ktlim to condense on given bin
                   klim=ktlim*frac
                   ! apply the limit
                   ! only if necessary
                   if (k(j,jesp).gt.klim) k(j,jesp)=klim
                endif
             enddo
          endif
       endif
    enddo

  end subroutine ssh_KLIMIT


  subroutine ssh_HPLFLIM(alfa,qih,N_size_loc,init_bulk_gas,&
       ce_kernal_coef_i,Kelvin_effect,surface_equilibrium_conc,ce_kernel)
    !------------------------------------------------------------------------
    !
    !     -- DESCRIPTION
    !     This subroutine computes the flux limitation for the
    !     condensation/evaporation flux. The algorithm is based
    !     on the limitation of the aerosol acidity rate.
    !
    !     The details may be found in the PhD Work of Edouard Debry,
    !     Chapter 10 (section 10.1.5) or in the reference:
    !     Pilinis et al: MADM, a new multicomponent aerosol dynamic model
    !     Aerosol Science and Technology 32, 482:502, 2000.
    !------------------------------------------------------------------------
    !
    !     -- INPUT VARIABLES
    !
    !     alfa: percentage of H+ allowed to c/e(0.1)
    !     qih: int H+ conc (µg)
    !     N_size_loc: size of vectors following below
    !     init_bulk_gas: bulk gas conc (µg.m-3)
    !     surface_equilibrium_conc: surface equilibrium concentration of aerosol species(µg/m^3)
    !     Kelvin_effect:kelvin coef (adim)
    !     ce_kernal_coef_i: c/e kernel coefficient          ([m3.s-1]).
    !
    !     -- OUTPU VARIABLES
    !
    !     ce_kernel: modified c/e kernel (µg/m^3)
    !------------------------------------------------------------------------
    implicit none
    integer N_size_loc
    double precision init_bulk_gas(N_size_loc),ce_kernal_coef_i(N_size_loc)
    double precision Kelvin_effect(N_size_loc),surface_equilibrium_conc(N_size_loc)
    double precision ce_kernel(N_size_loc),qih,alfa

    integer jesp,s
    double precision maa(N_size_loc),mkercd(N_size_loc)
    double precision cfa,cfb,cfc,cc
    double precision mih,melec,mlim,q

    !     Compute mol fluxes

    do s=1,nesp_isorropia
       jesp=isorropia_species(s)
       ! maa(*) in m3.mol.s-1.µg-1
       maa(jesp)= ce_kernal_coef_i(jesp)&   ! m3.s-1
            /molecular_weight_aer(jesp)        ! µg.mol-1
       ! mkercd(*) in mol.s-1
       mkercd(jesp)= ce_kernel(jesp)& ! µg.s-1
            /molecular_weight_aer(jesp)        ! µg.mol-1
    end do

    !     H+ limitation

    mih=qih/imw(IH)           ! mol of H+ in aerosol


    ! maximum of mih variation tolerated
    mlim=mih*alfa             ! mol.s-1 alfa=(0.1)

    ! electroneutrality  ! mol.s-1
    melec= 2.D0*mkercd(ESO4)+mkercd(ENO3)&
         +mkercd(ECl)-mkercd(ENH4)

    ! correction factor default value
    cc=0.D0
    ! correction calculation
    if (DABS(melec).gt.mlim) then
       ! we give to mlim the sign of melec
       mlim=DSIGN(mlim,melec)! returns the value of mlim with the sign of melec
       !to judge + or -
       ! cfa,cfb,cfc are coefficients of
       ! 2nd order eq : cfa*cc^2+cfb*cc+cfc=0
       ! satisfied by the correction factor cc
       cfa= maa(ENO3)*surface_equilibrium_conc(ENO3)*Kelvin_effect(ENO3)&
            +maa(ECl)*surface_equilibrium_conc(ECl)*Kelvin_effect(ECl)
       cfb=mlim-2.D0*maa(ESO4)*init_bulk_gas(ESO4)& ! mol.s-1
            -maa(ENO3)*init_bulk_gas(ENO3)&
            -maa(ECl) *init_bulk_gas(ECl)&
            +maa(ENH4)*init_bulk_gas(ENH4)
       cfc=-maa(ENH4)*surface_equilibrium_conc(ENH4)*Kelvin_effect(ENH4)
       ! one can note cfa>=0 and cfc<=0
       ! then there always exist a + but
       ! possibly zero root
       ! root computation
       if (cfa.gt.0.D0) then
          if (cfb*cfb-4.D0*cfa*cfc.lt.0.D0) then
             WRITE(6,*)'(hplflim.f): sqrt(<0)'
             WRITE(6,*)cfb,cfa,cfc,mlim,'Time'
             WRITE(6,*)'ASO4',ce_kernal_coef_i(ESO4),'ANH4',ce_kernal_coef_i(ENH4),'ANO3',&
                  ce_kernal_coef_i(ENO3),'ACl',ce_kernal_coef_i(ECl)
             WRITE(6,*)'SO4',mkercd(ESO4),'NH4',mkercd(ENH4),&
                  'NO3',mkercd(ENO3),'Cl',mkercd(ECl)
             WRITE(6,*)'SO4',maa(ESO4),init_bulk_gas(ESO4)
             WRITE(6,*)'NH4',maa(ENH4),init_bulk_gas(ENH4),surface_equilibrium_conc(ENH4)
             WRITE(6,*)'NO3',maa(ENO3),init_bulk_gas(ENO3),surface_equilibrium_conc(ENO3)
             WRITE(6,*)'Cl',maa(ECl),init_bulk_gas(ECl),surface_equilibrium_conc(ECl)
             STOP
          endif
          q=-5.D-01*(cfb+DSIGN(1.D0,cfb)&
               *DSQRT(cfb*cfb-4.D0*cfa*cfc))
          cc=DMAX1(q/cfa,cfc/q) ! we select the + root
       else
          if (cfb.ne.0.D0) cc=-cfc/cfb
       endif
    endif

    !     A correction is done if only upper calculation
    !     root has changed cc to a strictly positive
    !     value, otherwise it is considered as non stiff cases

    if (cc.gt.0.D0) then
       surface_equilibrium_conc(ENH4)=surface_equilibrium_conc(ENH4)/cc
       surface_equilibrium_conc(ENO3)=surface_equilibrium_conc(ENO3)*cc
       surface_equilibrium_conc(ECl) =surface_equilibrium_conc(ECl)*cc
       do jesp=ENH4,ECl
	  ce_kernel(jesp)= ce_kernal_coef_i(jesp)*(init_bulk_gas(jesp)-surface_equilibrium_conc(jesp)*Kelvin_effect(jesp))
       end do
    endif

  end subroutine ssh_HPLFLIM

  subroutine ssh_equi_const(TEMP, XK3, XK4, XK6, XK8, XK9, XK10)
    implicit none
    double precision TEMP, T0, T0T, COEF
    double precision XK3, XK4, XK6, XK8, XK9, XK10

    XK3  = 1.971D6   ! HCL(g)           <==> H(aq)     + CL(aq)
    XK4  = 2.511e6   ! HNO3(g)          <==> H(aq)     + NO3(aq) ! ISORR
    XK6  = 1.086D-16 ! NH4CL(s)         <==> NH3(g)    + HCL(g)
    XK8  = 37.661D0  ! NACL(s)          <==> NA(aq)    + CL(aq)
    XK9  = 11.977D0  ! NANO3(s)         <==> NA(aq)    + NO3(aq)
    XK10 = 5.746D-17 ! NH4NO3(s)        <==> NH3(g)    + HNO3(g) ! ISORR

    IF (INT(TEMP) .NE. 298) THEN   ! FOR T != 298K or 298.15K
       T0  = 298.15D0
       T0T = T0/TEMP
       COEF= 1.0+LOG(T0T)-T0T

       XK3 = XK3 *EXP( 30.20*(T0T-1.0) + 19.910*COEF)
       XK4 = XK4 *EXP( 29.17*(T0T-1.0) + 16.830*COEF) !ISORR
       XK6 = XK6 *EXP(-71.00*(T0T-1.0) +  2.400*COEF)
       XK8 = XK8 *EXP( -1.56*(T0T-1.0) + 16.900*COEF)
       XK9 = XK9 *EXP( -8.22*(T0T-1.0) + 16.010*COEF)
       XK10= XK10*EXP(-74.38*(T0T-1.0) +  6.120*COEF) ! ISORR
    ENDIF

  end subroutine ssh_equi_const

  subroutine ssh_DRYIN(Temp,qinti,N_size_loc,init_bulk_gas,ce_kernal_coef_i,&
       Kelvin_effect,surface_equilibrium_conc,ce_kernel)
    !------------------------------------------------------------------------
    !
    !     -- DESCRIPTION
    !     This subroutine computes the aerosol surface gas-phase
    !     concentration (through equilibrium) for dry aerosols.
    !
    !     The algorithms are detailed in Chapter 10 (section 10.1.6) of
    !     the PhD work of Edouard Debry.
    !     See also the reference:
    !     Pilinis et al: MADM, a new multicomponent aerosol dynamic model
    !     Aerosol Science and Technology 32, 482:502, 2000.
    !------------------------------------------------------------------------
    !
    !     -- INPUT VARIABLES
    !
    !     qih: int H+ conc (µg)
    !     N_size_loc: size of vectors following below
    !     init_bulk_gas: bulk gas conc (µg.m-3)
    !     surface_equilibrium_conc: surface equilibrium concentration of aerosol species(µg/m^3)
    !     Kelvin_effect:kelvin coef (adim)
    !     ce_kernal_coef_i: c/e kernel coefficient          ([m3.s-1]).
    !
    !     -- OUTPU VARIABLES
    !
    !     ce_kernel: modified c/e kernel (µg/m^3)
    !------------------------------------------------------------------------
    ! ISORROPIA commons needed by this routine,
    ! directly taken from isrpia.inc.
    implicit none
    double precision XK3,XK4,XK6,XK8,XK9,XK10
    ! COMMON /EQUK/ XK1,XK2,XK3,XK4,XK5,XK6,XK7,XK8,&
    !     	  XK9,XK10,XK11,XK12,XK13,XK14,&
    !     	  XKW,XK21,XK22,XK31,XK32,XK41,XK42 ! YK
!!!$OMP THREADPRIVATE(/EQUK/)
    integer N_size_loc
    double precision qinti(N_inside_aer),Kelvin_effect(N_size_loc)
    double precision ce_kernal_coef_i(N_size_loc),init_bulk_gas(N_size_loc)
    double precision surface_equilibrium_conc(N_size_loc),ce_kernel(N_size_loc)
    double precision Temp
    LOGICAL leq1,leq2,leq3,lr47,lr56
    integer icase,jesp,s
    double precision rgas1,maa(N_size_loc)
    double precision rk1,rk2,rk3
    double precision sat1,sat2,sat3a,sat3b
    double precision mkercd(N_size_loc),msat
    double precision cfa,cfb,cfc

    !     Initialization:
    !     1 stands for nh4no3 equilibrium
    !     2 stands for nh4cl equilibrium
    !     3 stands for nacl/nano3 equilibrium
    !     47 for nacl and nh4cl reactions
    !     56 for nano3 and nh4no3 reactions

    leq1=.false.
    leq2=.false.
    leq3=.false.
    lr47=.false.
    lr56=.false.

    rgas1=ATM/(RGAS*Temp)

    call ssh_equi_const(Temp, XK3, XK4, XK6, XK8, XK9, XK10)

    rk1= XK10*rgas1*rgas1*molecular_weight_aer(ENH4)&
         *molecular_weight_aer(ENO3)*Kelvin_effect(ENH4)*Kelvin_effect(ENO3)
    ! rk1 in (µg.m-3)2

    rk2=XK6*rgas1*rgas1*molecular_weight_aer(ENH4)&
         *molecular_weight_aer(ECl)*Kelvin_effect(ENH4)*Kelvin_effect(ECl)
    ! rk2 in (µg.m-3)2

    rk3=XK4*XK8/(XK3*XK9) ! adim

    sat1=init_bulk_gas(ENH4)*init_bulk_gas(ENO3) ! (µg.m-3)2
    sat2=init_bulk_gas(ENH4)*init_bulk_gas(ECl) ! (µg.m-3)2

    sat3a=init_bulk_gas(ECl)          ! µg.m-3
    sat3b=rk3*init_bulk_gas(ENO3)     ! µg.m-3

    do s=2,nesp_isorropia
       jesp=isorropia_species(s)
       ! maa(*) in m3.mol.s-1.µg-1
       maa(jesp)= ce_kernal_coef_i(jesp)&   ! m3.s-1
            /molecular_weight_aer(jesp)        ! µg.mol-1

       ! mkercd(*) in mol.s-1
       mkercd(jesp)= ce_kernel(jesp)& ! µg.s-1
            /molecular_weight_aer(jesp)        ! µg.mol-1
    end do

    msat=2.D0*mkercd(ESO4)-mkercd(ENH4) ! mol.s-1

    !     Determine which reaction is active

    if (qinti(SNH4NO3).gt.0.D0.OR.sat1.gt.rk1) then
       leq1=.true.
    endif

    if (qinti(SNH4Cl).gt.0.D0.OR.sat2.gt.rk2) then
       leq2=.true.
    endif

    if (qinti(SNaNO3).gt.0.D0.AND.qinti(SNaCl).gt.0.D0) then
       leq3=.true.
    endif

    if (qinti(SNaNO3).gt.0.D0.AND. sat3a.gt.sat3b) then
       leq3=.true.
    endif

    if (qinti(SNaCl).gt.0.D0.AND.sat3a.lt.sat3b) then
       leq3=.true.
    endif

    if (init_bulk_gas(ESO4).gt.0.D0) then
       if (qinti(SNaCl).gt.0.D0.OR.qinti(SNH4Cl).gt.0.D0) then
          lr47=.true.
       endif

       if (qinti(SNaNO3).gt.0.D0.OR.qinti(SNH4NO3).gt.0.D0) then
          lr56=.true.
       endif
    endif

    if (leq1.AND.leq2.AND.leq3) then
       PRINT *,'Warning from dryin.f: << solid : leq123 >>'
    endif

    !     Determine which case is relevant

    if (leq2.AND.leq3) then
       icase=1                ! R2 and R3 active
    elseif (leq1.AND.leq2) then
       icase=2                ! R1 and R2 active
    elseif (leq1.AND.leq3) then
       icase=3                ! R1 and R3 active
    elseif (leq1) then
       icase=4                ! only R1 active
    elseif (leq2) then
       icase=5                ! only R2 active
    elseif (leq3) then
       icase=6                ! only R3 active
    else                      ! no active equilibrium
       if (lr47) then
          icase=7             ! R4 or R7 active
       elseif (lr56) then
          icase=8             ! R5 or R6 active
       else                   ! nothing active

          if (msat.lt.0.D0) then
             icase=9          ! enough nh3 to neutralize so4
          else
             icase=10         ! not enough nh3 to neutralize so4
             ! in this case aerosol become acidic
          endif
       endif
    endif
    !     Solve each case
    ! icase 3 is not physical but used
    ! to determine the real icase
    if (icase.eq.3) then
       cfa=( maa(ENO3)*Kelvin_effect(ENO3)+maa(ECl)*Kelvin_effect(ECl)*rk3 )

       cfb=2.D0*mkercd(ESO4)+mkercd(ENO3)+mkercd(ECl)-mkercd(ENH4)

       cfc=rk1*maa(ENH4)*Kelvin_effect(ENH4)

       if (cfb*cfb+4.D0*cfa*cfc.le.0.D0) then
          WRITE(6,*)'(dryin.f): (1) sqrt(<0) '
          STOP
       endif
       surface_equilibrium_conc(ENO3)= (cfb+DSQRT(cfb*cfb+4.D0*cfa*cfc)) /(2.D0*cfa)
       !surface_equilibrium_conc(ECl)=rk3*surface_equilibrium_conc(ENO3)
       !surface_equilibrium_conc(ENH4)=rk1/surface_equilibrium_conc(ENO3)
       ce_kernel(ENO3)= ce_kernal_coef_i(ENO3) *( init_bulk_gas(ENO3)&
            -surface_equilibrium_conc(ENO3)*Kelvin_effect(ENO3) )
       ! test if no3 condenses
       if (ce_kernel(ENO3).gt.0.D0) then
          icase=2             ! if nh4no3 forms then real case=2
       else
          icase=1             ! if nacl forms then real case=1
       endif
    endif

    ! other cases
    if (icase.eq.1) then
       cfa=( maa(ENO3)*Kelvin_effect(ENO3)+maa(ECl)*Kelvin_effect(ECl)*rk3 )
       cfb=2.D0*mkercd(ESO4)+mkercd(ENO3)+mkercd(ECl)-mkercd(ENH4)
       cfc=rk2/rk3*maa(ENH4)*Kelvin_effect(ENH4)
       if (cfb*cfb+4.D0*cfa*cfc.le.0.D0) then
          WRITE(6,*)'(dryin.f): (2) sqrt(<0) '
          STOP
       endif
       surface_equilibrium_conc(ENO3)= (cfb+DSQRT(cfb*cfb+4.D0*cfa*cfc))/(2.D0*cfa)
       surface_equilibrium_conc(ECl)=rk3*surface_equilibrium_conc(ENO3)
       surface_equilibrium_conc(ENH4)=rk2/rk3/surface_equilibrium_conc(ENO3)
    elseif (icase.eq.2) then
       cfa=( maa(ENO3)*Kelvin_effect(ENO3)+maa(ECl)*Kelvin_effect(ECl)*rk2/rk1 )
       cfb=2.D0*mkercd(ESO4)+mkercd(ENO3)+mkercd(ECl)-mkercd(ENH4)
       cfc=rk1*maa(ENH4)*Kelvin_effect(ENH4)
       if (cfb*cfb+4.D0*cfa*cfc.le.0.D0) then
          WRITE(6,*)'(dryin.f): (3) sqrt(<0) '
          STOP
       endif
       surface_equilibrium_conc(ENO3)= (cfb+DSQRT(cfb*cfb+4.D0*cfa*cfc))/(2.D0*cfa)
       surface_equilibrium_conc(ENH4)=rk1/surface_equilibrium_conc(ENO3)
       surface_equilibrium_conc(ECl)=rk2/rk1*surface_equilibrium_conc(ENO3)
    elseif (icase.eq.4) then
       cfa=maa(ENO3)*Kelvin_effect(ENO3)
       cfb= 2.D0*mkercd(ESO4)+mkercd(ENO3)-mkercd(ENH4)
       cfc=rk1*maa(ENH4)*Kelvin_effect(ENH4)
       if (cfb*cfb+4.D0*cfa*cfc.le.0.D0) then
          WRITE(6,*)'(dryin.f): (4) sqrt(<0) '
          STOP
       endif
       surface_equilibrium_conc(ENO3)= (cfb+DSQRT(cfb*cfb+4.D0*cfa*cfc))/(2.D0*cfa)
       if (surface_equilibrium_conc(ENO3)==0.0.and.cfb<0) then       
          surface_equilibrium_conc(ENO3)=-cfc/cfb          
       endif
       surface_equilibrium_conc(ENH4)=rk1/surface_equilibrium_conc(ENO3)
       surface_equilibrium_conc(ECl)=init_bulk_gas(ECl)/Kelvin_effect(ECl)
    elseif (icase.eq.5) then
       cfa=maa(ENH4)*Kelvin_effect(ENH4)
       cfb= 2.D0*mkercd(ESO4)+mkercd(ECl)-mkercd(ENH4)
       cfc=maa(ECl)*rk2*Kelvin_effect(ECl)
       if (cfb*cfb+4.D0*cfa*cfc.le.0.D0) then
          WRITE(6,*)'(dryin.f): (5) sqrt(<0) '
          STOP
       endif
       surface_equilibrium_conc(ENH4)= (cfb+DSQRT(cfb*cfb+4.D0*cfa*cfc))/(2.D0*cfa)
       surface_equilibrium_conc(ECl)=rk2/surface_equilibrium_conc(ENH4)
       surface_equilibrium_conc(ENO3)=init_bulk_gas(ENO3)/Kelvin_effect(ENO3)
    elseif (icase.eq.6) then
       cfa=( maa(ENO3)*Kelvin_effect(ENO3)+rk3*maa(ECl)*Kelvin_effect(ECl) )
       cfb=2.D0*mkercd(ESO4)+mkercd(ECl)+mkercd(ENO3)
       surface_equilibrium_conc(ENO3)=cfb/cfa
       surface_equilibrium_conc(ECl)=rk3*surface_equilibrium_conc(ENO3)
       surface_equilibrium_conc(ENH4)=init_bulk_gas(ECl)/Kelvin_effect(ECl)
    elseif (icase.eq.7) then
       surface_equilibrium_conc(ENH4)=init_bulk_gas(ENH4)/Kelvin_effect(ENH4)
       surface_equilibrium_conc(ENO3)=init_bulk_gas(ENO3)/Kelvin_effect(ENO3)

       ce_kernel(ENH4)=0.D0
       ce_kernel(ENO3)=0.D0

       ce_kernel(ECl)=-molecular_weight_aer(ECl)*2.D0*mkercd(ESO4)
       surface_equilibrium_conc(ECl)= ( init_bulk_gas(ECl)-ce_kernel(ECl)&
            /ce_kernal_coef_i(ECl) )/Kelvin_effect(ECl)
    elseif (icase.eq.8) then
       surface_equilibrium_conc(ENH4)=init_bulk_gas(ENH4)/Kelvin_effect(ENH4)
       surface_equilibrium_conc(ECl)=init_bulk_gas(ECl)/Kelvin_effect(ECl)
       ce_kernel(ENH4)=0.D0
       ce_kernel(ECl)=0.D0
       ce_kernel(ENO3)=-molecular_weight_aer(ECl)*2.D0 *mkercd(ESO4)
       surface_equilibrium_conc(ENO3)= ( init_bulk_gas(ENO3)-ce_kernel(ENO3)&
            /ce_kernal_coef_i(ENO3) )/Kelvin_effect(ENO3)
    elseif (icase.eq.9) then
       surface_equilibrium_conc(ENO3)=init_bulk_gas(ENO3)/Kelvin_effect(ENO3)
       surface_equilibrium_conc(ECl)=init_bulk_gas(ECl)/Kelvin_effect(ECl)
       ce_kernel(ENO3)=0.D0
       ce_kernel(ECl)=0.D0
       ce_kernel(ENH4)= molecular_weight_aer(ENH4)*2.D0 *mkercd(ESO4)
       surface_equilibrium_conc(ENH4)= ( init_bulk_gas(ENH4)-ce_kernel(ENH4)&
            /ce_kernal_coef_i(ENH4) )/Kelvin_effect(ENH4)
    elseif (icase.eq.10) then
       surface_equilibrium_conc(ENO3)=init_bulk_gas(ENO3)/Kelvin_effect(ENO3)
       surface_equilibrium_conc(ECl)=init_bulk_gas(ECl)/Kelvin_effect(ECl)
       ce_kernel(ENO3)=0.D0
       ce_kernel(ECl)=0.D0
       !    no more electroneutrality in this case
    endif
    !     Giving out the kernel for case <=6
    if (icase.LE.6) then
       do s=3,nesp_isorropia
          jesp=isorropia_species(s)
          ce_kernel(jesp)= ce_kernal_coef_i(jesp)*( init_bulk_gas(jesp)-&
               surface_equilibrium_conc(jesp)*Kelvin_effect(jesp) )
       end do
    endif
  end subroutine ssh_DRYIN


  subroutine ssh_isoropia_drv(nesp_aer,&
       aero, gas, organion, watorg, ionic, proton, lwc, rh, Temp, &
       liquid)

    !------------------------------------------------------------------------
    !
    !     -- DESCRIPTION
    !
    !     This subroutine computes the equilibrium between inorganic aerosols
    !     and gas-phase (forward mode), taking in account organic liquid
    !     water content and organic ions.
    !     It calls ISORROPIA by Nenes et al.
    !
    !------------------------------------------------------------------------
    !
    !     -- INPUT VARIABLES
    !
    !     ORGANION: organic ions ([\mu mol.m^-3]).
    !     WATORG: organic liquid water content ([\mu g.m^-3]).
    !     Relative_Humidity: relative humidity 0< <1 ([]).
    !     Temperature: temperature ([Kelvin]).
    !
    !     -- INPUT/OUTPUT VARIABLES
    !
    !     AERO: aerosol bulk concentration ([\mu g.m^-3]).
    !     GAS: gas concentration ([\mu g.m^-3]).
    !
    !     -- OUTPUT VARIABLES
    !
    !     PROTON: hydronium ion concentration ([\mu g.m^-3]).
    !     LWC: total liquid water content ([\mu g.m^-3]).
    !------------------------------------------------------------------------

    IMPLICIT NONE

    integer nesp_aer
    double precision aero(nesp_aer), gas(nesp_aer)
    double precision organion, watorg, proton
    double precision lwc, rh, Temp
    double precision wi(N_inorganic),w(N_inorganic),gas2(3),cntrl(2), other(6)
    double precision liquid(N_liquid),solid(N_solid)
    double precision organion2, watorg2, ionic ,gammaH
    integer i,idx !,j
    !     mol neg charge in mol.m-3 */
    organion2 = organion * 1.D-6
    liquid=0.d0
    solid=0.d0
    w=0.d0
    gas2=0.d0
    other=0.d0
    !     organic water content converted from
    !     microg/m3 (aec output) to kg/m3 (isorropia input)
    watorg2 = watorg * 1.D-9
    cntrl(1) = 0.D0
    cntrl(2) = MTSBL
    !     concentration in microg.m-3
    gas(ENa) = 0.D0

    do i=1,nesp_isorropia
       idx = isorropia_species(i)
       wi(i) = aero(idx) + gas(idx)
       wi(i) = wi(i) / molecular_weight_aer(idx) ! microg.m-3 / microg.mol-1 = mol.m-3
    enddo

    call SSH_ISOROPIA(wi, rh, Temp, cntrl, w, gas2,&
         liquid, solid, other, organion2, watorg2)

    !     Aqueous phase total liquid water content and pH (proton) concentration
    lwc = liquid(IH2O) * imw(IH2O) ! microg.m-3 
    ionic = other(5)
    !    gammaH = 10**(-0.511 * (298.0/Temperature)**1.5 * sqrt(ionic)/(1+sqrt(ionic)))
    !! YK gammaH will be calculated later in SOAP

    proton = liquid(IH) * imw(IH) !* gammaH  ! microg.m-3 but equivalent to micromol.m-3

    gas(isorropia_species(2)) = 0.D0
    do i=3,nesp_isorropia
       idx = isorropia_species(i)
       gas(idx) = gas2(i-2) * molecular_weight_aer(idx)
    enddo

    do i=1,nesp_isorropia
       idx = isorropia_species(i)
       if (i > 2) then
          aero(idx) = (w(i) - gas2(i-2)) * molecular_weight_aer(idx)
       else
          aero(idx) = w(i) * molecular_weight_aer(idx)
       endif
    enddo

  end subroutine ssh_isoropia_drv

End module cThermodynamics
