module mod_sshchem

  use aInitialization
  use mod_sshchemkinetic
  
  implicit none

contains   
                                                    
  SUBROUTINE ssh_chem ()                             
                                                                        
!-----------------------------------------------------------------------
!                                                                       
!     -- DESCRIPTION                                                    
!                                                                       
!     This routine currently can not be used for GENOA reductions.      
!     Please use twostep solver.                                        
!                                                                       
!-----------------------------------------------------------------------
                                                                        
      IMPLICIT NONE 
                                                                       
!     ! exit                                                            
      print*, "ssh_chem has not been updated for GENOA reduction." 
      print*, "Please use twostep solver (tag_twostep = 1)." 
      stop 
                                                                        
  END SUBROUTINE ssh_chem
                                                                        
                                                                                                                                                
  SUBROUTINE ssh_chem_twostep (current_time,delta_t_now,nstep)
                                                                        
!-----------------------------------------------------------------------
!                                                                       
!     -- DESCRIPTION                                                    
!                                                                       
!     This routine computes the chemistry with the                      
!     TWO-STEP time numerical solver. It is based on                    
!     the application of a Gauss-Seidel iteration                       
!     scheme to the two-step implicit backward                          
!     differentiation (BDF2) formula.                                   
!                                                                       
!     Reference: Verwer J. (1994). Gauss-Seidel                         
!       iteration for stiff odes from chemical kinetics.                
!       Journal on Scientific Computing, 15:1243-1250                   
!                                                                       
!-----------------------------------------------------------------------
                                                                        
                                                                        
      IMPLICIT NONE 
                                                                        
      INCLUDE 'CONST.INC' 
      INCLUDE 'CONST_A.INC' 
      INCLUDE 'paraero.inc'
      
      ! use humidity, temperature, pressure at this timestep
      
      ! no emission
      ! use: longitude, latitude                               
      DOUBLE PRECISION Zangzen !solar zenith angle
      DOUBLE PRECISION ssh_muzero,DLmuzero 
      EXTERNAL ssh_muzero 
                                                                      
      INTEGER INUM, IDENS ! use: with_fixed_density
      PARAMETER (INUM = 1)                                              
      DOUBLE PRECISION RHOA 
      double precision rho_dry(n_size) 
      DOUBLE PRECISION MSF(n_size) 
      DOUBLE PRECISION DSF(n_size) 
      DOUBLE PRECISION DBF((n_size/n_fracmax)+1) 
      INTEGER idx_bs(n_size) 
                                                                        
!     Parameters initialized for the two-step solver                    
      integer m,j,i1,Jsp,i,Jb, nstep
      DOUBLE PRECISION current_time, delta_t_now ! current time, use delta_t
      DOUBLE PRECISION tschem,tfchem,tstep,tstep_min ! use min_adaptive_time_step
      DOUBLE PRECISION :: wk,dtnsave,error_max,c,gam,ratloss
      DOUBLE PRECISION dun,dzero,alpha
      PARAMETER (alpha=5.d0) 
      PARAMETER (dun=1.d0) 
      PARAMETER (dzero=0.d0)
      DOUBLE PRECISION rtol,atol ! use adaptive_time_step_tolerance 
      DOUBLE PRECISION conci(n_gas),concii(n_gas) 

      ! concs use: conversionfactor

      ! ZC,DLRkf,DLRki,chpr0,chlo0,dw: y,w,rk,prod,loss,dw
      ! gas_yield,rcn_rate, kinetic_rate,chem_prod, chem_loss,drv_knt_rate
      
     ! n_gas, n_reaction: ns,nr
     ! DLconc: concentration_gas_all
     ! DLconc_aer: concentration_mass
     ! DLnumconc_aer: concentration_number

!     genoa for RO2-RO2 reaction and constant profile                   
      ! use: ncst_gas,nRO2_chem,iRO2,iRO2_cst,tag_RO2 
      ! use: cst_gas_index(ncst_gas), RO2index(nRO2_chem), cst_gas
      DOUBLE PRECISION RO2, ZC_cst(ncst_gas) 
                                                                        
!     genoa keep_gp: gas-particle partitioning                          
      ! used for genoa with large timestep           
      integer s, keep_gp 
      integer aerosol_species_interact(n_gas) 
      DOUBLE PRECISION toadd, conc_tot, ZCtot_save(n_gas) 
      DOUBLE PRECISION ratio_gas(n_gas),ZCtot(n_gas),DLconc_save(n_gas) 
                                                                        
      ! use: Nsps_rcn(:,:), photo_rcn(:,:), TB_rcn(:,:) 
      ! use: fall_rcn(:,:), extra_rcn(:,:) 
      ! use: index_RCT(:), index_PDT(:) 
      ! kinetics                                                        
      ! use: Arrhenius(:,:), fall_coeff(:,:), ratio_PDT(:) , photo_ratio(:,:,:) 

                                                                        
!     Aerosol density converted in microg / microm^3.                   
      RHOA = fixed_density * 1.D-09 
!     Aerosol discretization converted in microm.                       
      DO Jb=1,(n_size/n_fracmax)+1 
         DBF(Jb) = diam_bound(Jb) * 1.D06 
      ENDDO 
                                                                        
!     relations between bin idx and size idx                            
      DO Jb=1,n_size 
      idx_bs(Jb)=(Jb-1)/n_fracmax+1 
      ENDDO 
                                                                        
!     With real number concentration.                                   
      IF (INUM.EQ.1) THEN
         ! Compute aerosol density                                           
         rho_dry = RHOA
         IDENS = not(with_fixed_density)
         ! for varying density                   
         IF (IDENS.EQ.1) THEN 
            DO Jb=1,n_size 
               CALL SSH_COMPUTE_DENSITY(n_size,n_aerosol, n_aerosol, TINYM, &
     &              concentration_mass,                                         &
     &              mass_density,Jb,rho_dry(Jb))                    
            ENDDO 
         ENDIF
         DO Jb = 1, n_size
            conc_tot = 0.d0
            DO Jsp = 1, n_aerosol
               conc_tot = conc_tot + concentration_mass(Jb,Jsp)
            ENDDO
            
!     Compute mass and diameter of each section                         
            IF (concentration_number(Jb) .GT. 0.d0) THEN 
               MSF(Jb) = conc_tot/concentration_number(Jb) 
            ELSE 
               MSF(Jb) = 0.d0 
            ENDIF 
                                                                        
            if ((concentration_number(Jb).GT. TINYN .or.                       &
     &           conc_tot.GT.TINYM)                                     &
     &           .AND. IDENS .EQ. 1) then                               
               DSF(Jb) = (MSF(Jb)/cst_PI6/rho_dry(Jb))**cst_FRAC3 
                                                                        
            else !sz   
               DSF(Jb) = DSQRT(DBF(idx_bs(Jb))* DBF(idx_bs(Jb)+1)) 
            endif 
                                                                        
            if (DSF(Jb) .LT. DBF(idx_bs(Jb)) .or.                       &
     &          DSF(Jb) .GT. DBF(idx_bs(Jb)+1)) THEN                    
               DSF(Jb) =  DSQRT(DBF(idx_bs(Jb)) * DBF(idx_bs(Jb)+1)) 
            endif 
                                                                        
         ENDDO 
                                                                        
      ELSE 
         DO Jb = 1, n_size 
                                                                !sz     
            DSF(Jb) = DSQRT(DBF(idx_bs(Jb)) * DBF(idx_bs(Jb)+1)) 
            MSF(Jb) = RHOA * cst_pi6 * DSF(Jb)**3 
         ENDDO 
      ENDIF 

!     Projection.                                                       
!     Conversion mug/m3 to molecules/cm3.                               
      DO Jsp=1,n_gas 
         gas_yield(Jsp) = concentration_gas_all(Jsp)* conversionfactor(Jsp) 
         ! genoa keep_gp                                                
         DLconc_save(Jsp)=concentration_gas_all(Jsp) 
         ZCtot(Jsp)=gas_yield(Jsp) 
         ratio_gas(Jsp)=1.d0 
      ENDDO 
                                                                        
      ! genoa keep_gp                                                   
      if (keep_gp==1) then 
         do s = 1, n_aerosol 
            Jsp=aerosol_species_interact(s) 
            !print*,s,Jsp                                               
            if (Jsp>0.) then 
               conc_tot=0.0D0 
               DO Jb=1,n_size 
                  conc_tot = conc_tot + concentration_mass(Jb,s) 
               ENDDO 
                                                                        
               ZCtot(Jsp) = (concentration_gas_all(Jsp)+conc_tot)*conversionfactor(Jsp) 
               if ( ZCtot(Jsp)>0.d0) then 
                  ratio_gas(Jsp)= gas_yield(Jsp)/ZCtot(Jsp) 
               else 
                  ratio_gas(Jsp)= 1.d0 
               endif 
            endif 
         ENDDO 
      endif 
      !print*,maxval(ratio_gas)                                         
                                                                        
!!  input constant concentration genoa                                  
      if (ncst_gas.gt.0) then 
        do i1=1,ncst_gas 
           ZC_cst(i1)=cst_gas(i1,nstep)* conversionfactor(cst_gas_index(i1)) 
           gas_yield(cst_gas_index(i1))=ZC_cst(i1) 
           ! genoa keep_gp                                              
           ZCtot(cst_gas_index(i1))=ZC_cst(i1) 
           ZCtot_save(cst_gas_index(i1))=ZCtot(cst_gas_index(i1)) 
           ratio_gas(i1)=1.d0 
        enddo 
      endif
      ZCtot_save(:)=ZCtot(:) 
                                                                        
!     two-step solver starts                                            
      !tolerence                                                        
      rtol = adaptive_time_step_tolerance 
      atol = 1.d-6*rtol 
      
      !m nombre d'iteration                                         
      m=2 
      ! initial step tstep: dtn  
      tstep_min = min_adaptive_time_step                
      tstep=tstep_min 
      ! initial timestep                                      
      tschem=current_time
      ! total chemistry step                       
      tfchem=tschem+delta_t_now 
                                                                        
      ! init conc.                                                      
      do Jsp=1,n_gas ! genoa keep_gp                         
         concii(Jsp)=ZCtot(Jsp)
         conci(Jsp)=ZCtot(Jsp)
      enddo 
                                                                        
      !premier calcul de l'ordre 1                                      
      do j=1,m 
                                                                        
         !computes chpr and chlo which are the vectors with all the terms
         ! compute RO2                                                  
         RO2=0.d0 
         if (tag_RO2.ne.0) then 
           ! add background if tag = 2,3,4 (from ro2 file)              
           if (tag_RO2.ge.2 .and. iRO2_cst.ne.0) then                                 
             RO2 = RO2 + ZC_cst(iRO2_cst) 
           endif 
           ! add from primary VOCs (RO2 list) tag = 1,3                 
           if (tag_RO2 .eq. 1 .or. tag_RO2 .eq. 3) then
             do Jsp =1,nRO2_chem 
               RO2=RO2+gas_yield(RO2index(Jsp)) 
             enddo 
           endif 
         endif 
                                                                        
        ! kinetic rate                                                  
        ! Compute zenithal angles                                       
        DLmuzero=ssh_muzero(tschem,longitude,latitude) 
        Zangzen=dabs(DACOS(DLmuzero)*180.D0/PI)
        
          ! aerosol_formation = F ! to use T option, change file to chem

          CALL SSH_Kinetic(Zangzen,RO2)                   
                                                                        
            ! keep inorganic constant                                   
        if (ncst_gas.gt.0) then 
            do i1=1,ncst_gas 
               gas_yield(cst_gas_index(i1))=ZC_cst(i1)
               ZCtot(cst_gas_index(i1))=ZC_cst(i1) 
            enddo 
        endif 

        call ssh_rates()  
        call ssh_fexprod()                          
        call ssh_dratedc()
        call ssh_fexloss()
                                                                        
        do Jsp=1,n_gas 
          if (chem_prod(Jsp)>dzero.or.chem_loss(Jsp)>dzero) then 
            !init        
            ratloss=chem_loss(Jsp)*ratio_gas(Jsp) 
            ZCtot(Jsp)=(conci(Jsp)+tstep*chem_prod(Jsp))                    &
     &                 /(dun+tstep*ratloss)                             
            gas_yield(Jsp)=ratio_gas(Jsp)*ZCtot(Jsp) 
              ! clip                                                    
            if(gas_yield(Jsp)<dzero) gas_yield(Jsp)=dzero 
            if(ZCtot(Jsp)<dzero) ZCtot(Jsp)=dzero 
          endif 
        enddo 
      enddo 
                                                                        
      ! keep inorganic constant                                         
      if (ncst_gas.gt.0) then 
        do i1=1,ncst_gas 
           gas_yield(cst_gas_index(i1))=ZC_cst(i1)
           ZCtot(cst_gas_index(i1))=ZC_cst(i1) 
        enddo 
      endif 
                                                                        
      ! calculate error                                                 
      error_max=0.d0 
      do Jsp=1,n_gas ! genoa keep_gp         
         wk=atol+rtol*ZCtot(Jsp)*ratio_gas(Jsp) 
         error_max=max(error_max,                                       &
     &        ratio_gas(Jsp)*abs(ZCtot(Jsp)-conci(Jsp))/wk)             
      enddo 
                                                                        
      ! set the current time                                            
      tschem= tschem + tstep 
      ! save timestep                                                   
      dtnsave=tstep 
                                                                        
      ! assgin the conc. after 1st order calculation                    
      do Jsp=1,n_gas 
         conci(Jsp)=ZCtot(Jsp) 
      enddo 
                                                                        
      ! pour evider diviser par zero                                    
      if(error_max>dzero) then 
        tstep=max(tstep_min,max(dun/alpha,min(alpha,0.8d0/              &
     &             (error_max**0.5d0)))*tstep)                          
      else 
        tstep=alpha*tstep 
      endif 
                                                                        
      tstep=min(tstep,tfchem-tschem) 
      if (tstep.gt.dzero) then 
          c=dtnsave/tstep 
      else 
         c=1 
      endif 
      gam=(c+dun)/(c+2.d0) 
                                                                        
      !les calculs suivants de l'ordre 2                                
      do while (tschem<tfchem) 
        do j=1,m 
            !computes chpr and chlo which are the vectors with all the t
                                                                        
            ! compute RO2 used in chem                                  
            RO2=0.d0 
            if (tag_RO2.ne.0) then 
                ! add background tag = 2,3,4 (from ro2 file)            
                if (tag_RO2.ge.2                                        &
     &              .and.iRO2_cst.ne.0) then                            
                   RO2 = RO2 + ZC_cst(iRO2_cst) 
                endif 
                ! add from primary VOCs (RO2 list) tag = 1,3            
                if (tag_RO2 .eq. 1 .or. tag_RO2 .eq. 3                  &
     &              ) then                                              
                   do Jsp =1,nRO2_chem 
                     RO2=RO2+gas_yield(RO2index(Jsp)) 
                   enddo 
                endif 
            endif 
                                                                        
            ! kinetic rate                                              
            ! Compute zenithal angles                                   
            DLmuzero=ssh_muzero(tschem,longitude,latitude) 
            Zangzen=dabs(DACOS(DLmuzero)*180.D0/PI) 

          ! aerosol_formation = F                                       
          CALL SSH_Kinetic(Zangzen,RO2)                   
                                                                        
            ! keep inorganic constant                                   
            if (ncst_gas.gt.0) then 
               do i1=1,ncst_gas 
                  gas_yield(cst_gas_index(i1))=ZC_cst(i1)
                  ZCtot(cst_gas_index(i1))=ZC_cst(i1) 
               enddo 
            endif 
                                                                        
            call ssh_rates()  
            call ssh_fexprod()                          
            call ssh_dratedc()
            call ssh_fexloss() 
                                                                        
            do Jsp=1,n_gas 
               if (chem_prod(Jsp)>dzero.or.chem_loss(Jsp)>dzero) then 
                  ratloss=chem_loss(Jsp)*ratio_gas(Jsp) 
                  ! concc(Jsp)=((c+1)*(c+1)*conci(Jsp)-concii(Jsp))/(c*c
                  ZCtot(Jsp)=(((c+dun)*(c+dun)*conci(Jsp)-              &
     &                  concii(Jsp))/(c*c+2.d0*c)+gam*tstep*            &
     &                 chem_prod(Jsp))/(dun+gam*tstep*ratloss)              
                  gas_yield(Jsp)=ratio_gas(Jsp)*ZCtot(Jsp) 
                  ! clip                                                
                  if(gas_yield(Jsp)<dzero) gas_yield(Jsp)=dzero 
                  if(ZCtot(Jsp)<dzero) ZCtot(Jsp)=dzero 
               endif 
            enddo 
        enddo 
                                                                        
        ! keep inorganic constant                                       
        if (ncst_gas.gt.0) then 
            do i1=1,ncst_gas 
               gas_yield(cst_gas_index(i1))=ZC_cst(i1)
               ZCtot(cst_gas_index(i1))=ZC_cst(i1) 
            enddo 
        endif 
                                                                        
         error_max=0.d0 ! genoa keep_gp                                    
         do Jsp=1,n_gas 
            wk=atol+rtol*abs(ZCtot(Jsp))*ratio_gas(Jsp) 
            error_max=max(error_max,                                    &
     &           abs(2.0d0*(c*ratio_gas(Jsp)*ZCtot(Jsp)-                &
     &           (dun+c)*ratio_gas(Jsp)*conci(Jsp)+                     &
     &           ratio_gas(Jsp)*concii(Jsp))/                           &
     &           (c*(c+dun)*wk)))                                       
         enddo 
                                                                        
        do while (error_max>10.0d0 .and. tstep>tstep_min) 
            tstep=max(tstep_min,max(dun/alpha,min(alpha,                &
     &                0.8d0/(error_max**0.5d0)))*tstep)                 
            tstep=min(tstep,tfchem-tschem) 
            if (tstep.gt.dzero) then 
                c=dtnsave/tstep 
            else 
                c=1 
            endif 
            gam=(c+1)/(c+2.d0) 
            do Jsp=1,n_gas !FCo ! genoa keep_gp               
               ZCtot(Jsp)=conci(Jsp) 
               gas_yield(Jsp)=ratio_gas(Jsp)*ZCtot(Jsp) 
            enddo 
                                                                        
          do j=1,m 
            !computes chpr and chlo which are the vectors with all the t
                                                                        
               ! compute RO2 used in chem                               
               RO2=0.d0 
               if (tag_RO2.ne.0) then 
                   ! add background tag = 2,3,4 (from ro2 file)         
                   if (tag_RO2.ge.2.and.iRO2_cst.ne.0) then                          
                      RO2 = RO2 + ZC_cst(iRO2_cst) 
                   endif 
                   ! add from primary VOCs (RO2 list) tag = 1,3         
                   if (tag_RO2 .eq. 1 .or. tag_RO2 .eq. 3) then                                            
                      do Jsp =1,nRO2_chem 
                        RO2=RO2+gas_yield(RO2index(Jsp)) 
                      enddo 
                   endif 
               endif 
                                                                        
            ! kinetic rate                                              
            ! Compute zenithal angles                                   
            DLmuzero=ssh_muzero(tschem,longitude,latitude) 
            Zangzen=dabs(DACOS(DLmuzero)*180.D0/PI) 

          ! aerosol_formation = F                                       
          CALL SSH_Kinetic(Zangzen,RO2)                   
                                                                        
            ! keep inorganic constant                                   
            if (ncst_gas.gt.0) then 
                do i1=1,ncst_gas 
                   gas_yield(cst_gas_index(i1))=ZC_cst(i1)
                   ZCtot(cst_gas_index(i1))=ZC_cst(i1) 
                enddo 
            endif 
                                                                        
            call ssh_rates()  
            call ssh_fexprod()                          
            call ssh_dratedc()
            call ssh_fexloss() 
                                                                        
            ! genoa keep_gp                              
            do Jsp=1,n_gas 
             if (chem_prod(Jsp)>dzero.or.chem_loss(Jsp)>dzero) then 
                ratloss=chem_loss(Jsp)*ratio_gas(Jsp) 
                                                                        
                ! concc(Jsp)=((c+1)*(c+1)*conci(Jsp)-concii(Jsp))/(
                ZCtot(Jsp)=(((c+dun)*(c+dun)*conci(Jsp)-           &
     &                 concii(Jsp))/(c*c+2.d0*c)+gam*tstep*      &
     &                 chem_prod(Jsp))/(dun+gam*tstep*ratloss)           
                gas_yield(Jsp)=ratio_gas(Jsp)*ZCtot(Jsp) 
                ! clip                                             
                if(ZCtot(Jsp)<dzero) ZCtot(Jsp)=dzero 
                if(gas_yield(Jsp)<dzero) gas_yield(Jsp)=dzero 
              endif 
            enddo 
          enddo 
                                                                        
          ! keep inorganic constant                                     
          if (ncst_gas.gt.0) then 
            do i1=1,ncst_gas 
               gas_yield(cst_gas_index(i1))=ZC_cst(i1) 
               ZCtot(cst_gas_index(i1))=ZC_cst(i1) 
            enddo 
          endif 
                                                                        
          error_max=0.d0 
          ! genoa keep_gp                                   
          do Jsp=1,n_gas 
               wk=atol+rtol*abs(ZCtot(Jsp))*ratio_gas(Jsp) 
               error_max=max(error_max,                                 &
     &              abs(2.0d0*(c*ZCtot(Jsp)*ratio_gas(Jsp)-             &
     &              (dun+c)*ratio_gas(Jsp)*conci(Jsp)                   &
     &              +ratio_gas(Jsp)*concii(Jsp))/                       &
     &              (c*(c+dun)*wk)))                                    
          enddo 
        enddo 
        !save timestep                                                  
        dtnsave=tstep 
        ! update current time                                           
        tschem=tschem+tstep 
        ! update timestep                                               
        if(error_max>dzero) then 
            tstep=max(tstep_min,max(dun/alpha,min(alpha,                &
     &                0.8d0/(error_max**0.5d0)))*tstep)                 
        else 
            tstep=alpha*tstep 
        endif 
                                                                        
        tstep=min(tstep,tfchem-tschem) 
                                                                        
        if (tstep.gt.dzero) then 
            c=dtnsave/tstep 
        else 
           c=1 
        endif 
        gam=(c+dun)/(c+2.d0) 
                                                                        
        do Jsp=1,n_gas 
            concii(Jsp)=conci(Jsp)
            conci(Jsp)=ZCtot(Jsp) 
            ! keep inorganic constant                                   
            if (ncst_gas.gt.0) then 
               do i1=1,ncst_gas 
                  concii(cst_gas_index(i1))=ZC_cst(i1) 
                  conci(cst_gas_index(i1))=ZC_cst(i1) 
                  gas_yield(cst_gas_index(i1))=ZC_cst(i1) 
                  ZCtot(cst_gas_index(i1))=ZC_cst(i1) 
               enddo 
            endif 
        enddo 
        !print*,'      timestep',tstep,tschem                           
      enddo 
                                                                        
!     two-step solver end

!     Storage in the array of chemical concentrations.                  
                                                                        
       DO Jsp=1,n_gas 
          ! NAN detection algorithm used                        
          if (gas_yield(Jsp).ne.gas_yield(Jsp)) then 
             write(*,*) "From chem 0D function:" 
             write(*,*) Jsp, concentration_gas_all(Jsp),conversionfactor(Jsp),   &
&                          gas_yield(Jsp)/conversionfactor(Jsp)             
             stop 
          endif 
          ! Conversion molecules/cm3 to mug/m3.                         
          concentration_gas_all(Jsp) = gas_yield(Jsp)/conversionfactor(Jsp)
       ENDDO 
                                                                        
      if (keep_gp==1) then 
         do s = 1, n_aerosol 
            Jsp=aerosol_species_interact(s) 
            if (Jsp>0) then 
           !if(ratio_gas(Jsp)<1.d0) then                            
              conc_tot=0.0D0 
              DO Jb=1,n_size !*conversionfactor
                 conc_tot=conc_tot+concentration_mass(Jb,s) 
              ENDDO 
                                                                    
              toadd=(ZCtot(Jsp)-ZCtot_save(Jsp))/conversionfactor(Jsp) 
              if (DLconc_save(Jsp)+toadd>=0.d0) then 
                 concentration_gas_all(Jsp)=DLconc_save(Jsp)+toadd 
                 !Not enough mass in the gas phase, have 
              else 
!                     print*,toadd,DLconc_save(Jsp)+conc_tot,           
!     s                    ZCtot(Jsp)/conversionfactor(Jsp),              
!     s                    ZCtot_save(Jsp)/conversionfactor(Jsp)          
                 toadd=(conc_tot+toadd+DLconc_save(Jsp))/conc_tot 
                 concentration_gas_all(Jsp)=0.d0 
                 DO Jb=1,n_size 
                    concentration_mass(Jb,s) = concentration_mass(Jb,s)*toadd 
                 ENDDO 
              endif 
                                                                    
              conc_tot=0.d0 
              DO Jb=1,n_size 
                 conc_tot=conc_tot+concentration_mass(Jb,s) 
              ENDDO 
!                  print*,conc_tot+concentration_gas_all(Jsp)                          
!     s                 ,ZCtot(Jsp)/conversionfactor(Jsp)                 
           !endif                                                   
            endif 
         ENDDO 
      endif 
                                                                        
      ! keep inorganic constant                                         
      if (ncst_gas.gt.0) then 
         do i1=1,ncst_gas 
            concentration_gas_all(cst_gas_index(i1))=   &
     &                gas_yield(cst_gas_index(i1))/conversionfactor(cst_gas_index(i1))               
         enddo 
      endif 
                                                                        
      ! output RO2                                                      
      if (iRO2.ne.0.and.tag_RO2.ne.0) then 
            ! output only ro2 from list                                 
            concentration_gas_all(iRO2)=0.d0 
            do Jsp =1,nRO2_chem 
               concentration_gas_all(iRO2)= concentration_gas_all(iRO2)+ &
     &                  gas_yield(RO2index(Jsp))/conversionfactor(RO2index(Jsp))                   
            enddo 
      endif 
                                                                        
      DO Jb=1,n_size 
         DO Jsp=1,n_aerosol 
            if (concentration_mass(jb,jsp) .ne. concentration_mass(jb,jsp)) then 
               write(*,*) "From chem 0D function (aerosol conc.):" 
               write(*,*) jb,jsp,concentration_mass(jb,jsp) 
               stop 
            endif 
         enddo 
      enddo
  END SUBROUTINE ssh_chem_twostep 
END module mod_sshchem                                 
