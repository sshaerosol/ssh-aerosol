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
      
      ! use humidity, temperature, pressure at this timestep
      
      ! no emission
      ! use: longitude, latitude                               
      DOUBLE PRECISION Zangzen !solar zenith angle
      DOUBLE PRECISION ssh_muzero,DLmuzero 
      EXTERNAL ssh_muzero 

!     Parameters initialized for the two-step solver                    
      integer m,j,i1,Jsp,Jb,s, nstep
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
      ! use: ncst_gas,nRO2_chem,iRO2_cst,tag_RO2 
      ! use: cst_gas_index(ncst_gas), RO2index(nRO2_chem,2), cst_gas
      ! use: nRO2_group, RO2out_index
      DOUBLE PRECISION ZC_cst(ncst_gas), RO2s(nRO2_group)
                                                                        
!     genoa keep_gp: gas-particle partitioning                          
      ! used for genoa with large timestep           
      !keep_gp, aerosol_species_interact(n_gas) 
      DOUBLE PRECISION toadd, conc_tot, ZCtot_save(n_gas) 
      DOUBLE PRECISION ratio_gas(n_gas),ZCtot(n_gas),DLconc_save(n_gas) 
                                                                        
      ! use: photo_rcn(:,:), TB_rcn(:,:) 
      ! use: fall_rcn(:), extra_rcn(:) 
      ! use: index_RCT(:,:), index_PDT(:,:) 
      ! kinetics    
      ! SumMc, YlH2O                                                    
      ! use: Arrhenius(:,:), fall_coeff(:,:), ratio_PDT(:) , photo_ratio(:,:,:) 

!     Projection.                                                       
!     Conversion mug/m3 to molecules/cm3.
      gas_yield = 0.d0 !init

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
         if (tag_RO2.eq.1.or.tag_RO2.eq.3) then
            call ssh_compute_ro2(RO2s)
         else
            RO2s = 0d0
         endif
         ! add background RO2
         if (tag_RO2.eq.2.or.tag_RO2.eq.3) then
            RO2s(1) = RO2s(1) + ZC_cst(iRO2_cst) 
         endif
         
        ! kinetic rate                                                  
        ! Compute zenithal angles                                       
        DLmuzero=ssh_muzero(tschem,longitude,latitude) 
        Zangzen=dabs(DACOS(DLmuzero)*180.D0/PI)
        
          ! aerosol_formation = F ! to use T option, change file to chem

          CALL SSH_Kinetic(Zangzen,RO2s)                   
                                                                        
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
            if (tag_RO2.eq.1.or.tag_RO2.eq.3) then
                call ssh_compute_ro2(RO2s)
            else
                RO2s = 0d0
            endif
            ! add background RO2
            if (tag_RO2.eq.2.or.tag_RO2.eq.3) then
                RO2s(1) = RO2s(1) + ZC_cst(iRO2_cst) 
            endif
                                                                        
            ! kinetic rate                                              
            ! Compute zenithal angles                                   
            DLmuzero=ssh_muzero(tschem,longitude,latitude) 
            Zangzen=dabs(DACOS(DLmuzero)*180.D0/PI) 

          ! aerosol_formation = F                                       
          CALL SSH_Kinetic(Zangzen,RO2s)                   
                                                                        
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
            if (tag_RO2.eq.1.or.tag_RO2.eq.3) then
                call ssh_compute_ro2(RO2s)
            else
                RO2s = 0d0
            endif
            ! add background RO2
            if (tag_RO2.eq.2.or.tag_RO2.eq.3) then
                RO2s(1) = RO2s(1) + ZC_cst(iRO2_cst) 
            endif
                                                                        
            ! kinetic rate                                              
            ! Compute zenithal angles                                   
            DLmuzero=ssh_muzero(tschem,longitude,latitude) 
            Zangzen=dabs(DACOS(DLmuzero)*180.D0/PI) 

          ! aerosol_formation = F                                       
          CALL SSH_Kinetic(Zangzen,RO2s)                   
                                                                        
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

      DO Jb=1,n_size 
         DO Jsp=1,n_aerosol 
            if (concentration_mass(jb,jsp) .ne. concentration_mass(jb,jsp)) then 
               write(*,*) "From chem 0D function (aerosol conc.):" 
               write(*,*) jb,jsp,concentration_mass(jb,jsp) 
               stop 
            endif 
         enddo 
      enddo
      
      endif 
                                                                        
      ! keep inorganic constant                                         
      if (ncst_gas.gt.0) then 
         do i1=1,ncst_gas 
            concentration_gas_all(cst_gas_index(i1))= cst_gas(i1,nstep)           
         enddo 
      endif 
                                                                        
      ! output RO2s
      if (tag_RO2.ne.0) then
        ! init
        do Jb = 1, nRO2_group
            i1 = RO2out_index(Jb) ! output id
            concentration_gas_all(i1)=0.d0 
        enddo    
        
        do Jsp =1,nRO2_chem
          Jb = RO2index(Jsp,1)! isps
          i1 = RO2index(Jsp,2)! igroup
          concentration_gas_all(i1)= concentration_gas_all(i1)+ &
     &                           concentration_gas_all(Jb)                  
        enddo
      endif
                                                                        
  END SUBROUTINE ssh_chem_twostep

END module mod_sshchem                                 
