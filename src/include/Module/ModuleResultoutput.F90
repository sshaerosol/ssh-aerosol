!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------

Module Resultoutput
  use aInitialization
  use dPhysicalbalance
  USE netcdf
  
  implicit none

contains

!============================================================================================
! could be usefull for ssh_write_output_netcdf

  subroutine handle_err(status)
    integer, intent(in) :: status
    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      !stop "Stopped"
    endif
  end subroutine handle_err


!============================================================================================

   subroutine ssh_save_report()
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This subroutine write report file "report.txt", which records most
!     settings, physical conditions and options adopted in the simulation.
!
!     File "report.txt" is saved in the directory : output_directory/
!     provided by user (namelist.ssh).
!
!------------------------------------------------------------------------
!
!     -- OUTPUT 
!     "report.txt"
!
!------------------------------------------------------------------------

     implicit none
     integer :: j,stat
     logical :: file_exists
     character (len=80) :: cmd

     !Create output directory
     cmd = trim('mkdir -p '// trim(output_directory))
     call system(cmd)

    ! delete the old report file if it exists under the current saving directory
    inquire (file = trim(output_directory) // "/" // "report.txt", exist = file_exists)
    if (file_exists) then
            open(unit=10, file = trim(output_directory) // "/" // "report.txt", status='old', iostat=stat)
            if (stat == 0) close(10, status='delete')
    endif
    ! write the new report
    open(unit=10,file=trim(output_directory) // "/" // "report.txt", status="new")

	write(unit=10,FMT=*) '<<<< Meteorological setup >>>>'
        write(unit=10,FMT=*) 'location', latitude, 'N','	', longitude,'E','	','Temperature', Temperature, 'K'
        write(unit=10,FMT=*) 'Pressure', Pressure, 'Pa', '	','Specific Humidity', Humidity, '	',&
				'Cloud attenuation field', attenuation, '	relative humidity', relative_humidity
	write(unit=10,FMT=*)
	write(unit=10,FMT=*) '<<<< Simulation time setup >>>>'
        write(unit=10,FMT=*) 'Begining time (from Jan. 1st)', initial_time, 's','	',&
				'Simulation Time', final_time,'s','	','Initial Time Step', delta_t,'s',&
         			'	','Number of iterations:', nt
	write(unit=10,FMT=*)
	write(unit=10,FMT=*) '<<<< Inition condition >>>>'
	if (tag_init == 0) write(unit=10,FMT=*) &
			'Internally mixed aerosol species are provided for initial condition.'
	write(unit=10,FMT=*)  'Gas-phase conc. input file :', init_gas_conc_file
	write(unit=10,FMT=*)  'Particle conc. input file :', init_aero_conc_mass_file
	write(unit=10,FMT=*)  'N_sizebin', N_sizebin
        if (with_init_num .eq. 1) write(unit=10,FMT=*)  'Aerosol number conc. is read from file :',&
							 init_aero_conc_num_file
	if (with_init_num .eq. 0) then
		write(unit=10,FMT=*) ' Aerosol number conc. is estimated from mass and diameter.' 
		write(unit=10,FMT=*)  '====== concentration_number : ======'
		write(unit=10,FMT=*)  concentration_number
	end if
	write(unit=10,FMT=*)
	write(unit=10,FMT=*) '<<<< Mixing state >>>>'
	if (tag_external == 1) write(unit=10,FMT=*)  'simulation is mixing-state resolved.', &
					'N_frac', N_frac ,'	','N_groups', N_groups
	if (tag_external == 0) write(unit=10,FMT=*)  'simulation is internally mixed.' ,&
					'N_frac', N_frac ,'	','N_groups', N_groups
	write(unit=10,FMT=*) 'frac_bound', frac_bound
	write(unit=10,FMT=*)
	write(unit=10,FMT=*) '<<<< Species lists >>>>'
	write(unit=10,FMT=*)  'gas phase species file :', species_list_file
	write(unit=10,FMT=*)  'particle species file :', aerosol_species_list_file
	write(unit=10,FMT=*)
	if (tag_chem == 0) write(unit=10,FMT=*) '<<<< Without Gas-phase chemistry >>>>'
	if (tag_chem == 1) then
		write(unit=10,FMT=*) '<<<< Gas-phase chemistry >>>>'
		write(unit=10,FMT=*)  'with_heterogeneous', with_heterogeneous,'	','with_adaptive', with_adaptive,&
					'	','adaptive time step tolerance', adaptive_time_step_tolerance,&
					'	','min adaptive time step', min_adaptive_time_step
	end if
	write(unit=10,FMT=*)
	write(unit=10,FMT=*) '<<<< Emissions >>>>'
	if (tag_emis == 1) then
		write(unit=10,FMT=*)  'With internally-mixed emissions.'
		write(unit=10,FMT=*)  'Gas-phase conc. emission file :', emis_gas_file
		write(unit=10,FMT=*)  'Particle conc. emission file :', emis_aero_mass_file
		if (with_emis_num == 1) write(unit=10,FMT=*) &
			'Emitted aerosol number conc. is read from file :', emis_aero_num_file
		if (with_emis_num == 0) write(unit=10,FMT=*) &
			'Emitted aerosol number conc. is estimated from mass and diameter.'
	end if
        if (tag_emis == 0) write(unit=10,FMT=*) 'Without emission.'
	write(unit=10,FMT=*)
	write(unit=10,FMT=*) '<<<< Particle Dynamic >>>>'
	if (with_cond == 1) then 
		write(unit=10,FMT=*)  'With condensation', '	','Cut_dim', Cut_dim,'	', 'ISOAPDYN', ISOAPDYN
	else
		write(unit=10,FMT=*)  'Without condensation'
	end if
	if (with_coag == 1) then 
		write(unit=10,FMT=*)  'With coagulation', '	','coefficient file :', Coefficient_file
	else
		write(unit=10,FMT=*)  'Without coagulation'
	end if
	if (with_nucl == 1) then 
		write(unit=10,FMT=*)  'With nucleation', '	','nucl_model_binary', nucl_model_binary
		write(unit=10,FMT=*)  'With nucleation', '	','nucl_model_ternary', nucl_model_ternary
		write(unit=10,FMT=*)  'With nucleation', '	','nucl_model_hetero', nucl_model_hetero
	else
		write(unit=10,FMT=*)  'Without nucleation'
	end if
	write(unit=10,FMT=*)  'DTAEROMIN', DTAEROMIN, '	','redistribution_method', redistribution_method
	write(unit=10,FMT=*)  'with_fixed_density', with_fixed_density,'	', 'fixed_density', fixed_density, 'kg/m^3'
	write(unit=10,FMT=*)
	write(unit=10,FMT=*) 'output directory :', trim(output_directory),'/'

   CLOSE(10)

   end subroutine ssh_save_report


!============================================================================================


  subroutine ssh_init_output()
    
    implicit none
    integer :: i,j,s,b
    
    allocate(output_time(nt+1))                  
    allocate(output_gas(nt+1,N_gas))            
    allocate(output_numb(nt+1,N_size+1))         
    allocate(output_diam(nt+1,N_size))           
    allocate(output_aero(nt+1,N_aerosol,N_size)) 
    allocate(output_TM(nt+1,N_aerosol,3))        
    allocate(output_special(nt+1,8))
    allocate(output_pH(nt+1,N_size))      
    
    t_out = 0 
    mixing_nb = N_size / N_sizebin ! 1 if internal, >1 if external

    do b=1, N_sizebin+1
      if (diam_bound(b) .le. 1.0d0 )  ipm1=b-1  ! ipm1 = last size bin number of PM 1 
      if (diam_bound(b) .le. 2.5d0 )  ipm25=b-1 ! ipm25 = last size bin number of PM 2.5 
      if (diam_bound(b) .le. 10.0d0 ) ipm10=b-1 ! ipm10 = last size bin number of PM 10 
    enddo

    ! initialization
    do i=1, nt+1
      output_time(i) = 0.0d0
      do s=1, n_gas
        output_gas(i,s) =  0.0d0
      enddo
      do b=1, N_size+1
        output_numb(i,b) = 0.0d0
      enddo
      do b=1, N_size
        output_diam(i,b) = 0.0d0
        output_pH(i,b) = 7.0d0
      enddo
      do s=1, N_aerosol
        do j=1, 3
          output_TM(i,s,j) = 0.0d0
        enddo
        do b=1, N_size
          output_aero(i,s,b) = 0.0d0
        enddo
      enddo
      do j=1, 8
        output_special(i,j) = 0.0d0
      enddo
    enddo
    
    ! maybe removal of previous reseults files should be called here?
    
  end subroutine ssh_init_output


  ! add for nout_gas and nout_aero
  subroutine ssh_init_output_genoa()
    
    implicit none
    integer :: i,j,s,b
    
    ! timestep for output
    t_out = 0
    
    ! save gas species
    if (nout_gas .gt. 0) then                  
        allocate(output_gas(nt+1,nout_gas))
        output_gas = 0d0
    endif
    
    ! save aerosol species
    if (nout_aero .gt. 0) then                  
        allocate(output_aero(nt+1,nout_aero,N_size))
        output_aero = 0d0
    endif
    
  end subroutine ssh_init_output_genoa

!============================================================================================

   
  subroutine ssh_save_concentration()
  
    implicit none
    integer :: s, b, i, j, jesp
    
    t_out = t_out+1
    output_time(t_out) = initial_time + (t_out - 1) * delta_t

    do s = 1, n_gas
      output_gas(t_out,s) =  concentration_gas_all(s) 
    enddo

    do b = 1, N_size   
      output_numb(t_out,b)        = concentration_number(b)
      output_numb(t_out,N_size+1) = output_numb(t_out,N_size+1) + concentration_number(b)
      output_diam(t_out,b)        = cell_diam_av(b)
      if (lwc_Nsize(b)>0.d0.and.proton_Nsize(b)>0.d0) then          
        output_pH(t_out,b)=-log10(proton_Nsize(b)/lwc_Nsize(b)*1.0e3)   
      endif
    enddo
     
    do s = 1, N_aerosol
      total_aero_mass(s) = 0.d0 
      do b = 1, N_size
        concentration_mass_tmp(b,s) = 0.d0
      enddo
    enddo

    do b = 1, N_size
      do s = 1, N_aerosol_layers
        jesp = List_species(s)
        total_aero_mass(jesp) = total_aero_mass(jesp) + concentration_mass(b,s)
        concentration_mass_tmp(b,jesp) = concentration_mass_tmp(b,jesp) + concentration_mass(b,s)
      enddo
    enddo
    
    do s = 1, N_aerosol
      do b = 1, N_size
        output_aero(t_out,s,b) = concentration_mass_tmp(b,s)
      enddo
      output_TM(t_out,s,1) = total_aero_mass(s)
      if (aerosol_species_interact(s) .gt. 0) then
        output_TM(t_out,s,2) = concentration_gas_all(aerosol_species_interact(s))
      endif
      output_TM(t_out,s,3)  = output_TM(t_out,s,1) + output_TM(t_out,s,2)
    enddo
    

    ! save organic, inorganic, PBC, Dust, PM2.5, PM10 results over each time step
    do s=1, N_aerosol
      if (aerosol_type(s).EQ.1) then 
        output_special(t_out,1) = total_aero_mass(s)                       !Mineral Dust
      elseif (aerosol_type(s).EQ.2) then
        output_special(t_out,2) = total_aero_mass(s)                       !Black Carbon
      elseif (aerosol_type(s).EQ.9) then
        output_special(t_out,3) = total_aero_mass(s)                       !Water
      elseif (aerosol_type(s).EQ.3) then 
        output_special(t_out,4) = output_special(t_out,4) + total_aero_mass(s) !Inorganics
      elseif (aerosol_type(s).EQ.4) then 
        output_special(t_out,5) = output_special(t_out,5) + total_aero_mass(s) !Organics
      endif
    enddo
   
    do s=1, N_aerosol-1	! without water
      if (ipm1 .GE. 1) then ! PM1
        do b=1, ipm1*mixing_nb
          output_special(t_out,6) = output_special(t_out,6) + concentration_mass_tmp(b,s)
        enddo
      endif
      if (ipm25 .GE. 1) then ! PM2.5
        do b=1, ipm25*mixing_nb
          output_special(t_out,7) = output_special(t_out,7) + concentration_mass_tmp(b,s)
        enddo
      endif
      if (ipm10 .GE. 1) then ! PM10
        do b=1, ipm25*mixing_nb
          output_special(t_out,8) = output_special(t_out,8) + concentration_mass_tmp(b,s)
        enddo
      endif
    enddo
      

  end subroutine ssh_save_concentration
       
  subroutine ssh_save_concentration_genoa()
  
    implicit none
    integer :: s, b, i, j, jesp, t
    
    ! update timestep - default 0
    t_out = t_out + 1
    
    ! save gas species concs
    do i = 1, nout_gas
      s = output_gas_index(i)
      output_gas(t_out,i) = concentration_gas_all(s) 
    enddo
    
    ! update aerosol concs at this timestep
    concentration_mass_tmp = 0d0   
    do b = 1, N_size
      do s = 1, N_aerosol_layers
        jesp = List_species(s)
        concentration_mass_tmp(b,jesp) = concentration_mass_tmp(b,jesp) + concentration_mass(b,s)
      enddo
    enddo
    
    ! save aerosol species concs
    do i = 1, nout_aero
      s = output_aero_index(i)
      do b = 1, N_size
         output_aero(t_out,i,b) = concentration_mass_tmp(b, s)
      enddo
    enddo

  end subroutine ssh_save_concentration_genoa
  
  subroutine ssh_save_total_soa_genoa()

    implicit none
    integer :: s, b, i, j, jesp, t
    
    ! time step has been updated in ssh_save_concentration()
    
    ! save gas species concs
    do i = 1, nout_gas
      s = output_gas_index(i)
      output_gas(t_out,i) = concentration_gas_all(s) 
    enddo
    
    ! update aerosol concs at this timestep
    total_aero_mass = 0d0 
    do b = 1, N_size
      do s = 1, N_aerosol_layers
        jesp = List_species(s)
        total_aero_mass(jesp) = total_aero_mass(jesp) + concentration_mass(b,s)
      enddo
    enddo
    
  ! save species for error computation
  ! index in total_soa: 1: total SOA + nout_err + nout_soa
  
  ! compute total SOA - index 1
  total_soa(t_out,1) = 0d0
  do s = 1, N_aerosol ! remove water and no-organics
    !if (aerosol_species_name(s).eq.'PBiMT') cycle
    !if (aerosol_species_name(s).eq.'PSOAlP') cycle
    ! only keep organics from primary vocs
    if (aerosol_type(s).ne.4.or.Index_groups(s).le.1) cycle
    total_soa(t_out,1) = total_soa(t_out,1) + total_aero_mass(s)
  end do
      
    ! index 2 - nout_err + 1
    ! save species for error computation
    do i = 1, nout_err
        ! check gas
        s = output_err_index(1,i) ! index for gas
     if (s.ne.0) then
        total_soa(t_out,1+i) = concentration_gas_all(s)
     else ! check aerosol
        s = output_err_index(2,i) ! index for aerosol
        if (s.ne.0) then
          total_soa(t_out,1+i) = total_aero_mass(s)
        else ! not found index - stop
          print*, "Output_err_index not found in both gas and aerosol list."
          print*, "Check err_species_list in the namelist."
          print*,i, output_err_index(1,i),output_err_index(2,i)
          stop
        endif
     endif
    enddo
          
    ! save group concentrations for aerosols
    do i = 1, nout_soa
      ! compute total SOA
      total_soa(t_out,1+nout_err+i) = 0d0
      
      do s = 1, N_aerosol ! remove water and no-organics
        if (aerosol_type(s).ne.4) cycle ! only keep organics from primary vocs
        if (Index_groups(s).ne.i+1) cycle ! sum concs from same group - example
        
        total_soa(t_out,1+nout_err+i) = total_soa(t_out,1+nout_err+i) + total_aero_mass(s)
      end do
    end do

    ! round values if need
    if (iout_round) then
        do i = 1, nout_total
            if (total_soa(t_out,i).lt.1d7) then
                j = anint(total_soa(t_out,i)*1d6)
                total_soa(t_out,i) = j/1d6
            end if
        enddo
    end if
    
  end subroutine ssh_save_total_soa_genoa
  
!============================================================================================

  SUBROUTINE ssh_write_output()
    IMPLICIT NONE
      
      if (ssh_standalone) write(*,*) ""
      if (ssh_logger) write(logfile,*) ""
      
      if (output_type == 3) then
        if (ssh_standalone) write(*,*) "====  Write Netcdf output files  ====="
        if (ssh_logger) write(logfile,*) "====  Write Netcdf output files  ====="
        call ssh_write_output_netcdf()
      else if (output_type .ne. 0) then
        if (ssh_standalone) write(*,*) "====  Write bin/txt output files  ====="
        if (ssh_logger) write(logfile,*) "====  Write bin/txt output files  ====="
        call ssh_write_output_bintxt()
      end if
      
      if (ssh_standalone) write(*,*) ""
      if (ssh_logger) write(logfile,*) ""
      
  END SUBROUTINE ssh_write_output
  

!============================================================================================



  SUBROUTINE ssh_write_output_netcdf()
    USE netcdf
    IMPLICIT NONE
    
    !file ID
    INTEGER :: ncid
    
    !dimension IDs & vardim IDs
    INTEGER :: tid, tvid  !temps
    INTEGER :: pid, pvid  !layer
    INTEGER :: sid, svid  !size
    INTEGER :: lbid, hbid, emid, sbid  !lower and higher bound

    !variable IDs
    INTEGER :: concid(5000) !nb mx of var
  
    !local var
    INTEGER :: s,b,i,j,k,m,isizmix 
    INTEGER :: ncstat
    CHARACTER(20) :: out_aero(8) 
    DOUBLE PRECISION :: lb(N_size), hb(N_size)
    INTEGER :: sb(N_size), em(N_size)

    
    out_aero(1) = 'Dust'
    out_aero(2) = 'Black_Carbon'
    out_aero(3) = 'Water'
    out_aero(4) = 'Inorganic'
    out_aero(5) = 'Organic'
    out_aero(6) = 'PM1'
    out_aero(7) = 'PM2.5'
    out_aero(8) = 'PM10'
    
!   to add after "ncstat line" to have information on potential error: ";call handle_err(ncstat)"
    
!=======
! gas file
    ncstat = nf90_create(trim(output_directory) // "gas.nc", NF90_NETCDF4, ncid)
    ncstat = nf90_put_att(ncid, NF90_GLOBAL, 'Title', 'SSH-aerosol ouput file: concentrations of gaseous species')  
    ncstat = nf90_def_dim(ncid,"Time",nt+1,tid)
    ncstat = nf90_def_var(ncid,"Time", NF90_DOUBLE,  (/ tid /), tvid)
    ncstat = nf90_put_att(ncid, tvid, 'units', 'seconds')    
    do s = 1, n_gas
      ncstat = nf90_def_var(ncid,trim(species_name(s)), NF90_DOUBLE,  (/ tid /), concid(s))
      ncstat = nf90_put_att(ncid, concid(s), 'units', 'microgram_m-3')
    enddo
    
    ncstat = nf90_put_var(ncid, tvid, output_time(1:nt+1))
    do s = 1, n_gas
      ncstat = nf90_put_var(ncid, concid(s), output_gas(1:nt+1,s))
    enddo
    ncstat = nf90_close(ncid)
!=======
   
    
!=======
! number file
    ncstat = nf90_create(trim(output_directory) // "number.nc", NF90_NETCDF4, ncid)
    ncstat = nf90_put_att(ncid, NF90_GLOBAL, 'Title', 'SSH-aerosol ouput file: numbers of particles')  
    ncstat = nf90_def_dim(ncid,"Time",nt+1,tid)
    ncstat = nf90_def_var(ncid,"Time", NF90_DOUBLE,  (/ tid /), tvid)
    ncstat = nf90_put_att(ncid, tvid, 'units', 'seconds')        
    isizmix = 0
    do b=1, N_sizebin
      do m=1, mixing_nb
        isizmix = isizmix+1
        ncstat = nf90_def_var(ncid,"sizemixbin_" // trim(str(isizmix)), NF90_DOUBLE,  (/ tid /), concid(isizmix))
        ncstat = nf90_put_att(ncid, concid(isizmix), 'units', 'particles_m-3')
        ncstat = nf90_put_att(ncid, concid(isizmix), 'lower_boundary', diam_bound(b))
        ncstat = nf90_put_att(ncid, concid(isizmix), 'higher_boundary', diam_bound(b+1))
        ncstat = nf90_put_att(ncid, concid(isizmix), 'ext_mix', m)
      enddo
    enddo
    if (isizmix .ne. N_size) write(*,*) 'ERROR in Netcdf writing: isizmix and N_size are different'
    ncstat = nf90_def_var(ncid,"total", NF90_DOUBLE,  (/ tid /), concid(N_size+1))
    ncstat = nf90_put_att(ncid, concid(N_size+1), 'units', 'particles_m-3')
    ncstat = nf90_put_att(ncid, concid(N_size+1), 'lower_boundary', diam_bound(1))
    ncstat = nf90_put_att(ncid, concid(N_size+1), 'higher_boundary', diam_bound(N_sizebin+1))
    
    ncstat = nf90_put_var(ncid, tvid, output_time(1:nt+1))
    do b = 1, N_size+1
      ncstat = nf90_put_var(ncid, concid(b), output_numb(1:nt+1,b))
    enddo
    ncstat = nf90_close(ncid)
!=======
    
    
!=======
! diameter file
    ncstat = nf90_create(trim(output_directory) // "diameter.nc", NF90_NETCDF4, ncid)
    ncstat = nf90_put_att(ncid, NF90_GLOBAL, 'Title', 'SSH-aerosol ouput file: average diameters of particles')  
    ncstat = nf90_def_dim(ncid,"Time",nt+1,tid)
    ncstat = nf90_def_var(ncid,"Time", NF90_DOUBLE,  (/ tid /), tvid)
    ncstat = nf90_put_att(ncid, tvid, 'units', 'seconds')    
    isizmix = 0
    do b=1, N_sizebin
      do m=1, mixing_nb
        isizmix = isizmix+1
        ncstat = nf90_def_var(ncid,"sizebin_" // trim(str(isizmix)), NF90_DOUBLE,  (/ tid /), concid(isizmix))
        ncstat = nf90_put_att(ncid, concid(isizmix), 'units', 'micrometer')
        ncstat = nf90_put_att(ncid, concid(isizmix), 'lower_boundary', diam_bound(b))
        ncstat = nf90_put_att(ncid, concid(isizmix), 'higher_boundary', diam_bound(b+1))
        ncstat = nf90_put_att(ncid, concid(isizmix), 'ext_mix', m)
      enddo
    enddo
    if (isizmix .ne. N_size) write(*,*) 'ERROR in Netcdf writing: isizmix and N_size are different'
    
    ncstat = nf90_put_var(ncid, tvid, output_time(1:nt+1))
    do b = 1, N_size
      ncstat = nf90_put_var(ncid, concid(b), output_diam(1:nt+1,b))
    enddo
    ncstat = nf90_close(ncid)
!=======
   
    
!=======
! TM file    
    ncstat = nf90_create(trim(output_directory) // "TM.nc", NF90_NETCDF4, ncid)
    ncstat = nf90_put_att(ncid, NF90_GLOBAL, 'Title', 'SSH-aerosol ouput file: concentrations of condensable species')  
    ncstat = nf90_def_dim(ncid,"Time",nt+1,tid)
    ncstat = nf90_def_dim(ncid,"Phase",3,pid)
    ncstat = nf90_def_var(ncid,"Time", NF90_DOUBLE,  (/ tid /), tvid)
    ncstat = nf90_put_att(ncid, tvid, 'units', 'seconds')
    ncstat = nf90_def_var(ncid,"Phase", NF90_CHAR,  (/ pid /), pvid)
    ncstat = nf90_put_att(ncid, pvid, 'phase', '1=aer 2=gas 3=tot')
    do s = 1, N_aerosol
      ncstat = nf90_def_var(ncid,aerosol_species_name(s), NF90_DOUBLE,  (/ tid , pid /), concid(s))
      ncstat = nf90_put_att(ncid, concid(s), 'units', 'microgram_m-3')
    enddo

    ncstat = nf90_put_var(ncid, tvid, output_time(1:nt+1))
    do s = 1, N_aerosol
      ncstat = nf90_put_var(ncid, concid(s), output_TM(1:nt+1,s,1:3))
    enddo
    ncstat = nf90_close(ncid)
!=======


!=======
! aero file 
    ncstat = nf90_create(trim(output_directory) // "aero.nc", NF90_NETCDF4, ncid)
    ncstat = nf90_put_att(ncid, NF90_GLOBAL, 'Title', 'SSH-aerosol ouput file: concentrations of aerosol species')  
    ncstat = nf90_def_dim(ncid,"Time",nt+1,tid)
    ncstat = nf90_def_dim(ncid,"Size_mix",N_size,sid)
    ncstat = nf90_def_var(ncid,"Time", NF90_DOUBLE,  (/ tid /), tvid)
    ncstat = nf90_put_att(ncid, tvid, 'units', 'seconds')
    
    ncstat = nf90_def_var(ncid,"lower_boundary", NF90_DOUBLE,  (/ sid /), lbid)
    ncstat = nf90_put_att(ncid, lbid, 'units', 'micrometers')
    ncstat = nf90_def_var(ncid,"higher_boundary", NF90_DOUBLE,  (/ sid /), hbid)
    ncstat = nf90_put_att(ncid, hbid, 'units', 'micrometers')
    ncstat = nf90_def_var(ncid,"Sizebin", NF90_INT,  (/ sid /), sbid)
    ncstat = nf90_put_att(ncid, sbid, 'title', 'Sizebin indice')
    ncstat = nf90_def_var(ncid,"ext_mix", NF90_INT,  (/ sid /), emid)      
    ncstat = nf90_put_att(ncid, emid, 'title', 'external mixing indice')
    
    isizmix = 0
    do b=1, N_sizebin
      do m=1, mixing_nb
        isizmix = isizmix+1
        lb(isizmix) = diam_bound(b)
        hb(isizmix) = diam_bound(b+1)
        sb(isizmix) = b
        em(isizmix) = m
      enddo
    enddo         
    if (isizmix .ne. N_size) write(*,*) 'ERROR in Netcdf writing: isizmix and N_size are different'
    
    do s = 1, N_aerosol
      ncstat = nf90_def_var(ncid,aerosol_species_name(s), NF90_DOUBLE,  (/ tid , sid /), concid(s))
      ncstat = nf90_put_att(ncid, concid(s), 'units', 'microgram_m-3')
    enddo
    do s = 1, 8
      ncstat = nf90_def_var(ncid,out_aero(s), NF90_DOUBLE,  (/ tid /), concid(N_aerosol+s))
      ncstat = nf90_put_att(ncid, concid(N_aerosol+s), 'units', 'microgram_m-3')
    enddo
    ncstat = nf90_def_var(ncid,"pH", NF90_DOUBLE,  (/ tid , sid /), concid(N_aerosol+8+1))
    
    ncstat = nf90_put_var(ncid, lbid, lb)
    ncstat = nf90_put_var(ncid, hbid, hb)
    ncstat = nf90_put_var(ncid, sbid, sb)
    ncstat = nf90_put_var(ncid, emid, em)
    
    ncstat = nf90_put_var(ncid, tvid, output_time(1:nt+1))
    do s = 1, N_aerosol
      ncstat = nf90_put_var(ncid, concid(s), output_aero(1:nt+1,s,1:N_size))
    enddo
    do s = 1, 8
      ncstat = nf90_put_var(ncid, concid(N_aerosol+s), output_special(1:nt+1,s))
    enddo
    ncstat = nf90_put_var(ncid, concid(N_aerosol+8+1), output_pH(1:nt+1,1:N_size))
    ncstat = nf90_close(ncid)
!=======
  
  END SUBROUTINE ssh_write_output_netcdf  
  

!============================================================================================


  SUBROUTINE ssh_write_output_bintxt()
!------------------------------------------------------------------------
!     -- OUTPUTS 
!     Mass concentrations of each gas-phase species:
!     >>> output_directory/gas/species name/".txt"(".bin")
!
!     Mass concentrations of each aerosol species in each grid cell:
!     >>> output_directory/aero/aerosol species name_i/".txt"(".bin"), i = 1,2,3...N_size
!
!     Mass concentrations of organic, inorganic, Black_Carbon, Dust, PM2.5, PM10:
!     >>> output_directory/aero/name/".txt"(".bin")
!
!     Total mass concentrations of each aerosol species:
!     >>> output_directory/TM/aerosol species name_TM/".txt"(".bin")
!
!     Total mass concentrations of each aerosol species + precursor:
!     >>> output_directory/TM/aerosol species name_precursor name_TM/".txt"(".bin")
!
!     Number concentrations of each grid cell :
!     >>> output_directory/number/NUMBER_i/".txt"(".bin"), i = 1,2,3...N_size
!
!     Total number concentrations :
!     >>> output_directory/number/TNUM/".txt"(".bin"), i = 1,2,3...N_size
!
!     Average diameter of each grid cell :
!     >>> output_directory/diameter/DIAMETER_i/".txt"(".bin"), i = 1,2,3...N_size
!
!     -- unitS : 
!     mass concentration [ug/m3] 
!     number concentration [#/m3]
!     diameter [um]
!------------------------------------------------------------------------  
    integer :: stat, s, b, t
    logical :: file_exists
    character (len=100) output_filename
    character (len=80) :: cmd
    character (len=10) :: out_dir(5) 
    
    character(20) :: out_aero(8) 
    character(4) :: out_type(2) = (/".txt",".bin"/) ! 1: text, 2: binary
    
    out_dir(1) = "/number/"
    out_dir(2) = "/gas/"	
    out_dir(3) = "/aero/"
    out_dir(4) = "/TM/"
    out_dir(5) = "/diameter/"
    
    out_aero(1) = 'Dust'
    out_aero(2) = 'Black_Carbon'
    out_aero(3) = 'Water'
    out_aero(4) = 'Inorganic'
    out_aero(5) = 'Organic.5'
    out_aero(6) = 'PM1'
    out_aero(7) = 'PM2.5'
    out_aero(8) = 'PM10'

    ! Create directory if it does not exist.
    do s = 1, 5
      cmd = trim('mkdir -p '// trim(output_directory) // out_dir(s))
      call system(cmd)
    end do


    ! gas phase
    do s = 1, n_gas
      output_filename = trim(output_directory) // "/gas/" // trim(species_name(s)) // trim(out_type(output_type))
      ! Remove if output file exist
      inquire (file = output_filename, exist = file_exists)
      if (file_exists) then
        open(unit=100, file = output_filename, status='old', iostat=stat)
        if (stat == 0) close(100, status='delete')
      endif
      ! create new file 
      open(unit=100,file=output_filename, status="new")
      do t=1, nt+1
        write(100,*) output_gas(t,s)                  ! concentrations only
!        write(100,*) output_time(t),output_gas(t,s)   ! time + concentrations
      enddo
      close(100)
    enddo

    ! aerosols
    do s = 1, N_aerosol
      do b = 1, N_size
        output_filename = trim(output_directory) // "/aero/" // trim(aerosol_species_name(s)) &
                  // "_" // trim(str(b)) // trim(out_type(output_type))
        ! Remove if output files exist
        inquire (file = output_filename, exist = file_exists)
        if (file_exists) then
          open(unit=100, file = output_filename, status='old', iostat=stat)
          if (stat == 0) close(100, status='delete')
        endif
        ! create new file 
        open(unit=100,file=output_filename, status="new")
        do t=1, nt+1
          write(100,*) output_aero(t,s,b)                  ! concentrations only
!          write(100,*) output_time(t),output_aero(t,s,b)   ! time + concentrations
        enddo
        close(100)
      enddo
    enddo

    ! number
    do b = 1, N_size
      output_filename = trim(output_directory) // "/number/NUMBER_"// trim(str(b)) // trim(out_type(output_type))
      ! Remove if output files exist
      inquire (file = output_filename, exist = file_exists)
      if (file_exists) then
        open(unit=100, file = output_filename, status='old', iostat=stat)
        if (stat == 0) close(100, status='delete')
      endif
      ! create new file 
      open(unit=100,file=output_filename, status="new")
      do t=1, nt+1
        write(100,*) output_numb(t,b)                  ! concentrations only
!        write(100,*) output_time(t),output_numb(t,b)   ! time + concentrations
      enddo
      close(100)
    end do

    ! total number
    output_filename = trim(output_directory) // "/number/TNUM"//trim(out_type(output_type))    
    ! Remove if output files exist
    inquire (file = output_filename, exist = file_exists)
    if (file_exists) then
      open(unit=100, file = output_filename, status='old', iostat=stat)
      if (stat == 0) close(100, status='delete')
    endif   
    ! create new file  
    open(unit=100,file=output_filename, status="new")
    do t=1, nt+1
      write(100,*) output_numb(t,N_size+1)                  ! concentrations only
!      write(100,*) output_time(t),output_numb(t,N_size+1)   ! time + concentrations
    enddo
    close(100)
          
    ! total_aero_mass for aerosols
    do s = 1, N_aerosol
      output_filename = trim(output_directory) // "/TM/" // trim(aerosol_species_name(s)) &
                         //'_TM'//trim(out_type(output_type)) 
      ! Remove if output files exist
      inquire (file = output_filename, exist = file_exists)
      if (file_exists) then
        open(unit=100, file = output_filename, status='old', iostat=stat)
        if (stat == 0) close(100, status='delete')
      endif
      ! create new file  
      open(unit=100,file=output_filename, status="new")
      do t=1, nt+1
        write(100,*) output_TM(t,s,1)                  ! concentrations only
!        write(100,*) output_time(t),output_TM(t,s,1)   ! time + concentrations
      enddo
      close(100)

      if (aerosol_species_interact(s) .gt. 0) then
         output_filename = trim(output_directory) // "/TM/" // trim(aerosol_species_name(s))&
           //'_'//trim(species_name(aerosol_species_interact(s)))//'_TM'// trim(out_type(output_type))
        ! Remove if output files exist
        inquire (file = output_filename, exist = file_exists)
        if (file_exists) then
          open(unit=100, file = output_filename, status='old', iostat=stat)
          if (stat == 0) close(100, status='delete')
        endif
        ! create new file  
        open(unit=100,file=output_filename, status="new")
        do t=1, nt+1
          write(100,*) output_TM(t,s,3)                  ! concentrations only
!          write(100,*) output_time(t),output_TM(t,s,3)   ! time + concentrations
        enddo 
        close(100)
      endif
    enddo

   ! cell_diam_av for aerosols
   do b = 1, N_size
     output_filename = trim(output_directory) // "/diameter/DIAMETER_"// &
                       trim(str(b)) // trim(out_type(output_type))
     ! Remove if output files exist
     inquire (file = output_filename, exist = file_exists)
     if (file_exists) then
       open(unit=100, file = output_filename, status='old', iostat=stat)
       if (stat == 0) close(100, status='delete')
     endif
    ! create new file       
     open(unit=100,file=output_filename, status="new")
     do t=1, nt+1
       write(100,*) output_diam(t,b)                  ! concentrations only
!       write(100,*) output_time(t),output_diam(t,b)   ! time + concentrations
     enddo
     close(100)
   end do

   ! aero specials
    do s=1, 8
      output_filename = trim(output_directory) // "/aero/" // trim(out_aero(s))//trim(out_type(output_type))
      ! Remove if output file exist
      inquire (file = output_filename, exist = file_exists)
      if (file_exists) then
        open(unit=100, file = output_filename, status='old', iostat=stat)
        if (stat == 0) close(100, status='delete')
      endif
      ! create new file  
      open(unit=100,file=output_filename, status="new")
      do t=1, nt+1
        write(100,*) output_special(t,s)                  ! concentrations only
!        write(100,*) output_time(t),output_special(t,s)   ! time + concentrations
      enddo
      close(100)
    enddo
    
    ! pH
    do b=1, N_size
      output_filename = trim(output_directory) // "/aero/pH_" // trim(str(b)) // trim(out_type(output_type))
      ! Remove if output files exist
      inquire (file = output_filename, exist = file_exists)
      if (file_exists) then
        open(unit=100, file = output_filename, status='old', iostat=stat)
        if (stat == 0) close(100, status='delete')
      endif
      ! create new file 
      open(unit=100,file=output_filename, status="new")
      do t=1, nt+1
        write(100,*) output_pH(t,b)                  ! concentrations only
!        write(100,*) output_time(t),output_pH(t,b)   ! time + concentrations
      enddo
      close(100)
    enddo
    

  END SUBROUTINE ssh_write_output_bintxt



  SUBROUTINE ssh_write_output_genoa()

    integer :: stat, s, b, t
    logical :: file_exists
    character (len=200) output_filename
    character (len=200) :: cmd
    character(4) :: out_type(2) = (/".txt",".bin"/) ! 1: text, 2: binary
        
    ! add output for gas
    if (nout_gas.gt.0.) then
    
      ! create folder
      cmd = trim('mkdir -p '// trim(output_directory) //"/gas/")
      call system(cmd)

     do s = 1, nout_gas
        output_filename = trim(output_directory) // "/gas/" &
           //trim(output_gas_species(s))// trim(out_type(output_type))
        
        ! Remove if output file exist
        inquire (file = output_filename, exist = file_exists)
        if (file_exists) then
          open(unit=100, file = output_filename, status='old', iostat=stat)
          if (stat == 0) close(100, status='delete')
          endif
          ! create new file 
          open(unit=100,file=output_filename, status="new")
          do t=1, nt+1
             write(100,*) output_gas(t,s)
          enddo
          close(100)
      enddo
    endif

    ! output for aerosols
    if (nout_aero.gt.0) then
      ! Create directory if it does not exist.
      cmd = trim('mkdir -p '// trim(output_directory) //"/aero/")
      call system(cmd)

      ! aerosols
      do b = 1, N_size
        do s = 1, nout_aero
            output_filename = trim(output_directory) // "/aero/"&
            //trim(output_aero_species(s)) &
            //"_"//trim(str(b)) // trim(out_type(output_type))
            ! Remove if output files exist
            inquire (file = output_filename, exist = file_exists)
            if (file_exists) then
                open(unit=100, file = output_filename, status='old', iostat=stat)
                if (stat == 0) close(100, status='delete')
            endif
            ! create new file 
            open(unit=100,file=output_filename, status="new")
            do t=1, nt+1
              write(100,*) output_aero(t,s,b) ! concentrations only
            enddo
            close(100)
        enddo
      enddo
    endif
    
    ! write concs.txt
    if (output_type .ne. 0) then
      ! Create directory if it does not exist.
      cmd = trim('mkdir -p '// trim(output_directory))
      call system(cmd)
    
      ! error computation - at least save total soa concs
      output_filename = trim(output_directory) // "/concs" // trim(out_type(output_type))
      ! Remove if output files exist
      inquire (file = output_filename, exist = file_exists)
      if (file_exists) then
          open(unit=100, file = output_filename, status='old', iostat=stat)
          if (stat == 0) close(100, status='delete')
      endif
      ! create new file
      open(unit=100,file=output_filename, status="new")
      if (iout_round) then
        do t=1, nt+1
          ! output total SOA
          if (total_soa(t,1).lt.1d7) then ! round
              write(100,'(999E15.6)') (total_soa(t,s), s=1, nout_total)
          else ! no round
              write(100,*) (total_soa(t,s), s=1, nout_total)
          end if
        enddo
      else ! no round up
        do t=1, nt+1
          write(100,*) (total_soa(t,s), s=1, nout_total)
        enddo
      endif
      close(100)
    endif
    
  END SUBROUTINE ssh_write_output_genoa


!============================================================================================


  character(len=20) function str(k)
    !   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
  end function str
  
  end module Resultoutput
