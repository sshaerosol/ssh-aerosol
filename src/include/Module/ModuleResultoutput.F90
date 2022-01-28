!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------

Module Resultoutput
  use aInitialization
  use dPhysicalbalance

  implicit none

! out_aero : array of file names; outpout time variation of organic, inorganic, PBC, Dust, PM2.5, PM10 results
  character(20), save :: out_aero(6) 
  character(4), save :: out_type(2) = (/".txt",".bin"/) ! 1: text, 2: binary

contains

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
		write(unit=10,FMT=*)  'With nucleation', '	','nucl_model', nucl_model
	else
		write(unit=10,FMT=*)  'Without nucleation'
	end if
	write(unit=10,FMT=*)  'DTAEROMIN', DTAEROMIN, '	','redistribution_method', redistribution_method
	write(unit=10,FMT=*)  'Method',dynamic_solver, '	 with_oligomerization', with_oligomerization
	write(unit=10,FMT=*)  'with_fixed_density', with_fixed_density,'	', 'fixed_density', fixed_density, 'kg/m^3'
	write(unit=10,FMT=*)
	write(unit=10,FMT=*) 'output directory :', trim(output_directory),'/'

   CLOSE(10)

   end subroutine ssh_save_report


  subroutine ssh_save_concentration()
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This subroutine records simulation results over each time step. 
!
!------------------------------------------------------------------------
!
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
!
!------------------------------------------------------------------------

    integer :: s, b, i, j, jesp
    !real (kind = 4) :: conc_save
    double precision :: conc_save
    character (len=200) output_filename
    integer :: out_aero_num(6) = (/2,1,7,33,0,33/)  ! last index of required aero species list; 0 for PM2.5

    ! **** output_directory/gas/
    ! save gas concentration results over each time step
    do s = 1, n_gas
          output_filename = trim(output_directory) // "/gas/" // trim(species_name(s)) // trim(out_type(output_type))
          open(unit=100,file=output_filename, status='old', position = "append")
	       write(100,*) concentration_gas_all(s) 
          close(100)
    enddo

     ! **** output_directory/number/
     ! save number concentration results over each time step
     do b = 1, N_size
          output_filename = trim(output_directory) // "/number/NUMBER_" // trim(str(b)) // trim(out_type(output_type))
          open(unit=100,file=output_filename, status="old", position = "append")
               write(100,*) concentration_number(b)
          close(100)
     end do

     ! save total number concentration of all particles
	conc_save = 0.d0
        do b = 1, N_size
	   conc_save = conc_save + concentration_number(b)
	end do

        output_filename = trim(output_directory) // "/number/TNUM" //  trim(out_type(output_type))  
        open(unit=100,file=output_filename, status="old", position = "append")
             write(100,*) conc_save 
        close(100)

     ! **** output_directory/diameter/
     ! save dry diacell_diam_av results over each time step
     do b = 1, N_size
        output_filename = trim(output_directory) // "/diameter/DIAMETER_" // trim(str(b)) &
                          //trim(out_type(output_type))
        open(unit=100,file=output_filename, status="old", position = "append")
             write(100,*) cell_diam_av(b)
        close(100)
     end do

    ! **** output_directory/TM/
    ! re-new total_aero_mass(N_aerosol)
        total_aero_mass = 0.d0
	do s = 1, N_aerosol_layers
           jesp = List_species(s)
	   do j = 1,N_size
              total_aero_mass(jesp) = total_aero_mass(jesp) + concentration_mass(j,s)
	   enddo
	end do

    ! save total mass for each aerosol species over each time step
    do s = 1, N_aerosol
       output_filename = trim(output_directory) // "/TM/" // trim(aerosol_species_name(s))&
                         //'_TM'// trim(out_type(output_type))  
       open(unit=100,file=output_filename, status='old', position = "append")
            write(100,*) total_aero_mass(s) 
       close(100)
       if (aerosol_species_interact(s) .gt. 0) then
          output_filename = trim(output_directory) // "/TM/" // trim(aerosol_species_name(s))&
            //'_'//trim(species_name(aerosol_species_interact(s)))//'_TM'// trim(out_type(output_type))
          open(unit=100,file=output_filename, status='old', position = "append")
               conc_save = total_aero_mass(s) + concentration_gas_all(aerosol_species_interact(s))
               write(100,*) conc_save
          close(100)
       end if
    enddo

    do s = 1, N_aerosol
       do b = 1, N_size
          concentration_mass_tmp(b ,s) = 0.d0
       enddo
    enddo
    
    do b = 1, N_size
       do s = 1, N_aerosol_layers
          jesp = List_species(s)
          concentration_mass_tmp(b ,jesp) = concentration_mass_tmp(b ,jesp) + concentration_mass(b ,s)
       enddo
    enddo

    ! **** output_directory/aero/
    ! save aerosol concentration results over each time step
    do s = 1, N_aerosol
       do b = 1, N_size
          output_filename = trim(output_directory) // "/aero/" // trim(aerosol_species_name(s)) &  
                  // "_" // trim(str(b)) // trim(out_type(output_type))
          open(unit=100,file=output_filename, status="old", position = "append")
               write(100,*) concentration_mass_tmp(b, s)
          close(100)
       end do
    end do

    ! save organic, inorganic, PBC, Dust, PM2.5, PM10 results over each time step
    ! generate files' names
    do s = 1, 6
        conc_save = 0.0
        output_filename = trim(output_directory) // "/aero/" // trim(out_aero(s))//trim(out_type(output_type))
        open(unit=100,file=output_filename, status="old", position = "append") ! index (/1,2,7,33,0,33/) 
	     if (s .eq. 5) then	! PM2.5
		do b = 1, N_sizebin
		     if (diam_bound(b) .ge. 2.5d0 ) exit
		     j = b	! j = last size bin number of PM 2.5 
		enddo
		do b = 1, N_aerosol-1	! without water
		     do i = 1, j
			conc_save = conc_save + concentration_mass_tmp(i, b)
		     enddo
		enddo
	     else if (out_aero_num(s) .le. 2) then ! PBC and Dust, not PM2.5
		conc_save = total_aero_mass(out_aero_num(s))
	     else if (out_aero_num(s) .gt. 2) then ! inorganic (index 3-7) organic(index 8-33) PM10 (index 1-33)
		do b = out_aero_num(s-1) + 1, out_aero_num(s)
			conc_save = conc_save + total_aero_mass(b)
		enddo
	     end if
             write(100,*) conc_save
        close(100)
    end do

  end subroutine ssh_save_concentration



  subroutine ssh_init_output_conc()
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This subroutine initiailize output files, which should be called 
!     before save_concentration()
!
!------------------------------------------------------------------------

    integer :: stat, s, b
    logical :: file_exists
    character (len=200) output_filename
    character (len=200) :: cmd
    character (len=10) :: out_dir(5) 
    out_dir(1) = "/number/"
    out_dir(2) = "/gas/"	
    out_dir(3) = "/aero/"
    out_dir(4) = "/TM/"
    out_dir(5) = "/diameter/"
    ! update out_aero, outpout names 
    out_aero(1) = 'Black_Carbon'
    out_aero(2) = 'Dust'
    out_aero(3) = 'Inorganic'
    out_aero(4) = 'Organic'
    out_aero(5) = 'PM2.5'
    out_aero(6) = 'PM10'

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
          ! creative new empty file 
          open(unit=100,file=output_filename, status="new")
          close(100)
    enddo


    ! organic, inorganic, PBC, Dust, PM2.5, PM10
    do s = 1, 6
          output_filename = trim(output_directory) // "/aero/" // trim(out_aero(s))//trim(out_type(output_type))
          ! Remove if output file exist
          inquire (file = output_filename, exist = file_exists)
          if (file_exists) then
             open(unit=100, file = output_filename, status='old', iostat=stat)
             if (stat == 0) close(100, status='delete')
          endif
          ! creative new empty file 
          open(unit=100,file=output_filename, status="new")
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
          ! creative new empty file 
          open(unit=100,file=output_filename, status="new")
          close(100)
       end do
    end do

    ! number
    do b = 1, N_size
       output_filename = trim(output_directory) // "/number/NUMBER_"// trim(str(b)) // trim(out_type(output_type))
       ! Remove if output files exist
       inquire (file = output_filename, exist = file_exists)
       if (file_exists) then
          open(unit=100, file = output_filename, status='old', iostat=stat)
          if (stat == 0) close(100, status='delete')
       endif
       ! creative new empty file 
       open(unit=100,file=output_filename, status="new")
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
    open(unit=100,file=output_filename, status="new")
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
       open(unit=100,file=output_filename, status="new")
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
          open(unit=100,file=output_filename, status="new")
          close(100)
       end if
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
      open(unit=100,file=output_filename, status="new")
      close(100)
   end do


  end subroutine ssh_init_output_conc


  subroutine ssh_delete_empty_file()
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This subroutine delete output file in where all the values are zero.
!
!------------------------------------------------------------------------
    integer :: stat, s, b
    real :: conc_value
    logical :: file_exists
    character (len=200) output_filename
    character (len=200) :: cmd
    character (len=10) :: out_dir(5) 
   ! gas phase
    do s = 1, n_gas
       output_filename = trim(output_directory) // "/gas/" // trim(species_name(s)) // trim(out_type(output_type))
       open(unit=100, file = output_filename, status='old', iostat=stat)
           conc_value = 0.0
           do while(stat .eq. 0)
              read(100, *,iostat=stat) conc_value
              if (conc_value .gt. 0.0) exit
           end do
       if (conc_value .eq. 0.0) close(100, status='delete') ! if all values are zero, delete file
       if (conc_value .ne. 0.0) close(100)
    enddo


    do s = 1, N_aerosol
    ! total_aero_mass for aerosols
       output_filename = trim(output_directory) // "/TM/" // trim(aerosol_species_name(s)) &
                         //'_TM'//trim(out_type(output_type)) 
       open(unit=100, file = output_filename, status='old', iostat=stat)
           conc_value = 0.0
           do while(stat .eq. 0)
              read(100, *,iostat=stat) conc_value
              if (conc_value .gt. 0.0) exit
           end do
       if (conc_value .eq. 0.0) close(100, status='delete')
       if (conc_value .ne. 0.0) close(100)

    ! total_mass for aerosols + precursor
       if (aerosol_species_interact(s) .gt. 0) then
          output_filename = trim(output_directory) // "/TM/" // trim(aerosol_species_name(s))&
            //'_'//trim(species_name(aerosol_species_interact(s)))//'_TM'// trim(out_type(output_type))
            open(unit=100, file = output_filename, status='old', iostat=stat)
                conc_value = 0.0
                do while(stat .eq. 0)
                   read(100, *,iostat=stat) conc_value
                   if (conc_value .gt. 0.0) exit
                 end do
            if (conc_value .eq. 0.0) close(100, status='delete')
            if (conc_value .ne. 0.0) close(100)
       end if

    ! aerosols mass conc. of each cell
       do b = 1, N_size
	  output_filename = trim(output_directory) // "/aero/" // trim(aerosol_species_name(s)) &
                  // "_" // trim(str(b)) // trim(out_type(output_type))
          open(unit=100, file = output_filename, status='old', iostat=stat)
              conc_value = 0.0
              do while(stat .eq. 0)
                 read(100, *,iostat=stat) conc_value
                 if (conc_value .gt. 0.0) exit
              end do
          if (conc_value .eq. 0.0) close(100, status='delete')
          if (conc_value .ne. 0.0) close(100)
       end do
    enddo

    ! number
    do b = 1, N_size
       output_filename = trim(output_directory) // "/number/NUMBER_"// trim(str(b)) // trim(out_type(output_type))
       open(unit=100, file = output_filename, status='old', iostat=stat)
          conc_value = 0.0
          do while(stat .eq. 0)
             read(100, *,iostat=stat) conc_value
             if (conc_value .gt. 0.0) exit
          end do
       if (conc_value .eq. 0.0) close(100, status='delete')
       if (conc_value .ne. 0.0) close(100)
    enddo

  end subroutine ssh_delete_empty_file

  character(len=20) function str(k)
    !   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
  end function str
  
  end module Resultoutput
