!!-----------------------------------------------------------------------
!!     Copyright (C) 2003-2014, ENPC - INRIA - EDF R&D
!!     Author(s): Shupeng Zhu
!!
!!     This file is part of the Size Composition Resolved Aerosol Model (SCRAM), a
!!     component of the air quality modeling system Polyphemus.
!
!!-----------------------------------------------------------------------

Module Resultoutput
  use aInitialization
  use dPhysicalbalance

  implicit none

contains

   subroutine save_report()
   implicit none
	integer :: j
   OPEN(UNIT=10,FILE=trim(output_directory) // "/" // "report.txt")
	write(unit=10,FMT=*),'<<<< Meteorological setup >>>>'
        write(unit=10,FMT=*),'location', latitude, 'N','	', longitude,'E','	','Temperature', Temperature, 'K'
        write(unit=10,FMT=*),'Pressure', Pressure, 'Pa', '	','Specific Humidity', Humidity, '	',&
				'Cloud attenuation field', attenuation, '	relative humidity', relative_humidity
	write(unit=10,FMT=*)
	write(unit=10,FMT=*),'<<<< Simulation time setup >>>>'
        write(unit=10,FMT=*),'Begining time (from Jan. 1st)', initial_time, 's','	',&
				'Simulation Time', final_time,'s','	','Initial Time Step', delta_t,'s',&
         			'	','Number of iterations:', nt
	write(unit=10,FMT=*)
	write(unit=10,FMT=*),'<<<< Inition condition >>>>'

	if (tag_init == 0) write(unit=10,FMT=*),&
			'Internally mixed aerosol species are provided for initial condition.'


	write(unit=10,FMT=*), 'Gas-phase conc. input file :', init_gas_conc_file
	write(unit=10,FMT=*), 'Particle conc. input file :', init_aero_conc_mass_file
        if (with_init_num .eq. 1) write(unit=10,FMT=*), 'Aerosol number conc. is read from file :',&
							 init_aero_conc_num_file
	if (with_init_num .eq. 0) then
		write(unit=10,FMT=*),' Aerosol number conc. is estimated from mass and diameter.' 
	end if
	write(unit=10,FMT=*), '====== concentration_number : ======'
	write(unit=10,FMT=*), concentration_number
	write(unit=10,FMT=*)
	write(unit=10,FMT=*), 'N_sizebin', N_sizebin
	write(unit=10,FMT=*), '===== diam_bounds : ====='
	write(unit=10,FMT=*), diam_bound
	write(unit=10,FMT=*)
	write(unit=10,FMT=*), '===== cell_diam_av : ====='
	write(unit=10,FMT=*), cell_diam_av 
	write(unit=10,FMT=*)
	write(unit=10,FMT=*), '===== wet_volume : ====='
	write(unit=10,FMT=*), wet_volume
	write(unit=10,FMT=*)
	write(unit=10,FMT=*), '===== concentration_mass : ====='
	write(unit=10,FMT=*), aerosol_species_name(4)
	write(unit=10,FMT=*), (concentration_mass(j,4), j = 1, N_size)
	write(unit=10,FMT=*)
	write(unit=10,FMT=*), aerosol_species_name(34)
	write(unit=10,FMT=*), (concentration_mass(j,34), j = 1, N_size)

           
	write(unit=10,FMT=*)
	write(unit=10,FMT=*),'<<<< Mixing state >>>>'
	if (tag_external == 1) write(unit=10,FMT=*), 'simulation is mixing-state resolved.', &
					'N_frac', N_frac ,'	','N_groups', N_groups
	if (tag_external == 0) write(unit=10,FMT=*), 'simulation is internally mixed.' ,&
					'N_frac', N_frac ,'	','N_groups', N_groups
	write(unit=10,FMT=*),'frac_bound', frac_bound
	write(unit=10,FMT=*)
	write(unit=10,FMT=*),'<<<< Species lists >>>>'
	write(unit=10,FMT=*), 'gas phase species file :', species_list_file
	write(unit=10,FMT=*), 'particle species file :', aerosol_species_list_file
	write(unit=10,FMT=*)
	if (tag_chem == 0) write(unit=10,FMT=*),'<<<< Without Gas-phase chemistry >>>>'
	if (tag_chem == 1) then
		write(unit=10,FMT=*),'<<<< Gas-phase chemistry >>>>'
		write(unit=10,FMT=*), 'with_heterogeneous', with_heterogeneous,'	','with_adaptive', with_adaptive,&
					'	','adaptive time step tolerance', adaptive_time_step_tolerance,&
					'	','min adaptive time step', min_adaptive_time_step
	end if

	write(unit=10,FMT=*)
	write(unit=10,FMT=*),'<<<< Emissions >>>>'
	if (tag_emis == 1) then
		write(unit=10,FMT=*), 'With internally-mixed emissions.'
		write(unit=10,FMT=*), 'Gas-phase conc. emission file :', emis_gas_file
		write(unit=10,FMT=*), 'Particle conc. emission file :', emis_aero_mass_file
		if (with_emis_num == 1) write(unit=10,FMT=*),&
			'Emitted aerosol number conc. is read from file :', emis_aero_num_file
		if (with_emis_num == 0) write(unit=10,FMT=*),&
			'Emitted aerosol number conc. is estimated from mass and diameter.'
	end if
        if (tag_emis == 0) write(unit=10,FMT=*),'Without emission.'

	write(unit=10,FMT=*)
	write(unit=10,FMT=*),'<<<< Particle Dynamic >>>>'
	if (with_cond == 1) then 
		write(unit=10,FMT=*), 'With condensation', '	','Cut_dim', Cut_dim,'	', 'ISOAPDYN', ISOAPDYN
	else
		write(unit=10,FMT=*), 'Without condensation'
	end if
	if (with_coag == 1) then 
		write(unit=10,FMT=*), 'With coagulation', '	','coefficient file :', Coefficient_file
	else
		write(unit=10,FMT=*), 'Without coagulation'
	end if
	if (with_nucl == 1) then 
		write(unit=10,FMT=*), 'With nucleation', '	','nucl_model', nucl_model
	else
		write(unit=10,FMT=*), 'Without nucleation'
	end if
	write(unit=10,FMT=*), 'DTAEROMIN', DTAEROMIN, '	','redistribution_method', redistribution_method
	write(unit=10,FMT=*), 'Method',dynamic_solver, '	 with_oligomerization', with_oligomerization
	write(unit=10,FMT=*), 'with_fixed_density', with_fixed_density,'	', 'fixed_density', fixed_density, 'kg/m^3'
	if (wet_diam_estimation == 1) write(unit=10,FMT=*), 'initial wet diameter is computed by Gerber.'
	if (wet_diam_estimation == 0) write(unit=10,FMT=*), 'initial wet diameter is computed by isorropia.'

	write(unit=10,FMT=*)
	write(unit=10,FMT=*),'output directory :', trim(output_directory),'/'

   CLOSE(10)

   end subroutine save_report




  subroutine save_concentration()

    integer :: s, b
    real (kind = 4) :: conc_save
    character (len=100) output_filename


    ! save gas concentration results over each time step
    do s = 1, n_gas
       if (output_type == 1) then
          output_filename = trim(output_directory) // "/gas/" // trim(species_name(s)) // ".txt"     !modify here

          OPEN(UNIT=100,FILE=output_filename, status='old', position = "append")

          conc_save = concentration_gas_all(s) 
          write(100,*) conc_save

       else if (output_type == 2) then
          output_filename = trim(output_directory) // "/gas/" // trim(species_name(s)) // ".bin"

          OPEN(UNIT=100,FILE=output_filename, status='old', form='unformatted', access='stream', position = 'append')

          conc_save = concentration_gas_all(s) 
          write(100) conc_save

       endif

       close(100)
    enddo

    ! save aerosol concentration results over each time step
    do s = 1, N_aerosol
       do b = 1, N_size
          if (output_type == 1) then
             output_filename = trim(output_directory) // "/aero/" // trim(aerosol_species_name(s)) &     !modify here
                  // "_" // trim(str(b)) // ".txt"
             OPEN(UNIT=100,FILE=output_filename, status="old", position = "append")
             conc_save = concentration_mass(b, s)
             write(100,*) concentration_mass(b, s)

          else if (output_type == 2) then
        
             output_filename = trim(output_directory) // "/aero/" // trim(aerosol_species_name(s)) &
                  // "_" // trim(str(b)) // ".bin"
             OPEN(UNIT=100,FILE=output_filename, status="old", form='unformatted', &
                  access='stream', position = 'append')
             conc_save = concentration_mass(b, s)
             write(100) conc_save

          end if

          close(100)
       end do
    end do

     ! save number concentration results over each time step
     do b = 1, N_size

          if (output_type == 1) then

             output_filename = trim(output_directory) // "/number/NUMBER_" // trim(str(b)) // ".txt"     !modify here
             OPEN(UNIT=100,FILE=output_filename, status="old", position = "append")
             conc_save = concentration_number(b)
             write(100,*) concentration_number(b)! conc_save

          else if (output_type == 2) then
        
             output_filename = trim(output_directory) // "/number/NUMBER_"// trim(str(b)) // ".bin"
             OPEN(UNIT=100,FILE=output_filename, status="old", form='unformatted', &
                  access='stream', position = 'append')
             conc_save = concentration_number(b)
             write(100) conc_save

          end if

          close(100)
       end do

	total_number = 0.d0
        do b = 1, N_size
		total_number = total_number + concentration_number(b)
	end do
 	if (output_type == 1) then

             output_filename = trim(output_directory) // "/number/TNUM.txt"     
             OPEN(UNIT=100,FILE=output_filename, status="old", position = "append")
             conc_save = total_number
             write(100,*) total_number !conc_save

          else if (output_type == 2) then
        
             output_filename = trim(output_directory) // "/number/TNUM.bin"
             OPEN(UNIT=100,FILE=output_filename, status="old", form='unformatted', &
                  access='stream', position = 'append')
             conc_save = total_number
             write(100) conc_save
	end if


    ! save total mass for each aerosol species over each time step
    do s = 1, N_aerosol
       if (output_type == 1) then
          output_filename = trim(output_directory) // "/TM/" // trim(aerosol_species_name(s)) //'_TM'// ".txt"  

          OPEN(UNIT=100,FILE=output_filename, status='old', position = "append")

          conc_save = total_mass(s) 
          write(100,*) conc_save

       else if (output_type == 2) then
          output_filename = trim(output_directory) // "/TM/" // trim(aerosol_species_name(s)) //'_TM'// ".bin"  

          OPEN(UNIT=100,FILE=output_filename, status='old', form='unformatted', access='stream', position = 'append')

          conc_save = total_mass(s) 
          write(100) conc_save

       endif

       close(100)
    enddo


     ! save cell_diam_av results over each time step
     do b = 1, N_size

        if (output_type == 1) then

             output_filename = trim(output_directory) // "/diameter/DIAMETER_" // trim(str(b)) // ".txt"     !modify here
             OPEN(UNIT=100,FILE=output_filename, status="old", position = "append")
             write(100,*) cell_diam_av(b)! conc_save

          else if (output_type == 2) then
        
             output_filename = trim(output_directory) // "/diameter/DIAMETER_"// trim(str(b)) // ".bin"
             OPEN(UNIT=100,FILE=output_filename, status="old", form='unformatted', &
                  access='stream', position = 'append')
             write(100) cell_diam_av(b)

          end if

          close(100)
       end do

  end subroutine save_concentration



  subroutine save_values()

    integer :: s, b
    real (kind = 4) :: conc_save
    character (len=100) output_filename

! save cell_diam_av
       if (output_type == 1) then
          output_filename = trim(output_directory) // "/cell_diam_av.txt"   

          OPEN(UNIT=100,FILE=output_filename, status='old', position = "append")
	  write(100,*) 'cell_diam_av'
          do s = 1, N_size
		write(100,*) cell_diam_av(s)
	  end do

       else if (output_type == 2) then
          output_filename = trim(output_directory) // "/cell_diam_av.bin"    

          OPEN(UNIT=100,FILE=output_filename, status='old', form='unformatted', access='stream', position = 'append')

          write(100) cell_diam_av

       endif

       close(100)



! save wet_volume

       if (output_type == 1) then
          output_filename = trim(output_directory) // "/wet_volume.txt"  

          OPEN(UNIT=100,FILE=output_filename, status='old', position = "append")

          write(100,*) wet_volume

       else if (output_type == 2) then
          output_filename = trim(output_directory) // "/wet_volume.bin"   

          OPEN(UNIT=100,FILE=output_filename, status='old', form='unformatted', access='stream', position = 'append')

          write(100) wet_volume

       endif

       close(100)





  end subroutine save_values




  ! ! ! ! ! ! ! ! ! ! ! Initiailize output files.
  subroutine init_output_conc()

    integer :: stat, s, b
    logical :: file_exists
    character (len=100) output_filename
    character (len=80) :: cmd
    character (len=10) :: out_dir(5)
    out_dir(1) = "/number/"
    out_dir(2) = "/gas/"	
    out_dir(3) = "/aero/"
    out_dir(4) = "/TM/"
    out_dir(5) = "/diameter/"
    ! Create directory if it does not exist.
    do s = 1, 5
       	cmd = trim('mkdir -p '// trim(output_directory) // out_dir(s))
       	call system(cmd)
    end do

    ! gas
    do s = 1, n_gas

       if (output_type == 1) then
          output_filename = trim(output_directory) // "/gas/" // trim(species_name(s)) // ".txt"   !modify here

          ! Remove if output files exist
          inquire (file = output_filename, exist = file_exists)
          if (file_exists) then
             open(unit=100, file = output_filename, status='old', iostat=stat)
             if (stat == 0) close(100, status='delete')
          endif

          OPEN(UNIT=100,FILE=output_filename, status="new")

       else if (output_type == 2) then

          output_filename = trim(output_directory) // "/gas/" // trim(species_name(s)) // ".bin"

          ! Remove if output files exist
          inquire (file = output_filename, exist = file_exists)
          if (file_exists) then
             open(unit=100, file = output_filename, status='old', iostat=stat)
             if (stat == 0) close(100, status='delete')
          endif
       
          OPEN(UNIT=100,FILE=output_filename, status="new", form='unformatted')
       end if

       close(100)
    enddo
    !aerosols
    do s = 1, N_aerosol
       do b = 1, N_size

          if (output_type == 1) then
          
             output_filename = trim(output_directory) // "/aero/" // trim(aerosol_species_name(s)) &     !modify here
                  // "_" // trim(str(b)) // ".txt"

             ! Remove if output files exist
             inquire (file = output_filename, exist = file_exists)
             if (file_exists) then
                open(unit=100, file = output_filename, status='old', iostat=stat)
                if (stat == 0) close(100, status='delete')
             endif
          
             OPEN(UNIT=100,FILE=output_filename, status="new")

          else if (output_type == 2) then
          
             output_filename = trim(output_directory) // "/aero/" // trim(aerosol_species_name(s)) &
                  // "_" // trim(str(b)) // ".bin"

             ! Remove if output files exist
             inquire (file = output_filename, exist = file_exists)
             if (file_exists) then
                open(unit=100, file = output_filename, status='old', iostat=stat)
                if (stat == 0) close(100, status='delete')
             endif
          
             OPEN(UNIT=100,FILE=output_filename, status="new", form='unformatted')
          end if

          close(100)
       end do
    end do

    ! number
        do b = 1, N_size

          if (output_type == 1) then
             output_filename = trim(output_directory) // "/number/NUMBER_"// trim(str(b)) // ".txt"     !modify here
             ! Remove if output files exist
             inquire (file = output_filename, exist = file_exists)
             if (file_exists) then
                open(unit=100, file = output_filename, status='old', iostat=stat)
                if (stat == 0) close(100, status='delete')
             endif
          
             OPEN(UNIT=100,FILE=output_filename, status="new")

          else if (output_type == 2) then
          
             output_filename = trim(output_directory) // "/number/NUMBER_"// trim(str(b)) // ".bin"
             ! Remove if output files exist
             inquire (file = output_filename, exist = file_exists)
             if (file_exists) then
                open(unit=100, file = output_filename, status='old', iostat=stat)
                if (stat == 0) close(100, status='delete')
             endif
          
             OPEN(UNIT=100,FILE=output_filename, status="new", form='unformatted')
          end if

          close(100)
       end do


  ! total number

          if (output_type == 1) then
             output_filename = trim(output_directory) // "/number/TNUM.txt"     
             ! Remove if output files exist
             inquire (file = output_filename, exist = file_exists)
             if (file_exists) then
                open(unit=100, file = output_filename, status='old', iostat=stat)
                if (stat == 0) close(100, status='delete')
             endif
          
             OPEN(UNIT=100,FILE=output_filename, status="new")

          else if (output_type == 2) then
          
             output_filename = trim(output_directory) // "/number/TNUM.bin"
             ! Remove if output files exist
             inquire (file = output_filename, exist = file_exists)
             if (file_exists) then
                open(unit=100, file = output_filename, status='old', iostat=stat)
                if (stat == 0) close(100, status='delete')
             endif
          
             OPEN(UNIT=100,FILE=output_filename, status="new", form='unformatted')
          end if

   ! total_mass for aerosols
    do s = 1, N_aerosol

       if (output_type == 1) then
          output_filename = trim(output_directory) // "/TM/" // trim(aerosol_species_name(s)) //'_TM'// ".txt"  

          ! Remove if output files exist
          inquire (file = output_filename, exist = file_exists)
          if (file_exists) then
             open(unit=100, file = output_filename, status='old', iostat=stat)
             if (stat == 0) close(100, status='delete')
          endif

          OPEN(UNIT=100,FILE=output_filename, status="new")

       else if (output_type == 2) then
          output_filename = trim(output_directory) // "/TM/" // trim(aerosol_species_name(s)) //'_TM'// ".bin"  

          ! Remove if output files exist
          inquire (file = output_filename, exist = file_exists)
          if (file_exists) then
             open(unit=100, file = output_filename, status='old', iostat=stat)
             if (stat == 0) close(100, status='delete')
          endif
       
          OPEN(UNIT=100,FILE=output_filename, status="new", form='unformatted')
       end if

       close(100)
    enddo

    
   ! cell_diam_av for aerosols
       ! number
        do b = 1, N_size

          
           if (output_type == 1) then

               
             output_filename = trim(output_directory) // "/diameter/DIAMETER_"// trim(str(b)) // ".txt"     !modify here
             ! Remove if output files exist
             inquire (file = output_filename, exist = file_exists)
             if (file_exists) then
                open(unit=100, file = output_filename, status='old', iostat=stat)
                if (stat == 0) close(100, status='delete')
             endif
          
             OPEN(UNIT=100,FILE=output_filename, status="new")

          else if (output_type == 2) then
          
             output_filename = trim(output_directory) // "/diameter/DIAMETER_"// trim(str(b)) // ".bin"
             ! Remove if output files exist
             inquire (file = output_filename, exist = file_exists)
             if (file_exists) then
                open(unit=100, file = output_filename, status='old', iostat=stat)
                if (stat == 0) close(100, status='delete')
             endif
          
             OPEN(UNIT=100,FILE=output_filename, status="new", form='unformatted')
          end if

          close(100)
       end do

	!report
	!output_filename = trim(output_directory) // "/" // "report.txt"
	!inquire (file = output_filename, exist = file_exists)
	!if (file_exists) then
        !        open(unit=100, file = output_filename, status='old', iostat=stat)
        !        if (stat == 0) close(100, status='delete')
        !endif
        !OPEN(UNIT=100,FILE=output_filename, status="new")

  end subroutine init_output_conc


  ! Initiailize output files.
  subroutine init_output_vl()

    integer :: stat, s, b
    logical :: file_exists
    character (len=100) output_filename

   ! cell_diam_av(n_size) for aerosols

       if (output_type == 1) then
          output_filename = trim(output_directory) // "/cell_diam_av.txt"  

          ! Remove if output files exist
          inquire (file = output_filename, exist = file_exists)
          if (file_exists) then
             open(unit=100, file = output_filename, status='old', iostat=stat)
             if (stat == 0) close(100, status='delete')
          endif

          OPEN(UNIT=100,FILE=output_filename, status="new")

       else if (output_type == 2) then
          output_filename = trim(output_directory) // "/cell_diam_av.bin"  

          ! Remove if output files exist
          inquire (file = output_filename, exist = file_exists)
          if (file_exists) then
             open(unit=100, file = output_filename, status='old', iostat=stat)
             if (stat == 0) close(100, status='delete')
          endif
       
          OPEN(UNIT=100,FILE=output_filename, status="new", form='unformatted')
       end if

       close(100)

   ! wet_volume(n_size) for aerosols

       if (output_type == 1) then
          output_filename = trim(output_directory) // "/wet_volume.txt"  

          ! Remove if output files exist
          inquire (file = output_filename, exist = file_exists)
          if (file_exists) then
             open(unit=100, file = output_filename, status='old', iostat=stat)
             if (stat == 0) close(100, status='delete')
          endif

          OPEN(UNIT=100,FILE=output_filename, status="new")

       else if (output_type == 2) then
          output_filename = trim(output_directory) // "/wet_volume.bin"  

          ! Remove if output files exist
          inquire (file = output_filename, exist = file_exists)
          if (file_exists) then
             open(unit=100, file = output_filename, status='old', iostat=stat)
             if (stat == 0) close(100, status='delete')
          endif
       
          OPEN(UNIT=100,FILE=output_filename, status="new", form='unformatted')
       end if

       close(100)




  end subroutine init_output_vl




  character(len=20) function str(k)
    !   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
  end function str
  
  end module Resultoutput
