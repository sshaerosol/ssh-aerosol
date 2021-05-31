!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------
!!
!!     -- DESCRIPTION
!!    This module contains methods read coefficient repartition data
!!
!!-----------------------------------------------------------------------
Module bCoefficientRepartition
  use aInitialization
  use netcdf
  use omp_lib
  implicit none

  ! Parameter and variable definitions
  type :: ptr_to_real_array
     integer n
     double precision, dimension(:), pointer :: arr
  end type ptr_to_real_array

  type :: ptr_to_integer_array
     integer n
     integer, dimension(:), pointer :: arr
  end type ptr_to_integer_array

  ! Repartition coefficients.
  type(ptr_to_real_array), dimension(:), allocatable :: repartition_coefficient
  ! Index of repartition coefficients.
  type(ptr_to_integer_array), dimension(:), allocatable :: index1_repartition_coefficient
  type(ptr_to_integer_array), dimension(:), allocatable :: index2_repartition_coefficient

  ! Monte Carlo method.
  integer (kind=8) :: number_monte_carlo

  ! MPI
  integer :: rank, rank_size, ierr, tag = 2, master = 0
  integer, parameter :: size_buf = 100000


contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CHANGE!!!!!!!!!!
  SUBROUTINE ssh_ComputeCoefficientRepartition()
    
    implicit none

    !!integer (kind=8), parameter :: Nmc = 10000
    integer, parameter :: Nalloc_tmp = 1000000 
    real, dimension(:), allocatable :: repartition_coefficient_tmp
    real, dimension(:), allocatable :: index1_repartition_coefficient_tmp
    real, dimension(:), allocatable :: index2_repartition_coefficient_tmp
    double precision, dimension(:), allocatable :: random_vector
    double precision :: dsum,f1,f2,m1,m2,m12,fact,f12
    integer :: i1,i2,l1,l2,j1,j2,n,s,i,j,k,l,g
    integer :: ki
    integer, dimension(40) :: seed

    if (ssh_standalone) write(*,*) "Compute coefficient repartition"
    if (ssh_logger) write(logfile,*) "Compute coefficient repartition"

    ! Allocate repartition coefficient
    allocate(repartition_coefficient(N_size))
    allocate(index1_repartition_coefficient(N_size))
    allocate(index2_repartition_coefficient(N_size))
    allocate(repartition_coefficient_tmp(Nalloc_tmp))
    allocate(index1_repartition_coefficient_tmp(Nalloc_tmp))
    allocate(index2_repartition_coefficient_tmp(Nalloc_tmp))
    allocate(random_vector(2+2*N_groups))

    ! set seed for random number generator
    seed(1:40) = 47281
    call random_seed(put=seed)
 
    !!allocate(concentration_index_iv(N_sizebin, N_fracmax))
    !!allocate(discretization_composition(N_fracmax, N_groups, 2))
    ! Compute repartition coefficient
    do j = 1,N_size
       repartition_coefficient_tmp = 0.0
       k = concentration_index(j, 1)
       i = concentration_index(j, 2)         
       ! ki = j !!!ki =  concentration_index_iv(k,i)
       l=1
	do j1 = 1,N_size
          do j2 = 1,N_size
             dsum = 0.0
!$omp parallel do reduction(+:dsum) private(random_vector,l1,l2,i1,i2,m1,m2,m12,fact,g,f1,f2,f12)
             do n = 1, Nmc
                call random_number(random_vector)
                l1 = concentration_index(j1, 1)
                l2 = concentration_index(j2, 1)
                
                i1 = concentration_index(j1, 2)
                i2 = concentration_index(j2, 2)
                
                m1 = discretization_mass(l1) +random_vector(1) &
                     * (discretization_mass(l1+1)-discretization_mass(l1))
                m2 = discretization_mass(l2) + random_vector(2) &
                     * (discretization_mass(l2+1)-discretization_mass(l2))
              
                m12 = m1 + m2
                fact = 0.0           
                if(k /= N_sizebin) then
                   if (m12 >= discretization_mass(k) &
                        .AND. m12 <= discretization_mass(k+1)) then
                      fact = 1.0      
                   endif
                endif
          
                if (k == N_sizebin) then
                   if (discretization_mass(N_sizebin) <= m12) then
                      fact = 1.0 
                   endif
                endif
              
                do s = 1,N_groups
                    g = s !!Index_groups(s)           ! g : group number!
                   f1 = discretization_composition(i1, g,1) + random_vector(2 + s) &
                        *(discretization_composition(i1, g,2) - discretization_composition(i1, g, 1))
                   f2 = discretization_composition(i2, g, 1) + random_vector(2 +N_groups + s) &
                        * (discretization_composition(i2, g, 2) - discretization_composition(i2, g, 1))
                   
                   f12 = (m1 * f1 + m2 * f2) / m12

                   if (f12 <discretization_composition(i, g, 1) &
                        .OR. f12 > discretization_composition(i,g , 2)) then
                      fact = 0.0
                   endif
                enddo
                dsum = dsum + fact
             enddo
!$omp end parallel do
          
             dsum = dsum / dble(Nmc)
           

             if (dsum > 0.0) then
                repartition_coefficient_tmp(l) = dsum
             
                index1_repartition_coefficient_tmp(l) = j1
                index2_repartition_coefficient_tmp(l) = j2
           
                l = l + 1
             endif
          enddo
       enddo

       repartition_coefficient(j)%n = 0
       do l = 1,Nalloc_tmp
          if (repartition_coefficient_tmp(l) > 0.0) then
             repartition_coefficient(j)%n = repartition_coefficient(j)%n + 1
          endif
          
     enddo
       index1_repartition_coefficient(j)%n = repartition_coefficient(j)%n
       index2_repartition_coefficient(j)%n = repartition_coefficient(j)%n

       ! YK to avoid zero Ncoef
       if (repartition_coefficient(j)%n == 0) then
          repartition_coefficient(j)%n = 1
          repartition_coefficient_tmp(1) = 0.0
          index1_repartition_coefficient_tmp(1) = j
          index2_repartition_coefficient_tmp(1) = j
       endif

       allocate(repartition_coefficient(j)%arr(repartition_coefficient(j)%n))
       allocate(index1_repartition_coefficient(j)%arr(repartition_coefficient(j)%n))
       allocate(index2_repartition_coefficient(j)%arr(repartition_coefficient(j)%n))
          
       do l = 1,repartition_coefficient(j)%n
          repartition_coefficient(j)%arr(l) = repartition_coefficient_tmp(l)
          index1_repartition_coefficient(j)%arr(l) = index1_repartition_coefficient_tmp(l)
          index2_repartition_coefficient(j)%arr(l) = index2_repartition_coefficient_tmp(l)
       enddo

    enddo

  end SUBROUTINE  ssh_ComputeCoefficientRepartition

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Read repartition coefficients and their indexes from file.
  ! Size and composition discretization must correspond.
  subroutine ssh_ReadCoefficientRepartition(file, tag_file)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine reads coefficient repartition data
!     based on its directory and format (*.nc; *.txt; *.bin)
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     file: name of coefficient repartition file
!     tag_file: tag of file format 1=NetCDF 2=TXT 3=BIN
!
!------------------------------------------------------------------------   
    implicit none 

    integer, intent(in) :: tag_file
    double precision:: TempCoef
    integer::j, l, i, Ncoef, varid, dimid, ncid, i1, i2 , Ncfix
    character (len=100), intent(in) :: file
    ! This is the name of the data file we will read.
    character (len = *), parameter :: NCP_NAME = "Ncompute"
    double precision::sum_coef(N_size,N_size)
    type(ptr_to_real_array), dimension(N_size) :: repartition_coefficient_nc
    type(ptr_to_integer_array), dimension(N_size) :: index1_repartition_coefficient_nc
    type(ptr_to_integer_array), dimension(N_size) :: index2_repartition_coefficient_nc

       !TMP
       character( len = 5 ) :: cTemp
       character( len = 5 ) :: Temp
       character (len = 11):: DIM_NAME
       character (len = 11):: VAR_NAME

    ! NetCDF input
    if (tag_file .eq. 1) then
       if (ssh_standalone) write(*,*) " Read repartition coefficients and their indexes from NetCDF file."
       if (ssh_logger) write(logfile,*) " Read repartition coefficients and their indexes from NetCDF file."
       ! Open the file.
       call ssh_check( nf90_open(file, NF90_NOWRITE, ncid) )
       ! Allocate memory.
       allocate(index1_repartition_coefficient(N_size))
       allocate(index2_repartition_coefficient(N_size))
       allocate(repartition_coefficient(N_size))
      !OPEN(12,file='coeff_reff.txt')
      do i = 1, N_size
             write( cTemp,'(i4)' ) (i-1)
             Temp=trim(adjustl(cTemp ))
             DIM_NAME = 'Ncoef_'//Temp
              ! Get the varid of the data variable, based on its name.
              call ssh_check( nf90_inq_dimid(ncid, DIM_NAME, dimid) )
              ! Read the data.
              call ssh_check( nf90_inquire_dimension(ncid, dimid, len= Ncoef) )

	      ! Print the variables.
	      index1_repartition_coefficient_nc(i)%n = Ncoef
	      index2_repartition_coefficient_nc(i)%n = Ncoef
	      repartition_coefficient_nc(i)%n= Ncoef
              allocate(index1_repartition_coefficient_nc(i)%arr(Ncoef))
              allocate(index2_repartition_coefficient_nc(i)%arr(Ncoef))
              allocate(repartition_coefficient_nc(i)%arr(Ncoef))
                  VAR_NAME= 'index1_'//Temp
		  ! Get the varid of the data variable, based on its name.
		  call ssh_check( nf90_inq_varid(ncid, VAR_NAME, varid) )
		  ! Read the data.
		  call ssh_check( nf90_get_var(ncid, varid, index1_repartition_coefficient_nc(i)%arr))

                  VAR_NAME= 'index2_'//Temp
		  ! Get the varid of the data variable, based on its name.
		  call ssh_check( nf90_inq_varid(ncid, VAR_NAME, varid) )
		  ! Read the data.
		  call ssh_check( nf90_get_var(ncid, varid, index2_repartition_coefficient_nc(i)%arr) )
                  VAR_NAME= 'coef_'//Temp
		  ! Get the varid of the data variable, based on its name.
		  call ssh_check( nf90_inq_varid(ncid, VAR_NAME, varid) )
		  ! Read the data.
		  call ssh_check( nf90_get_var(ncid, varid, repartition_coefficient_nc(i)%arr) )
		  !Problem: In NetCDF index of grid start from 0; However, in bin and fortran index of grid start from 1
                  Ncfix=0
                  do l=1,repartition_coefficient_nc(i)%n
                  index1_repartition_coefficient_nc(i)%arr(l)=index1_repartition_coefficient_nc(i)%arr(l)+1
                  i1=index1_repartition_coefficient_nc(i)%arr(l)
                  index2_repartition_coefficient_nc(i)%arr(l)=index2_repartition_coefficient_nc(i)%arr(l)+1
                  i2=index2_repartition_coefficient_nc(i)%arr(l)
                  enddo
      end do
 
   do j = 1, N_size
         Ncfix=repartition_coefficient_nc(j)%n
         Ncoef=Ncfix
         do l=1,repartition_coefficient_nc(j)%n! ckeck symmetric case
	    i1=index1_repartition_coefficient_nc(j)%arr(l)
	    i2=index2_repartition_coefficient_nc(j)%arr(l)
	    if(i1/=i2) then
		Ncoef=Ncoef+1
	    endif
        enddo
        repartition_coefficient(j)%n=Ncoef
        index1_repartition_coefficient(j)%n=Ncoef
        index2_repartition_coefficient(j)%n=Ncoef
        allocate(repartition_coefficient(j)%arr(Ncoef))
        allocate(index1_repartition_coefficient(j)%arr(Ncoef))
        allocate(index2_repartition_coefficient(j)%arr(Ncoef))
        do l=1,repartition_coefficient_nc(j)%n! Reassign the coefficient and add missing symmetric coefficients
	    i1=index1_repartition_coefficient_nc(j)%arr(l)
	    index1_repartition_coefficient(j)%arr(l)=i1
	    i2=index2_repartition_coefficient_nc(j)%arr(l)
	    index2_repartition_coefficient(j)%arr(l)=i2
            Tempcoef=repartition_coefficient_nc(j)%arr(l)
	    repartition_coefficient(j)%arr(l)=Tempcoef
	    if(i1/=i2) then
	      Ncfix=Ncfix+1
	      index1_repartition_coefficient(j)%arr(Ncfix)=i2!set symmetric coefficient
	      index2_repartition_coefficient(j)%arr(Ncfix)=i1
	      repartition_coefficient(j)%arr(Ncfix)=Tempcoef
	      if(Ncfix>Ncoef) then
		  print*,"Error: CoefficientRepartition Aarry out of bounds...."
		  STOP
	      endif
	    endif
	  enddo
    end do
    sum_coef=0.d0
    do i = 1, N_size
          do l=1,repartition_coefficient(i)%n
             i1=index1_repartition_coefficient(i)%arr(l)
             i2=index2_repartition_coefficient(i)%arr(l)
             sum_coef(i1,i2)=sum_coef(i1,i2)+repartition_coefficient(i)%arr(l)
           enddo
    enddo
    
    ! TXT input
    elseif (tag_file.eq.2) then
    if (ssh_standalone) write(*,*) " Read repartition coefficients and their indexes from TXT file."
    if (ssh_logger) write(logfile,*) " Read repartition coefficients and their indexes from TXT file."
    open(11, file = file , status = "old")  
       allocate(index1_repartition_coefficient(N_size))
       allocate(index2_repartition_coefficient(N_size))
       allocate(repartition_coefficient(N_size))

    do j = 1,N_size
       read(11,*)repartition_coefficient_nc(j)%n
       allocate(repartition_coefficient_nc(j)%arr(repartition_coefficient_nc(j)%n))
       allocate(index1_repartition_coefficient_nc(j)%arr(repartition_coefficient_nc(j)%n))
       allocate(index2_repartition_coefficient_nc(j)%arr(repartition_coefficient_nc(j)%n))

       do l=1,repartition_coefficient_nc(j)%n 
          read(11,*)i1
          index1_repartition_coefficient_nc(j)%arr(l)=i1+1
          read(11,*)i2
          index2_repartition_coefficient_nc(j)%arr(l)=i2+1
          read(11,*)repartition_coefficient_nc(j)%arr(l)
       enddo
    enddo
    !reversed and ignore symmetric case
     do j = 1, N_size
         Ncfix=repartition_coefficient_nc(j)%n
         Ncoef=Ncfix
         do l=1,repartition_coefficient_nc(j)%n! ckeck symmetric case
	    i1=index1_repartition_coefficient_nc(j)%arr(l)
	    i2=index2_repartition_coefficient_nc(j)%arr(l)
	    if(i1/=i2) then
		Ncoef=Ncoef+1
	    endif
        enddo
        repartition_coefficient(j)%n=Ncoef
        index1_repartition_coefficient(j)%n=Ncoef
        index2_repartition_coefficient(j)%n=Ncoef
        allocate(repartition_coefficient(j)%arr(Ncoef))
        allocate(index1_repartition_coefficient(j)%arr(Ncoef))
        allocate(index2_repartition_coefficient(j)%arr(Ncoef))
        do l=1,repartition_coefficient_nc(j)%n! Reassign the coefficient and add missing symmetric coefficients
	    i1=index1_repartition_coefficient_nc(j)%arr(l)
	    index1_repartition_coefficient(j)%arr(l)=i1
	    i2=index2_repartition_coefficient_nc(j)%arr(l)
	    index2_repartition_coefficient(j)%arr(l)=i2
            Tempcoef=repartition_coefficient_nc(j)%arr(l)
	    repartition_coefficient(j)%arr(l)=Tempcoef
	    if(i1/=i2) then
	      Ncfix=Ncfix+1
	      index1_repartition_coefficient(j)%arr(Ncfix)=i2!set symmetric coefficient
	      index2_repartition_coefficient(j)%arr(Ncfix)=i1
	      repartition_coefficient(j)%arr(Ncfix)=Tempcoef
	      if(Ncfix>Ncoef) then
		  print*,"Error: CoefficientRepartition Aarry out of bounds...."
		  STOP
	      endif
	    endif
	  enddo
    end do
    

      do i = 1, N_size
	  TempCoef=0.d0
          do l=1,repartition_coefficient(i)%n
             i1=index1_repartition_coefficient(i)%arr(l)
             i2=index2_repartition_coefficient(i)%arr(l)
             TempCoef=TempCoef+repartition_coefficient(i)%arr(l)
           enddo
           if (ssh_standalone) write(*,*) 'coef sum',i,TempCoef
           if (ssh_logger) write(logfile,*) 'coef sum',i,TempCoef
    enddo
    close(11)

    ! BIN input
    else
    if (ssh_standalone) write(*,*) " Read repartition coefficients and their indexes from BIN file."
    if (ssh_logger) write(logfile,*) " Read repartition coefficients and their indexes from BIN file."
    open(11, file = file ,form="unformatted" , status = "old")  

    allocate(repartition_coefficient(N_size))
    allocate(index1_repartition_coefficient(N_size))
    allocate(index2_repartition_coefficient(N_size))
    
    do j = 1,N_size
       read(11) repartition_coefficient(j)%n
       index1_repartition_coefficient(j)%n = repartition_coefficient(j)%n
       index2_repartition_coefficient(j)%n = repartition_coefficient(j)%n
       allocate(repartition_coefficient(j)%arr(repartition_coefficient(j)%n))
       allocate(index1_repartition_coefficient(j)%arr(repartition_coefficient(j)%n))
       allocate(index2_repartition_coefficient(j)%arr(repartition_coefficient(j)%n))

       do l=1,repartition_coefficient(j)%n 
          read(11)index1_repartition_coefficient(j)%arr(l)
          i1=index1_repartition_coefficient(j)%arr(l)
          read(11)index2_repartition_coefficient(j)%arr(l)
          i2=index2_repartition_coefficient(j)%arr(l)
          read(11)repartition_coefficient(j)%arr(l)
       enddo
    enddo
    close(11)
    endif
    
  end subroutine ssh_ReadCoefficientRepartition

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Write repartition coefficients and their indexes from file.
  ! Size and composition discretization must correspond.
  recursive subroutine ssh_WriteCoefficientRepartition(file, tag_file)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine writes coefficient repartition data
!     based on its directory and format (*.nc; *.txt; *.bin)
!   
!------------------------------------------------------------------------
!      
!     -- INPUT VARIABLES
!      
!     file: name of coefficient repartition file
!     tag_file: tag of file format 1=NetCDF 2=TXT 3=BIN
!      
!------------------------------------------------------------------------   
    implicit none

    ! Inputs
    integer, intent(in) :: tag_file
    character (len=100), intent(in) :: file

    ! Local variables
    integer, parameter :: write_unit = 11
    integer :: j, l
    logical :: is_file

    integer :: dimid, ncid
    integer :: index1_varid, index2_varid, coef_varid
    character (len = 11) :: dim_name, var_name
    character (len = 5) :: ctemp, temp
 
    ! NetCDF output
    if (tag_file.eq.1) then

       ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
       ! overwrite this file, if it already exists.
       call ssh_check( nf90_create(file, NF90_CLOBBER, ncid) )

       ! Define the dimensions. NetCDF will hand back an ID for each. 
       call ssh_check( nf90_def_dim(ncid, "Nmc", Nmc, dimid) )
       call ssh_check( nf90_def_dim(ncid, "Nsize", N_size, dimid) )
       call ssh_check( nf90_def_dim(ncid, "Nb", N_sizebin, dimid) )
       call ssh_check( nf90_def_dim(ncid, "Nc", N_fracmax, dimid) )
       do j = 1, N_size
          write(ctemp, '(i4)') (j-1)
          temp = trim(adjustl(ctemp))
          dim_name = 'Ncoef_'//temp

          call ssh_check( nf90_def_dim(ncid, dim_name, repartition_coefficient(j)%n, dimid) )

          ! Define the variables
          var_name = 'index1_'//temp
          call ssh_check( nf90_def_var(ncid, var_name, NF90_INT, dimid, index1_varid) )
          var_name = 'index2_'//temp
          call ssh_check( nf90_def_var(ncid, var_name, NF90_INT, dimid, index2_varid) )
          var_name = 'coef_'//temp
          call ssh_check( nf90_def_var(ncid, var_name, NF90_DOUBLE, dimid, coef_varid) )

          ! End define mode. This tells netCDF we are done defining metadata.
          call ssh_check( nf90_enddef(ncid) )

          ! Index must decrease by 1 because it increases by 1 in 
          ! ReadCoefficientRepartition()
          ! coagulation-coef program generates NetCDF file.
          ! Index starts from 1 in that generated file.
          do l = 1, repartition_coefficient(j)%n
             index1_repartition_coefficient(j)%arr(l) = &
                  index1_repartition_coefficient(j)%arr(l) - 1
             index2_repartition_coefficient(j)%arr(l) = &
                  index2_repartition_coefficient(j)%arr(l) - 1
          enddo

          ! Write 'index1' data
          call ssh_check( nf90_put_var(ncid, index1_varid, index1_repartition_coefficient(j)%arr) )
          ! Write 'index2' data
          call ssh_check( nf90_put_var(ncid, index2_varid, index2_repartition_coefficient(j)%arr) )
          ! Write 'coefficient' data
          call ssh_check( nf90_put_var(ncid, coef_varid, repartition_coefficient(j)%arr) )
          ! Change into define mode
          call ssh_check( nf90_redef(ncid) )

          ! Index must increase again by 1 after writing to the file.
          do l = 1, repartition_coefficient(j)%n
             index1_repartition_coefficient(j)%arr(l) = &
                  index1_repartition_coefficient(j)%arr(l) + 1
             index2_repartition_coefficient(j)%arr(l) = &
                  index2_repartition_coefficient(j)%arr(l) + 1
          enddo

       enddo
       
       ! Close the file.
       call ssh_check( nf90_close(ncid) )

       ! If we got this far, everything worked as expected. 
      if (ssh_standalone) write(*,*) "*** SUCCESS writing file "//file
      if (ssh_logger) write(logfile,*) "*** SUCCESS writing file "//file

    ! TXT
    else if (tag_file.eq.2) then
 
      ! Check if there is a file
      inquire(file = file, exist = is_file)
      ! If there is a file, it is replaced
      if (is_file) then
        open(write_unit, file = file , status = "replace")
      ! Otherwise it is created
      else
        open(write_unit, file = file , status = "new")
      endif

      ! Write data
      !   index is 0-based
      do j = 1, N_size
        write(write_unit,*) repartition_coefficient(j)%n
        do l = 1, repartition_coefficient(j)%n 
          write(write_unit,*) index1_repartition_coefficient(j)%arr(l) - 1
          write(write_unit,*) index2_repartition_coefficient(j)%arr(l) - 1
          write(write_unit,*) repartition_coefficient(j)%arr(l)
        enddo
      enddo

      ! Close file
      close(write_unit)

    ! BIN
    else if (tag_file.eq.3) then

      ! Check if there is a file
      inquire(file = file, exist = is_file)
      ! If there is a file, it is replaced
      if (is_file) then
        open(write_unit, file = file, form = "unformatted", status = "replace")
      ! Otherwise it is created
      else
        open(write_unit, file = file, form = "unformatted", status = "new")
      endif

      ! Write data
      !   index is 1-based
      do j = 1, N_size
        write(write_unit) repartition_coefficient(j)%n
        do l = 1, repartition_coefficient(j)%n
          write(write_unit) index1_repartition_coefficient(j)%arr(l)
          write(write_unit) index2_repartition_coefficient(j)%arr(l)
          write(write_unit) repartition_coefficient(j)%arr(l)
        enddo
      enddo

      ! Close file
      close(write_unit)

    ! Unknown fomat
    else

      if (ssh_standalone) write(*,*) "Warning: Unknown writer format. Fallback to TXT."
      if (ssh_logger) write(logfile,*) "Warning: Unknown writer format. Fallback to TXT."
      call ssh_WriteCoefficientRepartition(file, 2)

    endif

  end subroutine ssh_WriteCoefficientRepartition

  subroutine ssh_check(status)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine reads the netCDF file based on status
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     status: variables to be read from netCDF files
!
!------------------------------------------------------------------------ 
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
       if (ssh_standalone) write(*,*) "Error: NetCDF file. ", &
            trim(nf90_strerror(status))
       if (ssh_logger) write(logfile,*) "Error: NetCDF file. ", &
            trim(nf90_strerror(status))
       stop "NetCDF"
    end if
  end subroutine ssh_check
  
  subroutine ssh_check_repart_coeff()
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine checks the quality of coefficient repartition data
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!------------------------------------------------------------------------
    implicit none
    integer::k,i,j,l
    double precision:: sumcr(N_size,N_size)
    
    sumcr=0.d0
    do k=1,N_size 
       do l=1,repartition_coefficient(k)%n    
          i=index1_repartition_coefficient(k)%arr(l)! index of grid 1
          j=index2_repartition_coefficient(k)%arr(l)! index of grid 2
	  sumcr(i,j)=sumcr(i,j)+repartition_coefficient(k)%arr(l)
	  !print*,k,"=",i,"+",j,repartition_coefficient(k)%arr(l)*100.d0,"%"
      enddo
    enddo
    !do i=1,N_size
      !do j=1,N_size
	!print*,i,j,sumcr(i,j)
      !enddo
    !enddo
    
    do k=1,N_size 
       do l=1,repartition_coefficient(k)%n    
          i=index1_repartition_coefficient(k)%arr(l)! index of grid 1
          j=index2_repartition_coefficient(k)%arr(l)! index of grid 2
	  if(sumcr(i,j).ne.1.d0) then
! 	    print*,'repartition_coefficient not right'
! 	    stop
	    if(repartition_coefficient(k)%arr(l)*sumcr(i,j).gt.0.d0) then
	      repartition_coefficient(k)%arr(l)=repartition_coefficient(k)%arr(l)/sumcr(i,j)
	    endif
	  endif
      enddo
    enddo    
  end subroutine

  subroutine ssh_DeallocateCoefficientRepartition
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine Deallocate Coefficient Repartition arrays
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!------------------------------------------------------------------------  
    implicit none

    integer :: j

    do j=1,N_size
       deallocate (repartition_coefficient(j)%arr)
       deallocate (index1_repartition_coefficient(j)%arr)
       deallocate (index2_repartition_coefficient(j)%arr)
    end do

    deallocate (repartition_coefficient)
    deallocate (index1_repartition_coefficient)
    deallocate (index2_repartition_coefficient)
  end subroutine ssh_DeallocateCoefficientRepartition

end module bCoefficientRepartition
