!!-----------------------------------------------------------------------
!!     Copyright (C) 2019, CEREA- ENPC - EDF R&D - INERIS
!!     Author(s): Youngseob Kim
!!
!!     This file is part of the Size Composition Resolved Aerosol Model (SCRAM), a
!!     component of the SSH-aerosol model.
!!
!!     SSH-aerosol is a free software; you can redistribute it and/or modify
!!     it under the terms of the GNU General Public License as published
!!     by the Free Software Foundation; either version 2 of the License,
!!     or (at your option) any later version.
!!
!!     SSH-aerosol is distributed in the hope that it will be useful, but
!!     WITHOUT ANY WARRANTY; without even the implied warranty of
!!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
!!     General Public License for more details.
!!
!!-----------------------------------------------------------------------
!!
!!     -- DESCRIPTION
!!    
!!-----------------------------------------------------------------------

MODULE mod_photolysis
  use aInitialization
  
  implicit none

  real(4), dimension(:, :, :), allocatable :: file_rates_real
  double precision, dimension(:, :, :), allocatable :: file_rates

  integer :: photolysis_date_min = 0 ! in seconds
  integer :: photolysis_delta_t = 1 ! in days
  
  double precision, dimension(:), allocatable :: time_angle_photolysis
  double precision, dimension(:), allocatable :: latitude_photolysis
  double precision, dimension(:), allocatable :: altitude_photolysis
  
contains
 
  subroutine read_photolysis()
    !------------------------------------------------------------------------
    !
    !     -- DESCRIPTION
    !     
    !     This subroutine is based on InitPhotolysis of Polair3DChemistry.cxx
    !------------------------------------------------------------------------
    !
    !     -- INPUT VARIABLES
    !
    !
    !------------------------------------------------------------------------
    
    implicit none
    integer :: j, k, r, t
    character (len=80) :: dir_photolysis 

    double precision :: time_angle
    integer :: angle_in, j_in, k_in
    double precision :: alpha_angle, alpha_y, alpha_z 
    double precision :: one_alpha_angle, one_alpha_y, one_alpha_z
    integer :: day, nb_days
    integer :: mypos

    ! Get time record to read binary fiels.
    ! Photolysis data at this record are be used.
    day = (current_time - photolysis_date_min) / 86400 / &
         photolysis_delta_t + 1
    
    ! Allocate arrays for photolysis
    call allocate_photolysis()    

    ! Interpolation
    ! Along z
    k_in = 1
    alpha_z = 0.d0

    ! Along y (latitude).
    j_in = int((latitude - latitude_min) / delta_latitude) + 1
    alpha_y = (latitude - latitude_min - dble(j_in - 1) * delta_latitude) &
         / delta_latitude

    ! Time angle.
    time_angle = dble(current_time) / 3600.d0 - 12.d0 + longitude / 15.d0
    nb_days = int(time_angle / 24.d0)
    time_angle = abs(time_angle - 24.d0 * dble(nb_days))
    if (time_angle > 12.d0) then
       time_angle = 24.d0 - time_angle
    end if
    
    angle_in = int((time_angle - time_angle_min) / delta_time_angle) + 1
    alpha_angle = (time_angle - time_angle_min - &
         dble(angle_in - 1) * delta_time_angle) / &
         delta_time_angle

    one_alpha_angle = 1.d0 - alpha_angle
    one_alpha_y = 1.d0 - alpha_y
    one_alpha_z = 1.d0 - alpha_z
    
    ! Read photolysis rate
    do r = 1, n_photolysis

       ! Read from binary files.
       open(unit = 50, file = trim(photolysis_dir)//"/"//&
            trim(photolysis_name(r))//".bin", &
            form = 'unformatted', access = 'stream', status = 'old')
       mypos = 1 + 4 * (n_time_angle * n_latitude * n_altitude) * &
            (day - 1)
       read(unit = 50, pos = mypos) &
            (((file_rates_real(t, j, k), k = 1, n_altitude), &
            j = 1, n_latitude), t = 1, n_time_angle)
       close(50)
       
       ! Conversion from real to double precision.
       do k = 1, n_altitude
          do j = 1, n_latitude
             do t = 1, n_time_angle
                file_rates(t, j, k) = file_rates_real(t, j, k)
             end do
          end do
       end do
       
       if (angle_in >= n_time_angle) then
          photolysis_rate(r) = 0.d0
       else
          photolysis_rate(r) = &
               one_alpha_z * one_alpha_y * one_alpha_angle &
               * file_rates(angle_in, j_in, k_in) &
               + alpha_z * one_alpha_y * one_alpha_angle &
               * file_rates(angle_in, j_in, k_in + 1) &
               + one_alpha_z * alpha_y * one_alpha_angle &
               * file_rates(angle_in, j_in + 1, k_in) &
               + alpha_z * alpha_y * one_alpha_angle &
               * file_rates(angle_in, j_in + 1, k_in + 1) &
               + one_alpha_z * one_alpha_y * alpha_angle &
               * file_rates(angle_in + 1, j_in, k_in) &
               + alpha_z * one_alpha_y * alpha_angle &
               * file_rates(angle_in + 1, j_in, k_in + 1) &
               + one_alpha_z * alpha_y * alpha_angle &
               * file_rates(angle_in + 1, j_in + 1, k_in) &
               + alpha_z * alpha_y * alpha_angle & 
               * file_rates(angle_in + 1, j_in + 1, k_in + 1)
       end if
    end do
    
    call deallocate_photolysis()
    
  end subroutine read_photolysis


  subroutine allocate_photolysis()

    implicit none

    integer :: i
    
    ! Set time angle for photolysis rate
    time_angle_min = 0.d0
    delta_time_angle = 1.d0
    n_time_angle = 9
    allocate(time_angle_photolysis(n_time_angle))
    do i = 1, n_time_angle    
       time_angle_photolysis(i) = time_angle_min + i * delta_time_angle
    end do

    ! Set latitude for photolysis rate
    latitude_min = 0.d0
    delta_latitude = 10.d0
    n_latitude = 10
    allocate(latitude_photolysis(n_latitude))
    do i = 1, n_latitude    
       latitude_photolysis(i) = latitude_min + i * delta_latitude
    end do

    ! Set altitude for photolysis rate
    allocate(altitude_photolysis(n_altitude))
    do i = 1, n_altitude
       altitude_photolysis(i) = altitude_photolysis_input(i)
    end do
        
    allocate(file_rates(n_time_angle, n_latitude, n_altitude))
    allocate(file_rates_real(n_time_angle, n_latitude, n_altitude))
    
  end subroutine allocate_photolysis

  subroutine deallocate_photolysis()

    if (allocated(time_angle_photolysis)) deallocate(time_angle_photolysis)
    if (allocated(latitude_photolysis)) deallocate(latitude_photolysis)
    if (allocated(altitude_photolysis)) deallocate(altitude_photolysis)
    if (allocated(file_rates))  deallocate(file_rates)
    if (allocated(file_rates_real))  deallocate(file_rates_real)
    
  end subroutine deallocate_photolysis
  
  
end MODULE mod_photolysis
