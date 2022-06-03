!!-----------------------------------------------------------------------
!!     Copyright (C) 2020 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------
!!
!!     -- DESCRIPTION
!!    This module reads meteorological data.
!!-----------------------------------------------------------------------

module mod_meteo

  use aInitialization

  implicit none

contains

  subroutine ssh_read_meteo()

    implicit none

    integer :: i, j, count, ierr, nt_meteo, ind_t
    double precision, dimension(:), allocatable :: t_meteo
    double precision, dimension(:), allocatable :: temperature_in    
    double precision, dimension(:), allocatable :: pressure_in
    double precision, dimension(:), allocatable :: rh_in, sh_in
    double precision :: ctime, ratio, rh_tmp, sh_tmp
    
    ! Read meteo data
    if (imeteo) then
       open(unit = 10, file = meteo_file, status = "old")

       count = 0
       ierr = 0
       do while(ierr .eq. 0)
          read(10, *, iostat=ierr)
          if (ierr == 0) count = count + 1
       end do
       nt_meteo = count - 1 ! excluding the header

       allocate(t_meteo(nt_meteo))
       allocate(temperature_in(nt_meteo))
       allocate(pressure_in(nt_meteo))
       allocate(rh_in(nt_meteo))
       allocate(sh_in(nt_meteo))
       
       rewind 10
       read(10, *) ! Read a header line (#)
       do i = 1, nt_meteo
          read(10, *) t_meteo(i), temperature_in(i), pressure_in(i), rh_in(i), sh_in(i)
       enddo
       close(10)

       ! Interpolation
       do i = 1, nt
          if (nt_meteo == 1) then
             temperature_array(i) = temperature_in(1)
             pressure_array(i) = pressure_in(1)
             call ssh_compute_sh(rh_in(1), temperature_in(1), &
                  pressure_in(1), sh_in(1))
             humidity_array(i) = sh_in(1)
             relative_humidity_array(i) = rh_in(1)
          else
             if (t_meteo(1) > initial_time .or. t_meteo(nt_meteo) < final_time) then
                write(*,*) "Error: missing meteorological data."
                stop
             endif
             ctime = initial_time + (i - 1) * delta_t
             ind_t = 0
             do j = 1, nt_meteo - 1
                if (ctime >= t_meteo(j) .and. ctime < t_meteo(j + 1)) then
                   ind_t = j
                endif
             enddo
             if (ind_t == 0) then
                write(*,*) "Error: missing meteorological data."
                stop
             endif
             
             ! Interpolation ratio
             if (t_meteo(ind_t + 1) == t_meteo(ind_t)) then
                write(*,*) "Error: same times for meteo data.", &
                     t_meteo(ind_t + 1), t_meteo(ind_t)
             endif
             ratio = (ctime - t_meteo(ind_t)) / (t_meteo(ind_t + 1) - t_meteo(ind_t))
             temperature_array(i) = ratio * temperature_in(ind_t + 1) + &
                  (1 - ratio) * temperature_in(ind_t)
             pressure_array(i) = ratio * pressure_in(ind_t + 1) + &
                  (1 - ratio) * pressure_in(ind_t)
             rh_tmp = ratio * rh_in(ind_t + 1) + &
                  (1 - ratio) * rh_in(ind_t)
             sh_tmp = ratio * sh_in(ind_t + 1) + &
                  (1 - ratio) * sh_in(ind_t)
             call ssh_compute_sh(rh_tmp, temperature_array(i), &
                  pressure_array(i), sh_tmp)
             relative_humidity_array(i) = rh_tmp
             humidity_array(i) = sh_tmp
             
          endif
       enddo
       
    else
       ! Use constant values given in namelist.input
       call ssh_compute_sh(relative_humidity, temperature, pressure, humidity)
       do i = 1, nt
          temperature_array(i) = temperature
          pressure_array(i) = Pressure
          humidity_array(i) = humidity
          relative_humidity_array(i) = relative_humidity
       enddo
    endif
    
  end subroutine ssh_read_meteo

  ! 
  subroutine ssh_compute_sh(rh, temp, pres, sh)

    implicit none

    double precision :: rh, temp, pres, sh
    
    if  (rh .gt. 0d0 ) then
       if (      rh.lt.Threshold_RH_inf &
            .or. rh.gt.Threshold_RH_sup) then
          if (ssh_standalone) write(*,*) 'Warning : clipping relative humidity.'
          if (ssh_logger) write(logfile,*) 'Warning : clipping relative humidity.'
          rh = DMIN1(DMAX1(rh, Threshold_RH_inf), Threshold_RH_sup)
       endif
       call ssh_compute_psat_sh(rh, temp, pres, pressure_sat, sh)
    else
       call ssh_compute_psat_rh(sh, temp, pres, pressure_sat, rh)
       if (      rh.lt.Threshold_RH_inf &
            .or. rh.gt.Threshold_RH_sup) then
          if (ssh_standalone) write(*,*) 'Warning : clipping relative humidity.'
          if (ssh_logger) write(logfile,*) 'Warning : clipping relative humidity.'
          rh = DMIN1(DMAX1(rh, Threshold_RH_inf), Threshold_RH_sup)
          call ssh_compute_psat_sh(rh, temp, pres, pressure_sat, sh)
       endif
    end if
    
  end subroutine ssh_compute_sh
  
end module mod_meteo
