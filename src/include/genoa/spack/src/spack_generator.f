C------------------------------------------------------------------------
C Copyright (C) 2015, ENPC
C     Author(s):Bruno Sportisse, Sylvain Doré
C
C This file is part of the air quality modeling system Polyphemus.
C
C Polyphemus is developed in the INRIA project-team CLIME and in
C the ENPC - EDF R&D joint laboratory CEREA.
C
C Polyphemus is free software; you can redistribute it and/or modify it under
C the terms of the GNU General Public License as published by the Free
C Software Foundation; either version 2 of the License, or (at your option)
C any later version.
C
C Polyphemus is distributed in the hope that it will be useful, but WITHOUT
C ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
C FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
C more details.
C
C For more information, visit the Polyphemus web site:
C      http://cerea.enpc.fr/polyphemus/
C------------------------------------------------------------------------

      program commande
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Master code for SPACK:
C     Simplified Preprocessor for Atmospheric Chemical Kinetics.
C
C     The input files are:
C     - a file describing the mechanism in a symbolic way,
C     - a file of chemical species.
C
C     The output files are F77 routines describing the chemical
C     production terms.
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     -- INPUT/OUTPUT VARIABLES
C
C     -- OUTPUT VARIABLES
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Bruno Sportisse, CEREA, 2003.
C
C------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      include 'spack/parametre'
      include 'spack/ficcom'
      common/comprec/bpsave(4,nrmax)

      dimension y0(nespmax)

      character(len=256) :: config_file

      ! Parse Command line
      call get_command_argument(1, config_file)
      call get_command_argument(2, filespecies)
      call get_command_argument(3, filemeca)

!zhizhao      print *, 'Species file: ', filespecies(1:LEN_TRIM(filespecies))
!zhizhao      print *, 'Reactions file: ', filemeca(1:LEN_TRIM(filemeca))
!zhizhao      print *, 'Configuration file: ',
!zhizhao     &                            config_file(1:LEN_TRIM(config_file))

      ! Parse Config file
      call ssh_read_config_file(config_file)
      if (LEN_TRIM(mechanism_name).eq.0) then
         print *, 'ERROR: A mechanism name must be specified.'
         print *, 'Please specify the mechanism name in ', config_file
         stop 1
      endif

      ! Initialization of data
      ! NEQ has the right dimension.
      call ssh_lectdata(y0,neq,indicaq)

      do i=1,nrmax
         bpsave(1,i) =bp(1,i)
      enddo
      end


      ! This function was taken and adapted from:
      ! http://jblevins.org/log/control-file
      subroutine ssh_read_config_file(config_file)
          implicit none

          character(len=256) :: config_file

          ! Input related variables
          character(len=100) :: buffer, label
          integer :: pos_separator
          integer, parameter :: fh = 15
          integer :: ios = 0
          integer :: line = 0

          ! Control file variables
          character *16 mechanism_name
          common mechanism_name

          character *24 function_suffix
          common function_suffix

          logical aerosol_formation
          common aerosol_formation

          mechanism_name=""
          function_suffix=""
          aerosol_formation=.FALSE.

          open(fh, file=config_file)
          ! ios is negative if an end of record condition is encountered or if
          ! an endfile condition was detected.  It is positive if an error was
          ! detected.  ios is zero otherwise.

          do while (ios == 0)
             read(fh, '(A)', iostat=ios) buffer
             if (ios == 0) then
                line = line + 1
                if (LEN_TRIM(buffer).eq.0) then
                  cycle
                endif
                if (buffer(1:1).eq.'#') then
                  cycle
                endif

                pos_separator = scan(buffer, "=")
                label = buffer(1:pos_separator-1)
                buffer = buffer(pos_separator+1:)
                select case (label)
                case ('mechanism_name')
                   read(buffer, *, iostat=ios) mechanism_name
!zhizhao                   print *, '  |-- Read mechanism_name: ',
!zhizhao     &                       mechanism_name
                case ('function_suffix')
                   if (LEN_TRIM(buffer).ne.0) then
                       read(buffer, *, iostat=ios) function_suffix
                   endif
!zhizhao                   print *, '  |-- Read function_suffix: ',
!zhizhao     &                       function_suffix
                   if (LEN_TRIM(function_suffix).ne.0) then
                     function_suffix = "_"//function_suffix
                   end if
                case ('aerosol_formation')
                   read(buffer, *, iostat=ios) aerosol_formation
!zhizhao                   print *, '  |-- Read aerosol_formation: ',
!zhizhao     &                       aerosol_formation
                case default
!zhizhao                   print *, '  |-- /!\ Skipping invalid label at line',
!zhizhao     &                       line
                end select
             endif
          end do
          return
          end
