!!-----------------------------------------------------------------------
!!     Copyright (C) 2012-2018, ENPC - EDF R&D - INERIS
!!     Author(s): Hilel Dergaoui, Edouard Debry, Karine Sartelet
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
!!    This module contains methods to compute particle coagulation rate
!!-----------------------------------------------------------------------
Module gCoagulation
  use dPhysicalbalance
  use bCoefficientRepartition
  use aInitialization
  use omp_lib, only: omp_get_thread_num  
  implicit none

contains

  subroutine Gain (distribution, coagulation_rate_gain,c_number)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine computes the gain rate of coagulation
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     c_number: aerosol number concentration(#/m^3)
!     distribution:number/mass distribution before the coagulation(~/m^3)
!
!     -- OUTPUT VARIABLES
!
!     coagulation_rate_gain: the gain rate of coagulation(~/m^3/s)
!------------------------------------------------------------------------   
    implicit none
    integer::k,i,j,l
    double precision::distribution(N_size)
    double precision::coagulation_rate_gain(N_size)
    double precision ::c_number(N_size)
    double precision :: gain_term
    
    coagulation_rate_gain = 0.d0
    

    do k=1,N_size 
      gain_term=0.d0
      do l=1,repartition_coefficient(k)%n
	!check all possible repartition_coefficient combination of grid 1 and grid 2 into grid k
	i=index1_repartition_coefficient(k)%arr(l)! index of grid 1
	j=index2_repartition_coefficient(k)%arr(l)! index of grid 2
! 	if(IsNaN(kernel_coagulation(i,j)) ) then
! 	  print*, "kernel_coagulation(",i,j,")=",kernel_coagulation(i,j)!The problem is here NaN
! 	  kernel_coagulation(i,j)=0.d0
! 	endif
	if(c_number(i).gt.0.d0) then
	  gain_term=gain_term+kernel_coagulation(j,i)*repartition_coefficient(k)%arr(l)&
	      *c_number(i)*distribution(j)
! 	    if(IsNaN(gain_term*0.d0)) then
! 	      print*,'Error! IsNaN(gain_term)',k,i,j,gain_term
! 	      print*,kernel_coagulation(j,i),repartition_coefficient(k)%arr(l),&
! 	      c_number(i),distribution(j)
! 	      stop
! 	    endif
	endif
      enddo
       
      coagulation_rate_gain(k) =gain_term
    enddo
  end subroutine Gain
  
  subroutine Loss(distribution, coagulation_rate_loss,c_number)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine computes the loss rate of coagulation
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     c_number: aerosol number concentration(#/m^3)
!     distribution:number/mass distribution before the coagulation(~/m^3)
!
!     -- OUTPUT VARIABLES
!
!     coagulation_rate_loss: the loss rate of coagulation(~/m^3/s)
!------------------------------------------------------------------------   
    implicit none 
    integer::j,i
    double precision:: distribution(N_size)
    double precision:: coagulation_rate_loss(N_size)
    double precision :: loss_term
    double precision ::c_number(N_size)
    
    coagulation_rate_loss= 0.d0
    
    do j=1,N_size 
      loss_term=0.d0
      do i=1,N_size
	if(c_number(i).gt.0.d0) then
! 	if(IsNaN(kernel_coagulation(j,i)) ) then
! 	  print*, "kernel_coagulation(",j,i,")=",kernel_coagulation(j,i)
! 	  kernel_coagulation(j,i)=0.d0
! 	endif
	!loss by coagulation between grid k and i
	  loss_term= loss_term + kernel_coagulation(j,i) &
		* c_number(i)
	endif
      enddo
!       if(IsNaN(loss_term*0.d0)) then
! 	print*,'Error! IsNaN(gain_term)',i,j
! 	loss_term=0.d0
! 	stop
!       endif
      coagulation_rate_loss(j) =loss_term*distribution(j)
    enddo
     
  end subroutine loss
 
  subroutine Rate(rate_number,rate_mass,c_number,c_mass)
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!     This subroutine computes coagulation rate for number and each species
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     c_number: aerosol number concentration(#/m^3)
!     c_mass: aerosol mass concentration(µg/m^3)
!
!     -- OUTPUT VARIABLES
!
!     rate_number: the coagulation rate of number concentration(#/m^3/s)
!     rate_mass: the loss rate of mass concentration(µg/m^3/s)
!------------------------------------------------------------------------   
    implicit none
    integer :: j,i,jesp
    double precision::distribution(N_size)
    double precision::coagulation_rate_loss(N_size)
    double precision::coagulation_rate_gain(N_size)
    double precision ::c_number(N_size)
    double precision ::c_mass(N_size,N_aerosol)
    double precision ::rate_number(N_size)
    double precision ::rate_mass(N_size,N_aerosol)
    double precision ::total_mass_tmp
	
    do i=1,(N_aerosol-1)!loop by species
      jesp=List_species(i)
      do j = 1,N_size! Reassigned distribution by mass of each species
	distribution(j) = c_mass(j,jesp)
      enddo
      call  Gain (distribution,coagulation_rate_gain,c_number)
      call  Loss (distribution,coagulation_rate_loss,c_number)
      !call  Ratebalance(coagulation_rate_gain,coagulation_rate_loss)
      do j = 1,N_size
	  rate_mass(j,jesp) = coagulation_rate_gain(j) - coagulation_rate_loss(j)!
	  if(IsNaN(rate_mass(j,jesp)*0.d0)) print*,"infinity/NaN mass",j,jesp,rate_mass(j,jesp)
      enddo
    enddo

    if(with_number.eq.1) then
      do  j = 1,N_size
	distribution(j) = c_number(j)
      enddo
      call  Gain (distribution, coagulation_rate_gain,c_number)
      call  Loss (distribution, coagulation_rate_loss,c_number)

      do j = 1,N_size
	rate_number(j) =0.5d0*coagulation_rate_gain(j) - coagulation_rate_loss(j)!
	if(IsNaN(rate_number(j)*0.d0)) print*,"infinity/NaN numb",j,rate_number(j)
      enddo

    else
      do  j = 1,N_size
	total_mass_tmp=0.d0
	do i=1,(N_aerosol-1)
	  total_mass_tmp=total_mass_tmp+rate_mass(j,jesp)
	enddo
	rate_number(j) = total_mass_tmp/cell_diam_av(j)
      enddo
    endif    
    
  end subroutine Rate
  
!   subroutine Ratebalance(rate_gain,rate_loss)
!     implicit none
!     integer :: j,jesp
!     double precision::rate_loss(N_size)
!     double precision::rate_gain(N_size)
!     double precision ::total_g,total_l,d_rate
!       total_g=0.d0
!       total_l=0.d0
!       do j = 1,N_size
! 	  total_g=total_g+rate_gain(j)
! 	  total_l=total_l+rate_loss(j)
!       enddo
!       d_rate=total_g-total_l
!       if(d_rate.ne.0.d0) rate_gain(N_size)=rate_gain(N_size)-d_rate
! 
!   end subroutine Ratebalance
 
end module gCoagulation
