!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
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

  subroutine ssh_Gain (distribution, coagulation_rate_gain,c_number)
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
    integer::k,i,j,l,jesp
    double precision::distribution(N_size,N_aerosol_layers+1)
    double precision::coagulation_rate_gain(N_size,N_aerosol_layers+1)
    double precision ::c_number(N_size)
    double precision :: gain_term
    
    coagulation_rate_gain = 0.d0
    
    do k=1,N_size
      gain_term = 0.d0
      do l=1,repartition_coefficient(k)%n
	!check all possible repartition_coefficient combination of grid 1 and grid 2 into grid k
	i=index1_repartition_coefficient(k)%arr(l)! index of grid 1
	j=index2_repartition_coefficient(k)%arr(l)! index of grid 2
	if(c_number(i).gt.0.d0) then
	  gain_term=kernel_coagulation(j,i)*repartition_coefficient(k)%arr(l)&
	      *c_number(i)
        else
          gain_term = 0.d0
	endif
        do jesp=1,N_aerosol_layers+1
          coagulation_rate_gain(k,jesp) = coagulation_rate_gain(k,jesp) + gain_term*distribution(j,jesp)
        enddo
      enddo
       
    enddo
  end subroutine ssh_Gain
  
  subroutine ssh_Loss(distribution, coagulation_rate_loss,c_number)
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
    integer::j,i,jesp
    double precision:: distribution(N_size,N_aerosol_layers+1)
    double precision:: coagulation_rate_loss(N_size,N_aerosol_layers+1)
    double precision :: loss_term
    double precision ::c_number(N_size)
    
    coagulation_rate_loss= 0.d0
    
    do j=1,N_size 
      loss_term=0.d0
      do i=1,N_size
	if(c_number(i).gt.0.d0) then
	!loss by coagulation between grid k and i
	  loss_term= loss_term + kernel_coagulation(j,i) &
		* c_number(i)
	endif
      enddo
      do jesp=1,N_aerosol_layers+1
          coagulation_rate_loss(j,jesp) =loss_term*distribution(j,jesp)
      enddo
    enddo
     
  end subroutine ssh_loss
 
  subroutine ssh_Rate(rate_number,rate_mass,c_number,c_mass)
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
    double precision::coagulation_rate_loss(N_size,N_aerosol_layers+1)
    double precision::coagulation_rate_gain(N_size,N_aerosol_layers+1)
    double precision ::c_number(N_size)
    double precision ::c_mass(N_size,N_aerosol_layers)
    double precision ::distribution(N_size,N_aerosol_layers+1)
    double precision ::rate_number(N_size)
    double precision ::rate_mass(N_size,N_aerosol_layers)
    double precision ::total_mass_tmp

    distribution = 0.0
    do jesp=1,N_aerosol_layers !loop by species
      do j = 1,N_size
         distribution(j,jesp) = c_mass(j,jesp)
      enddo
    enddo
    do  j = 1,N_size
	distribution(j,N_aerosol_layers+1) = c_number(j)
    enddo
    call  ssh_Gain (distribution,coagulation_rate_gain,c_number)
    call  ssh_Loss (distribution,coagulation_rate_loss,c_number)
    do jesp=1,N_aerosol_layers!loop by species
      do j = 1,N_size
	  rate_mass(j,jesp) = rate_mass(j,jesp) + coagulation_rate_gain(j,jesp) - coagulation_rate_loss(j,jesp)
      enddo
    enddo

    do j = 1,N_size
	rate_number(j) = rate_number(j) + 0.5d0*coagulation_rate_gain(j,N_aerosol_layers+1) - &
                                          coagulation_rate_loss(j,N_aerosol_layers+1)
    enddo

    
  end subroutine ssh_Rate
  
end module gCoagulation
