C-----------------------------------------------------------------------
C     Copyright (C) 2019 CEREA (ENPC) - INERIS
C     SSH-aerosol is distributed under the GNU General Public License v3
C-----------------------------------------------------------------------

c     Function: compute_coagulation_free_transition
c
c     Computes coagulation kernels for monodispersed
c     aerosols in the free transition regime.
c     2005/3/23: cleaning (Bruno Sportisse).
c
c     Parameters:
c     dp - aerosol diameter (µm).
c     vmp - mean particle velocity (m/s).
c     stick - sticking probability 0< <1 ().
c     cdifp - diffusion coefficient (m^2/s).
c     deltap - particle Knudsen number (µm).
c
c     Returns:
c     kercg - coagulation kernel (m^3/s).
      subroutine ssh_compute_coagulation_free_transition(dp, cdifp, 
     $     deltap, vmp, stick, kercg)

      double precision pi
      parameter(pi = 3.141592653589d0)

      double precision dp, cdifp, deltap
      double precision vmp, stick, kercg

      double precision beta, dpp

      dpp = dp * 1.d-06         ! convert µm to m

      beta = 1.d0 /
     &     ( dp / (dp + deltap) ! adim
     &     + 8.d0 * stick       ! adim
     &     * cdifp              ! m2.s - 1
     &     / vmp                ! m.s - 1
     &     / dpp )              ! m

      kercg = 8.d0 * pi         ! adim
     &     * cdifp              ! m2.s - 1
     &     * dpp                ! m
     &     * beta               ! adim

      end


c     Function: compute_coagulation_free_molecular
c
c     Computes coagulation kernels for monodispersed
c     aerosols in the free molecular regime.
c     2005/3/23: cleaning (Bruno Sportisse).
c
c     Parameters:
c     dp - aerosol diameter (µm).
c     vmp - mean particle velocity (m/s).
c     stick - sticking probability 0< <1 ().
c
c     Returns:
c     kercg - coagulation kernel (m^3/s).
      subroutine ssh_compute_coagulation_free_molecular(dp, vmp, stick,
     $     kercg)

      double precision pi
      parameter(pi = 3.141592653589d0)

      double precision dp, vmp, stick, kercg

      double precision dpp

      dpp = dp * 1.d-06         ! convert µm to m

      kercg = pi                ! adim
     &     * dpp * dpp          ! m2
     &     * vmp                ! m.s - 1
     &     * stick              ! adim

      end


c     Function: compute_coagulation_continuous
c
c     Computes coagulation kernels for monodispersed
c     aerosols in the continuous regime.
c     2005/3/23: cleaning (Bruno Sportisse).
c
c     Parameters:
c     dp - aerosol diameter (µm).
c     cdifp - diffusion coefficient (m^2/s).
c
c     Returns:
c     kercg - coagulation kernel (m^3/s).
      subroutine ssh_compute_coagulation_continuous(dp, cdifp, kercg)

      double precision pi
      parameter(pi = 3.141592653589d0)

      double precision dp, cdifp, kercg

      double precision dpp

      dpp = dp * 1.d-06         ! convert µm to m

      kercg = 8.d0 * pi         ! adim
     &     * cdifp              ! m2.s - 1
     &     * dpp                ! m

      end
