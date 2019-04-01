      SUBROUTINE compute_tab(wavelength_tab,wavelength_tab_size)

      IMPLICIT NONE


      integer i,j,k,l,w
      integer Nre, Nimag, Ndiam, Nlegendre
      PARAMETER(Nre=88,Nimag=100,Ndiam=100)
      PARAMETER(Nlegendre=8)
      double precision index_re_class(Nre+1),
     $     index_imag_class(Nimag+1)
      double precision diameters(Ndiam+1)

      double precision pi
      PARAMETER (pi=3.1416)

      integer NPL
      PARAMETER(NPL=20000)
      double precision AL1(NPL)
      double precision Reff,Alb

c     Change lambda for an other wavelength.
      double precision lambda0
      integer lambda0_nm
c      PARAMETER (lambda0=3.00D-1)

      character*3 filen

      integer wavelength_tab_size
      double precision wavelength_tab(wavelength_tab_size)

      double precision qext, qsca

      real Qsca_save(Nre+1,Nimag+1,Ndiam+1),
     $     Qext_save(Nre+1,Nimag+1,Ndiam+1),
     $     Phase_Function_save(Nre+1,Nimag+1,Ndiam+1,Nlegendre)

      integer NUC, NUC_grid

      real hx

      double precision Mierayon
      double precision  Mieind_re
      double precision  Mieind_im
      double precision BB

      NUC_grid = 3
      NUC = 2

c     WRITE GRID TABULATION

      OPEN(NUC_grid,FILE='Mish_Grid_Mie.dat')

      do i=1,Nre+1
         index_re_class(i) = 1.11D0 + (i-1)*1D-2
      enddo
      do i=1,Nimag+1
         index_imag_class(i) = (i-1) * 4.4D-3
      enddo

      hx = log(20/1D-2) / Ndiam
      do i=1,Ndiam+1
         diameters(i) = 1D-2 * exp((i-1)*hx)
      enddo

      WRITE(NUC_grid,*) '#index_re_class= '
      WRITE(NUC_grid,*) (index_re_class(i),i=1,Nre+1)
      WRITE(NUC_grid,*) '#index_imag_class= '
      WRITE(NUC_grid,*) (index_imag_class(i),i=1,Nimag+1)
      WRITE(NUC_grid,*) '#diameters= '
      WRITE(NUC_grid,*) (diameters(i),i=1,Ndiam+1)


      do w=1,wavelength_tab_size

         lambda0 = wavelength_tab(w)
         lambda0_nm = lambda0*1.D3
         write(filen,'(I3)') lambda0_nm

      OPEN(NUC,FILE='Mish_efficiency_factors_tab_' // filen // '.dat')

         do i=1,Nre+1
            do j=1,Nimag+1
               do k=1,Ndiam+1
                  Mierayon= diameters(k)
                  Mieind_re=index_re_class(i)
                  Mieind_im=index_imag_class(j)

                  BB=0.d0
                  call MieMish(lambda0,4,Mierayon,BB,
     $                 Mieind_re,Mieind_im,Qext,Qsca,
     $                 Reff,Alb,AL1)

                  Qsca_save(i,j,k)=Qsca
                  Qext_save(i,j,k)=Qext
                  do l=1,Nlegendre
                     Phase_Function_save(i,j,k,l)=AL1(l)
                  enddo
               enddo
            enddo
         enddo

         WRITE(NUC,*) '# Qsca'
         WRITE(NUC,*) (((Qsca_save(i,j,k),k=1,Ndiam+1),j=1,Nimag+1)
     $        ,i=1,Nre+1)
         WRITE(NUC,*) '# Qext'
         WRITE(NUC,*) (((Qext_save(i,j,k),k=1,Ndiam+1),j=1,Nimag+1)
     $        ,i=1,Nre+1)
         WRITE(NUC,*) '# Phase function'
         WRITE(NUC,*) ((((Phase_Function_save(i,j,k,l),l=1,Nlegendre),
     $        k=1,Ndiam+1), j=1,Nimag+1),i=1,Nre+1)

      enddo

      END
