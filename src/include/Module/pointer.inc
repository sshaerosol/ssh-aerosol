!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------

!C*************************************************
!C     
!C     LIST OF POINTERS                           
!C     
!C*************************************************
      INTEGER EMD,EBC
      INTEGER ENa,ESO4,ENH4,ENO3,ECl&
          ,EBiA2D,EBiA1D,EBiA0D,EAGLY&
          ,EAMGLY,EBiMT,EBiPER,EAnBlP,EAnBmP&
          ,EBiBlP,EBiBmP,EBiNGA,EBiNIT3,EBiNIT&
          ,EBiMGA,EBiDER,ESOAhP,EPOAlP,EPOAmP&
          ,EAnCLP,ESOAlP,ESOAmP,EPOAhP,EH2O&
          ,EMonomer,EDimer,EBiA3D,EorgNIT,EH2O_layers
      INTEGER IH,INa,INH4,ICl,ISO4,IHSO4,INO3&
          ,IH2O,INH3,IHCl,IHNO3,IOH
      INTEGER SNaNO3,SNH4NO3,SNACl,SNH4Cl,SLC&
          ,SNa2SO4,SNH42S4,SNaHSO4,SNH4HS4
!C*****
      ! PARAMETER ( EMD =1,&
      !      EBC =2,&
      !      ENa =3,&
      !      ESO4=4,&
      !      ENH4=5,&
      !      ENO3=6,&
      !      ECl =7,&
      !      EBiA2D=8,&
      !      EBiA1D=9,&
      !      EBiA0D=10,&
      !      EAGLY=11,&
      !      EAMGLY=12,&
      !      EBiMT=13,&
      !      EBiPER=14,&
      !      EBiDER=15,&
      !      EBiMGA=16,&
      !      EAnBlP=17,&
      !      EAnBmP= 18,&
      !      EBiBlP= 19,&
      !      EBiBmP= 20,&
      !      EBiNGA= 21,&
      !      EBiNIT3=22,&
      !      EBiNIT=23,&
      !      EAnCLP=24,&
      !      EMonomer=25,&
      !      EDimer=26,&
      !      EBiA3D=27,&
      !      EorgNIT=28,&
      !      ESOAlP=29,&
      !      ESOAmP=30,&
      !      ESOAhP=31,&
      !      EPOAlP= 32,&
      !      EPOAmP= 33,&
      !      EPOAhP= 34,&
      !      EH2O=35)
      
      PARAMETER (IH   =1,&
          INa  =2,&
          INH4 =3,&
          ICl  =4,&
          ISO4 =5,&
          IHSO4=6,&
          INO3 =7,&
          IH2O =8,&
          INH3 =9,&
          IHCl =10,&
          IHNO3=11,&
          IOH  =12)
      
      PARAMETER (SNaNO3 =13,&
          SNH4NO3=14,&
          SNACl  =15,&
          SNH4Cl =16,&
          SNa2SO4=17,&
          SNH42S4=18,&
          SNaHSO4=19,&
          SNH4HS4=20,&
          SLC    =21)

!C*************************************************
