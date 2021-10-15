!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------

!C*************************************************
!C     
!C     LIST OF POINTERS                           
!C     
!C*************************************************
      INTEGER IH,INa,INH4,ICl,ISO4,IHSO4,INO3&
          ,IH2O,INH3,IHCl,IHNO3,IOH
      INTEGER SNaNO3,SNH4NO3,SNACl,SNH4Cl,SLC&
          ,SNa2SO4,SNH42S4,SNaHSO4,SNH4HS4
!C*****
      INTEGER ENa,ESO4,ENH4,ENO3,ECl,&
          EH2O,EH2O_layers
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
