C=======================================================================
C
C *** ISORROPIA CODE
C *** INCLUDE FILE 'ISRPIA.INC'
C *** THIS FILE CONTAINS THE DECLARATIONS OF THE GLOBAL CONSTANTS
C     AND VARIABLES. 
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NCOMP=5,NIONS=7,NGASAQ=3,NSLDS=9,NPAIR=13,NZSR=100,
     &           NERRMX=25)
C
C *** INPUT VARIABLES **************************************************
C
      INTEGER METSTBL
      COMMON /SSHINPT/ W(NCOMP), WAER(NCOMP), TEMP, RH, IPROB, METSTBL,
     &              NADJ

C
C *** WATER ACTIVITIES OF PURE SALT SOLUTIONS **************************
C
      COMMON /SSHZSR / AWAS(NZSR), AWSS(NZSR), AWAC(NZSR), AWSC(NZSR),
     &              AWAN(NZSR), AWSN(NZSR), AWSB(NZSR), AWAB(NZSR),
     &              AWSA(NZSR), AWLC(NZSR)

C
C *** DELIQUESCENCE RELATIVE HUMIDITIES ********************************
C
      INTEGER WFTYP
      COMMON /SSHDRH / DRH2SO4,  DRNH42S4, DRNAHSO4, DRNACL,   DRNANO3, 
     &              DRNA2SO4, DRNH4HS4, DRLC,     DRNH4NO3, DRNH4CL

      COMMON /SSHMDRH/ DRMLCAB,  DRMLCAS,  DRMASAN,  DRMG1,    DRMG2,
     &              DRMG3,    DRMH1,    DRMH2,    DRMI1,    DRMI2,
     &              DRMI3,    DRMQ1,    DRMR1,    DRMR2,    DRMR3,
     &              DRMR4,    DRMR5,    DRMR6,    DRMR7,    DRMR8,
     &              DRMR9,    DRMR10,   DRMR11,   DRMR12,   DRMR13,
     &              WFTYP

C
C *** VARIABLES FOR LIQUID AEROSOL PHASE *******************************
C
      DOUBLE PRECISION MOLAL, MOLALR, M0
      REAL IONIC
      LOGICAL CALAOU, CALAIN, FRST, DRYF
      COMMON /SSHIONS/ MOLAL(NIONS), MOLALR(NPAIR), GAMA(NPAIR),
     &              ZZ(NPAIR),
     &              Z(NIONS),     GAMOU(NPAIR),  GAMIN(NPAIR),M0(NPAIR),
     &              GASAQ(NGASAQ),
     &              EPSACT,       COH,           CHNO3,       CHCL,         
     &              WATER,        IONIC,         IACALC,      
     &              FRST,         CALAIN,        CALAOU,    DRYF

C
C *** VARIABLES FOR SOLID AEROSOL PHASE ********************************
C
      COMMON /SSHSALT/ CH2SO4,  CNH42S4, CNH4HS4, CNACL,   CNA2SO4, 
     &              CNANO3,  CNH4NO3, CNH4CL,  CNAHSO4, CLC

C
C *** VARIABLES FOR GAS PHASE ******************************************
C
      COMMON /SSHGAS / GNH3, GHNO3, GHCL 

C
C *** EQUILIBRIUM CONSTANTS ********************************************
C
      COMMON /SSHEQUK/ XK1, XK2, XK3, XK4, XK5, XK6, 
     &              XK7, XK8, XK9, XK10,
     &              XK11,XK12,XK13,XK14,XKW, 
     &              XK21,XK22,XK31,XK32,XK41,
     &              XK42

C
C *** MOLECULAR WEIGHTS ************************************************
C
      DOUBLE PRECISION IMW
      COMMON /SSHOTHR/ R, IMW(NIONS), WMW(NCOMP), SMW(NPAIR)

C
C *** SOLUTION/INFO VARIABLES ******************************************
C
      CHARACTER SCASE*15
      COMMON /SSHCASE/ SULRATW, SULRAT, SODRAT, SCASE

C
      COMMON /SSHSOLN/ EPS, MAXIT, NSWEEP, NDIV, ICLACT

C
C *** ERROR SYSTEM *****************************************************
C
      CHARACTER ERRMSG*40
      INTEGER   ERRSTK, NOFER   
      LOGICAL   STKOFL   
      COMMON /SSHEROR/ STKOFL, NOFER, ERRSTK(NERRMX), ERRMSG(NERRMX)

C
C *** GENERIC VARIABLES ************************************************
C
      CHARACTER VERSION*15
      COMMON /SSHCGEN/ GREAT, TINY, TINY2, ZERO, ONE, VERSION

C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS  
C
      DOUBLE PRECISION ORGANION, WATORG 
      COMMON /SSHORG/ ORGANION, WATORG
C
C *** END OF INCLUDE FILE **********************************************
C
