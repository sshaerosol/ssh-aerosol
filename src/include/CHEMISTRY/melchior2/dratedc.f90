!------------------------------------------------------------------------
!     This file was automatically generated by SPACK.
!------------------------------------------------------------------------

subroutine ssh_dratedc                        ( &
    ns,nr,rk,y,dw)
 
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This routine computes the derivative of reaction  rates.
!     This routine is automatically generated by SPACK.
!     Mechanism: MELCHIOR2_poa.reacti
!     Species: _ciMELCHIOR2_poa.spe
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     Ns: chemical species number.
!     NR: reaction number.
!     RK: kinetic rates.
!     Y: chemical concentrations.
!
!     -- INPUT/OUTPUT VARIABLES
!
!     -- OUTPUT VARIABLES
!
!     DW: derivative of reaction rates wrt Y.
!
!------------------------------------------------------------------------
!
!     -- REMARKS
!
!------------------------------------------------------------------------
!
!     -- MODIFICATIONS
!
!------------------------------------------------------------------------
!
!     -- AUTHOR(S)
!
!     SPACK.
!
!------------------------------------------------------------------------
 
  implicit none
 
  integer nr,ns
  double precision rk(nr),y(ns),dw(nr,ns)
 
 
 
  dw(   1,  73) =  rk(   1) * Y(  88)
  dw(   1,  88) =  rk(   1) * Y(  73)
  dw(   2,  73) =  rk(   2) * Y(  87)
  dw(   2,  87) =  rk(   2) * Y(  73)
  dw(   3,  73) =  rk(   3) * Y(  74)
  dw(   3,  74) =  rk(   3) * Y(  73)
  dw(   4,  73) =  rk(   4) * Y(  83)
  dw(   4,  83) =  rk(   4) * Y(  73)
  dw(   5,  88) =  rk(   5) * Y(  83)
  dw(   5,  83) =  rk(   5) * Y(  88)
  dw(   6,  87) =  rk(   6) * Y(  74)
  dw(   6,  74) =  rk(   6) * Y(  87)
  dw(   7,  83) =  rk(   7) * Y(  74)
  dw(   7,  74) =  rk(   7) * Y(  83)
  dw(   8,  35) =  rk(   8) * Y(  74)
  dw(   8,  74) =  rk(   8) * Y(  35)
  dw(   9,  69) =  rk(   9) * Y(  74)
  dw(   9,  74) =  rk(   9) * Y(  69)
  dw(  10,  68) =  rk(  10) * Y(  74)
  dw(  10,  74) =  rk(  10) * Y(  68)
  dw(  11,  83) =  rk(  11) * Y(  83)
  dw(  11,  83) =  rk(  11) * Y(  83)
  dw(  12,  83) =  rk(  12) * Y(  83)
  dw(  12,  83) =  rk(  12) * Y(  83)
  dw(  13,  84) =  rk(  13) * Y(  83)
  dw(  13,  83) =  rk(  13) * Y(  84)
  dw(  14,  84) =  rk(  14) * Y(  35)
  dw(  14,  35) =  rk(  14) * Y(  84)
  dw(  15,  84) =  rk(  15) * Y(  88)
  dw(  15,  88) =  rk(  15) * Y(  84)
  dw(  16,  87) =  rk(  16) * Y(  84)
  dw(  16,  84) =  rk(  16) * Y(  87)
  dw(  17,  84) =  rk(  17) * Y(  87)
  dw(  17,  87) =  rk(  17) * Y(  84)
  dw(  18,  20) =  rk(  18)
  dw(  19,  20) =  rk(  19)
  dw(  20,  20) =  rk(  20)
  dw(  21,  88) =  rk(  21) * Y(  74)
  dw(  21,  74) =  rk(  21) * Y(  88)
  dw(  22,  55) =  rk(  22) * Y(  74)
  dw(  22,  74) =  rk(  22) * Y(  55)
  dw(  23,  88) =  rk(  23) * Y(  88)
  dw(  23,  88) =  rk(  23) * Y(  88)
  dw(  24,  87) =  rk(  24)
  dw(  25,  21) =  rk(  25) * Y(  81)
  dw(  25,  81) =  rk(  25) * Y(  21)
  dw(  26,  21) =  rk(  26) * Y(  74)
  dw(  26,  74) =  rk(  26) * Y(  21)
  dw(  27,  77) =  rk(  27) * Y(  88)
  dw(  27,  88) =  rk(  27) * Y(  77)
  dw(  28,  77) =  rk(  28) * Y(  83)
  dw(  28,  83) =  rk(  28) * Y(  77)
  dw(  29,  77) =  rk(  29) * Y(  84)
  dw(  29,  84) =  rk(  29) * Y(  77)
  dw(  30,  81) =  rk(  30) * Y(  77)
  dw(  30,  77) =  rk(  30) * Y(  81)
  dw(  31,  86) =  rk(  31) * Y(  77)
  dw(  31,  77) =  rk(  31) * Y(  86)
  dw(  32,  29) =  rk(  32) * Y(  74)
  dw(  32,  74) =  rk(  32) * Y(  29)
  dw(  33,  78) =  rk(  33) * Y(  88)
  dw(  33,  88) =  rk(  33) * Y(  78)
  dw(  34,  78) =  rk(  34) * Y(  83)
  dw(  34,  83) =  rk(  34) * Y(  78)
  dw(  35,  78) =  rk(  35) * Y(  78)
  dw(  35,  78) =  rk(  35) * Y(  78)
  dw(  36,  78) =  rk(  36) * Y(  84)
  dw(  36,  84) =  rk(  36) * Y(  78)
  dw(  37,  81) =  rk(  37) * Y(  78)
  dw(  37,  78) =  rk(  37) * Y(  81)
  dw(  38,  86) =  rk(  38) * Y(  78)
  dw(  38,  78) =  rk(  38) * Y(  86)
  dw(  39,  30) =  rk(  39) * Y(  74)
  dw(  39,  74) =  rk(  39) * Y(  30)
  dw(  40,  40) =  rk(  40) * Y(  74)
  dw(  40,  74) =  rk(  40) * Y(  40)
  dw(  41,  85) =  rk(  41) * Y(  88)
  dw(  41,  88) =  rk(  41) * Y(  85)
  dw(  42,  85) =  rk(  42) * Y(  84)
  dw(  42,  84) =  rk(  42) * Y(  85)
  dw(  43,  85) =  rk(  43) * Y(  83)
  dw(  43,  83) =  rk(  43) * Y(  85)
  dw(  44,  32) =  rk(  44) * Y(  74)
  dw(  44,  74) =  rk(  44) * Y(  32)
  dw(  45,  85) =  rk(  45) * Y(  81)
  dw(  45,  81) =  rk(  45) * Y(  85)
  dw(  46,  85) =  rk(  46) * Y(  86)
  dw(  46,  86) =  rk(  46) * Y(  85)
  dw(  47,  85) =  rk(  47) * Y(  87)
  dw(  47,  87) =  rk(  47) * Y(  85)
  dw(  48,  34) =  rk(  48)
  dw(  49,  34) =  rk(  49) * Y(  74)
  dw(  49,  74) =  rk(  49) * Y(  34)
  dw(  50,  59) =  rk(  50) * Y(  84)
  dw(  50,  84) =  rk(  50) * Y(  59)
  dw(  51,  41) =  rk(  51) * Y(  84)
  dw(  51,  84) =  rk(  51) * Y(  41)
  dw(  52,  62) =  rk(  52) * Y(  74)
  dw(  52,  74) =  rk(  52) * Y(  62)
  dw(  53,  39) =  rk(  53) * Y(  84)
  dw(  53,  84) =  rk(  53) * Y(  39)
  dw(  54,  65) =  rk(  54) * Y(  88)
  dw(  54,  88) =  rk(  54) * Y(  65)
  dw(  55,  65) =  rk(  55) * Y(  84)
  dw(  55,  84) =  rk(  55) * Y(  65)
  dw(  56,  65) =  rk(  56) * Y(  83)
  dw(  56,  83) =  rk(  56) * Y(  65)
  dw(  57,  42) =  rk(  57) * Y(  74)
  dw(  57,  74) =  rk(  57) * Y(  42)
  dw(  58,  19) =  rk(  58) * Y(  74)
  dw(  58,  74) =  rk(  58) * Y(  19)
  dw(  59,  59) =  rk(  59) * Y(  74)
  dw(  59,  74) =  rk(  59) * Y(  59)
  dw(  60,  41) =  rk(  60) * Y(  74)
  dw(  60,  74) =  rk(  60) * Y(  41)
  dw(  61,  80) =  rk(  61) * Y(  74)
  dw(  61,  74) =  rk(  61) * Y(  80)
  dw(  62,  64) =  rk(  62) * Y(  74)
  dw(  62,  74) =  rk(  62) * Y(  64)
  dw(  63,  15) =  rk(  63) * Y(  74)
  dw(  63,  74) =  rk(  63) * Y(  15)
  dw(  64,  63) =  rk(  64) * Y(  74)
  dw(  64,  74) =  rk(  64) * Y(  63)
  dw(  65,  16) =  rk(  65) * Y(  74)
  dw(  65,  74) =  rk(  65) * Y(  16)
  dw(  66,  54) =  rk(  66) * Y(  74)
  dw(  66,  74) =  rk(  66) * Y(  54)
  dw(  67,  53) =  rk(  67) * Y(  74)
  dw(  67,  74) =  rk(  67) * Y(  53)
  dw(  68,  56) =  rk(  68) * Y(  74)
  dw(  68,  74) =  rk(  68) * Y(  56)
  dw(  69,  33) =  rk(  69) * Y(  74)
  dw(  69,  74) =  rk(  69) * Y(  33)
  dw(  70,  56) =  rk(  70) * Y(  74)
  dw(  70,  74) =  rk(  70) * Y(  56)
  dw(  71,  80) =  rk(  71) * Y(  84)
  dw(  71,  84) =  rk(  71) * Y(  80)
  dw(  72,  64) =  rk(  72) * Y(  84)
  dw(  72,  84) =  rk(  72) * Y(  64)
  dw(  73,  81) =  rk(  73) * Y(  84)
  dw(  73,  84) =  rk(  73) * Y(  81)
  dw(  74,  86) =  rk(  74) * Y(  84)
  dw(  74,  84) =  rk(  74) * Y(  86)
  dw(  75,  59) =  rk(  75) * Y(  73)
  dw(  75,  73) =  rk(  75) * Y(  59)
  dw(  76,  41) =  rk(  76) * Y(  73)
  dw(  76,  73) =  rk(  76) * Y(  41)
  dw(  77,  39) =  rk(  77) * Y(  73)
  dw(  77,  73) =  rk(  77) * Y(  39)
  dw(  78,  40) =  rk(  78) * Y(  73)
  dw(  78,  73) =  rk(  78) * Y(  40)
  dw(  79,  53) =  rk(  79) * Y(  73)
  dw(  79,  73) =  rk(  79) * Y(  53)
  dw(  80,  81) =  rk(  80) * Y(  88)
  dw(  80,  88) =  rk(  80) * Y(  81)
  dw(  81,  86) =  rk(  81) * Y(  88)
  dw(  81,  88) =  rk(  81) * Y(  86)
  dw(  82,  81) =  rk(  82) * Y(  83)
  dw(  82,  83) =  rk(  82) * Y(  81)
  dw(  83,  86) =  rk(  83) * Y(  83)
  dw(  83,  83) =  rk(  83) * Y(  86)
  dw(  84,  81) =  rk(  84) * Y(  81)
  dw(  84,  81) =  rk(  84) * Y(  81)
  dw(  85,  86) =  rk(  85) * Y(  81)
  dw(  85,  81) =  rk(  85) * Y(  86)
  dw(  86,  86) =  rk(  86) * Y(  86)
  dw(  86,  86) =  rk(  86) * Y(  86)
  dw(  87,  86) =  rk(  87) * Y(  87)
  dw(  87,  87) =  rk(  87) * Y(  86)
  dw(  88,  31) =  rk(  88)
  dw(  89,  31) =  rk(  89) * Y(  74)
  dw(  89,  74) =  rk(  89) * Y(  31)
  dw(  90,  44) =  rk(  90) * Y(  74)
  dw(  90,  74) =  rk(  90) * Y(  44)
  dw(  91,  73) =  rk(  91)
  dw(  92,  87) =  rk(  92)
  dw(  93,  84) =  rk(  93)
  dw(  94,  84) =  rk(  94)
  dw(  95,  35) =  rk(  95)
  dw(  96,  69) =  rk(  96)
  dw(  97,  55) =  rk(  97)
  dw(  98,  80) =  rk(  98)
  dw(  99,  80) =  rk(  99)
  dw( 100,  64) =  rk( 100)
  dw( 101,  26) =  rk( 101)
  dw( 102,  54) =  rk( 102)
  dw( 103,  16) =  rk( 103)
  dw( 104,  15) =  rk( 104)
  dw( 105,  63) =  rk( 105)
  dw( 106,  20) =  rk( 106)
  dw( 107,  56) =  rk( 107)
  dw( 108,  33) =  rk( 108)
  dw( 109,  31) =  rk( 109)
  dw( 110,  32) =  rk( 110)
  dw( 111,  30) =  rk( 111)
  dw( 112,  29) =  rk( 112)
  dw( 113,  24) =  rk( 113)
  dw( 114,  36) =  rk( 114) * Y(  84)
  dw( 114,  84) =  rk( 114) * Y(  36)
  dw( 115,  45) =  rk( 115) * Y(  84)
  dw( 115,  84) =  rk( 115) * Y(  45)
  dw( 116,  57) =  rk( 116) * Y(  84)
  dw( 116,  84) =  rk( 116) * Y(  57)
  dw( 117,  37) =  rk( 117) * Y(  84)
  dw( 117,  84) =  rk( 117) * Y(  37)
  dw( 118,  36) =  rk( 118) * Y(  74)
  dw( 118,  74) =  rk( 118) * Y(  36)
  dw( 119,  45) =  rk( 119) * Y(  74)
  dw( 119,  74) =  rk( 119) * Y(  45)
  dw( 120,  57) =  rk( 120) * Y(  74)
  dw( 120,  74) =  rk( 120) * Y(  57)
  dw( 121,  37) =  rk( 121) * Y(  74)
  dw( 121,  74) =  rk( 121) * Y(  37)
  dw( 122,   4) =  rk( 122) * Y(  74)
  dw( 122,  74) =  rk( 122) * Y(   4)
  dw( 123,  36) =  rk( 123) * Y(  73)
  dw( 123,  73) =  rk( 123) * Y(  36)
  dw( 124,  45) =  rk( 124) * Y(  73)
  dw( 124,  73) =  rk( 124) * Y(  45)
  dw( 125,  57) =  rk( 125) * Y(  73)
  dw( 125,  73) =  rk( 125) * Y(  57)
  dw( 126,  37) =  rk( 126) * Y(  73)
  dw( 126,  73) =  rk( 126) * Y(  37)
  dw( 127,  39) =  rk( 127) * Y(  74)
  dw( 127,  74) =  rk( 127) * Y(  39)
  dw( 128,  39) =  rk( 128) * Y(  84)
  dw( 128,  84) =  rk( 128) * Y(  39)
  dw( 129,  82) =  rk( 129) * Y(  83)
  dw( 129,  83) =  rk( 129) * Y(  82)
  dw( 130,  48) =  rk( 130) * Y(  74)
  dw( 130,  74) =  rk( 130) * Y(  48)
  dw( 131,  82) =  rk( 131) * Y(  86)
  dw( 131,  86) =  rk( 131) * Y(  82)
  dw( 132,  82) =  rk( 132) * Y(  81)
  dw( 132,  81) =  rk( 132) * Y(  82)
  dw( 133,  82) =  rk( 133) * Y(  88)
  dw( 133,  88) =  rk( 133) * Y(  82)
  dw( 134,  82) =  rk( 134) * Y(  84)
  dw( 134,  84) =  rk( 134) * Y(  82)
  dw( 135,  75) =  rk( 135) * Y(  74)
  dw( 135,  74) =  rk( 135) * Y(  75)
  dw( 136,  75) =  rk( 136) * Y(  84)
  dw( 136,  84) =  rk( 136) * Y(  75)
  dw( 137,  75) =  rk( 137) * Y(  73)
  dw( 137,  73) =  rk( 137) * Y(  75)
  dw( 138,  76) =  rk( 138) * Y(  88)
  dw( 138,  88) =  rk( 138) * Y(  76)
  dw( 139,  76) =  rk( 139) * Y(  83)
  dw( 139,  83) =  rk( 139) * Y(  76)
  dw( 140,  76) =  rk( 140) * Y(  81)
  dw( 140,  81) =  rk( 140) * Y(  76)
  dw( 141,  76) =  rk( 141) * Y(  87)
  dw( 141,  87) =  rk( 141) * Y(  76)
  dw( 142,  61) =  rk( 142)
  dw( 143,  49) =  rk( 143) * Y(  74)
  dw( 143,  74) =  rk( 143) * Y(  49)
  dw( 144,  61) =  rk( 144) * Y(  74)
  dw( 144,  74) =  rk( 144) * Y(  61)
  dw( 145,  61) =  rk( 145) * Y(  84)
  dw( 145,  84) =  rk( 145) * Y(  61)
  dw( 146,  60) =  rk( 146) * Y(  74)
  dw( 146,  74) =  rk( 146) * Y(  60)
  dw( 147,  60) =  rk( 147) * Y(  84)
  dw( 147,  84) =  rk( 147) * Y(  60)
  dw( 148,   7) =  rk( 148) * Y(  74)
  dw( 148,  74) =  rk( 148) * Y(   7)
  dw( 149,  18) =  rk( 149) * Y(  74)
  dw( 149,  74) =  rk( 149) * Y(  18)
  dw( 150,  14) =  rk( 150) * Y(  74)
  dw( 150,  74) =  rk( 150) * Y(  14)
  dw( 151,  71) =  rk( 151) * Y(  83)
  dw( 151,  83) =  rk( 151) * Y(  71)
  dw( 152,  71) =  rk( 152) * Y(  86)
  dw( 152,  86) =  rk( 152) * Y(  71)
  dw( 153,  71) =  rk( 153) * Y(  81)
  dw( 153,  81) =  rk( 153) * Y(  71)
  dw( 154,  71) =  rk( 154) * Y(  88)
  dw( 154,  88) =  rk( 154) * Y(  71)
  dw( 155,  71) =  rk( 155) * Y(  84)
  dw( 155,  84) =  rk( 155) * Y(  71)
  dw( 156,  72) =  rk( 156) * Y(  83)
  dw( 156,  83) =  rk( 156) * Y(  72)
  dw( 157,  72) =  rk( 157) * Y(  86)
  dw( 157,  86) =  rk( 157) * Y(  72)
  dw( 158,  72) =  rk( 158) * Y(  81)
  dw( 158,  81) =  rk( 158) * Y(  72)
  dw( 159,  72) =  rk( 159) * Y(  88)
  dw( 159,  88) =  rk( 159) * Y(  72)
  dw( 160,  72) =  rk( 160) * Y(  84)
  dw( 160,  84) =  rk( 160) * Y(  72)
  dw( 161,  17) =  rk( 161) * Y(  74)
  dw( 161,  74) =  rk( 161) * Y(  17)
  dw( 162,  83) =  rk( 162)
  dw( 163,  87) =  rk( 163)
  dw( 164,  84) =  rk( 164)
  dw( 165,  20) =  rk( 165)
  dw( 166,   9) =  rk( 166) * Y(  74)
  dw( 166,  74) =  rk( 166) * Y(   9)
  dw( 167,  10) =  rk( 167) * Y(  74)
  dw( 167,  74) =  rk( 167) * Y(  10)
  dw( 168,   8) =  rk( 168) * Y(  74)
  dw( 168,  74) =  rk( 168) * Y(   8)

return
end

