!------------------------------------------------------------------------
!     Copyright (C) 2001-2008, ENPC - INRIA - EDF R&D
!
!     This file is part of the air quality modeling system Polyphemus.
!
!     Polyphemus is developed in the INRIA - ENPC joint project-team
!     CLIME and in the ENPC - EDF R&D joint laboratory CEREA.
!
!     Polyphemus is free software; you can redistribute i and/or modify
!     it under the terms of the GNU General Public License as published
!     by the Free Software Foundation; either version 2 of the License,
!     or (at your option) any later version.
!
!     Polyphemus is distributed in the hope that it will be useful, but
!     WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
!     General Public License for more details.
!
!     For more information, visit the Polyphemus web site:
!     http://cerea.enpc.fr/polyphemus/
!------------------------------------------------------------------------

subroutine rates90                        ( &
    ns,nr,rk,y,w)
 
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This routine computes the reaction rates.
!     This routine is automatically generated by SPACK.
!     Mechanism: CB05_poa.reactions  
!     Species: _ciCB05_poa.species_
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
!     W: reaction rates.
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
  double precision rk(nr),y(ns),w(nr)
 
 
 
  w(  1) =  rk(  1) * Y( 90)
  w(  2) =  rk(  2) * Y( 82)
  w(  3) =  rk(  3) * Y( 87) * Y( 92)
  w(  4) =  rk(  4) * Y( 82) * Y( 90)
  w(  5) =  rk(  5) * Y( 82) * Y( 90)
  w(  6) =  rk(  6) * Y( 82) * Y( 92)
  w(  7) =  rk(  7) * Y( 90) * Y( 87)
  w(  8) =  rk(  8) * Y( 87)
  w(  9) =  rk(  9) * Y( 87)
  w( 10) =  rk( 10) * Y(  4)
  w( 11) =  rk( 11) * Y(  4)
  w( 12) =  rk( 12) * Y( 87) * Y( 83)
  w( 13) =  rk( 13) * Y( 87) * Y( 93)
  w( 14) =  rk( 14) * Y( 84)
  w( 15) =  rk( 15) * Y( 84)
  w( 16) =  rk( 16) * Y( 84) * Y( 92)
  w( 17) =  rk( 17) * Y( 84) * Y( 90)
  w( 18) =  rk( 18) * Y( 84) * Y( 90)
  w( 19) =  rk( 19) * Y( 25)
  w( 20) =  rk( 20) * Y( 25)
  w( 21) =  rk( 21) * Y( 92) * Y( 92)
  w( 22) =  rk( 22) * Y( 92) * Y( 90)
  w( 23) =  rk( 23) * Y( 83) * Y( 92)
  w( 24) =  rk( 24) * Y( 41)
  w( 25) =  rk( 25) * Y( 83) * Y( 41)
  w( 26) =  rk( 26) * Y( 41) * Y( 41)
  w( 27) =  rk( 27) * Y( 90) * Y( 83)
  w( 28) =  rk( 28) * Y( 83) * Y( 74)
  w( 29) =  rk( 29) * Y( 93) * Y( 92)
  w( 30) =  rk( 30) * Y( 93) * Y( 90)
  w( 31) =  rk( 31) * Y( 36)
  w( 32) =  rk( 32) * Y( 83) * Y( 36)
  w( 33) =  rk( 33) * Y( 93) * Y( 93)
  w( 34) =  rk( 34) * Y( 93) * Y( 93)
  w( 35) =  rk( 35) * Y( 40)
  w( 36) =  rk( 36) * Y( 83) * Y( 40)
  w( 37) =  rk( 37) * Y(  4)
  w( 38) =  rk( 38) * Y( 83)
  w( 39) =  rk( 39) * Y( 83) * Y( 82)
  w( 40) =  rk( 40) * Y( 83) * Y( 83)
  w( 41) =  rk( 41) * Y( 83) * Y( 83)
  w( 42) =  rk( 42) * Y( 83) * Y( 93)
  w( 43) =  rk( 43) * Y( 93) * Y( 82)
  w( 44) =  rk( 44) * Y( 40) * Y( 82)
  w( 45) =  rk( 45) * Y( 84) * Y( 82)
  w( 46) =  rk( 46) * Y( 84) * Y( 83)
  w( 47) =  rk( 47) * Y( 84) * Y( 93)
  w( 48) =  rk( 48) * Y( 84) * Y( 87)
  w( 49) =  rk( 49) * Y( 84) * Y( 84)
  w( 50) =  rk( 50) * Y( 36)
  w( 51) =  rk( 51) * Y( 74)
  w( 52) =  rk( 52) * Y( 25)
  w( 53) =  rk( 53) * Y( 91) * Y( 92)
  w( 54) =  rk( 54) * Y( 85) * Y( 92)
  w( 55) =  rk( 55) * Y( 91) * Y( 93)
  w( 56) =  rk( 56) * Y( 85) * Y( 93)
  w( 57) =  rk( 57) * Y( 91) * Y( 91)
  w( 58) =  rk( 58) * Y( 85) * Y( 85)
  w( 59) =  rk( 59) * Y( 91) * Y( 85)
  w( 60) =  rk( 60) * Y( 69) * Y( 83)
  w( 61) =  rk( 61) * Y( 69)
  w( 62) =  rk( 62) * Y( 61) * Y( 83)
  w( 63) =  rk( 63) * Y( 61)
  w( 64) =  rk( 64) * Y( 83) * Y( 78)
  w( 65) =  rk( 65) * Y( 83) * Y(  6)
  w( 66) =  rk( 66) * Y( 80) * Y( 92)
  w( 67) =  rk( 67) * Y( 80) * Y( 93)
  w( 68) =  rk( 68) * Y( 80) * Y( 80)
  w( 69) =  rk( 69) * Y( 44) * Y( 83)
  w( 70) =  rk( 70) * Y( 44)
  w( 71) =  rk( 71) * Y( 23) * Y( 83)
  w( 72) =  rk( 72) * Y( 81) * Y( 83)
  w( 73) =  rk( 73) * Y( 81)
  w( 74) =  rk( 74) * Y( 81)
  w( 75) =  rk( 75) * Y( 81) * Y( 82)
  w( 76) =  rk( 76) * Y( 81) * Y( 84)
  w( 77) =  rk( 77) * Y( 81) * Y( 93)
  w( 78) =  rk( 78) * Y( 43)
  w( 79) =  rk( 79) * Y( 43) * Y( 92)
  w( 80) =  rk( 80) * Y( 43) * Y( 93)
  w( 81) =  rk( 81) * Y( 56) * Y( 83)
  w( 82) =  rk( 82) * Y( 88) * Y( 82)
  w( 83) =  rk( 83) * Y( 88) * Y( 83)
  w( 84) =  rk( 84) * Y( 88) * Y( 84)
  w( 85) =  rk( 85) * Y( 88)
  w( 86) =  rk( 86) * Y( 89) * Y( 92)
  w( 87) =  rk( 87) * Y( 89) * Y( 90)
  w( 88) =  rk( 88) * Y( 31)
  w( 89) =  rk( 89) * Y( 31)
  w( 90) =  rk( 90) * Y( 89) * Y( 93)
  w( 91) =  rk( 91) * Y( 89) * Y( 80)
  w( 92) =  rk( 92) * Y( 89) * Y( 91)
  w( 93) =  rk( 93) * Y( 89) * Y( 89)
  w( 94) =  rk( 94) * Y( 62) * Y( 83)
  w( 95) =  rk( 95) * Y( 62)
  w( 96) =  rk( 96) * Y( 77) * Y( 83)
  w( 97) =  rk( 97) * Y( 71) * Y( 82)
  w( 98) =  rk( 98) * Y( 71) * Y( 83)
  w( 99) =  rk( 99) * Y( 71) * Y( 84)
  w(100) =  rk(100) * Y( 71)
  w(101) =  rk(101) * Y( 86) * Y( 92)
  w(102) =  rk(102) * Y( 86) * Y( 90)
  w(103) =  rk(103) * Y( 38)
  w(104) =  rk(104) * Y( 38)
  w(105) =  rk(105) * Y( 38) * Y( 83)
  w(106) =  rk(106) * Y( 86) * Y( 93)
  w(107) =  rk(107) * Y( 89) * Y( 80)
  w(108) =  rk(108) * Y( 89) * Y( 91)
  w(109) =  rk(109) * Y( 86) * Y( 86)
  w(110) =  rk(110) * Y( 86) * Y( 89)
  w(111) =  rk(111) * Y( 70) * Y( 83)
  w(112) =  rk(112) * Y( 39)
  w(113) =  rk(113) * Y( 39)
  w(114) =  rk(114) * Y( 39) * Y( 90)
  w(115) =  rk(115) * Y( 82) * Y( 58)
  w(116) =  rk(116) * Y( 83) * Y( 58)
  w(117) =  rk(117) * Y( 87) * Y( 58)
  w(118) =  rk(118) * Y( 84) * Y( 58)
  w(119) =  rk(119) * Y( 82) * Y( 55)
  w(120) =  rk(120) * Y( 83) * Y( 55)
  w(121) =  rk(121) * Y( 87) * Y( 55)
  w(122) =  rk(122) * Y( 84) * Y( 55)
  w(123) =  rk(123) * Y( 59) * Y( 82)
  w(124) =  rk(124) * Y( 59) * Y( 83)
  w(125) =  rk(125) * Y( 59) * Y( 87)
  w(126) =  rk(126) * Y( 59) * Y( 84)
  w(127) =  rk(127) * Y(  8) * Y( 83)
  w(128) =  rk(128) * Y( 24) * Y( 92)
  w(129) =  rk(129) * Y( 24)
  w(130) =  rk(130) * Y( 83) * Y( 63)
  w(131) =  rk(131) * Y( 63) * Y( 84)
  w(132) =  rk(132) * Y( 65) * Y( 90)
  w(133) =  rk(133) * Y( 65) * Y( 93)
  w(134) =  rk(134) * Y( 64)
  w(135) =  rk(135) * Y( 64) * Y( 83)
  w(136) =  rk(136) * Y( 64) * Y( 87)
  w(137) =  rk(137) * Y( 83) * Y(  9)
  w(138) =  rk(138) * Y( 83) * Y( 45)
  w(139) =  rk(139) * Y( 45)
  w(140) =  rk(140) * Y( 82) * Y( 67)
  w(141) =  rk(141) * Y( 83) * Y( 67)
  w(142) =  rk(142) * Y( 87) * Y( 67)
  w(143) =  rk(143) * Y( 84) * Y( 67)
  w(144) =  rk(144) * Y( 83) * Y( 68)
  w(145) =  rk(145) * Y( 87) * Y( 68)
  w(146) =  rk(146) * Y( 84) * Y( 68)
  w(147) =  rk(147) * Y( 68)
  w(148) =  rk(148) * Y( 10) * Y( 83)
  w(149) =  rk(149) * Y( 83) * Y(  5)
  w(150) =  rk(150) * Y( 83) * Y(  7)
  w(151) =  rk(151) * Y( 90) * Y( 67)
  w(152) =  rk(152) * Y( 93)
  w(153) =  rk(153) * Y( 90)
  w(154) =  rk(154) * Y( 84)
  w(155) =  rk(155) * Y( 25)
  w(156) =  rk(156) * Y( 46) * Y( 83)
  w(157) =  rk(157) * Y( 46) * Y( 87)
  w(158) =  rk(158) * Y( 46) * Y( 84)
  w(159) =  rk(159) * Y( 49) * Y( 83)
  w(160) =  rk(160) * Y( 49) * Y( 87)
  w(161) =  rk(161) * Y( 49) * Y( 84)
  w(162) =  rk(162) * Y( 79) * Y( 93)
  w(163) =  rk(163) * Y( 79) * Y( 80)
  w(164) =  rk(164) * Y( 79) * Y( 89)
  w(165) =  rk(165) * Y( 79) * Y( 92)
  w(166) =  rk(166) * Y( 79) * Y( 84)
  w(167) =  rk(167) * Y( 73) * Y( 93)
  w(168) =  rk(168) * Y( 73) * Y( 80)
  w(169) =  rk(169) * Y( 73) * Y( 89)
  w(170) =  rk(170) * Y( 73) * Y( 92)
  w(171) =  rk(171) * Y( 73) * Y( 84)
  w(172) =  rk(172) * Y( 75) * Y( 93)
  w(173) =  rk(173) * Y( 75) * Y( 80)
  w(174) =  rk(174) * Y( 75) * Y( 89)
  w(175) =  rk(175) * Y( 75) * Y( 84)
  w(176) =  rk(176) * Y( 75) * Y( 92)
  w(177) =  rk(177) * Y( 66) * Y( 83)
  w(178) =  rk(178) * Y( 66) * Y( 84)
  w(179) =  rk(179) * Y( 76) * Y( 92)
  w(180) =  rk(180) * Y( 76) * Y( 93)
  w(181) =  rk(181) * Y( 76) * Y( 80)
  w(182) =  rk(182) * Y( 76) * Y( 90)
  w(183) =  rk(183) * Y( 57)
  w(184) =  rk(184) * Y( 57) * Y( 83)
  w(185) =  rk(185) * Y( 57) * Y( 84)
  w(186) =  rk(186) * Y( 28)
  w(187) =  rk(187) * Y( 12) * Y( 83)
  w(188) =  rk(188) * Y( 60) * Y( 83)
  w(189) =  rk(189) * Y( 60) * Y( 84)
  w(190) =  rk(190) * Y( 60) * Y( 87)
  w(191) =  rk(191) * Y( 15) * Y( 83)
  w(192) =  rk(192) * Y( 16) * Y( 83)
  w(193) =  rk(193) * Y( 17) * Y( 83)
  w(194) =  rk(194) * Y( 46) * Y( 82)
  w(195) =  rk(195) * Y( 49) * Y( 82)
  w(196) =  rk(196) * Y( 60) * Y( 82)
  w(197) =  rk(197) * Y( 47)
  w(198) =  rk(198) * Y( 21)
  w(199) =  rk(199) * Y( 22)
  w(200) =  rk(200) * Y( 32)
  w(201) =  rk(201) * Y( 48) * Y( 47)
  w(202) =  rk(202) * Y( 48) * Y( 93)
  w(203) =  rk(203) * Y( 48) * Y( 92)
  w(204) =  rk(204) * Y( 32) * Y( 48)
  w(205) =  rk(205) * Y( 32) * Y( 93)
  w(206) =  rk(206) * Y( 50) * Y( 83)

return
end

