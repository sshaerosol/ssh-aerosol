
      subroutine write_common_header_90(nwrite)
      integer nwrite

      write(nwrite,200)
      write(nwrite,101)
 101  format('!',5x,'Copyright (C) 2001-2008, ENPC - INRIA - EDF R&D')
      write(nwrite,150)
      write(nwrite,102)
 102  format('!',5x,'This file is part of the air quality ',
     2      'modeling system Polyphemus.')
      write(nwrite,150)
      write(nwrite,103)
 103  format('!',5x,'Polyphemus is developed in the INRIA - ENPC ',
     2      'joint project-team')
      write(nwrite,104)
 104  format('!',5x,'CLIME and in the ENPC - EDF R&D joint ',
     2      'laboratory CEREA.')
      write(nwrite,150)
      write(nwrite,105)
 105  format('!',5x,'Polyphemus is free software; you can ',
     2      'redistribute i and/or modify')
      write(nwrite,106)
 106  format('!',5x,'it under the terms of the GNU General ',
     2     'Public License as published')
      write(nwrite,107)
 107  format('!',5x,'by the Free Software Foundation; either ',
     2     'version 2 of the License,')
      write(nwrite,108)
 108  format('!',5x,'or (at your option) any later version.')
      write(nwrite,150)
      write(nwrite,109)
 109  format('!',5x,'Polyphemus is distributed in the hope ',
     2     'that it will be useful, but')
      write(nwrite,110)
 110  format('!',5x,'WITHOUT ANY WARRANTY; without even ',
     2     'the implied warranty of')
      write(nwrite,111)
 111  format('!',5x,'MERCHANTABILITY or FITNESS FOR A ',
     2     'PARTICULAR PURPOSE. See the GNU')
      write(nwrite,112)
 112  format('!',5x,'General Public License for more details.')
      write(nwrite,150)
      write(nwrite,113)
 113  format('!',5x,'For more information, visit ',
     2     'the Polyphemus web site:')
      write(nwrite,114)
 114  format('!',5x,'http://cerea.enpc.fr/polyphemus/')
      write(nwrite,200)
      write(nwrite,300)


 150  format('!',6x,a65)
 200  format('!-----------------------------------',
     2     '-------------------------------------')
 300  format('')

      RETURN
      END
