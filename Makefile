#-*-makefile-*-

include		./Makefile.hdr

LDFLAGS	=	 -L${NETCDFLIB} 

.PHONY:		src clean

RM_LIST =	*.a *.o *.e *.mod core 

src:		
		( $(CD) src ; $(MAKE) all )

clean:
		$(RM) $(RM_LIST)
		( $(CD) src    ; $(MAKE) clean )
