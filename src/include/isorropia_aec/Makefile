#-*-makefile-*-
# iso

include		../../../Makefile.hdr

RM_LIST =	*.a *.o *.e *.mod core 


OBJS =		isocom.o isofwd.o isorev.o   

all:		libiso.a

libiso.a:	$(OBJS)
	        $(RM) libiso.a
		$(AR) r libiso.a $(OBJS)

%.o:		%.f
		$(REALFC) -o $@ $(F77FLAGS) -c $<
clean:
		$(RM) $(RM_LIST)


