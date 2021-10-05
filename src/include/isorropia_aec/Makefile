#-*-makefile-*-
# iso

include		../../../Makefile.hdr

RM_LIST =	*.a *.o *.e *.mod core 


OBJS =		isocom.o isofwd.o isorev.o   

all:		libssh-iso.a

libssh-iso.a:	$(OBJS)
	        $(RM) libssh-iso.a
		$(AR) r libssh-iso.a $(OBJS)

%.o:		%.f
		$(REALFC) -o $@ $(F77FLAGS) -c $<
clean:
		$(RM) $(RM_LIST)


