#-*-makefile-*-
# iso

include		../../../../../Makefile.hdr

RM_LIST =	*.a *.o *.e *.mod core 


OBJS =		cinet.o  generator.o  gestion.o  initcinet.o  initialisation.o  lectci.o  lectcinet.o  write_common_header.o  write_end.o  write_header.o


all:		libspack.a

libspack.a:	$(OBJS)
	        $(RM) libspack.a
		$(AR) r libspack.a $(OBJS)

%.o:		%.f
		$(REALFC) -o $@ $(F77FLAGS) -c $<
clean:
		$(RM) $(RM_LIST)


