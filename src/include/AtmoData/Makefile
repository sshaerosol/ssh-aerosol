#-*-makefile-*-
# iso

include		../../../Makefile.hdr

RM_LIST =	*.a *.o *.e *.mod core 


OBJS =		Aerosol.o  compute_bidisperse_coagulation.o  compute_binary_nucleation.o  compute_monodisperse_coagulation.o  compute_ternary_nucleation.o  compute_ternary_nucleation_merikanto.o  Meteorology.o

all:		libssh-atmo.a

libssh-atmo.a:	$(OBJS)
	        $(RM) lssh-atmo.a
		$(AR) r libssh-atmo.a $(OBJS)

%.o:		%.f
		$(REALFC) -o $@ $(F77FLAGS) -c $<
clean:
		$(RM) $(RM_LIST)


