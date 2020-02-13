#-*-makefile-*-
# iso

include		../../../Makefile.hdr

RM_LIST =	*.a *.o *.e *.mod core 


OBJS =		Aerosol.o  compute_bidisperse_coagulation.o  compute_binary_nucleation.o  compute_monodisperse_coagulation.o  compute_ternary_nucleation.o  compute_ternary_nucleation_merikanto.o  Meteorology.o

all:		libatmo.a

libatmo.a:	$(OBJS)
	        $(RM) latmo.a
		$(AR) r libatmo.a $(OBJS)

%.o:		%.f
		$(REALFC) -o $@ $(F77FLAGS) -c $<
clean:
		$(RM) $(RM_LIST)


