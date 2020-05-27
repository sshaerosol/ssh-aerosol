#-*-makefile-*-
# iso

include		../../../Makefile.hdr

RM_LIST =	*.a *.o *.e *.mod core 

OBJS =		compute_density.o  euler_coupled.o  euler_coupled_norm.o  euler_mass.o  euler_number.o  euler_number_norm.o  hemen.o  hemen_norm.o  locate.o  moving_diam.o  redistribution.o  test_mass_nb.o

all:		libssh-rdb.a

libssh-rdb.a:	$(OBJS)
	        $(RM) libssh-rdb.a
		$(AR) r libssh-rdb.a $(OBJS)

%.o:		%.f
		$(REALFC) -o $@ $(F77FLAGS) -c $<

%.o:		%.f90 
		$(MPIFC) -I.. -I. -o $@ $(FFLAGS) -I../INC/ -c $<

clean:
		$(RM) $(RM_LIST)


