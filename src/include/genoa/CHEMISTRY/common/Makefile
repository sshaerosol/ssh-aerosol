#-*-makefile-*-
# iso

include		../../../../Makefile.hdr

RM_LIST =	*.a *.o *.e *.mod core 


OBJS =		angzen.edf.o  chem.o  hetrxn.o  dimensions.o  dratedc.o  fexprod.o  fexloss.o  kinetic.o  rates.o


all:		libssh-chem.a

libssh-chem.a:	$(OBJS)
	        $(RM) libssh-chem.a
		$(AR) r libssh-chem.a $(OBJS)

%.o:		%.f
		$(REALFC) -o $@ $(F77FLAGS) -I../../INC/ -c $<

%.o:            %.f90
		$(MPIFC) -I.. -I. -o $@ $(FFLAGS) -I../INC/ -c $<

clean:
		$(RM) $(RM_LIST)


