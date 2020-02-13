#-*-makefile-*-
# iso

include		../../../../Makefile.hdr

RM_LIST =	*.a *.o *.e *.mod core 


OBJS =		angzen.edf.o  chem.o  hetrxn.o  roschem.o  solvlin.o dimensions.o  dratedc.o  fexchem.o  jacdchemdc.o  kinetic.o  LU_decompose.o  LU_solve.o  LU_solve_tr.o  rates.o

all:		libchem.a

libchem.a:	$(OBJS)
	        $(RM) libchem.a
		$(AR) r libchem.a $(OBJS)

%.o:		%.f
		$(REALFC) -o $@ $(F77FLAGS) -I../../INC/ -c $<

%.o:            %.f90
		$(MPIFC) -I.. -I. -o $@ $(FFLAGS) -I../INC/ -c $<

clean:
		$(RM) $(RM_LIST)


