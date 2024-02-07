#-*-makefile-*-
# iso

include		../../../Makefile.hdr

RM_LIST =	*.a *.o *.e *.mod core 

LDFLAGS =       -L${NETCDFLIBS} $(CDFLIBS) -L${HDF5LIB} -I $(NETCDFINC)

OBJS =		init.o ModuleThermodynamics.o repartition.o balance.o redis.o condens.o coagu.o m_soap.o bulk.o congret.o result.o discret.o ModuleAdaptstep.o ModuleEmissions.o ModulePhotolysis.o compute_general.o ModuleMeteo.o

all:		libssh-modules.a

libssh-modules.a:	$(OBJS) 
	        $(RM) libssh-modules.a
		$(AR) r libssh-modules.a $(OBJS)

init.o:         ModuleInitialization.F90 
		$(MPIFC) -I.. -I. -o $@ $(FFLAGS) -I../INC/ -c $<

repartition.o:  ModuleCoefficientRepartition.F90 
		$(MPIFC) -I.. -I. -o $@ $(FFLAGS) -I../INC/ $(LDFLAGS) -c $<

balance.o:      ModulePhysicalbalance.F90 
		$(MPIFC) -I.. -I. -o $@ $(FFLAGS) -I../INC/ -c $<

redis.o:        ModuleRedistribution.F90 
		$(MPIFC) -I.. -I. -o $@ $(FFLAGS) -I../INC/ -c $<

condens.o:      ModuleCondensation.F90 
		$(MPIFC) -I.. -I. -o $@ $(FFLAGS) -I../INC/ -c $<

coagu.o:        ModuleCoagulation.F90 
		$(MPIFC) -I.. -I. -o $@ $(FFLAGS) -I../INC/ -c $<

m_soap.o:       ModuleSOAP.F90 
		$(MPIFC) -I.. -I. -o $@ $(FFLAGS) -I../INC/ -c $<

bulk.o:         ModuleBulkequibrium.F90 
		$(MPIFC) -I.. -I. -o $@ $(FFLAGS) -I../INC/ -c $<

congret.o:      ModuleCongregation.F90 
		$(MPIFC) -I.. -I. -o $@ $(FFLAGS) -I../INC/ -c $<

result.o:       ModuleResultoutput.F90  
		$(MPIFC) -I.. -I. -o $@ $(FFLAGS) -I../INC/ -c $<

discret.o:      ModuleDiscretization.F90   
		$(MPIFC) -I.. -I. -o $@ $(FFLAGS) -I../INC/ -c $<

%.o:		%.f
		$(REALFC) -o $@ $(F77FLAGS) -c $<

%.o:		%.F90 init.o
		$(MPIFC) -I.. -I. -o $@ $(FFLAGS) -I../INC/ -c $<

clean:
		$(RM) $(RM_LIST)


