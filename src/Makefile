#-*-makefile-*-

# Cuisine et dependances ...

include		../Makefile.hdr

RM_LIST =	*.a *.o *.e *.mod core 

LDFLAGS =	-L. -lssh-modules -lssh-iso -lssh-atmo -lssh-chem -lssh-rdb -lssh-soap -lssh-talos -lnetcdff \
			-L${NETCDFLIBS} $(CDFLIBS) -L${HDF5LIB} -I $(BLITZLIB) -I $(NETCDFINC)

.PHONY:		all libssh-modules main libssh-iso libssh-atmo libssh-chem libssh-rdb libssh-soap libssh-talos

all:		libssh-aerosol.a 
#ssh-aerosol

libssh-aerosol.a:	ModuleSSH-aerosol.o
	        $(RM) libssh-aerosol.a
		$(AR) r libssh-aerosol.a ModuleSSH-aerosol.o

libssh-modules:		
		( $(CD) include/Module/ && $(MAKE) )
		$(LN) include/Module/libssh-modules.a .

libssh-iso:		
		( $(CD) include/isorropia_aec/ && $(MAKE) )
		$(LN) include/isorropia_aec/libssh-iso.a .

libssh-atmo:		
		( $(CD) include/AtmoData/ && $(MAKE) )
		$(LN) include/AtmoData/libssh-atmo.a .

libssh-chem:		
		( $(CD) include/CHEMISTRY/common/ && $(MAKE) )
		$(LN) include/CHEMISTRY/common/libssh-chem.a .

libssh-rdb:		
		( $(CD) include/RDB/ && $(MAKE) )
		$(LN) include/RDB/libssh-rdb.a .

libssh-soap:		
		( $(CD) include/SOAP && $(MAKE) )
		$(LN) include/SOAP/libssh-soap.a .

libssh-talos:		
		( $(CD) include/Talos && $(MAKE) )
		$(LN) include/Talos/libssh-talos.a .

ssh-aerosol.o:	ssh-aerosol.f90 
		$(MPIFC) -I.. -I. -o $@ $(FFLAGS) -Iinclude/Module/ $(LDFLAGS) -c $<


#ssh-aerosol:	libssh-modules libssh-iso libssh-atmo libssh-chem libssh-rdb libssh-soap libssh-talos main ssh-aerosol.o
#		$(MPIFC) -Iinclude/Module/ -o ssh-aerosol $(FFLAGS) \
#		ssh-aerosol.o $(LDFLAGS)


ModuleSSH-aerosol.o:	ModuleSSH-aerosol.f90 libssh-modules libssh-iso libssh-atmo libssh-chem libssh-rdb libssh-soap libssh-talos main 
		$(MPIFC) -Iinclude/Module/ -o ModuleSSH-aerosol.o $(FFLAGS) \
		-c ModuleSSH-aerosol.f90 $(LDFLAGS)

clean:
	$(RM) $(RM_LIST)
	( $(CD) include/isorropia_aec     ; $(MAKE) clean )
	( $(CD) include/AtmoData     ; $(MAKE) clean )
	( $(CD) include/CHEMISTRY/common     ; $(MAKE) clean )	
	( $(CD) include/Module     ; $(MAKE) clean )
	( $(CD) include/RDB     ; $(MAKE) clean )	
	( $(CD) include/SOAP     ; $(MAKE) clean )

