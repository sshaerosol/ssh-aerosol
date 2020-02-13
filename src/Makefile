#-*-makefile-*-

# Cuisine et dependances ...

include		../Makefile.hdr

RM_LIST =	*.a *.o *.e *.mod core 

LDFLAGS =	-L. -lmodules -liso -latmo -lchem -lrdb -lsoap -ltalos -lnetcdff \
			-L${NETCDFLIB} $(CDFLIBS) -L${HDF5LIB} -I $(BLITZINC)


.PHONY:		all libmodules main libiso libatmo libchem librdb libsoap libtalos

all:		ssh-aerosol libssh-aerosol.a

libssh-aerosol.a:	ModuleSSH-aerosol.o
	        $(RM) libssh-aerosol.a
		$(AR) r libssh-aerosol.a ModuleSSH-aerosol.o

libmodules:		
		( $(CD) include/Module/ && $(MAKE) )
		$(LN) include/Module/libmodules.a .

libiso:		
		( $(CD) include/isorropia_aec/ && $(MAKE) )
		$(LN) include/isorropia_aec/libiso.a .

libatmo:		
		( $(CD) include/AtmoData/ && $(MAKE) )
		$(LN) include/AtmoData/libatmo.a .

libchem:		
		( $(CD) include/CHEMISTRY/common/ && $(MAKE) )
		$(LN) include/CHEMISTRY/common/libchem.a .

librdb:		
		( $(CD) include/RDB/ && $(MAKE) )
		$(LN) include/RDB/librdb.a .

libsoap:		
		( $(CD) include/SOAP && $(MAKE) )
		$(LN) include/SOAP/libsoap.a .

libtalos:		
		( $(CD) include/Talos && $(MAKE) )
		$(LN) include/Talos/libtalos.a .

ssh-aerosol.o:	ssh-aerosol.f90 
		$(MPIFC) -I.. -I. -o $@ $(FFLAGS) -Iinclude/Module/ -c $<


ssh-aerosol:	libmodules libiso libatmo libchem librdb libsoap libtalos main ssh-aerosol.o
		$(MPIFC) -o ssh-aerosol $(FFLAGS) \
		ssh-aerosol.o $(LDFLAGS)

ModuleSSH-aerosol.o:	ModuleSSH-aerosol.f90 libmodules libiso libatmo libchem librdb libsoap libtalos main 
		$(MPIFC) -Iinclude/Module/ -o ModuleSSH-aerosol.o $(FFLAGS) \
		-c ModuleSSH-aerosol.f90 $(LDFLAGS)

clean:
	$(RM) $(RM_LIST)
	( $(CD) include/isorropia_aec     ; $(MAKE) clean )
	( $(CD) include/AtmoData     ; $(MAKE) clean )
	( $(CD) include/CHEMISTRY/common     ; $(MAKE) clean )	
	( $(CD) include/Module     ; $(MAKE) clean )
	( $(CD) include/RDB     ; $(MAKE) clean )	


