#-*-makefile-*-

# Cuisine et dependances ...

include		../Makefile.hdr

RM_LIST =	*.a *.o *.e *.mod core 

SSHPATH =       ../src/

LDFLAGS =	-L. -lssh-aerosol -lmodules -liso -latmo -lchem -lrdb -lsoap -ltalos -lnetcdff -I${SSHPATH} -L${SSHPATH}  \
			-L${NETCDFLIB} $(CDFLIBS) -L${HDF5LIB}

all:		launch

launch:		../src/libmodules.a ../src/libiso.a ../src/libatmo.a ../src/libchem.a ../src/librdb.a ../src/libsoap.a ../src/libtalos.a ../src/libssh-aerosol.a launch.f90  
		$(MPIFC) -o launch $(FFLAGS) launch.f90 ${LDFLAGS}

clean:
	$(RM) launch	


