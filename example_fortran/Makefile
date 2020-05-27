#-*-makefile-*-

# Cuisine et dependances ...

include		../Makefile.hdr

RM_LIST =	*.a *.o *.e *.mod core 

SSHPATH =       ../src/

LDFLAGS =	-L. -lssh-aerosol -lssh-modules -lssh-iso -lssh-atmo -lssh-chem -lssh-rdb -lssh-soap -lssh-talos -lnetcdff -I${SSHPATH} -L${SSHPATH}  \
			-L${NETCDFLIB} $(CDFLIBS) -L${HDF5LIB}

all:		launch

launch:		../src/libssh-modules.a ../src/libssh-iso.a ../src/libssh-atmo.a ../src/libssh-chem.a ../src/libssh-rdb.a ../src/libssh-soap.a ../src/libssh-talos.a ../src/libssh-aerosol.a launch.f90  
		$(MPIFC) -o launch $(FFLAGS) launch.f90 ${LDFLAGS}

clean:
	$(RM) launch	


