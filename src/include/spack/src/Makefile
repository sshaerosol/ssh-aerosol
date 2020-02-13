include		../../../../Makefile.hdr

RM_LIST =	*.a *.o *.e *.mod core spack_generator 
#call_soap: soap.o unifac.o isofwd.o isorev.o isocom.o
#	ifort soap.o unifac.o isofwd.o isorev.o isocom.o -o call_soap -lstdc++


OBJS =		soap.o

all:	LU_generator spack_input find_perm spack_generator

libspack.a:		
	( $(CD) spack/ && $(MAKE) )
	$(LN) spack/libspack.a .

LU_generator: LU_generator.cpp
	icc -o LU_generator LU_generator.cpp -I /ccc/cont004/home/ineris/couvidaf/blitz-0.10/ -I ../../SeldonData/ -I ../../Talos/

spack_input: spack_input.cpp
	icc -o spack_input spack_input.cpp -I /ccc/cont004/home/ineris/couvidaf/blitz-0.10/ -I ../../SeldonData/ -I ../../Talos/

find_perm: find_perm.cpp
	icc -o find_perm find_perm.cpp -I /ccc/cont004/home/ineris/couvidaf/blitz-0.10/ -I ../../SeldonData/ -I ../../Talos/

spack_generator.o:	spack_generator.f
		$(REALFC) -I. -o spack_generator.o $(F77FLAGS) -c $<

spack_generator:	libspack.a spack_generator.o 
		$(REALFC) -o spack_generator $(F77FLAGS) spack_generator.o -lspack -Lspack/

clean:
	(cd spack ; rm *.o libspack.a)
	$(RM) $(RM_LIST)

