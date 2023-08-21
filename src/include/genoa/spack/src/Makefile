include         ../../../../Makefile.hdr

RM_LIST =	*.a *.o *.e *.mod core spack_generator LU_generator spack_input spack_generator find_perm
#call_soap: soap.o unifac.o isofwd.o isorev.o isocom.o
#	ifort soap.o unifac.o isofwd.o isorev.o isocom.o -o call_soap -lstdc++


OBJS =		LU_generator spack_generator

all:	LU_generator spack_input find_perm spack_generator

libspack.a:
	( $(CD) spack/ && $(MAKE) )
	$(LN) spack/libspack.a .

LU_generator: LU_generator.cpp
	$(REALCC) -o LU_generator LU_generator.cpp -I $(BLITZLIB) -I ../../SeldonData/ -I ../../Talos/

spack_input: spack_input.cpp
	$(REALCC) -o spack_input spack_input.cpp -I $(BLITZLIB) -I ../../SeldonData/ -I ../../Talos/

find_perm: find_perm.cpp
	$(REALCC) -o find_perm find_perm.cpp -I $(BLITZLIB) -I ../../SeldonData/ -I ../../Talos/

spack_generator:	spack_generator.o libspack.a
		$(REALFC) -o spack_generator $(F77FLAGS) spack_generator.o -lspack -Lspack/

spack_generator.o:	spack_generator.f
		$(REALFC) -I. -o spack_generator.o $(F77FLAGS) -c $<

clean:
	(rm $(RM_LIST) ; cd spack ; rm *.o libspack.a)
	$(RM) $(RM_LIST)

