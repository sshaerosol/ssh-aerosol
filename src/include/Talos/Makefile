include		../../../Makefile.hdr

RM_LIST =	*.a *.o *.e *.mod core 
#call_soap: soap.o unifac.o isofwd.o isorev.o isocom.o
#	ifort soap.o unifac.o isofwd.o isorev.o isocom.o -o call_soap -lstdc++


OBJS =		Talos.o 

all:	libtalos.a 

libtalos.a:	$(OBJS)
	        $(RM) libtalos.a
		$(AR) r libtalos.a $(OBJS)

Talos.o: Talos.cpp
	$(REALCC) -o Talos.o -c Talos.cpp -I $(BLITZLIB)

clean:
	$(RM) $(RM_LIST)
