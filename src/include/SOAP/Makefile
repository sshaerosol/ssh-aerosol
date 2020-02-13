include		../../../Makefile.hdr

RM_LIST =	*.a *.o *.e *.mod core 
#call_soap: soap.o unifac.o isofwd.o isorev.o isocom.o
#	ifort soap.o unifac.o isofwd.o isorev.o isocom.o -o call_soap -lstdc++


OBJS =		soap.o 

all:	libsoap.a 

libsoap.a:	$(OBJS)
	        $(RM) libsoap.a
		$(AR) r libsoap.a $(OBJS)

soap.o: soap.cpp
	$(REALCC) -o soap.o -c soap.cpp -I $(BLITZLIB)

clean:
	$(RM) $(RM_LIST)
