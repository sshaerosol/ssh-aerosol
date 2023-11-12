include		${SSH_dir}/Makefile.hdr

RM_LIST =	*.a *.o *.e *.mod core 
#call_soap: soap.o unifac.o isofwd.o isorev.o isocom.o
#	ifort soap.o unifac.o isofwd.o isorev.o isocom.o -o call_soap -lstdc++


OBJS =		soap.o 

all:	libssh-soap.a 

libssh-soap.a:	$(OBJS)
	        $(RM) libssh-soap.a
		$(AR) r libssh-soap.a $(OBJS)

soap.o: soap.cpp
	$(CC) -o soap.o $(CCFLAGS) -c soap.cpp -I $(BLITZINC)

clean:
	$(RM) $(RM_LIST)
