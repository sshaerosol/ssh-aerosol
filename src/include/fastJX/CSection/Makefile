#
#  Makefile for Fast-J Averaged X-section Code
#
OBJ1 = wn10.f
OBJ2 = add_xsect.f
OBJ3 = add_FJX.f
FCR  = ifort
OPT  = -static -O2 -r8

wn10 : $(OBJ1)
	$(FCR) $(OPT) -o wn10 $(OBJ1)

add_xsect : $(OBJ2)
	$(FCR) $(OPT) -o add_xsect $(OBJ2)

add_FJX : $(OBJ3)
	$(FCR) $(OPT) -o add_FJX $(OBJ3)

clean :
	rm wn10 add_xsect add_FJX core


$(OBJ3) : cmn_bins.f
