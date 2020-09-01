FC = gfortran
FFLAGS = -lstdc++ -g

all: iceIV ice4R_ortho

iceIV: iceIV.f90  ../../../install/lib/libmbx.so
	$(FC) iceIV.f90 ../../../install/lib/libmbx.so $(FFLAGS) -o iceIV

ice4R_ortho: ice4R_ortho.f90 ../../../install/lib/libmbx.so
	$(FC) ice4R_ortho.f90 ../../../install/lib/libmbx.so $(FFLAGS) -o ice4R_ortho
	
clean: 
	rm -f iceIV ice4R_ortho

.PHONY: all clean
