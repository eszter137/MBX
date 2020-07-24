FC = gfortran
FFLAGS = -lstdc++ -g

all: test_pbc_s2x3x2 test_pbc_ortho

test_pbc_s2x3x2: test_pbc_s2x3x2.f90  ../../../install/lib/libmbx.so
	$(FC) test_pbc_s2x3x2.f90 ../../../install/lib/libmbx.so $(FFLAGS) -o test_pbc_s2x3x2

test_pbc_ortho: test_pbc_ortho.f90 ../../../install/lib/libmbx.so
	$(FC) test_pbc_ortho.f90 ../../../install/lib/libmbx.so $(FFLAGS) -o test_pbc_ortho
	
clean: 
	rm -f test_pbc_s2x3x2 test_pbc_ortho

.PHONY: all clean
