FC = gfortran
FFLAGS = -lstdc++ -g

all: test_pbc_s2x3x2

test_pbc_s2x3x2: test_pbc_s2x3x2.f90  ../../../install/lib/libmbx.so
	$(FC) test_pbc_s2x3x2.f90 ../../../install/lib/libmbx.so $(FFLAGS) -o test_pbc_s2x3x2

clean: 
	rm -f test_pbc_s2x3x2

.PHONY: all clean
