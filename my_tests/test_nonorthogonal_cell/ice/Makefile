FC = gfortran
FFLAGS = -lstdc++ -g

all: test_pbc

test_pbc: test_pbc.f90 ../../../install/lib/libmbx.so
	$(FC) test_pbc.f90 ../../../install/lib/libmbx.so $(FFLAGS) -o test_pbc


clean: 
	rm -f test_pbc

.PHONY: all clean
