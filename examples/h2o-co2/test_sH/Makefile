FC = gfortran
FFLAGS = -lstdc++ -g

all: test_pbc test_pbc2 test_no3b test_no3b2

test_pbc: test_pbc.f90 ../../../install/lib/libmbx.so
	$(FC) test_pbc.f90 ../../../install/lib/libmbx.so $(FFLAGS) -o test_pbc

test_pbc2: test_pbc2.f90 ../../../install/lib/libmbx.so
	$(FC) test_pbc2.f90 ../../../install/lib/libmbx.so $(FFLAGS) -o test_pbc2

test_no3b: test_no3b.f90 ../../../install/lib/libmbx.so
	$(FC) test_no3b.f90 ../../../install/lib/libmbx.so $(FFLAGS) -o test_no3b

test_no3b2: test_no3b2.f90 ../../../install/lib/libmbx.so
	$(FC) test_no3b2.f90 ../../../install/lib/libmbx.so $(FFLAGS) -o test_no3b2

clean: 
	rm -f test_pbc test_pbc2 test_no3b test_no3b2

.PHONY: all clean
