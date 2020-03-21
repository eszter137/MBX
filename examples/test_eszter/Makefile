FC = gfortran
FFLAGS = -lstdc++ -g

all: test_gas_phase test_pbc test_pbc_v2 test_pbc_v3 test_pbc_v4 test_pbc_v5

test_pbc: test_pbc.f90 ../../install/lib/libmbx.so
	$(FC) test_pbc.f90 ../../install/lib/libmbx.so $(FFLAGS) -o test_pbc

test_pbc_v2: test_pbc_v2.f90 ../../install/lib/libmbx.so
	$(FC) test_pbc_v2.f90 ../../install/lib/libmbx.so $(FFLAGS) -o test_pbc_v2

test_pbc_v3: test_pbc_v3.f90 ../../install/lib/libmbx.so
	$(FC) test_pbc_v3.f90 ../../install/lib/libmbx.so $(FFLAGS) -o test_pbc_v3

test_pbc_v4: test_pbc_v4.f90 ../../install/lib/libmbx.so
	$(FC) test_pbc_v4.f90 ../../install/lib/libmbx.so $(FFLAGS) -o test_pbc_v4

test_pbc_v5: test_pbc_v5.f90 ../../install/lib/libmbx.so
	$(FC) test_pbc_v5.f90 ../../install/lib/libmbx.so $(FFLAGS) -o test_pbc_v5

test_gas_phase: test_gas_phase.f90 ../../install/lib/libmbx.so
	$(FC) test_gas_phase.f90 ../../install/lib/libmbx.so $(FFLAGS) -o test_gas_phase

clean: 
	rm -f test_gas_phase test_pbc test_pbc_v2 test_pbc_v3 test_pbc_v4 test_pbc_v5

.PHONY: all clean
