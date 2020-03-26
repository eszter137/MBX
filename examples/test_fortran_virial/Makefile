FC = gfortran
FFLAGS = -lstdc++ -g

all: test_pbc_v2 test_pbc_v2_no_mww test_pbc_v2_no_www_mww test_pbc_v2_no_pips test_pbc_v2_no_grad_call  test_pbc_v2_no_mww_no_grad_call test_pbc_v2_no_www_mww_no_grad_call


test_pbc_v2: test_pbc_v2.f90 ../../install/lib/libmbx.so
	$(FC) test_pbc_v2.f90 ../../install/lib/libmbx.so $(FFLAGS) -o test_pbc_v2
test_pbc_v2_no_mww: test_pbc_v2_no_mww.f90 ../../install/lib/libmbx.so
	$(FC) test_pbc_v2_no_mww.f90 ../../install/lib/libmbx.so $(FFLAGS) -o test_pbc_v2_no_mww
test_pbc_v2_no_www_mww: test_pbc_v2_no_www_mww.f90 ../../install/lib/libmbx.so
	$(FC) test_pbc_v2_no_www_mww.f90 ../../install/lib/libmbx.so $(FFLAGS) -o test_pbc_v2_no_www_mww
test_pbc_v2_no_pips: test_pbc_v2_no_pips.f90 ../../install/lib/libmbx.so
	$(FC) test_pbc_v2_no_pips.f90 ../../install/lib/libmbx.so $(FFLAGS) -o test_pbc_v2_no_pips
test_pbc_v2_no_grad_call: test_pbc_v2_no_grad_call.f90 ../../install/lib/libmbx.so
	$(FC) test_pbc_v2_no_grad_call.f90 ../../install/lib/libmbx.so $(FFLAGS) -o test_pbc_v2_no_grad_call
test_pbc_v2_no_mww_no_grad_call: test_pbc_v2_no_mww_no_grad_call.f90 ../../install/lib/libmbx.so
	$(FC) test_pbc_v2_no_mww_no_grad_call.f90 ../../install/lib/libmbx.so $(FFLAGS) -o test_pbc_v2_no_mww_no_grad_call
test_pbc_v2_no_www_mww_no_grad_call: test_pbc_v2_no_www_mww_no_grad_call.f90 ../../install/lib/libmbx.so
	$(FC) test_pbc_v2_no_www_mww_no_grad_call.f90 ../../install/lib/libmbx.so $(FFLAGS) -o test_pbc_v2_no_www_mww_no_grad_call


clean: 
	rm -f  test_pbc_v2 test_pbc_v2_no_mww test_pbc_v2_no_www_mww test_pbc_v2_no_pips test_pbc_v2_no_grad_call test_pbc_v2_no_mww_no_grad_call test_pbc_v2_no_www_mww_no_grad_call



.PHONY: all clean
