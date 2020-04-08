FC = gfortran
FFLAGS = -lstdc++ -g

all: test_pbc_v2 test_pbc_v2_no_co2w test_pbc_v2_no_ww_www_co2ww


test_pbc_v2: test_pbc_v2.f90 ../../install/lib/libmbx.so
	$(FC) test_pbc_v2.f90 ../../install/lib/libmbx.so $(FFLAGS) -o test_pbc_v2
test_pbc_v2_no_co2w: test_pbc_v2_no_co2w.f90 ../../install/lib/libmbx.so
	$(FC) test_pbc_v2_no_co2w.f90 ../../install/lib/libmbx.so $(FFLAGS) -o test_pbc_v2_no_co2w
test_pbc_v2_no_ww_www_co2ww: test_pbc_v2_no_ww_www_co2ww.f90 ../../install/lib/libmbx.so
	$(FC) test_pbc_v2_no_ww_www_co2ww.f90 ../../install/lib/libmbx.so $(FFLAGS) -o test_pbc_v2_no_ww_www_co2ww

clean: 
	rm -f  test_pbc_v2 test_pbc_v2_no_co2w test_pbc_v2_no_ww_www_co2ww 



.PHONY: all clean
