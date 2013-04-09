F90 = gfortran
FFLAGS = -O0 -g

all : reg_tests.exe analytic_nonlin.exe analytic_sys.exe

reg_tests.exe : reg_tests.f90 IRKSolve.f90
	$(F90) $(FFLAGS) $^ -o $@

analytic_nonlin.exe : analytic_nonlin.f90 IRKSolve.f90
	$(F90) $(FFLAGS) $^ -o $@

analytic_sys.exe : analytic_sys.f90 IRKSolve.f90
	$(F90) $(FFLAGS) $^ -o $@


clean :
	\rm -f *.o *~ *.exe
