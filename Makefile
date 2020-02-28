# -*- Makefile -*-

FC      = gfortran
OPTIM   = -g -C

FITSDIR = -L/usr/lib -lcfitsio
LAPACK  = -L/usr/lib -llapack -lblas
HEALPIX = -L/usr/local/src/Healpix_3.50/lib -lhealpix
HEALINC = -I/usr/local/src/Healpix_3.50/include
OUTPUT  = fit_dust

OBJS    = dust_fit_v03.o

fit_ame: $(OBJS)
	$(FC) $(OBJS) $(HEALPIX) $(FITSDIR) -fopenmp -o $(OUTPUT)

# Compilation stage
%.o : %.f90
	$(FC) $(OPTIM) $(HEALINC) $(LAPACK) $(CFITSIO) -c $<

# Cleaning command
.PHONY: clean
clean:
	rm *.o *~ fit_dust