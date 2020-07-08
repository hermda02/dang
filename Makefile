# -*- Makefile -*-

FC      = gfortran
OPTIM   = -g -C -O3

FITSDIR = -L/usr/lib -lcfitsio
LAPACK  = -L/usr/lib -llapack -lblas
HEALPIX = -L/usr/local/src/Healpix_3.50/lib -lhealpix
HEALINC = -I/usr/local/src/Healpix_3.50/include
OUTPUT  = dang

OBJS    = init_mod.o hashtbl.o param_mod.o linalg_mod.o foreground_mod.o data_mod.o dang.o

dang: $(OBJS)
	$(FC) $(OBJS) $(HEALPIX) $(FITSDIR) $(LAPACK) -fopenmp -o $(OUTPUT)

# Dependencies
linalg_mod.o           : init_mod.o
data_mod.o             : init_mod.o
foreground_mod.o       : init_mod.o
param_mod.o            : init_mod.o hashtbl.o
dang.o : init_mod.o param_mod.o linalg_mod.o data_mod.o foreground_mod.o

# Compilation stage
%.o : %.f90
	$(FC) $(OPTIM) $(HEALINC) $(LAPACK) $(CFITSIO) -c $<

# Cleaning command
.PHONY: clean
clean:
	rm *.o *.mod *~ dang
