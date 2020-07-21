# -*- Makefile -*-

FC      = ifort
OPTIM   = -g -c

MPICHINC = -I/astro/local/opt/Intel/compilers_and_libraries_2018.3.222/linux/mpi/intel64/include
#MPICH    = -L/astro/local/opt/Intel/compilers_and_libraries_2018.3.222/linux/mpi/intel64/lib -lmpi
FITSDIR  = -L/mn/stornext/u3/hke/local/lib -lcfitsio
LAPACK   = -L/mn/stornext/u3/hke/local/lib -llapack -lblas
HEALPIX  = -L/mn/stornext/u3/hke/local/lib -lhealpix
HEALINC  = -I/mn/stornext/u3/hke/local/include
OUTPUT   = dang

OBJS    = init_mod.o utility_mod.o hashtbl.o param_mod.o linalg_mod.o data_mod.o dang.o

dang: $(OBJS)
	$(FC) $(OBJS) $(HEALPIX) $(FITSDIR) $(LAPACK) -qopenmp -o $(OUTPUT)

# Dependencies
linalg_mod.o           : init_mod.o utility_mod.o
data_mod.o             : init_mod.o utility_mod.o
#foreground_mod.o       : init_mod.o utility_mod.o
param_mod.o            : init_mod.o hashtbl.o utility_mod.o
dang.o : init_mod.o utility_mod.o param_mod.o linalg_mod.o data_mod.o #foreground_mod.o

# Compilation stage
%.o : %.f90
	$(FC) -fpp $(OPTIM) $(HEALINC) $(LAPACK) $(CFITSIO) -qopenmp -parallel -c $<

# Cleaning command
.PHONY: clean
clean:
	rm *.o *.mod *~ dang
