# New Makefile attempt for super-dope optimized sickness

F90       = gfortran
MPF90     = mpifort
F90FLAGS  = -g -C -traceback -03

#LAPACK
MKLPATH         = $(MKLROOT)
LAPACK_INCLUDE  =
LAPACK_LINK     = -shared-intel -Wl,-rpath,$(MKLPATH) -L$(MKLPATH) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread

#CFITSIO
CFITSIO_INCLUDE =
CFITSIO_LINK    = -L/usr/lib -lcfitsio

# HEALPIX
HEALPIX         = /usr/local/src/Healpix_3.70
HEALPIX_INCLUDE = -I$(HEALPIX)/include
HEALPIX_LINK    = -L$(HEALPIX)/lib -lhealpix

#Combine them
F90COMP         = $(F90FLAGS) $(LAPACK_INCLUDE) $(CFITSIO_INCLUDE) $(HEALPIX_INCLUDE)
LINK            = $(HEALPIX_LINK) $(CFITSIO_LINK) $(LAPACK_LINK)

# Executable
all : dang

dang :
	cd src; $(MAKE) 

# Compilation stage
%.o : %.f90
	$(MPF90) -fpp $(F90COMP) -fopenmp -parallel -c $<

# Cleaning command
.PHONY: clean
clean :
	@cd src; $(MAKE) clean
