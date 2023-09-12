# New Makefile attempt for super-dope optimized sickness

# GNU Fortran
F90       = gfortran
MPF90     = mpif90
F90FLAGS  = -g -C 

LOCAL=/home/daniher/.local/

LAPACK_LINK     = -L/usr/lib -llapack
BLAS_LINK       = -L/usr/lib -lblas

HEALPIX         = /home/daniher/HEALPix/Healpix_3.82

# General variables:

#CFITSIO
CFITSIO_INCLUDE =
CFITSIO_LINK    = -L$(LOCAL)/lib -lcfitsio

# HEALPIX
HEALPIX_INCLUDE = -I$(HEALPIX)/include
HEALPIX_LINK    = -L$(HEALPIX)/lib -lhealpix

#Combine them
F90COMP         = $(F90FLAGS) $(LAPACK_INCLUDE) $(CFITSIO_INCLUDE) $(HEALPIX_INCLUDE)
LINK            = $(HEALPIX_LINK) $(CFITSIO_LINK) $(LAPACK_LINK) $(BLAS_LINK)

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
