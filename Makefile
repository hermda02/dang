# New Makefile attempt for super-dope optimized sickness

# Load variables from the config file
include config/config.gnu_laptop

export F90COMP := $(F90FLAGS) $(LAPACK_INCLUDE) $(CFITSIO_INCLUDE) $(HEALPIX_INCLUDE)
export LINK    := $(HEALPIX_LINK) $(CFITSIO_LINK) $(LAPACK_LINK) $(BLAS_LINK)

# Executable
all : dang

test :
	python3 -m unittest discover -s tests -v

fortran-test :
	$(MAKE) -C tests/fortran run

check : test fortran-test

dang :
	cd src; $(MAKE) 

# Compilation stage
%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

# Cleaning command
.PHONY: clean test fortran-test check
clean :
	@cd src; $(MAKE) clean
