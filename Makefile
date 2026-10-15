# New Makefile attempt for super-dope optimized sickness

# Load variables from the config file
include config/config.gnu_laptop

export F90COMP := $(F90FLAGS) $(LAPACK_INCLUDE) $(CFITSIO_INCLUDE) $(HEALPIX_INCLUDE)
export LINK    := $(HEALPIX_LINK) $(CFITSIO_LINK) $(LAPACK_LINK) $(BLAS_LINK)

# Executable
all : dang

dang :
	cd src; $(MAKE) 

# Compilation stage
%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

# Cleaning command
.PHONY: clean stub-cg-check stub-cg-polflag-regression
clean :
	@cd src; $(MAKE) clean

stub-cg-check:
	./scripts/stub_build_polarization.sh

stub-cg-polflag-regression:
	./scripts/stub_regression_cg_polflags.sh
