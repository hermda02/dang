# New Makefile attempt for super-dope optimized sickness

# Load variables from the config file
include config/config.gnu_desktop

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
.PHONY: clean
clean :
	@cd src; $(MAKE) clean
