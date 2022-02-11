# New Makefile attempt for super-dope optimized sickness

F90       = gfortran
MPF90     = mpif90
F90FLAGS  = -g -C

LOCAL=/usr/local/src

#LAPACK
MKLPATH         = $(MKLROOT)
LAPACK_INCLUDE  =
LAPACK_LINK     = -shared-intel -Wl,-rpath,$(MKLPATH) -L$(MKLPATH)  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread

#CFITSIO
CFITSIO_INCLUDE =
CFITSIO_LINK    = -L$(LOCAL)/lib -lcfitsio

# HEALPIX
HEALPIX         = /usr/local/src/Healpix_3.70
HEALPIX_INCLUDE = -I$(HEALPIX)/include
HEALPIX_LINK    = -L$(HEALPIX)/lib -lhealpix

#Combine them
F90COMP         = $(F90FLAGS) $(LAPACK_INCLUDE) $(CFITSIO_INCLUDE) $(HEALPIX_INCLUDE)
LINK            = $(HEALPIX_LINK) $(CFITSIO_LINK) $(LAPACK_LINK)
OBJS            = math_tools.o locate_mod.o spline_1D_mod.o dang_util_mod.o hashtbl.o dang_param_mod.o dang_bp_mod.o  dang_linalg_mod.o dang_component_mod.o dang_data_mod.o dang_sample_mod.o dang_swap_mod.o dang.o
OUTPUT          = dang

# Executable
dang: $(OBJS)
	$(MPF90) $(OBJS)  -o $(OUTPUT) $(LINK)

# Dependencies
spline_1D_mod.o        : math_tools.o locate_mod.o
dang_param_mod.o       : dang_util_mod.o hashtbl.o
dang_component_mod.o   : dang_util_mod.o dang_param_mod.o dang_bp_mod.o 
dang_bp_mod.o          : dang_util_mod.o dang_param_mod.o
dang_data_mod.o        : dang_util_mod.o dang_param_mod.o dang_bp_mod.o dang_component_mod.o
dang_linalg_mod.o      : spline_1D_mod.o dang_util_mod.o dang_param_mod.o dang_data_mod.o dang_component_mod.o
dang_sample_mod.o      : dang_util_mod.o dang_param_mod.o dang_bp_mod.o dang_data_mod.o dang_component_mod.o dang_linalg_mod.o  
dang_swap_mod.o        : dang_util_mod.o dang_param_mod.o dang_data_mod.o dang_component_mod.o dang_linalg_mod.o  
dang.o : dang_util_mod.o dang_param_mod.o dang_bp_mod.o dang_linalg_mod.o dang_sample_mod.o dang_data_mod.o dang_swap_mod.o

# Compilation stage
%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

# Cleaning command
.PHONY: clean
clean:
	rm *.o *.mod *~ dang
