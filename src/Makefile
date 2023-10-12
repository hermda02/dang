# Compile the dang objects into executable

OBJS            = math_tools.o locate_mod.o spline_1D_mod.o dang_util_mod.o hashtbl.o dang_param_mod.o dang_bp_mod.o  dang_linalg_mod.o dang_component_mod.o dang_data_mod.o dang_cg_mod.o dang_lnl_mod.o dang_sample_mod.o dang.o
OUTPUT          = dang

# Executable
dang: $(OBJS)
	$(MPF90) $(OBJS) $(F90FLAGS) -o $(OUTPUT) $(LINK) -fopenmp

# Dependencies
spline_1D_mod.o        : math_tools.o locate_mod.o
dang_param_mod.o       : dang_util_mod.o hashtbl.o
dang_component_mod.o   : dang_util_mod.o dang_param_mod.o dang_bp_mod.o 
dang_bp_mod.o          : dang_util_mod.o dang_param_mod.o
dang_data_mod.o        : dang_util_mod.o dang_param_mod.o dang_bp_mod.o dang_component_mod.o
dang_linalg_mod.o      : spline_1D_mod.o dang_util_mod.o dang_param_mod.o dang_data_mod.o dang_component_mod.o
dang_lnl_mod.o         : dang_util_mod.o dang_param_mod.o dang_bp_mod.o dang_data_mod.o dang_component_mod.o dang_linalg_mod.o dang_cg_mod.o  
dang_sample_mod.o      : dang_util_mod.o dang_param_mod.o dang_bp_mod.o dang_data_mod.o dang_component_mod.o dang_linalg_mod.o dang_cg_mod.o  dang_lnl_mod.o
dang_cg_mod.o          : dang_util_mod.o dang_param_mod.o dang_bp_mod.o dang_data_mod.o dang_component_mod.o 
dang.o : dang_util_mod.o dang_param_mod.o dang_bp_mod.o dang_linalg_mod.o dang_sample_mod.o dang_data_mod.o

# Compilation stage
%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

# Cleaning command
.PHONY: clean
clean:
	rm *.o *.mod *~ dang
