# Compile the dang objects into executable

OBJS            = math_tools.o \
                  locate_mod.o \
                  spline_1d_mod.o \
                  spline_2d_mod.o \
                  dang_util_mod.o \
                  hashtbl.o \
                  dang_param_mod.o \
                  dang_bp_mod.o \
                  dang_linalg_mod.o \
                  dang_component_mod.o \
                  dang_mixmat_mod.o \
                  dang_mixmat_1d_mod.o \
                  dang_mixmat_2d_mod.o \
                  dang_mbb_comp_mod.o \
                  dang_powlaw_comp_mod.o \
                  dang_freefree_comp_mod.o \
                  dang_cmb_comp_mod.o \
                  dang_T_cmb_comp_mod.o \
                  dang_lognormal_comp_mod.o \
                  dang_hi_fit_comp_mod.o \
                  dang_template_comp_mod.o \
                  dang_monopole_comp_mod.o \
                  dang_signal_mod.o \
                  dang_data_mod.o \
                  dang_cg_mod.o \
                  dang_lnl_mod.o \
                  dang_sample_mod.o \
                  dang.o
OUTPUT          = dang

# Executable
dang: $(OBJS)
	$(MPF90) $(OBJS) $(F90FLAGS) -o $(OUTPUT) $(LINK) -fopenmp

# Dependencies
spline_1D_mod.o              : math_tools.o locate_mod.o
dang_param_mod.o             : dang_util_mod.o hashtbl.o
dang_component_mod.o         : dang_util_mod.o dang_param_mod.o dang_bp_mod.o dang_mixmat_mod.o
dang_bp_mod.o                : dang_util_mod.o dang_param_mod.o
dang_mixmat_mod.o            : dang_util_mod.o dang_bp_mod.o
dang_mixmat_1d_mod.o         : dang_util_mod.o dang_bp_mod.o dang_component_mod.o dang_mixmat_mod.o spline_1d_mod.o
dang_mixmat_2d_mod.o         : dang_util_mod.o dang_bp_mod.o dang_component_mod.o dang_mixmat_mod.o spline_2d_mod.o
dang_powlaw_comp_mod.o       : dang_component_mod.o dang_mixmat_1d_mod.o
dang_freefree_comp_mod.o     : dang_component_mod.o dang_mixmat_1d_mod.o
dang_cmb_comp_mod.o          : dang_component_mod.o dang_mixmat_1d_mod.o
dang_T_cmb_comp_mod.o        : dang_component_mod.o dang_mixmat_1d_mod.o
dang_hi_fit_comp_mod.o       : dang_component_mod.o dang_mixmat_1d_mod.o
dang_mbb_comp_mod.o          : dang_component_mod.o dang_mixmat_2d_mod.o
dang_lognormal_comp_mod.o    : dang_component_mod.o dang_mixmat_2d_mod.o
dang_monopole_comp_mod.o     : dang_component_mod.o
dang_template_comp_mod.o     : dang_component_mod.o
dang_signal_mod.o            : dang_component_mod.o dang_param_mod.o dang_freefree_comp_mod.o dang_hi_fit_comp_mod.o dang_monopole_comp_mod.o dang_template_comp_mod.o dang_lognormal_comp_mod.o dang_mbb_comp_mod.o dang_powlaw_comp_mod.o dang_T_cmb_comp_mod.o dang_cmb_comp_mod.o 
dang_data_mod.o              : dang_util_mod.o dang_param_mod.o dang_bp_mod.o dang_component_mod.o dang_signal_mod.o
dang_linalg_mod.o            : spline_1D_mod.o dang_util_mod.o dang_param_mod.o dang_data_mod.o dang_component_mod.o dang_signal_mod.o
dang_lnl_mod.o               : dang_util_mod.o dang_param_mod.o dang_bp_mod.o dang_data_mod.o dang_component_mod.o dang_signal_mod.o dang_linalg_mod.o dang_cg_mod.o  
dang_sample_mod.o            : dang_util_mod.o dang_param_mod.o dang_bp_mod.o dang_data_mod.o dang_component_mod.o dang_signal_mod.o dang_linalg_mod.o dang_cg_mod.o  dang_lnl_mod.o
dang_cg_mod.o                : dang_util_mod.o dang_param_mod.o dang_bp_mod.o dang_data_mod.o dang_component_mod.o dang_signal_mod.o
dang.o                       : dang_util_mod.o dang_param_mod.o dang_bp_mod.o dang_linalg_mod.o dang_sample_mod.o dang_data_mod.o dang_signal_mod.o

# Compilation stage
%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

# Cleaning command
.PHONY: clean
clean:
	rm *.o *.mod *~ dang
