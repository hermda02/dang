FC = gfortran
FFLAGS ?= -O0 -g -Wall -Wextra -fcheck=all
LIBS ?= -llapack -lblas

ROOT := ../..
SRC := $(ROOT)/src
BUILD := $(ROOT)/build/fortran_tests

TEST_EXES := $(BUILD)/test_locate_mod.exe $(BUILD)/test_math_tools.exe $(BUILD)/test_math_tools_linear.exe $(BUILD)/test_math_tools_full.exe $(BUILD)/test_spline_1d_mod.exe $(BUILD)/test_spline_2d_mod.exe $(BUILD)/test_hashtbl.exe $(BUILD)/test_dang_util_mod.exe $(BUILD)/test_dang_linalg_mod.exe $(BUILD)/test_dang_mixmat_mod.exe

.PHONY: all run clean

all: $(TEST_EXES)

run: all
	$(BUILD)/test_locate_mod.exe
	$(BUILD)/test_math_tools.exe
	$(BUILD)/test_math_tools_linear.exe
	$(BUILD)/test_math_tools_full.exe
	$(BUILD)/test_spline_1d_mod.exe
	$(BUILD)/test_spline_2d_mod.exe
	$(BUILD)/test_hashtbl.exe
	$(BUILD)/test_dang_util_mod.exe
	$(BUILD)/test_dang_linalg_mod.exe
	$(BUILD)/test_dang_mixmat_mod.exe

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/healpix_types.o: healpix_types_stub.f90 | $(BUILD)
	$(FC) $(FFLAGS) -c $< -J$(BUILD) -o $@

$(BUILD)/external_lib_stubs.o: external_lib_stubs.f90 $(BUILD)/healpix_types.o
	$(FC) $(FFLAGS) -I$(BUILD) -c $< -J$(BUILD) -o $@

$(BUILD)/dang_model_stubs.o: dang_model_stubs.f90 $(BUILD)/healpix_types.o
	$(FC) $(FFLAGS) -I$(BUILD) -c $< -J$(BUILD) -o $@

$(BUILD)/locate_mod.o: $(SRC)/locate_mod.f90 $(BUILD)/healpix_types.o
	$(FC) $(FFLAGS) -I$(BUILD) -c $< -J$(BUILD) -o $@

$(BUILD)/math_tools.o: $(SRC)/math_tools.f90 $(BUILD)/healpix_types.o
	$(FC) $(FFLAGS) -I$(BUILD) -c $< -J$(BUILD) -o $@

$(BUILD)/spline_1d_mod.o: $(SRC)/spline_1d_mod.f90 $(BUILD)/locate_mod.o $(BUILD)/math_tools.o
	$(FC) $(FFLAGS) -I$(BUILD) -c $< -J$(BUILD) -o $@

$(BUILD)/spline_2d_mod.o: $(SRC)/spline_2d_mod.f90 $(BUILD)/spline_1d_mod.o $(BUILD)/math_tools.o
	$(FC) $(FFLAGS) -I$(BUILD) -c $< -J$(BUILD) -o $@

$(BUILD)/dang_util_mod.o: $(SRC)/dang_util_mod.f90 $(BUILD)/external_lib_stubs.o $(BUILD)/healpix_types.o
	$(FC) $(FFLAGS) -I$(BUILD) -c $< -J$(BUILD) -o $@

$(BUILD)/dang_linalg_mod.o: $(SRC)/dang_linalg_mod.f90 $(BUILD)/dang_util_mod.o $(BUILD)/dang_model_stubs.o
	$(FC) $(FFLAGS) -I$(BUILD) -c $< -J$(BUILD) -o $@

$(BUILD)/hashtbl.o: $(SRC)/hashtbl.f90 | $(BUILD)
	$(FC) $(FFLAGS) -c $< -J$(BUILD) -o $@

$(BUILD)/dang_mixmat_mod.o: $(SRC)/dang_mixmat_mod.f90 $(BUILD)/dang_util_mod.o $(BUILD)/dang_model_stubs.o
	$(FC) $(FFLAGS) -I$(BUILD) -c $< -J$(BUILD) -o $@

$(BUILD)/dang_mixmat_1d_mod.o: $(SRC)/dang_mixmat_1d_mod.f90 $(BUILD)/dang_mixmat_mod.o $(BUILD)/dang_model_stubs.o $(BUILD)/spline_1d_mod.o
	$(FC) $(FFLAGS) -I$(BUILD) -c $< -J$(BUILD) -o $@

$(BUILD)/dang_mixmat_2d_mod.o: $(SRC)/dang_mixmat_2d_mod.f90 $(BUILD)/dang_mixmat_mod.o $(BUILD)/dang_model_stubs.o $(BUILD)/spline_2d_mod.o
	$(FC) $(FFLAGS) -I$(BUILD) -c $< -J$(BUILD) -o $@

$(BUILD)/test_locate_mod.o: test_locate_mod.f90 $(BUILD)/locate_mod.o
	$(FC) $(FFLAGS) -I$(BUILD) -c $< -J$(BUILD) -o $@

$(BUILD)/test_math_tools.o: test_math_tools.f90 $(BUILD)/math_tools.o
	$(FC) $(FFLAGS) -I$(BUILD) -c $< -J$(BUILD) -o $@

$(BUILD)/test_math_tools_linear.o: test_math_tools_linear.f90 $(BUILD)/math_tools.o
	$(FC) $(FFLAGS) -I$(BUILD) -c $< -J$(BUILD) -o $@

$(BUILD)/test_math_tools_full.o: test_math_tools_full.f90 $(BUILD)/math_tools.o
	$(FC) $(FFLAGS) -I$(BUILD) -c $< -J$(BUILD) -o $@

$(BUILD)/test_spline_1d_mod.o: test_spline_1d_mod.f90 $(BUILD)/spline_1d_mod.o
	$(FC) $(FFLAGS) -I$(BUILD) -c $< -J$(BUILD) -o $@

$(BUILD)/test_spline_2d_mod.o: test_spline_2d_mod.f90 $(BUILD)/spline_2d_mod.o
	$(FC) $(FFLAGS) -I$(BUILD) -c $< -J$(BUILD) -o $@

$(BUILD)/test_hashtbl.o: test_hashtbl.f90 $(BUILD)/hashtbl.o
	$(FC) $(FFLAGS) -I$(BUILD) -c $< -J$(BUILD) -o $@

$(BUILD)/test_dang_util_mod.o: test_dang_util_mod.f90 $(BUILD)/dang_util_mod.o
	$(FC) $(FFLAGS) -I$(BUILD) -c $< -J$(BUILD) -o $@

$(BUILD)/test_dang_linalg_mod.o: test_dang_linalg_mod.f90 $(BUILD)/dang_linalg_mod.o
	$(FC) $(FFLAGS) -I$(BUILD) -c $< -J$(BUILD) -o $@

$(BUILD)/mixmat_test_types.o: mixmat_test_types.f90 $(BUILD)/dang_model_stubs.o
	$(FC) $(FFLAGS) -I$(BUILD) -c $< -J$(BUILD) -o $@

$(BUILD)/test_dang_mixmat_mod.o: test_dang_mixmat_mod.f90 $(BUILD)/mixmat_test_types.o $(BUILD)/dang_mixmat_1d_mod.o $(BUILD)/dang_mixmat_2d_mod.o
	$(FC) $(FFLAGS) -I$(BUILD) -c $< -J$(BUILD) -o $@

$(BUILD)/test_locate_mod.exe: $(BUILD)/healpix_types.o $(BUILD)/locate_mod.o $(BUILD)/test_locate_mod.o
	$(FC) $(FFLAGS) $^ -o $@

$(BUILD)/test_math_tools.exe: $(BUILD)/healpix_types.o $(BUILD)/math_tools.o $(BUILD)/test_math_tools.o
	$(FC) $(FFLAGS) $^ -o $@ $(LIBS)

$(BUILD)/test_math_tools_linear.exe: $(BUILD)/healpix_types.o $(BUILD)/math_tools.o $(BUILD)/test_math_tools_linear.o
	$(FC) $(FFLAGS) $^ -o $@ $(LIBS)

$(BUILD)/test_math_tools_full.exe: $(BUILD)/healpix_types.o $(BUILD)/math_tools.o $(BUILD)/test_math_tools_full.o
	$(FC) $(FFLAGS) $^ -o $@ $(LIBS)

$(BUILD)/test_spline_1d_mod.exe: $(BUILD)/healpix_types.o $(BUILD)/locate_mod.o $(BUILD)/math_tools.o $(BUILD)/spline_1d_mod.o $(BUILD)/test_spline_1d_mod.o
	$(FC) $(FFLAGS) $^ -o $@ $(LIBS)

$(BUILD)/test_spline_2d_mod.exe: $(BUILD)/healpix_types.o $(BUILD)/locate_mod.o $(BUILD)/math_tools.o $(BUILD)/spline_1d_mod.o $(BUILD)/spline_2d_mod.o $(BUILD)/test_spline_2d_mod.o
	$(FC) $(FFLAGS) $^ -o $@ $(LIBS)

$(BUILD)/test_hashtbl.exe: $(BUILD)/hashtbl.o $(BUILD)/test_hashtbl.o
	$(FC) $(FFLAGS) $^ -o $@

$(BUILD)/test_dang_util_mod.exe: $(BUILD)/healpix_types.o $(BUILD)/external_lib_stubs.o $(BUILD)/dang_util_mod.o $(BUILD)/test_dang_util_mod.o
	$(FC) $(FFLAGS) $^ -o $@

$(BUILD)/test_dang_linalg_mod.exe: $(BUILD)/healpix_types.o $(BUILD)/external_lib_stubs.o $(BUILD)/dang_model_stubs.o $(BUILD)/dang_util_mod.o $(BUILD)/dang_linalg_mod.o $(BUILD)/test_dang_linalg_mod.o
	$(FC) $(FFLAGS) $^ -o $@ $(LIBS)

$(BUILD)/test_dang_mixmat_mod.exe: $(BUILD)/healpix_types.o $(BUILD)/external_lib_stubs.o $(BUILD)/dang_model_stubs.o $(BUILD)/dang_util_mod.o $(BUILD)/locate_mod.o $(BUILD)/math_tools.o $(BUILD)/spline_1d_mod.o $(BUILD)/spline_2d_mod.o $(BUILD)/dang_mixmat_mod.o $(BUILD)/dang_mixmat_1d_mod.o $(BUILD)/dang_mixmat_2d_mod.o $(BUILD)/mixmat_test_types.o $(BUILD)/test_dang_mixmat_mod.o
	$(FC) $(FFLAGS) $^ -o $@ $(LIBS)

clean:
	rm -rf $(BUILD)
