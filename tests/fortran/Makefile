FC = gfortran
FFLAGS ?= -O0 -g -Wall -Wextra -fcheck=all
LIBS ?= -llapack -lblas

ROOT := ../..
SRC := $(ROOT)/src
BUILD := $(ROOT)/build/fortran_tests

TEST_EXES := $(BUILD)/test_locate_mod.exe $(BUILD)/test_math_tools.exe $(BUILD)/test_math_tools_linear.exe $(BUILD)/test_math_tools_full.exe $(BUILD)/test_spline_1d_mod.exe

.PHONY: all run clean

all: $(TEST_EXES)

run: all
	$(BUILD)/test_locate_mod.exe
	$(BUILD)/test_math_tools.exe
	$(BUILD)/test_math_tools_linear.exe
	$(BUILD)/test_math_tools_full.exe
	$(BUILD)/test_spline_1d_mod.exe

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/healpix_types.o: healpix_types_stub.f90 | $(BUILD)
	$(FC) $(FFLAGS) -c $< -J$(BUILD) -o $@

$(BUILD)/locate_mod.o: $(SRC)/locate_mod.f90 $(BUILD)/healpix_types.o
	$(FC) $(FFLAGS) -I$(BUILD) -c $< -J$(BUILD) -o $@

$(BUILD)/math_tools.o: $(SRC)/math_tools.f90 $(BUILD)/healpix_types.o
	$(FC) $(FFLAGS) -I$(BUILD) -c $< -J$(BUILD) -o $@

$(BUILD)/spline_1d_mod.o: $(SRC)/spline_1d_mod.f90 $(BUILD)/locate_mod.o $(BUILD)/math_tools.o
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

clean:
	rm -rf $(BUILD)
