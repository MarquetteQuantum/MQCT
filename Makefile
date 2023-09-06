#---------------------------------------------------------------------------------------------
#--------------                       Makefile user Guide                      ---------------
#---------------------------------------------------------------------------------------------
# This makefile is written by Bikramaditya Mandal July 09, 2023. 
# This makefile is intended to compile and link MQCT code with PES file(s).
# User should indicate appropriate direcotry for the PES and the file name for an executable.
# Make sure that all necessary modules are loaded.
# After appropriate changes in the makefile, user should clean once using $make clean
# then user should follow standard process of compiling by using $make.
# The binaries/executables should be present in the directory ./bin
#
# $make clean: removes all module files, links, old binaries/executables of MQCT & PES codes.
# $make: compiles the PES, MQCT code files, and stores executable file in the directory ./bin 


#---------------------------------------------------------------------------------------------
#--------------   This is where the user is expected to change accordingly     ---------------
#---------------------------------------------------------------------------------------------
pes_dir = ./PES_DIRECTORIES/PES_NH3_H2O
pes_files =  PES_NH3H2O.o PES_NH3_H2O.o isotopologue_coordinates_conversion.o lapack.o
exe_name = mqct_NH3H2O
src_dir = ./MQCT_SOURCE_CODES

# Compiler and flags:
FC = mpifort
FFLAGS = -g -traceback
LDFLAGS = -liomp5 -lmkl_intel_thread -lmkl_rt -lmkl_core -lmkl_intel_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_scalapack_lp64

# Dont forget to load manually all necessary compiler modules prior to using "make" comand.


#---------------------------------------------------------------------------------------------
#--------------                   DO NOT CHANGE ANYTHING BELOW                 ---------------
#---------------------------------------------------------------------------------------------
CURDIR = .
# all rules
all: pes mqct $(exe_name)

#PES Rules
$(info $(pes_files))
PES_OBJS := $(pes_files:%=$(pes_dir)/%)
pes: $(PES_OBJS)

#MQCT Source codes rules
src_files = head.o iotest.o propagator.o print_output.o asym_top.o vib_diatomic.o pot_comp.o matrix_ini.o mqct_eqtns.o pes_sys_type.o suppl_routines.o
SRC_OBJS := $(src_files:%=$(src_dir)/%)
mqct: $(SRC_OBJS)

$(exe_name):
	$(FC) $(PES_OBJS) $(SRC_OBJS) -o $@.exe $(FFLAGS) $(LDFLAGS)
	-mv *.mod $(src_dir)/
	-mv $@.exe ./bin/

%.o : %.f
	$(FC) -c $(FFLAGS) $(LDFLAGS) $< -o $@
	
%.o : %.f90
	$(FC) -c $(FFLAGS) $(LDFLAGS) $< -o $@

clean:
	rm -f $(pes_dir)/*.o
	rm -f $(src_dir)/*.o
	rm -f $(src_dir)/*.mod
	rm -f $(CURDIR)/*.o
	rm -f $(CURDIR)/*.mod
	rm -f $(CURDIR)/$(exe_name).exe
	rm -f $(CURDIR)/bin/$(exe_name).exe
