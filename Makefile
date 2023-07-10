# this makefile is written by Bikramaditya Mandal July 09, 2023
# this makefile is intended to compile the MQCT code with the PES file
# user should indicate appropriate direcotry and file name for PES and the executable
# make sure that the necessary modules are loaded
#---------------------------------------------------------------------------------------------

# Compiler and flags
FC = mpifort
FFLAGS = -g -traceback
LDFLAGS = -liomp5 -lmkl_intel_thread -lmkl_rt -lmkl_core -lmkl_intel_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_scalapack_lp64

CURDIR = .
pes_dir = ./PES_DIRECTORIES/PES_N2O
src_dir = ./MQCT_SOURCE_CODES

pes_files = PES_N2O.o
src_files = head.o iotest.o propagator.o print_output.o asym_top.o vib_diatomic.o pot_comp.o matrix_ini.o mqct_eqtns.o pes_sys_type.o suppl_routines.o

all: pes mqct mqct_N2O

pes: $(pes_dir)/$(pes_files)

mqct: $(src_dir)/$(src_files)

mqct_N2O:
	$(FC) $(pes_dir)/$(pes_files) $(src_dir)/$(src_files) -o $@.exe $(FFLAGS) $(LDFLAGS)
	mv *.mod $(src_dir)/
	mv *.exe $(src_dir)/
	ln -s $(src_dir)/$@.exe ./bin/

%.o : %.f
	$(FC) -c $(FFLAGS) $(LDFLAGS) $< -o $@

clean:
	rm -f $(pes_dir)/*.o
	rm -f $(src_dir)/*.o
	rm -f $(src_dir)/*.mod
	rm -f $(src_dir)/*.exe
	rm -f $(CURDIR)/bin/*
