# MQCT_2024


## Description:
A program named MQCT_2024 is developed for calculations of rotationally and vibrationally inelastic scattering of molecules using the mixed quantum/classical theory approach. Calculations of collisions between two general asymmetric top rotors are possible, which is a feature unavailable in other existing codes. Vibrational states of diatomic molecules can also be included in the basis to carry out calculations of ro-vibrational excitation and quenching. This version includes an important addition -- adiabatic trajectory method (AT-MQCT), in which the propagations of the equations of motion for classical and quantum parts of the system are decoupled. This approximate method is much faster, which permits to carry out calculations for larger molecular systems and at higher collision energies than it was possible before. The method is general and as such it can be applied to any molecule + molecule inelastic scattering problem.


## Compiling and running the code:
Compilation of the PES function and of MQCT code is done by using a `Makefile`. Our version of `Makefile` is provided as an example in the main directory. Users are expected to change variables in the Makefile as needed for their computer architecture, directory structure and PES function. There are two levels of parallelization in the code. At the first level, propagation of each trajectory can be done by multiple processors used as a group (e.g., all processors of a node) to compute right-hand sides of the classical and quantum equations of motion. This requires some minimal message passing. At the second level, propagation of different trajectories can be assigned to different groups, which requires virtually no message passing. 


## Preparing the input:
The input file for MQCT should have extension * *.* inp, and its name should be placed in the file `INPUT_NAME.inp`. This permits users to store multiple input files (e.g., for different molecules) in the program directory, but run actual calculations with one specific input file. Initial conditions for classical and quantum degrees of freedom in the system are defined separately in two blocks of the input file, called `$SYSTEM` and `$BASIS`. Potential energy surface describes interaction between quantum and classical parts of the system and is defined in the third block of the input file, called `$POTENTIAL`. In MQCT program there are ten system types from the simplest rigid-diatom + atom, to the most general case of two asymmetric-top rotor molecules. For complicated molecular systems with many states and large state-to-state transition matrix, it is recommended that the matrix is pre-computed in a separate run of MQCT (any number of processors can be used), and only after that the trajectories are propagated following parallelization strategy described above. Namely, the following optional keywords:

`SAVE_MTRX = YES, PROG_RUN = NO`

permit to compute transition matrix, save it into the file and stop (without doing the calculations of collision). Then the program is run again to read the transition matrix (computed previously) and perform massively parallel trajectory calculations using large number of processors. Keywords required for this are:

`READ_MTRX = YES, PROG_RUN = YES`

This approach is also convenient when multiple calculations are needed with different input parameters (such as collision energy, initial state, number of trajectories, time step, etc.) but with the same basis set, which determines the matrix size. The matrix must be computed only once, can be saved in the file and then reused later as many times as needed. The file name is `MTRX_UF.dat` for the binary form (unformatted) and is `MTRX.dat` for the formatted option of the matrix. Note that all intermediate data files created or used by in the code have extension * *.* dat.


## Understanding the output:
All output files have extension * *.* out. System setup is written into the file `USER_INPUT_CHECK.out` and should be checked by user for correctness. The file `STATES.out` (written if the option `PRNT_STATES=YES` is chosen) contains the list of all quantum states involved in calculations, including the channel number, the values of *j* and *m*, and the assigned quantum numbers. Major results are found in `CROSS_SECTIONS.out`. Other problem-specific output files are discussed in user manual.


## NOTE:
For description of other options, available for the expert calculations, please refer to the PDF file [MQCT USER MANUAL.pdf](./MQCT_USER_MANUAL.pdf) available in the root directory. For citing the code, please refer to the article associated with the code as [MQCT_2020](https://doi.org/10.1016/j.cpc.2020.107155) and [AT_MQCT](https://doi.org/10.1021/acs.jpclett.2c03328). 


## Program architecture:
The actual program includes 12 individual files listed and described below: 

**head.f**: Uses MPI to parallelize calculations and communicates between other files.

**propagator.f**: Propagates the collision trajectories for classical degrees of freedom.

**iotest.f**: Contains all global and local variables and the parser of input.

**mqct_eqtns.f**: Propagates the equations of motion for quantum degrees of freedom.

**print_output.f**: Prints the output files with the results of calculations.

**suppl_routines.f**: Contains several subroutines from Numerical Recipes.

**matrix_ini.f**: Computes elements of the state-to-state transition matrix.

**pot_comp.f**: Incorporates the potential energy surface for calculations.

**pes_sys_type.f**: Executes PES subroutines for various systems.

**basis_wave_functions.f**: Computes wave functions for diatoms and symmetric tops.

**asym_top.f**: Computes wave function for asymmetric top rotor systems.

**vib_diatomic.f**: Computes vibrational wave functions for diatomic systems.


