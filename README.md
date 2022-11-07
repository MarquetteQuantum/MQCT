# MQCT_2022


## Description:

A program named MQCT_2022 is developed for calculations of rotationally and vibrationally inelastic scattering of molecules using the mixed quantum/classical theory approach. Calculations of collisions between two general asymmetric top rotors are now possible, which is a feature unavailable in other existing codes. Vibrational states of diatomic molecules can also be included in the basis to carry out calculations of ro-vibrational excitation and quenching. This version includes an important addition -- adiabatic trajectory method (AT-MQCT), in which the propagations of the equations of motion for classical and quantum parts of the system are decoupled. This approximate method is much faster, which permits to carry out calculations for larger molecular systems and at higher collision energies than it was possible before. The method is general and as such it can be applied to any molecule + molecule inelastic scattering problem.



## Compiling and running the code:
User supplied subroutine for the PES should be compiled first, to create an object file, for example `PES_H2O+He.o`. It should be copied into (or linked to) the main program directory `/MQCT_2022`. The MQCT code itself is compiled independently from the PES to create the object file `head.o`, and then is linked with the desired PES to create an executable file.

Examples of this procedure are given in the files `./comp_MQCT` and `./link_ALL`. These can be executed as commands, after changing access: 

`chmod +x ./comp_PES ./comp_MQCT ./link_ALL`

`./comp_PES` 

`./comp_MQCT` 

`./link_ALL`

**Input file** for MQCT should have the extension * *.* inp, and its name should be placed in the file INPUT_NAME.inp. This permits user to store multiple input files (e.g., for different molecules) in the program directory, but run actual calculations with one specific input file. There are two general ways of running the code. In the straightforward approach, which is also the default, the program computes elements of the state-to-state transition matrix and then propagates trajectories for collisions, all in a single run. In the optional two-step approach, which we recommend following, the program is run first with small number of processors to compute transition matrix, save it into the file and stop (without doing the calculations of collision). This is done by indicating the following optional keywords:

SAVE_MTRX = YES, PROG_RUN = NO

Then the program is run again to read the transition matrix (computed previously) and perform massively parallel trajectory calculations using large number of processors. Keywords required for this are:

READ_MTRX = YES, PROG_RUN = YES

This approach is also convenient when multiple calculations are needed with different input parameters (such as collision energy, initial state, number of trajectories, time step, etc.) but with the same basis set, which determines the matrix size. Clearly, the matrix must be computed only once, can be saved in the file and then reused later as many times as needed. The file name is MTRX_UF.dat for the binary form (unformatted) and is MTRX.dat for the formatted option of the matrix. Note that all intermediate data files created or used by in the code have extension * *.* dat.

## Understanding the output:
All output files have extension * *.* out. System setup is written into the file USER_INPUT_CHECK.out and should be checked by user for correctness. The file STATES.out (written if the option PRNT_STATES=YES is chosen) contains the list of all quantum states involved in calculations, including the channel number, the values of *j<sub>12</sub>* and *m<sub>12</sub>*, and the assigned quantum numbers. Major results are found in CROSS_SECTIONS.out. Other problem-specific output files are discussed in the main text of the article. 

## NOTE:
For description of other options, available for the expert calculations, please refer to the PDF file [MQCT USER MANUAL.pdf](./write_filename.pdf) available in the root directory. For citing the code, please refer to the article associated with the code as [MQCT_2022.pdf](./write_filename.pdf) in the root directory or see it on [the journal webpage](https:insert_link).

## Program architecture:

The actual program includes 12 individual files: head.f, propagator.f, iotest.f, mqct_eqtns.f, print_output.f, suppl_routines.f, matrix_ini.f, pot_comp.f, pes_sys_type.f, basis_wave_functions.f, asym_top.f, and vib_diatomic.f. The description of each files is as follows.

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


