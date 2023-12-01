#!/bin/bash -l
#SBATCH -q debug  
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH -C cpu
#SBATCH -J H2OHe

module load PrgEnv-intel/8.3.3;

time srun -n 16 /pscratch/sd/b/bostan_d/example_io_tests/bin/mqct_CH3COOH_He.exe