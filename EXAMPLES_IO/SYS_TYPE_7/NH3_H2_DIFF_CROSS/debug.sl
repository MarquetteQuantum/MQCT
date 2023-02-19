#!/bin/bash -l
#SBATCH -q debug
#SBATCH -N 10
#SBATCH -t 00:30:00
#SBATCH -C knl
#SBATCH -J MQCT_H2O_He

time srun -n 640 ./mqct_exe
