#!/bin/bash

#SBATCH --job-name="C6H6_He_AT_MQCT_2nd_run"
#SBATCH --partition=batch
#SBATCH --time="unlimited"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=%x-%j.log

module load openmpi-geib-cuda10.2-intel/4.0.5 ;

time mpirun /mmfs1/home/4168mandalb/MQCT/mqct_C6H6_He
