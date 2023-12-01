#!/bin/bash

#SBATCH --job-name="MQCT"
#SBATCH --partition=batch
#SBATCH --time="unlimited"
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --output=%x-%j.log

module load openmpi-geib-cuda10.2-intel/4.0.5 


time mpirun ../../bin/mqct_ND3D2.exe
