#!/bin/bash

#SBATCH --job-name="water_water_AT_MQCT"
#SBATCH --partition=batch
#SBATCH --time="unlimited"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=%x-%j.log

module load openmpi-geib-cuda10.2-intel/4.0.5 ;

time mpirun /mmfs1/home/9319bostand/MQCT-2023_research/November_test/bin/mqct_H2OH2.exe
