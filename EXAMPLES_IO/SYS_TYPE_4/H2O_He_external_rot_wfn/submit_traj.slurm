#!/bin/bash

#SBATCH --job-name="MQCT"
#SBATCH --partition=batch
#SBATCH --time="unlimited"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=%x-%j.log

module load openmpi-geib-cuda10.2-intel/4.0.5 ;

time mpirun /mmfs1/home/4168mandalb/v20_21/mqct_v20_21_H2O_He


