#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tianze.peng@mail.mcgill.ca
#SBATCH --account=def-straub
#SBATCH --nodes=1
#SBATCH --job-name="Read Diags"
#SBATCH --output=results
#SBATCH --ntasks-per-node=01
#SBATCH --mem-per-cpu=0     # memory; default unit is megabytes
#SBATCH --time=0-01:45           # time (DD-HH:MM)
module load matlab/2018b
time matlab -nodisplay -nosplash -nodesktop -r "run('./postproc.m');exit;"> process.log

