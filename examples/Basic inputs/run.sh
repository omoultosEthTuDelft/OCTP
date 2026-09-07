#!/bin/bash
#SBATCH --job-name="PostProcess_example"
#SBATCH -p parallel-16
#SBATCH -n 16
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu=1G


srun lmp_avx -i simulation.in # computing with n cpu cores.
wait

exit 
