#!/bin/bash
#SBATCH --job-name="PostProcess_example"
#SBATCH -p parallel-12
#SBATCH -n 12
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu=1G

lmp=~/software/lammps/lam*22/src/ # getting the correct run file location
mpirun $lmp/lmp_mpi < simulation.in # computing with n cpu cores.
wait

exit 0

