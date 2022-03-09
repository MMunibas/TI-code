#!/bin/bash

#SBATCH --job-name=PhBrEq
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --partition=long

module load gcc/gcc-9.2.0-openmpi-4.0.2-ib

cd /data2/devereux/MDCM-symmetry/MDCM-dyna/PhBr/deltaG/equilibrate

#srun /opt/cluster/programs/charmm/developer/dev-release-dcm/build/cmake/charmm -i equil.inp -o equil.out
srun /opt/cluster/programs/charmm/developer/dev-release-dcm/build/cmake/charmm -i equil-gas.inp -o equil-gas.out
