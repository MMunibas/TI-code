#!/bin/bash

######################################################
# qsub arguments
######################################################

#SBATCH --job-name=test_stable
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --partition=epyc
#SBATCH --mem-per-cpu=1000
######################################################
# lib 
######################################################


cd PPP

echo "Job running on "`hostname`

module load gcc/gcc-9.2.0-openmpi-4.1.2-ib
CHARMM=/home/devereux/bin/dev-release-dcm/build/cmake/charmm
#module load gcc/gcc-9.2.0-openmpi-4.0.2-ib
#CHARMM=/opt/cluster/programs/charmm/developer/dev-release-dcm/build/cmake/charmm

srun $CHARMM < III.inp > III.log

# Read TI results from III.log
if [ -e III.log ]; then
  dG=$( grep "EPRTOT=" III.log | awk '{print $8}' | tail -1 )
  echo "$dG" > dG.dat

  # Read variance from III.log
  var=$( grep "DIFFLC=" III.log | awk '{print $10}' | tail -1 )
  echo "$var" > var.dat

fi

