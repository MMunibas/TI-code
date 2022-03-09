#!/bin/bash
#

#SBATCH --job-name=test_stable
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --partition=infinite 

module load gcc/gcc-9.2.0-openmpi-3.1.4

cd PPP

echo "Job running on "`hostname`

CHARMM=/home/devereux/charmm/c45a2/exec/charmm-cpu-v2
host=$(hostname | sed "s/node//g" | sed "s/.cluster//g")
if (( host > 48 )); then
  CHARMM=/home/devereux/charmm/c45a2/exec/charmm-cpu-v3
fi

mpirun --bind-to none -np 8 $CHARMM < lambda_III.inp > lambda_III.log
mpirun --bind-to none -np 1 $CHARMM < lambda_III_trajread.inp > lambda_III_trajread.log

# calculate mean electrostatic energy from III_trajread.log
if [ -e lambda_III_trajread.log ]; then
  total=0
  count=0
  for i in $( grep "ENER EXTERN>" lambda_III_trajread.log | awk '{print $4}' ); do
    total=$( echo "scale=3; $total+$i" | bc )
    ((count++))
  done
  mean=$( echo "scale=5; DDD * $total / $count" | bc )
  echo "$mean" | bc > dG.dat

  # calculate variance (should not exceed 0.5 or we may need to subdivide lambda windows)
  var=0
  for i in $( grep "ENER EXTERN>" lambda_III_trajread.log | awk '{print $4}' ); do
    tt=$( echo "scale=5; $i - $mean" | bc )
    tt=$( echo "scale=5; $tt * $tt" | bc )
    var=$( echo "scale=5; $var + $tt" | bc )
  done
  var=$( echo "scale=5; $var / $count" | bc )
  echo "$var" > var.dat

fi

