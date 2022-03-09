#!/bin/bash

WORKDIR=/data2/devereux/MDCM-symmetry/MDCM-dyna/PhBr/deltaG
REFDIR=$WORKDIR/ref
SCRIPTDIR=$WORKDIR/scripts
EQUDIR=$WORKDIR/equilibrate

rm -r lambda_*
rm -r not_converged

cd $WORKDIR/mdcm_gas/r2.05e0.32

module load gcc/gcc-9.2.0-openmpi-4.1.2-ib

$SCRIPTDIR/ti.sh -dcm $REFDIR/12charges.dcm -par $WORKDIR/mdcm_gas/r2.05e0.32/phbr.par -top $REFDIR/phbr.top -slu $REFDIR/step5_slu.pdb -gas -res $EQUDIR/equil-gas.res

