#!/bin/bash

WORKDIR=/data2/devereux/MDCM-symmetry/MDCM-dyna/PhBr/deltaG
REFDIR=$WORKDIR/ref
SCRIPTDIR=$WORKDIR/scripts
EQUDIR=$WORKDIR/equilibrate

cd $WORKDIR/mdcm_solv/pme/r2.05e0.32
rm -r lambda_*
rm -r not_converged

module load gcc/gcc-9.2.0-openmpi-4.1.2-ib

$SCRIPTDIR/ti.sh -dcm $REFDIR/12charges.dcm -par $WORKDIR/mdcm_solv/pme/r2.05e0.32/phbr.par -top $REFDIR/phbr.top -slu $REFDIR/step5_slu.pdb -solv $REFDIR/step5_slv.pdb -res $EQUDIR/equil.res

