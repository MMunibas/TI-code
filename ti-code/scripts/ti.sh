#!/bin/bash

# Mike Devereux. Jan 2017.
#
# Description
# ===========
#
# Script to run manual CHARMM thermodynamic integration for solvation free energy of a solute
# molecule in a solvent box for a given parameter set. 
#
# TI calculations consist of 4 terms (VDW gas-phase, VDW solution-phase, ELEC gas-phase and
# ELEC solution-phase). In VDW gas-phase the electrostatics are turned off and slow-growth
# TI with CHARMM's "PERT" routines and PSSP soft-core option are used. The same approach is
# taken for VDW in solution but for the solute only (i.e. standard solvent parameters are used
# throughout).
#
# For electrostatic terms, two approaches are possible. For CHARMM or MDCM charges the
# standard CHARMM PERT routines can be used with slow growth. For multipoles or optionally for
# MDCM charges to allow more direct comparison with multipole results, delta G is evaluated by
# running a trajectory at a given Lambda, then reading in the CHARMM dcd file and evaluating
# the energy at Lambda=1 for each step and using the mean Lambda=1 energy as dG for that window
# (see Eq. 14, Bereau et al., JCTC, 2013, 9, 5450-5459). This means that TI takes place in 2
# stages per window, an initial simulation step run with scaled MTP parameters and a
# subsequent analysis step that reads in the trajectory and evaluates energies with unscaled
# (lamda=1) parameters.
#
# Usage
# =====
#
# The different options are described by calling the ti.sh script with no arguments. Examples
# are given at the end of this section.
#
# Before calling ti.sh, you need to create an equilibrated, solvated system with a restart file
# that can be used as a starting point for TI calculations. This is to avoid running long
# equilibration steps for every lambda window unnecessarily.
#
# You will also need to provide 2 PDB files, one for the solute and one for the solvent with 
# the solute removed (i.e. the solvent with a cavity where the solute should be). This way
# the 2 PDBs can be used to run both gas-phase and solvated systems, with some minor modification
# it would be possible to adapt the code to accept PDBs for the gas-phase molecule and whole
# solvated system including solute instead.
#
# Once These files have been prepared, you need to modify the user-defined variables in this
# script (at the start of the code below), and adapt the template files in the 
# scripts/templates folder for your own purposes. The ".sh" submission scripts should be
# adapted for the cluster environment you wish to use (job scheduler options, environment
# variables etc.). Note that "tmpl-mtpl.sh" handles the 2-step TI approach described above
# for multipolar simulations (without "PERT" in CHARMM), while "tmpl-vdw.sh" handles both VDW
# and electrostatic simulations that use CHARMM's PERT module with slow-growth.
#
# Once all files have been prepared, call the script once for each TI energy component (VDW,
# electrostatic, gas-phase and solution-phase). The script will launch runs for each TI
# window, check for crashes and resubmit, and upon job completion it will check whether any 
# windows need to be subdivided, and if so it will move the unconverged data to a new folder,
# split the TI window and submit 2 new jobs.
#
# For slow-growth with PERT the variance printed to the output can be used to check the need
# to subdivide. With the 2-step approach for MTPs the only rigorous approach is to re-run with
# smaller windows to check whether results remain the same for each windowsand have thus
# converged.
#
# For PERT slow-growth runs, it is best-practise (for puslished results) to repeat calculations
# in the reverse direction, a corresponding flag "-back" is provided for this purpose.
#
# Results are written to "deltaG.dat" at the end of each run, note that the sign of each term
# will depend on your definition of the lambda=0 and lambda=1 states
#
# Examples
# ========
#
# VDW gas-phase:
# ti.sh -vdw  $WORKDIR/phbr.par -top  $REFDIR/phbr.top -slu  $REFDIR/step5_slu.pdb -gas -res $EQUDIR/equil-gas.res
#
# DCM gas-phase:
# ti.sh -dcm $REFDIR/12charges.dcm -par $WORKDIR/phbr.par -top $REFDIR/phbr.top -slu $REFDIR/step5_slu.pdb -gas -res $EQUDIR/equil-gas.res
#
# VDW solution-phase:
# ti.sh -vdw $WORKDIR/phbr.par -top $REFDIR/phbr.top -slu $REFDIR/step5_slu.pdb -solv $REFDIR/step5_slv.pdb -res $EQUDIR/equil.res
#
# DCM solution-phase:
# ti.sh -dcm $REFDIR/12charges.dcm -par $WORKDIR/phbr.par -top $REFDIR/phbr.top -slu $REFDIR/step5_slu.pdb -solv $REFDIR/step5_slv.pdb -res $EQUDIR/equil.res
#

#---------------------------------------------------------------------------------------------

# Set some variables:

# User-defined variables:

D_LAMBDA=0.05 # default size of lambda window
NSTEPS=150000  # default number of integration steps per lambda window
NEQUIL=50000 # number of dynamics equilibration steps (not used for T.I.)
TMPL=/data2/devereux/MDCM-symmetry/MDCM-dyna/PhBr/deltaG/scripts
#CHMM=/opt/cluster/programs/charmm/developer/dev-release-dcm/build/cmake/charmm
CHMM=/home/devereux/bin/dev-release-dcm/build/cmake/charmm
RES=/home/devereux/mdcm-phf-ti/mdcm-ti/ref/step5.res
LSTART=0.0 # start lambda value
LSTOP=1.0 # end lambda value


# Fixed variables (program defaults)

GAS=1     # gas-phase run by default
SOLV=0    # solution-phase available as option
mtpl=0    # flag for multipole run
dcm=0     # DCM approach with PERT module (slow growth)
dcm2=0    # DCM run with tristan windowed integration
vdw=0
BACK=0    # Flag to turn on "backwards" slow growth run as convergence check

# flag to see whether this is a "master" or subdivided lambda run
SUBDIV=0

printf "\e[1;32m\n\n\n\tSTARTING THERMODYNAMIC INTEGRATION SCRIPT AT `date`\e[0m\n\n"

# get cwd:
WORK_DIR=`pwd`
# get folder of ti.sh script:
pushd `dirname $0` > /dev/null
SCRIPT_DIR=`pwd -P`
popd > /dev/null


######################################################################################################
#
#                                     FUNCTIONS
#

function monitor_job {
  printf "\n\n\n"
  while true; do
    sleep 1m
    nrun=0
    for i in $(seq 1 $NWIN); do
      squeue -u $USER |  tail -n +2 | awk '{print $1}' | grep ${jobid[$i]} > /dev/null && ((nrun++))
    done
    tim=`date +%H:%M`
    day=`date +%d/%m/%y`
    if [ $nrun -eq 0 ]; then
      printf "\n\n\n \e[1;32mAll jobs finished at $tim on $day, checking output files...\e[0m\n\n\n"
      break
    fi
#    echo -ne "$nrun / $NWIN Jobs Still Running (last checked $tim)     \r"
    for i in {1..30}; do
      sleep 1
#      echo -ne "$nrun / $NWIN Jobs Still Running (last checked $tim) .        \r"
      sleep 1
#      echo -ne "$nrun / $NWIN Jobs Still Running (last checked $tim) ..       \r"
      sleep 1
#      echo -ne "$nrun / $NWIN Jobs Still Running (last checked $tim) ...      \r"
      sleep 1
#      echo -ne "$nrun / $NWIN Jobs Still Running (last checked $tim) ....     \r"
      sleep 1
#      echo -ne "$nrun / $NWIN Jobs Still Running (last checked $tim) .....    \r"
    done
  done
  echo -ne '\n'
}


function phelp {
  printf "\n\n Usage: ti.sh [-mtp [punchfile] -top [topfile] -par [parfile] -slu [pdbfile]]\n"
  printf "              [-dcm [chgfile] -top [topfile] -par [parfile] -slu [pdbfile] -res [restartfile]]\n"
  printf "              [-dcm2 [chgfile] -top [topfile] -par [parfile] -slu [pdbfile] -res [restartfile]]\n"
  printf "              [-vdw [parfile] -top [topfile] -slu [pdbfile] [-solv [pdbfile] -qsol [chgfile]] -res [restartfile]]\n\n"
  printf " -dcm:\t\t request a DCM run with slow-growth (PERT module) integration\n"
  printf " -dcm2:\t\t request a DCM run with windowed integration\n"
  printf " -mtp:\t\t request an MTP run\n"
  printf " -vdw:\t\t request a vdw run, add -qsol to specify lpun or dcm file for solvent molecules\n\n\n"
  printf " -qsol:\t\t optionally specify lpun or dcm file for solvent molecules for vdw run\n\n\n"
  printf " -gas:\t\t gas-phase run (default)\n"
  printf " -solv [pdbfile]:\t\t solution-phase run, pure solvent PDB must be given\n\n\n"
  printf " -lstart:\t\t start lambda value (default=0)\n"
  printf " -lstop:\t\t end lambda value (default=1)\n"
  printf " -dlam:\t\t lambda step size (default=0.1)\n\n\n"
  printf " -back:\t\t backwards slow growth run\n\n\n"
  exit
}

function resub() {
  cd lambda_$LAMBDA
  jobid[$i]=`sbatch lambda_$LAMBDA.sh | grep "Submitted batch" | awk '{print $4}'`
  printf " Job ID "${jobid[$i]}"\n\n\n"
  cd ..
}

###############################################################################################


# Parse arguments

if [ $# -eq 0 ]; then
  phelp
fi

while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -mtp)
    mtpl=1
    PUNFILE=$2
    TFILE=`basename $2`
    ARGS="$ARGS -mtp $PUNFILE"
    if [ -z $PUNFILE ]; then
      printf "\e[1;31mError: Punch file must be specified for mtp run. Type ti.sh -h for help.\e[0m\n\n"
      exit
    fi
    if ! [ -e $PUNFILE ]; then
      printf "\e[1;31mError: Punch file $PUNFILE not found in path `pwd`\e[0m\n\n"
      exit
    fi
    if [ $dcm -eq 1 ] || [ $vdw -eq 1 ] || [ $dcm2 -eq 1 ]; then
      printf "\e[1;31mError: options DCM, VDW and MTP are mutually exclusive\e[0m\n\n"
      exit
    fi
    if [ $SUBDIV -eq "0" ]; then
      printf "\n\n Multipole (mtp) run selected\n\n"
    fi
    shift # past argument
    ;;
    -dcm)
    dcm=1
    DCMFILE=$2
    TFILE=`basename $2`
    ARGS="$ARGS -dcm $DCMFILE"
    if [ -z $DCMFILE ]; then
      printf "\e[1;31mError: Charge file must be specified for DCM run. Type ti.sh -h for help.\e[0m\n\n"
      exit
    fi
    if ! [ -e $DCMFILE ]; then
      printf "\e[1;31mError: DCM charge file $DCMFILE not found in path `pwd`\e[0m\n\n"
      exit
    fi
    if [ $mtpl -eq 1 ] || [ $vdw -eq 1 ] || [ $dcm2 -eq 1 ]; then   
      printf "\e[1;31mError: options DCM, VDW and MTP are mutually exclusive\e[0m\n\n"
      exit
    fi
    if [ $SUBDIV -eq "0" ]; then
      printf "\n\n DCM run with slow growth (PERT module) integration selected\n\n"
    fi
    shift # past argument
    ;;
    -qsol)
    DCMFILE=$2
    TFILE=`basename $2`
    ARGS="$ARGS -qsol $DCMFILE"
    if [ -z $DCMFILE ]; then
      printf "\e[1;31mError: Charge file not specified with -qsol. Type ti.sh -h for help.\e[0m\n\n"
      exit
    fi
    if ! [ -e $DCMFILE ]; then
      printf "\e[1;31mError: charge file $DCMFILE not found in path `pwd`\e[0m\n\n"
      exit
    fi
    if [ $SUBDIV -eq "0" ]; then
      printf "\n\n Charge file for VDW run solvent molecules selected\n\n"
    fi
    shift # past argument
    ;;
    -dcm2)
    dcm2=1
    DCMFILE=$2
    TFILE=`basename $2`
    ARGS="$ARGS -dcm2 $DCMFILE"
    if [ -z $DCMFILE ]; then
      printf "\e[1;31mError: Charge file must be specified for DCM run. Type ti.sh -h for help.\e[0m\n\n"
      exit
    fi
    if ! [ -e $DCMFILE ]; then
      printf "\e[1;31mError: DCM charge file $DCMFILE not found in path `pwd`\e[0m\n\n"
      exit
    fi
    if [ $mtpl -eq 1 ] || [ $vdw -eq 1 ] || [ $dcm -eq 1 ]; then
      printf "\e[1;31mError: options DCM, VDW and MTP are mutually exclusive\e[0m\n\n"
      exit
    fi
    if [ $SUBDIV -eq "0" ]; then
      printf "\n\n DCM run with windowed integration selected\n\n"
    fi
    shift # past argument
    ;;
    -vdw)
    vdw=1
    PARAM=$2
    TFILE=`basename $2`
    ARGS="$ARGS -vdw $PARAM"
    if [ -z $PARAM ]; then
      printf "\e[1;31mError: Parameter file must be specified for vdw run. Type ti.sh -h for help.\e[0m\n\n"
      exit
    fi
    if ! [ -e $DCMFILE ]; then
      printf "\e[1;31mError: DCM charge file $DCMFILE not found in path `pwd`\e[0m\n\n"
      exit
    fi
    if [ $mtpl -eq 1 ] || [ $dcm -eq 1 ] || [ $dcm2 -eq 1 ]; then
      printf "\e[1;31mError: options DCM, VDW and MTP are mutually exclusive\e[0m\n\n"
      exit
    fi
    if [ $SUBDIV -eq "0" ]; then
      printf "\n\n Van der Waals (vdw) run selected\n\n"
    fi
    shift # past argument
    ;;
    -res)
    RES=$2
    TFILE=`basename $2`
    ARGS="$ARGS -res $RES"
    if [ -z $RES ]; then
      printf "\e[1;31mError: Restart file must be specified with -res. Type ti.sh -h for help.\e[0m\n\n"
      exit
    fi
    if [ $SUBDIV -eq "0" ]; then
      printf "\n\n Restart file $RES specified\n\n"
    fi
    shift # past argument
    ;;
    -top)
    TOPOL="$2"  # topology file for solute
    TFILE=`basename $2`
    ARGS="$ARGS -top $TOPOL"
    if [ -z $TOPOL ]; then
      printf "\e[1;31mError: Topology file must be specified. Type ti.sh -h for help.\e[0m\n\n"
      exit
    fi
    if ! [ -e $TOPOL ]; then
      printf "\e[1;31mError: Topology file $TOPOL not found in path `pwd`\e[0m\n\n"
      exit
    fi
    shift # past argument
    ;;
    -slu)
    SLU="$2"  # PDB file for solute
    TFILE=`basename $2`
    ARGS="$ARGS -slu $SLU"
    if [ -z $SLU ]; then
      printf "\e[1;31mError: Solute PDB file must be specified. Type ti.sh -h for help.\e[0m\n\n"
      exit
    fi
    if ! [ -e $SLU ]; then
      printf "\e[1;31mError: PDB file $SLU not found in path `pwd`\e[0m\n\n"
      exit
    fi
    shift # past argument
    ;;
    -par)
    PARAM="$2"  # parameter file for solute
    TFILE=`basename $2`
    ARGS="$ARGS -par $PARAM"
    if [ -z $PARAM ]; then
      printf "\e[1;31mError: Parameter file must be specified. Type ti.sh -h for help.\e[0m\n\n"
      exit
    fi
    if ! [ -e $PARAM ]; then
      printf "\e[1;31mError: Parameter file $PARAM not found in path `pwd`\e[0m\n\n"
      exit
    fi
    shift # past argument
    ;;
    -solv)
    SOLV=1 # solvated run selected
    GAS=0
    SLV="$2"  # parameter file for solute
    TFILE=`basename $2`
    ARGS="$ARGS -solv $SLV"
    printf " Solution phase run selected\n\n"
    if [ -z $SLV ]; then
      printf "\e[1;31mError: Solvent PDB file must be specified. Type ti.sh -h for help.\e[0m\n\n"
      exit
    fi
    if ! [ -e $SLV ]; then
      printf "\e[1;31mError: Solvent PDB file $SLV not found in path `pwd`\e[0m\n\n"
      exit
    fi
    shift # past argument
    ;;
    -gas)
    SOLV=0 
    GAS=1 # gas-phase run selected
    ARGS="$ARGS -gas"
    printf " Gas phase run selected\n\n"
    ;;
    -subdivide)
    SUBDIV=1
    ;;
    -back)
    BACK=1
    ARGS="$ARGS -back"
    printf " Backwards slow growth run selected\n\n"
    ;;
    -lstart)
    LSTART=$2
    if (( $(echo "$LSTART < 0.0" | bc -l) )); then
      printf "\e[1;31mError: LSTART (start lambda value) must be >= 0\e[0m\n\n"
    fi
    if (( $(echo "$LSTART > 1.0" | bc -l) )); then
      printf "\e[1;31mError: LSTART (start lambda value) must be <= 1\e[0m\n\n"
    fi
    printf " Start lambda value = $LSTART\n\n"
    shift # past argument
    ;;
    -lstop)
    LSTOP=$2
    if (( $(echo "$LSTOP < 0.0" | bc -l) )); then
      printf "\e[1;31mError: LSTOP (end lambda value) must be >= 0\e[0m\n\n"
    fi
    if (( $(echo "$LSTOP > 1.0" | bc -l) )); then
      printf "\e[1;31mError: LSTOP (end lambda value) must be <= 1\e[0m\n\n"
    fi
    printf " End lambda value = $LSTOP\n\n"
    shift # past argument
    ;;
    -dlam)
    D_LAMBDA=$2
    if (( $(echo "$D_LAMBDA < 0.0" | bc -l) )); then
      printf "\e[1;31mError: D_LAMBDA (lambda window size) must be >= 0\e[0m\n\n"
    fi
    if (( $(echo "$D_LAMBDA > 1.0" | bc -l) )); then
      printf "\e[1;31mError: D_LAMBDA (lambda window size) must be <= 1\e[0m\n\n"
    fi
    printf " lambda window size = $D_LAMBDA\n\n"
    shift # past argument
    ;;
    -h)
      phelp
    ;;
    *)
      printf "\e[1;31mError: unkown option $1. Type ti.sh -h for help.\e[0m\n\n"
      exit
    ;;
esac
shift # past argument or value
done

# Check consistency of arguments

if [ $dcm -eq 0 ] && [ $mtpl -eq 0 ] && [ $vdw -eq 0 ] && [ $dcm2 -eq 0 ]; then   
  printf "\e[1;31mError: either DCM, VDW or MTP must be selected. Run ti.sh -h for details.\e[0m\n\n"
  exit
fi

if [ -z $TOPOL ]; then
  printf "\e[1;31mError: Topology file not specified. Type ti.sh -h for help.\e[0m\n\n"
  exit
fi

if [ -z $SLU ]; then
  printf "\e[1;31mError: Solvent PDB file not specified. Type ti.sh -h for help.\e[0m\n\n"
  exit
fi

if [ -z $PARAM ]; then
  printf "\e[1;31mError: Parameter file not specified. Type ti.sh -h for help.\e[0m\n\n"
  exit
fi



# check that lambda window makes sense:
LRANGE=`bc <<< "scale=7; $LSTOP - $LSTART"`
TMP=`bc <<< "scale = 7; $LRANGE / $D_LAMBDA"`
NWIN=`bc <<< "scale = 0; $LRANGE / $D_LAMBDA"`
TMP2=`bc <<< "scale = 7; $TMP / $NWIN"`

if [ `bc <<< "$TMP2 == 1.0"` = 0 ]; then
  printf "\e[1;31mError: D_LAMBDA must form NWIN complete windows between $LSTART and $LSTOP (currently $TMP)\e[0m\n\n"
  exit 0
fi

printf " $NWIN lambda windows will be run\n\n"


# check that sufficient time steps were requested
TMP=`bc <<< "scale = 1; $NSTEPS - $NEQUIL"`
if [ $TMP -lt 100000 ]; then
  printf "\e[1;31m Severe Warning: short run ($TMP production steps) requested. Please ensure job was equilibrated \n beforehand, plus check that NEQUIL is small enough to leave sufficient production data. \n Modify NSTEPS variable in ti.sh to increase simulation length, NEQUIL to reduce equilibration steps.\e[0m\n\n"
fi

# check that NEQUIL makes sense
if [ $NEQUIL -lt 1 ] || [ $NEQUIL -ge $NSTEPS ]; then
  printf "\e[1;31m Error: NEQUIL parameter in ti.sh is $NEQUIL. It must be greater than 0 and less than NSTEPS (currently $NSTEPS).\e[0m\n\n"
  exit 0
fi


# check that templates folder exists in directory "$TMPL":
if ! [ -e $TMPL/templates ]; then
  printf "\e[1;31mError: templates/ folder containing CHARMM template files does not exist in directory\e[0m\n\n"
  exit 0
else
  if [ $SUBDIV -eq "0" ]; then
    printf " Templates folder exists at $TMPL\n\n"
  fi
fi


# check that no folders named lambda_* already exist:
if [ $SUBDIV -eq "0" ]; then
  if ( ls -1p | grep "lambda_" | grep "/" > /dev/null); then
    printf "\e[1;31mFolders named lambda_xxxx already exist in this directory. Please remove them and restart\e[0m\n\n"
    exit 0
  fi
fi

# Create some folders:

if [ $SUBDIV -eq "0" ]; then
  printf " Creating lambda_* subfolders in directory `pwd`\n\n"
fi

LAMBDA1=`bc <<< "scale=7; $D_LAMBDA / 2"`

for i in $(seq 1 $NWIN); do

  LAMBDA=`bc <<< "scale = 7; $LSTART + $D_LAMBDA * $i - $LAMBDA1"`
  mkdir lambda_$LAMBDA

done

if [ $SUBDIV -eq "0" ]; then
  printf " Creating common/ folder for common CHARMM files\n\n"
  if [ -e common ]; then
    printf "\e[1;33m Warning: common/ folder containing CHARMM common files already exists in this directory, overwriting...\e[0m\n\n"
  fi
  mkdir -p common || exit 0

  # copy CHARMM topology files for run to common folder
  #cp $TMPL/templates/*.top common || exit 0

  # copy solute topology file to common folder
  TFILE=`basename $TOPOL`
  diff $TOPOL common/$TFILE >& /dev/null || cp $TOPOL common 

  # copy solvent CHARMM parameter files for run to common folder
  cp $TMPL/templates/*.par common || exit 0
  cp $RES common || exit 0
  # copy solute parameter file to common folder
  TFILE=`basename $PARAM`
  diff $PARAM common/$PARAM >& /dev/null || cp $PARAM common

  # copy PDB files for run to common folder
  #cp $TMPL/templates/*.pdb common || exit 0

  # copy solute PDB file to common folder
  TFILE=`basename $SLU`
  diff $SLU common/$SLU >& /dev/null || cp $SLU common

  # copy solvent PDB file to common folder
  if [ $SOLV -eq 1 ]; then
    TFILE=`basename $SLV`
    diff $SLV common/$SLV >& /dev/null || cp $SLV common
  fi

  # copy stream file to common folder
  #cp $TMPL/templates/*.str common || exit 0

  # copy zeroed topology file to common folder
  #if [ $SUBDIV -eq "0" ]; then
  printf " Creating zeroed topology file in common folder\n\n"
  #fi
  $TMPL/zero_topology.pl $TOPOL $SLU; mv zeroed.top common || exit 0
fi

TOPOL=`basename $TOPOL`
SLU=`basename $SLU`
RES=`basename $RES`
if [ $SOLV -eq 1 ]; then
  SLV=`basename $SLV`
fi
PARAM=`basename $PARAM`


##### MULTIPOLE JOBS

if [ $mtpl -eq 1 ]; then  
  if [ $BACK -eq "1" ]; then
    printf "\e[1;33m Warning: Backwards slow growth selected for MTPL run without slow-growth, which\n doesn't make sense. Running forwards T.I. as usual.\e[0m\n\n"
  fi
  # copy lpun files to common folder
  [ -e $TMPL/templates/*.lpun ] && cp $TMPL/templates/*.lpun common # lpun files for surrounding moieties
  cp $PUNFILE common/lambda_1.0.lpun || exit 0 # solute lpun file
  printf " Writing scaled input files and submission scripts\n\n"
  # write input file for this lambda:
  if [ $SOLV -eq 1 ]; then
    tmplfile=tmpl-solv-mtpl
  else
    tmplfile=tmpl-gas-mtpl
  fi
  for i in $(seq 1 $NWIN); do
    TMP=`bc <<< "scale = 7; $NSTEPS - $NEQUIL"`
    LAMBDA=`bc <<< "scale = 7; $LSTART + $D_LAMBDA * $i - $LAMBDA1"`
    sed "s/LLL/$LAMBDA/g" $TMPL/templates/$tmplfile.inp > lambda_$LAMBDA/lambda_$LAMBDA.inp
    sed -i "s/SSS/$TMP/g" lambda_$LAMBDA/lambda_$LAMBDA.inp
    sed -i "s:RRR:$PARAM:g" lambda_$LAMBDA/lambda_$LAMBDA.inp
    sed -i "s:EEE:$NEQUIL:g" lambda_$LAMBDA/lambda_$LAMBDA.inp
    sed -i "s:TTT:$SLU:g" lambda_$LAMBDA/lambda_$LAMBDA.inp
    sed -i "s:UUU:$SLV:g" lambda_$LAMBDA/lambda_$LAMBDA.inp
    sed "s/LLL/$LAMBDA/g" "$TMPL/templates/"$tmplfile"-trajread.inp" > lambda_$LAMBDA/lambda_$LAMBDA"_trajread.inp"
    NREAD=$( bc <<< "scale=0; $TMP / 100" ) # no. steps stored in trajectory file (assumes write every 100)
    sed "s/NNN/$NREAD/g" "$TMPL/templates/loop.str" > common/loop.str
  # create submission script
    sed "s/NNN/ti$LAMBDA/g" $TMPL/templates/tmpl-mtpl.sh > lambda_$LAMBDA/lambda_$LAMBDA.sh
    sed -i "s:CCC:$CHMM:g" lambda_$LAMBDA/lambda_$LAMBDA.sh
    sed -i "s:PPP:`pwd`/lambda_$LAMBDA/:g" lambda_$LAMBDA/lambda_$LAMBDA.sh
    sed -i "s:III:$LAMBDA:g" lambda_$LAMBDA/lambda_$LAMBDA.sh
    sed -i "s:DDD:$D_LAMBDA:g" lambda_$LAMBDA/lambda_$LAMBDA.sh
    chmod ug+x lambda_$LAMBDA/lambda_$LAMBDA.sh

  # Run charmm job
    printf " \e[1;32mSubmitting job lambda_$LAMBDA.inp\e[0m\n"
    cd lambda_$LAMBDA
    jobid[$i]=`sbatch lambda_$LAMBDA.sh | grep "Submitted batch" | awk '{print $4}'`
    printf " Job ID "${jobid[$i]}"\n\n\n"
    cd ..
  done
fi


##### DCM JOBS WITH TRISTAN INTEGRATION

if [ $dcm2 -eq 1 ]; then
  if [ $BACK -eq "1" ]; then
    printf "\e[1;33m Warning: Backwards slow growth selected for DCM run without slow-growth, which\n doesn't make sense. Running forwards T.I. as usual.\e[0m\n\n"
  fi

  # copy lpun files to common folder
  [ -e $TMPL/templates/*.dcm ] && cp $TMPL/templates/*.dcm common # lpun files for surrounding moieties
  cp $DCMFILE common/lambda_1.0.dcm || exit 0 # solute lpun file
  if [ $SUBDIV -eq "0" ]; then
    printf " Writing scaled input files and submission scripts\n\n"
  fi
  # write input file for this lambda:
  if [ $SOLV -eq 1 ]; then
    tmplfile=tmpl-solv-dcm2
  else
    tmplfile=tmpl-gas-dcm2
  fi
  for i in $(seq 1 $NWIN); do
    TMP=`bc <<< "scale = 7; $NSTEPS - $NEQUIL"`
    LAMBDA=`bc <<< "scale = 7; $LSTART + $D_LAMBDA * $i - $LAMBDA1"`
    $TMPL/scale_dcm.pl common/lambda_1.0.dcm $LAMBDA > lambda_$LAMBDA/lambda_$LAMBDA.dcm
    sed "s/LLL/$LAMBDA/g" $TMPL/templates/$tmplfile.inp > lambda_$LAMBDA/lambda_$LAMBDA.inp
    sed -i "s/SSS/$TMP/g" lambda_$LAMBDA/lambda_$LAMBDA.inp
    sed -i "s:RRR:$PARAM:g" lambda_$LAMBDA/lambda_$LAMBDA.inp
    sed -i "s:TTT:$SLU:g" lambda_$LAMBDA/lambda_$LAMBDA.inp
    sed -i "s:UUU:$SLV:g" lambda_$LAMBDA/lambda_$LAMBDA.inp
    sed -i "s:EEE:$NEQUIL:g" lambda_$LAMBDA/lambda_$LAMBDA.inp
    sed -i "s:ZZZ:$RES:g" lambda_$LAMBDA/lambda_$LAMBDA.inp
    sed "s/LLL/$LAMBDA/g" "$TMPL/templates/"$tmplfile"-trajread.inp" > lambda_$LAMBDA/lambda_$LAMBDA"_trajread.inp"
    sed -i "s:RRR:$PARAM:g" lambda_$LAMBDA/lambda_$LAMBDA"_trajread.inp"
    sed -i "s:TTT:$SLU:g" lambda_$LAMBDA/lambda_$LAMBDA"_trajread.inp"
    sed -i "s:UUU:$SLV:g" lambda_$LAMBDA/lambda_$LAMBDA"_trajread.inp"
    NREAD=$( bc <<< "scale=0; $TMP / 100" ) # no. steps stored in trajectory file (assumes write every 100)
    sed "s/NNN/$NREAD/g" "$TMPL/templates/loop.str" > common/loop.str
  # create submission script
    sed "s/NNN/ti$LAMBDA/g" $TMPL/templates/tmpl-mtpl.sh > lambda_$LAMBDA/lambda_$LAMBDA.sh
    sed -i "s:CCC:$CHMM:g" lambda_$LAMBDA/lambda_$LAMBDA.sh
    sed -i "s:PPP:`pwd`/lambda_$LAMBDA/:g" lambda_$LAMBDA/lambda_$LAMBDA.sh
    sed -i "s:III:$LAMBDA:g" lambda_$LAMBDA/lambda_$LAMBDA.sh
    sed -i "s:DDD:$D_LAMBDA:g" lambda_$LAMBDA/lambda_$LAMBDA.sh
    chmod ug+x lambda_$LAMBDA/lambda_$LAMBDA.sh

  # Run charmm job
    printf " \e[1;32mSubmitting job lambda_$LAMBDA.inp\e[0m\n"
    cd lambda_$LAMBDA
    jobid[$i]=`sbatch lambda_$LAMBDA.sh | grep "Submitted batch" | awk '{print $4}'`
    printf " Job ID "${jobid[$i]}"\n\n\n"
    cd ..
  done
fi

###### DCM JOBS WITH SLOW GROWTH (PERT MODULE)

if [ $dcm -eq 1 ]; then
  if [ $SUBDIV -eq "0" ]; then
    # copy dcm charge files to common folder
    [ -e $TMPL/templates/*.dcm ] && cp $TMPL/templates/*.dcm common # DCM files for surrounding moieties
    cp $DCMFILE common/lambda_1.0.dcm || exit 0 # solute DCM file
  fi
  printf " Writing scaled input files and submission scripts\n\n"
  # write input file for this lambda:
  if [ $SOLV -eq 1 ]; then
    tmplfile=tmpl-solv-dcm
  else
  tmplfile=tmpl-gas-dcm
  fi
  if [ $SUBDIV -eq "0" ]; then
    # create zeroed DCM charge file
    printf " Creating zeroed DCM charge file as state LAMBDA=0\n\n"
    $TMPL/zero_dcm.pl $DCMFILE common/$SLU > common/lambda_0.0.dcm
  fi
  for i in $(seq 1 $NWIN); do
    LAMBDA=`bc <<< "scale = 7; $LSTART + $D_LAMBDA * $i - $LAMBDA1"`
    if [ $BACK -eq "0" ]; then
      TIVARS=`$TMPL/lambda-pars.pl $LAMBDA $D_LAMBDA $NSTEPS $NEQUIL "f" | sed -n 1p`
      TIVARS2=`$TMPL/lambda-pars.pl $LAMBDA $D_LAMBDA $NSTEPS $NEQUIL "f" | sed -n 2p`
    else
      TIVARS=`$TMPL/lambda-pars.pl $LAMBDA $D_LAMBDA $NSTEPS $NEQUIL "b" | sed -n 1p`
      TIVARS2=`$TMPL/lambda-pars.pl $LAMBDA $D_LAMBDA $NSTEPS $NEQUIL "b" | sed -n 2p`
    fi
    sed "s/LLL/$TIVARS\n$TIVARS2/g" $TMPL/templates/$tmplfile.inp > lambda_$LAMBDA/lambda_$LAMBDA.inp
    sed -i "s/SSS/$NSTEPS/g" lambda_$LAMBDA/lambda_$LAMBDA.inp
    sed -i "s:RRR:$PARAM:g" lambda_$LAMBDA/lambda_$LAMBDA.inp
    sed -i "s:TTT:$SLU:g" lambda_$LAMBDA/lambda_$LAMBDA.inp
    sed -i "s:UUU:$SLV:g" lambda_$LAMBDA/lambda_$LAMBDA.inp
    sed -i "s:ZZZ:$RES:g" lambda_$LAMBDA/lambda_$LAMBDA.inp
  # create submission script
    sed "s/NNN/ti$LAMBDA/g" $TMPL/templates/tmpl-vdw.sh > lambda_$LAMBDA/lambda_$LAMBDA.sh
    sed -i "s:CCC:$CHMM:g" lambda_$LAMBDA/lambda_$LAMBDA.sh
    sed -i "s:PPP:`pwd`/lambda_$LAMBDA/:g" lambda_$LAMBDA/lambda_$LAMBDA.sh
    sed -i "s:III:lambda_$LAMBDA:g" lambda_$LAMBDA/lambda_$LAMBDA.sh
    chmod ug+x lambda_$LAMBDA/lambda_$LAMBDA.sh

  # Run charmm job
    printf " \e[1;32mSubmitting job lambda_$LAMBDA.inp\e[0m\n"
    cd lambda_$LAMBDA
    jobid[$i]=`sbatch lambda_$LAMBDA.sh | grep "Submitted batch" | awk '{print $4}'`
    printf " Job ID "${jobid[$i]}"\n\n\n"
    cd ..
  done
fi

###### VDW JOBS

if [ $vdw -eq 1 ]; then
  if [ $SUBDIV -eq "0" ]; then
    # copy lpun files to common folder
    [ -e $TMPL/templates/*.lpun ] && echo "copying MTP file to common"
    [ -e $TMPL/templates/*.lpun ] && cp $TMPL/templates/*.lpun common # lpun files for surrounding moieties
    printf " Writing scaled input files and submission scripts\n\n"
    if [ ! -z $DCMFILE ] && [ -e $DCMFILE ]; then
      cp $DCMFILE common # lpun files for surrounding moieties
      echo "copying DCM file $DCMFILE"
    fi
  fi
  # write input file for this lambda:
  if [ $SOLV -eq 1 ]; then
    tmplfile=tmpl-solv-vdw
  else
    tmplfile=tmpl-gas-vdw
  fi
  for i in $(seq 1 $NWIN); do
    LAMBDA=`bc <<< "scale = 7; $LSTART + $D_LAMBDA * $i - $LAMBDA1"`
    if [ $BACK -eq "0" ]; then
      TIVARS=`$TMPL/lambda-pars.pl $LAMBDA $D_LAMBDA $NSTEPS $NEQUIL "f" | sed -n 1p`
      TIVARS2=`$TMPL/lambda-pars.pl $LAMBDA $D_LAMBDA $NSTEPS $NEQUIL "f" | sed -n 2p`
    else
      TIVARS=`$TMPL/lambda-pars.pl $LAMBDA $D_LAMBDA $NSTEPS $NEQUIL "b" | sed -n 1p`
      TIVARS2=`$TMPL/lambda-pars.pl $LAMBDA $D_LAMBDA $NSTEPS $NEQUIL "b" | sed -n 2p`
    fi
    if [ -z "$TIVARS" ]; then
      printf "\n\n\n \e[1;31m Encountered a fatal problem with the lambda-pars.pl script (null output). Exiting.\e[0m\n\n"
      exit
    fi 
    if [ -z "$TIVARS2" ]; then
      printf "\n\n\n \e[1;31m Encountered a fatal problem with the lambda-pars.pl script (null output). Exiting.\e[0m\n\n"
      exit 
    fi
    sed "s/LLL/$TIVARS\n PSSP $TIVARS2/g" $TMPL/templates/$tmplfile.inp > lambda_$LAMBDA/lambda_$LAMBDA.inp
    sed -i "s/SSS/$NSTEPS/g" lambda_$LAMBDA/lambda_$LAMBDA.inp
    sed -i "s:RRR:$PARAM:g" lambda_$LAMBDA/lambda_$LAMBDA.inp
    sed -i "s:TTT:$SLU:g" lambda_$LAMBDA/lambda_$LAMBDA.inp
    sed -i "s:UUU:$SLV:g" lambda_$LAMBDA/lambda_$LAMBDA.inp
    sed -i "s:DDD:$DCMFILE:g" lambda_$LAMBDA/lambda_$LAMBDA.inp
    sed -i "s:ZZZ:$RES:g" lambda_$LAMBDA/lambda_$LAMBDA.inp
  # create submission script
    sed "s/NNN/ti$LAMBDA/g" $TMPL/templates/tmpl-vdw.sh > lambda_$LAMBDA/lambda_$LAMBDA.sh
    sed -i "s:CCC:$CHMM:g" lambda_$LAMBDA/lambda_$LAMBDA.sh
    sed -i "s:PPP:`pwd`/lambda_$LAMBDA/:g" lambda_$LAMBDA/lambda_$LAMBDA.sh
    sed -i "s:III:lambda_$LAMBDA:g" lambda_$LAMBDA/lambda_$LAMBDA.sh
    chmod ug+x lambda_$LAMBDA/lambda_$LAMBDA.sh

  # Run charmm job
    printf " \e[1;32mSubmitting job lambda_$LAMBDA.inp\e[0m\n"
    cd lambda_$LAMBDA
    jobid[$i]=`sbatch lambda_$LAMBDA.sh | grep "Submitted batch" | awk '{print $4}'`
    printf " Job ID "${jobid[$i]}"\n\n\n"
    cd ..
  done
fi

printf "\n\n\n Submitted jobs at `date`\n\n"

monitor_job

printf " Checking to see whether all jobs finished successfully\n\n"


# Check to see whether jobs completed successfully
for j in 1 2; do
    CRASH=0
    for i in $(seq 1 $NWIN); do
      LAMBDA=`bc <<< "scale = 7; $LSTART + $D_LAMBDA * $i - $LAMBDA1"`
      # Did trajectory run?
      if ! [ -e lambda_$LAMBDA/lambda_$LAMBDA.log ]; then
        printf "\n\n\n \e[1;31m Job lambda_$LAMBDA.inp has no log file.\e[0m\n\n"
        CRASH=1
        if [ $j -eq 2 ]; then
          printf "\n\n\n \e[1;31m Exiting script, please check for problems with job submission.\e[0m\n\n\n"
          exit
        fi
        LOG=lambda_$LAMBDA.log
        resub
        continue 
      fi
      # Did analysis job run? (MTPL runs only)
      if [ $mtpl -eq 1 ]; then
        if ! [ -e lambda_$LAMBDA/lambda_$LAMBDA"_trajread.log" ]; then
          printf "\n\n\n \e[1;31m Job lambda_"$LAMBDA"_trajread.inp has no log file.\e[0m\n\n"
          CRASH=1
          if [ $j -eq 2 ]; then
            printf "\n\n\n \e[1;31m Exiting script, please check for problems with job submission.\e[0m\n\n\n"
            exit
          fi
          LOG=lambda_$LAMBDA"_trajread.log"
          resub
          continue
        fi
      fi
      # Did trajectory finish correctly?
      if ! ( grep " NORMAL TERMINATION" lambda_$LAMBDA/lambda_$LAMBDA.log > /dev/null ); then
        if [ $j -eq 2 ]; then
          printf "\n\n\n \e[1;31m Job lambda_$LAMBDA.log crashed for second time. Exiting script, please examine CHARMM output file lambda_$LAMBDA.log for more information.\e[0m\n\n\n"
          exit 
        fi
        CRASH=1
        printf "\n\n\n \e[1;31m Job lambda_$LAMBDA.log appears to have crashed. Moving log file to lambda_$LAMBDA/lambda_$LAMBDA.log.crash and resubmitting\e[0m\n\n\n"
        LOG=lambda_$LAMBDA.log
        mv lambda_$LAMBDA/$LOG lambda_$LAMBDA/$LOG.crash
        resub
        continue
      else
        printf " Job lambda_$LAMBDA.inp finished successfully\n"
      fi
      # Did analysis job finish correctly? (MTPL runs only)
      if [ $mtpl -eq 1 ] || [ $dcm2 -eq 1 ] ; then
        if ! ( grep " NORMAL TERMINATION" lambda_$LAMBDA/lambda_$LAMBDA"_trajread.log" > /dev/null ); then
          if [ $j -eq 2 ]; then
            printf "\n\n\n \e[1;31m Job lambda_"$LAMBDA"_trajread.log crashed for second time. Exiting script, please examine CHARMM output file lambda_"$LAMBDA"_trajread.log for more information.\e[0m\n\n\n"
            exit 
          fi
          CRASH=1
          printf "\n\n\n \e[1;31m Job lambda_"$LAMBDA"_trajread.log appears to have crashed. Moving log file to lambda_$LAMBDA/lambda_"$LAMBDA"_trajread.log.crash and resubmitting\e[0m\n\n\n"
          LOG=lambda_$LAMBDA"_trajread.log"
          mv lambda_$LAMBDA/$LOG lambda_$LAMBDA/$LOG.crash
          resub
          continue
        else
          printf " Job lambda_"$LAMBDA"_trajread.inp finished successfully\n"
        fi
      fi
    done
    if [ $CRASH -eq 0 ]; then
      break
    fi
    monitor_job
done

printf "\n\n \e[1;32mAll jobs completed successfully\e[0m\n\n\n"

# If we reach here, all jobs should have finished successfully (with "NORMAL TERMINATION"), so we can
# check whether runs converged and subdivide lambda windows as necessary

[ -e delta_G_elec.dat ] && rm delta_G_elec.dat
for i in $(seq 1 $NWIN); do
  LAMBDA=`bc <<< "scale = 7; $LSTART + $D_LAMBDA * $i - $LAMBDA1"`
  dGi=`cat lambda_$LAMBDA/dG.dat`
  var=`cat lambda_$LAMBDA/var.dat`
  printf " Lambda $LAMBDA:  delta_G = $dGi  variance = $var\n"

  if [ $mtpl -eq 0 ] && [ $dcm2 -eq 0 ]; then
    if (( $(echo "$var > 0.5" |bc -l) )); then
      printf "\e[1;33m Warning: Variance $var for lambda=$LAMBDA is greater than 0.5. Reducing the lambda window size.\e[0m\n\n"
      mkdir -p not_converged
      mv lambda_$LAMBDA "not_converged/lambda_"$LAMBDA"_d"$D_LAMBDA
      TLAMBDA=`bc <<< "scale=7; $D_LAMBDA / 2.0"`
      TSTART_LAMBDA=`bc <<< "scale=7; $LAMBDA - $TLAMBDA"`
      TSTOP_LAMBDA=`bc <<< "scale=7; $LAMBDA + $TLAMBDA"`
      echo "Launching new run with LSTART=$TSTART_LAMBDA, LSTOP=$TSTOP_LAMBDA, delta Lambda=$TLAMBDA"
      echo "$SCRIPT_DIR/ti.sh $ARGS -lstart $TSTART_LAMBDA -lstop $TSTOP_LAMBDA -dlam $TLAMBDA"
      $SCRIPT_DIR/ti.sh $ARGS -lstart $TSTART_LAMBDA -lstop $TSTOP_LAMBDA -dlam $TLAMBDA -subdivide &
    fi
  fi
done

# Now keep waiting until all subdivided window jobs are finished
wait
#for i in `pgrep -P "$$"`; do
#  if ps -p $i > /dev/null; then
#    wait $i
#  fi
#done

# Gather final results
if [ $SUBDIV -eq "0" ]; then
  dG=0
  printf " Gathering free energies for each window from dG.dat files\n\n\n"
  for i in lambda_*/; do
    dGi=`cat $i/dG.dat`
    dG=$( echo $dG+$dGi | bc )
  done

  printf "\n\n \e[1;32mTotal delta_G = $dG kcal/mol\e[0m\n\n\n"
  echo $dG > delta_G.dat

  printf " Exiting, final delta_G written to `pwd`/delta_G.dat\n"
  printf " Suggestions / bug reports to Michael.Devereux@unibas.ch\n\n"
fi
