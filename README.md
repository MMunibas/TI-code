# TI-code
Shell script to perform thermodynamic integration calculations with CHARMM
 Mike Devereux. Jan 2017.

 Description
 ===========

 Script to run manual CHARMM thermodynamic integration for solvation free energy of a solute
 molecule in a solvent box for a given parameter set. 

 TI calculations consist of 4 terms (VDW gas-phase, VDW solution-phase, ELEC gas-phase and
 ELEC solution-phase). In VDW gas-phase the electrostatics are turned off and slow-growth
 TI with CHARMM's "PERT" routines and PSSP soft-core option are used. The same approach is
 taken for VDW in solution but for the solute only (i.e. standard solvent parameters are used
 throughout).

 For electrostatic terms, two approaches are possible. For CHARMM or MDCM charges the
 standard CHARMM PERT routines can be used with slow growth. For multipoles or optionally for
 MDCM charges to allow more direct comparison with multipole results, delta G is evaluated by
 running a trajectory at a given Lambda, then reading in the CHARMM dcd file and evaluating
 the energy at Lambda=1 for each step and using the mean Lambda=1 energy as dG for that window
 (see Eq. 14, Bereau et al., JCTC, 2013, 9, 5450-5459). This means that TI takes place in 2
 stages per window, an initial simulation step run with scaled MTP parameters and a
 subsequent analysis step that reads in the trajectory and evaluates energies with unscaled
 (lamda=1) parameters.

 Usage
 =====

 The different options are described by calling the ti.sh script with no arguments. Examples
 are given at the end of this section.

 Before calling ti.sh, you need to create an equilibrated, solvated system with a restart file
 that can be used as a starting point for TI calculations. This is to avoid running long
 equilibration steps for every lambda window unnecessarily.

 You will also need to provide 2 PDB files, one for the solute and one for the solvent with 
 the solute removed (i.e. the solvent with a cavity where the solute should be). This way
 the 2 PDBs can be used to run both gas-phase and solvated systems, with some minor modification
 it would be possible to adapt the code to accept PDBs for the gas-phase molecule and whole
 solvated system including solute instead.

 Once These files have been prepared, you need to modify the user-defined variables in this
 script (at the start of the code below), and adapt the template files in the 
 scripts/templates folder for your own purposes. The ".sh" submission scripts should be
 adapted for the cluster environment you wish to use (job scheduler options, environment
 variables etc.). Note that "tmpl-mtpl.sh" handles the 2-step TI approach described above
 for multipolar simulations (without "PERT" in CHARMM), while "tmpl-vdw.sh" handles both VDW
 and electrostatic simulations that use CHARMM's PERT module with slow-growth.

 Once all files have been prepared, call the script once for each TI energy component (VDW,
 electrostatic, gas-phase and solution-phase). The script will launch runs for each TI
 window, check for crashes and resubmit, and upon job completion it will check whether any 
 windows need to be subdivided, and if so it will move the unconverged data to a new folder,
 split the TI window and submit 2 new jobs.

 For slow-growth with PERT the variance printed to the output can be used to check the need
 to subdivide. With the 2-step approach for MTPs the only rigorous approach is to re-run with
 smaller windows to check whether results remain the same for each windowsand have thus
 converged.

 For PERT slow-growth runs, it is best-practise (for puslished results) to repeat calculations
 in the reverse direction, a corresponding flag "-back" is provided for this purpose.

 Results are written to "deltaG.dat" at the end of each run, note that the sign of each term
 will depend on your definition of the lambda=0 and lambda=1 states

 Examples
 ========

 VDW gas-phase:
 ti.sh -vdw  $WORKDIR/phbr.par -top  $REFDIR/phbr.top -slu  $REFDIR/step5_slu.pdb -gas -res $EQUDIR/equil-gas.res

 DCM gas-phase:
 ti.sh -dcm $REFDIR/12charges.dcm -par $WORKDIR/phbr.par -top $REFDIR/phbr.top -slu $REFDIR/step5_slu.pdb -gas -res $EQUDIR/equil-gas.res

 VDW solution-phase:
 ti.sh -vdw $WORKDIR/phbr.par -top $REFDIR/phbr.top -slu $REFDIR/step5_slu.pdb -solv $REFDIR/step5_slv.pdb -res $EQUDIR/equil.res

 DCM solution-phase:
 ti.sh -dcm $REFDIR/12charges.dcm -par $WORKDIR/phbr.par -top $REFDIR/phbr.top -slu $REFDIR/step5_slu.pdb -solv $REFDIR/step5_slv.pdb -res $EQUDIR/equil.res
