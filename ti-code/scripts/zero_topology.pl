#!/usr/bin/perl

# script to take topology file and solute PDB as input and return the same topology with zeroed charges
# for solute residue only as output

use strict;
use File::Basename;

if(@ARGV+0 != 2){
  die "Usage: zero-topology.pl <topology-file> <solute PDB file>\n";
}
print " running zero-topology with Topology $ARGV[0] and PDB $ARGV[1]\n";
my $slu="NULL";
my $flag=0;

my $pdb=$ARGV[1];
open(PDB,"<$ARGV[1]");

while(<PDB>){

  chomp;
  my @a=split;

  if(lc $a[0] eq lc "ATOM"){
    if(lc $a[3] ne lc $slu && lc $slu ne lc "NULL") {
      die "zero_topology.pl: more than one residue type found in pdb file!\n(this should contain the solute molecule only...)\n";}
    $slu=$a[3];
    $flag=1;
  }

}

if($flag == 0){
  die "no residue found in PDB file $pdb\n";
}

open(INP,"<$ARGV[0]");
open(OUT,">zeroed.top");

my $flag=0;

while(<INP>){

  chomp;
  my @a=split;

  if(lc $a[0] eq "resi"){
    if(lc $slu eq lc $a[1]){
      $flag=1;
    }else{
      $flag=0;
    }
  }

  if($flag == 1 && lc $a[0] eq "atom"){
    printf OUT "ATOM   %4s %4s %7.4f\n",$a[1],$a[2],0.0;
  }else{
    print OUT "$_\n";
  }

}
