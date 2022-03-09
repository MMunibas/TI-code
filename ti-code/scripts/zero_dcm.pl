#!/usr/bin/perl

# script reads a PDB file for solute in TI calculation to determine solute residue name,
# then zeroes charges and polarizabilities for the solute in a DCM paramter file


use strict;

my $slu="NULL"; # solute residue name
my $nres=0; # no. residues
my $nchg=0; # no. DCM charges
my $npol=0; # no. polarizable sites
my $natm=0; # no. atoms in residue
my $nframe=0; # no. DCM axis frames in residue
my $line=1;
my $line2=1;
my $res=1;
my $flag=0;
my $str;


my $pdb=$ARGV[1];
open(PDB,"<$ARGV[1]");

while(<PDB>){

  chomp;
  my @a=split;

  if(lc $a[0] eq lc "ATOM"){
    if(lc $a[3] ne lc $slu && lc $slu ne lc "NULL") {
      die "zero_dcm.pl: more than one residue type found in pdb file!\n(this should contain the solute molecule only...)\n";}
    $slu=$a[3]; 
    $flag=1;
  }

}

if($flag == 0){
  die "no residue found in PDB file $pdb\n";
}

my $file=$ARGV[0];
open(INP,"<$ARGV[0]");

if(@ARGV+0 != 2){ die "usage: zero_dcm.pl <dcmfile> <pdbfile>\n"};

while(<INP>){

  chomp;
  print "$_\n"; # reads no. of residues or blank space after residues
  my @a=split; 

  if($line==1){
    $nres=$a[0];
  }

  if($line>1){

    for(my $n=0; $n<$nres; $n++){
      $str=readline(INP); # residue name
      my @b=split ' ', $str;
      my $this=$b[0];
      print "$str";
      $str=readline(INP); # no. frames
      @b=split ' ', $str;
      $nframe=$b[0];  # no. frames
      print "$str";
      for(my $m=0; $m<$nframe; $m++){  # loop over frames for this residue
        $str=readline(INP);  # atom indices for frame
        print "$str";
        @b=split ' ', $str;
        my $na=3; # default no. atoms in frame
        if($b[2] == 0) { $na=2; } # diatomic
        if($b[1] == 0) { $na=1; } # monatomic
        for(my $l=0; $l<$na; $l++){ # loop over atoms
          $str=readline(INP);
          print "$str";
          my @c=split ' ', $str;
          $nchg=$c[0];
          $npol=$c[1];
          for(my $k=0; $k<$nchg; $k++){ # loop over chgs
            $str=readline(INP);
            my @d=split ' ', $str;
            if(@d+0 < 4){ die "Could not parse coordinates, read $d[0] $d[1] $d[2],\nhas file syntax changed?\n"; }
            if(lc $this eq lc $slu){
              printf "     %12.10f      %12.10f      %12.10f       0.0000000000\n",$d[0],$d[1],$d[2];
            }else{
              print "$str";
            }
          }
          for(my $k=0; $k<$npol; $k++){ # loop over polarizabilities
            $str=readline(INP);
            my @d=split ' ', $str;
            if(@d+0 < 1){ die "Could not parse polarizability, read $d[0],\nhas file syntax changed?\n"; }
            if(lc $this eq lc $slu){
              print "  0.0000000\n";
            }else{
              print "$str";
            }
          }
        }
      }
      $str=readline(INP); # blank line
      print "$str";
    }
  }
  $line++;

}
