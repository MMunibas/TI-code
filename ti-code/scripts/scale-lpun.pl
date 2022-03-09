#!/usr/bin/perl

# M DEVEREUX, Script to write lpun file with scaled multipoles (factor sqrt lambda as
# calc is then lambda * V = (sqrt_lam * MTPL1) * (sqrt_lam * MTPL2) / R^n

use strict;

if(@ARGV+0 != 2){
  die "usage: scale-lpun.pl <lpun_file> <lambda>\n"
}

my $lambda=$ARGV[1];
my $sqrt_lam=sqrt($lambda);
my $rank=0;
my $flag=0;
my $l=0;

open(INP,"<$ARGV[0]");

while(<INP>){

  chomp;
  my $line=$_;
  my @a=split;

  if($flag==2){ # line containing multipoles
    for(my $n=0; $n<@a+0; $n++){ # scale multipoles
      printf "%13.10f ",$a[$n]*$sqrt_lam;
    }
    print "\n";
    $l++;
    if($l>$rank){ # end of multipole block for this atom
      $flag=0;
      $l=0;
      $rank=0;
    }
  }else{
    print "$line\n";
  }

  if($flag==1){ # line containing atom indices
    $flag=2;
  }

  if($a[5] eq "Rank"){
    $rank=$a[6];
    $flag=1;
  }

}
