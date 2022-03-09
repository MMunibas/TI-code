#!/usr/bin/perl

use strict;

my $nres=0;
my $nchg=0;
my $natm=0;
my $nframe=0;
my $line=1;
my $line2=1;
my $res=1;
my $flag=0;
my $str;

my $file=$ARGV[0];
open(INP,"<$ARGV[0]");

if(@ARGV+0 != 2){ die "usage: scale_dcm.pl <dcmfile> <lambda>\n"};

my $lambda=$ARGV[1];
my $sqrt_lam=sqrt($lambda);

while(<INP>){

  chomp;
  print "$_\n";
  my @a=split;

  if($line==1){
    $nres=$a[0];
  }

  if($line>2){

    for(my $n=0; $n<$nres; $n++){
      $str=readline(INP);
      print "$str";
      my @b=split ' ', $str;
      $nframe=$b[0];  # no. frames
      for(my $m=0; $m<$nframe; $m++){  # loop over frames for this residue
        $str=readline(INP);  # atom indices for frame
        print "$str";
        for(my $l=0; $l<3; $l++){ # loop over atoms
          $str=readline(INP);
          print "$str";
          my @c=split ' ', $str;
          $nchg=$c[0];
          for(my $k=0; $k<$nchg; $k++){ # loop over chgs
            $str=readline(INP);
            my @d=split ' ', $str;
            if(@d+0 < 4){ die "Could not parse coordinates, read $d[0] $d[1] $d[2],\nhas file syntax changed?\n"; }
            printf "     %12.10f      %12.10f      %12.10f       %13.10f\n",$d[0],$d[1],$d[2],$d[3]*$sqrt_lam;
          }
        }
      }

    }
    $str=readline(INP);  # blank line after each residue
    print "$str";
  }

  $line++;

}
