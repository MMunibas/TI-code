#!/usr/bin/perl

use strict;

# parameter:

if(@ARGV+0 != 5){
  die "usage: lambda-pars.pl <lambda> <window_size> <nsteps> <nequil> <f/b>\n";
}

my $LAMBDA=$ARGV[0];

my $delta=$ARGV[1];

my $PSTART=$ARGV[3];

my $nsteps=$ARGV[2];

my $back=$ARGV[4];

my $minsteps=$PSTART+100000;

my $halfd=$delta/2.0;

my $LSTART = $LAMBDA-$halfd;
my $LSTOP = $LAMBDA+$halfd;
my $PSTOP = $nsteps;
my $LINCR=$delta;

if($back eq "b"){
  printf " LSTART %8.6f LAMBDA %8.6f  LSTOP %8.6f  PSTART %7i -
  PSTOP  %7i  PSLOW LINCR %8.6f -\n",$LSTOP,$LAMBDA,$LSTART,$PSTART,$PSTOP,$LINCR;
}else{
  printf " LSTART %8.6f LAMBDA %8.6f  LSTOP %8.6f  PSTART %7i -
  PSTOP  %7i  PSLOW LINCR %8.6f -\n",$LSTART,$LAMBDA,$LSTOP,$PSTART,$PSTOP,$LINCR;
}
