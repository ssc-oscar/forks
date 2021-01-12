#!/usr/bin/perl

use warnings;
use strict;
use File::Temp qw/ :POSIX /;

my $maxAPI = 0;
$maxAPI = $ARGV[1] if defined $ARGV[1];
my $maxP = 0;
$maxP = $ARGV[2] if defined $ARGV[2];

my $t0 = time();
print STDERR "starting at $t0\n";
open A, "|gzip>$ARGV[0].names";
open B, "|gzip>$ARGV[0].versions";
open C, "|gzip>$ARGV[0].degree";
my (%f2num);
my $n = 0;
my $i = 0;

my %deg;


sub getId {
  my $v1 = $_[0];
  if (!defined $f2num{$v1}){
    $f2num{$v1} = $i+0;
    print A "$v1\n";
    $i++;
  }
  return $f2num{$v1};
}

while(<STDIN>){
  chop();
  my ($p, @vs) = split(/\;/, $_, -1);
  my @vs1 = (); 
  $p = "PRJ_$p";
  next if $maxAPI > 0 && $#vs > $maxAPI;
  my $id0 = getId ($p);
  for my $v1 (@vs){
    $v1 =~ s/^\s*//; $v1 =~ s/\s*$//;
    $v1 =~ s/ as .*$//;
    $v1 =~ s/^\.*//;
    $v1 =~ s|^/*||;
    $v1 =~ s|File.expand_path\('||;
    $v1 =~ s|'||g;
    # get last component of the path, perhaps?
    push @vs1, getId($v1);
  }
  my $pid = $id0;
  for my $id1 (sort @vs1){ 
    if ($pid ne $id1){
      print B "$id0 $id1\n";
      $deg{$id0}++;
      $deg{$id1}++;
    }
    $pid = $id1;
    #$edges{$id1}{$id0}++;
  }
  $n++;
  if (!($n%100000000)) { print STDERR "$n lines and $i nodes in $ARGV[0] done\n";}
}

for my $no (0..($i-1)){
  my $dg = 0;
  $dg = $deg{$no} if defined $deg{$no}; 
  print C "$no;$dg\n";
}

my $t1 = time();
print STDERR "finished at $t1 over ".($t1-$t0).". Has $i nodes and $n edges\n"; 
