#!/usr/bin/perl

use warnings;
use strict;
use File::Temp qw/ :POSIX /;


my $t0 = time();
print STDERR "starting at $t0\n";
open A, "|gzip>$ARGV[0].names";
my (%f2num);
my $n = 0;
my $i = 0;
my @degree;
my @degreeW;
my %edges;

my %lp;
open L, "zcat largePFS2b.gz|";
while(<L>){
  chop ();
  my ($p, $c) = split(/\;/, $_, -1);
  $lp{$p}++ if $c > 100000;
}

while(<STDIN>){
  chop();
  my ($w, @vs) = split(/\;/, $_, -1);
  for my $v1 (@vs){
    if (!defined $f2num{$v1}){
      $f2num{$v1} = $i+0;
      print A "$v1\n";
      $i++;
    }
  }
  my $v0 = shift @vs;
  my $id0 = $f2num{$v0};
  for my $v1 (@vs){
    my $id1 = $f2num{$v1};
    my $key = "$id0 $id1";
    $edges{$key} += $w+0;
    $degreeW[$id0] += $w;
    $degreeW[$id1] += $w;
    $degree[$id0]++;
    $degree[$id1]++;
  }
  $n ++;
  if (!($n%100000000)) { print STDERR "$n lines and $i nodes in $ARGV[0] done\n";}
}

open B, "|gzip>$ARGV[0].versions";
open C, "|gzip>$ARGV[0].weights";
$n = 0;
while (my ($ee, $w) = each %edges){
  print B "$ee\n";
  print C "$w\n";
  $n++;
}
open D, "|gzip>$ARGV[0].degree";
for my $no (0..($i-1)){
  my $a = 0;
  my $b = 0;
  ($a, $b) = ($degree[$no], $degreeW[$no]) if defined $degree[$no]; 
  print D "$a $b\n";
}
#print B "".($i-1)." ".($i-1)."\n";#ensure a complete list of vertices
my $t1 = time();
print STDERR "finished at $t1 over ".($t1-$t0).". Has $i nodes and $n edges\n"; 
