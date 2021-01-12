#!/usr/bin/perl

use warnings;
use strict;
use File::Temp qw/ :POSIX /;


my $t0 = time();
print STDERR "starting at $t0\n";
open A, "|gzip>$ARGV[0].names";
open B, "|gzip>$ARGV[0].versions";
open C, "|gzip>$ARGV[0].weights";
my (%f2num);
my $n = 0;
my $i = 0;
my @degree;
my @degreeW;

while(<STDIN>){
  chop();
  my ($w, @vs) = split(/\;/, $_, -1);
  my @vs1 = (); 
  my $cl = "cl$n";# use to identify groups of projects sharing these commits: echo p | getValues -f p2c | cut -d\; -f2 | lsort 1G -u > $p.cs; 
  push @vs, $cl;
  for my $v1 (@vs){
    push @vs1, $v1;
    if (!defined $f2num{$v1}){
      $f2num{$v1} = $i+0;
      print A "$v1\n";
      $i++;
    }
  }
  my $v0 = shift @vs1;
  for my $v1 (@vs1){
    my $id0 = $f2num{$v0};
    my $id1 = $f2num{$v1};
    print B "$id0 $id1\n";
    $degreeW[$id0]+=$w;
    $degreeW[$id1]+=$w;
    $degree[$id0]++;
    $degree[$id1]++;
    print C "$w\n";
  } 
  $n ++;
  if (!($n%100000000)) { print STDERR "$n lines and $i nodes in $ARGV[0] done\n";}
}
open C, "|gzip>$ARGV[0].degree";
for my $no (0..($i-1)){
  print C "$degree[$no];$degreeW[$no]\n";
}
#print B "".($i-1)." ".($i-1)."\n";#ensure a complete list of vertices
my $t1 = time();
print STDERR "finished at $t1 over ".($t1-$t0)."\n"; 
