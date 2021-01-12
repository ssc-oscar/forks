use warnings;
use strict;
use File::Temp qw/ :POSIX /;

my %mp;
my $i = 0;
open A, "|gzip > $ARGV[0].versions";
open B, "|gzip > $ARGV[0].names";
open C, "zcat $ARGV[0]|";
while(<C>){
  chop();
  my ($a, $b) = split (/\;/, $_, -1);
  if (!defined $mp{$a}){
    $mp{$a} = $i;
    print B "$a\n";
    $i++;
  } 
}

open C, "zcat $ARGV[0]|";
while(<C>){
  chop();
  my ($a, $b) = split (/\;/, $_, -1);
  if (!defined $mp{$b}){
    print STDERR "$a $b $mp{$a} $mp{$b}\n";
    exit();
  }
  print A "$mp{$a} $mp{$b}\n";
}

