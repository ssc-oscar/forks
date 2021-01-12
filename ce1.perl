use warnings;
use strict;
use File::Temp qw/ :POSIX /;

my %mp;
my $i = 0;
my %c;

while(<STDIN>){
  chop();
  my ($a, $b) = split (/\;/, $_, -1);
  if (!defined $mp{$a}){
    $mp{$a} = $i;
    #print B "$a\n";
    $i++;
  } 
  $c{$b}{$mp{$a}}++;
}

while(my ($k, $v) = each %c){
  for my $i (keys %$v){
    print " $i";
  }
  print "\n";
}

