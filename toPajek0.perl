use strict;
use warnings;

#*Vertices N
#1 "1"
#*Edges
#1 2 1.0

my %n;
my %e;
while(<STDIN>){
  chop ();
  my ($a, $t, $w) = split(/ /);
  ($a, $t) = sort ($a, $t);
  next if $a eq $t;
  $e{"$a $t"} += $w;
  $n{$a}++;
  $n{$t}++;
}

my @ns = keys %n;
my %ri;
print "*Vertices ".($#ns+1)."\n";
for my $i (1..($#ns+1)){
  print "$i ".'"'.$ns[$i-1].'"'."\n";
  $ri{$ns[$i-1]} = $i;
}
print "*Edges\n";
while (my ($k, $v) = each %e){
  my ($f, $t) = split(/ /, $k);
  print "$ri{$f} $ri{$t} $v\n";
}
