use strict;
use warnings;

#*Vertices N
#1 "1"
#*Edges
#1 2 1.0
#P2AFullS.nA2A.2000
open A, "zcat $ARGV[0]|";
my $id = 0;
my @a;
while(<A>){
  chop ();
  push @a, $_;
}

my %n;
my %e;
open A, "$ARGV[1]";
while(<A>){
  chop ();
  my ($f, $t, $w) = split(/ /);
  ($f, $t) = sort ($f, $t);
  next if $f eq $t;
  $e{"$f $t"} += $w;
  $n{$f}++; $n{$t}++;
}

my @ns = keys %n;
my %ri;
print "*Vertices ".($#ns+1)."\n";
for my $i (1..($#ns+1)){
  my $nn = $a[$ns[$i-1]];
  $nn =~ s|\s+| |g;
  $nn =~ s|[^A-Za-z0-9]| |g;
  $nn =~ s|^ ||;
  $nn =~ s| $||;
  $nn = substr ($nn, 0, 100) if length($nn) > 100;
  print "$i ".'"'.$nn.'"'."\n";
  $ri{$ns[$i-1]} = $i;
}
print "*Edges\n";
while (my ($k, $v) = each %e){
  my ($f, $t) = split(/ /, $k);
  print "$ri{$f} $ri{$t} $v\n";
}
