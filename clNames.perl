use strict;
use warnings;

my %rank;
open A, "zcat P2AFullS.nA2APw.2000.crank.map|"; 
while (<A>){
  chop();
  my ($f, $t)=split(/;/, $_);
  $rank{$f}=$t;
};
  
open A, "zcat P2AFullS.nA2AP.2000.names|";
my $i=0;
my %l;
my $marked = 0;
while(<A>){
  chop();
  if (defined $rank{$_}){
    $l{$i} = $_;
    $marked ++;
  }
  $i++;
};
print STDERR "marked=$marked\n";

$i=0;
my %cl;
my %cn;
$marked = 0;
while(<STDIN>){
  chop();
  my ($c,$w) = split(/;/);
  if (defined $l{$i}){
    $marked ++;
    $cn{$c}{$l{$i}}++;
  }
  $i++;
};
print STDERR "marked=$marked\n";

for my $c (keys %cn){
  my @as=keys %{$cn{$c}};
  #next if $#as < 10;
  for my $a (@as){
    print "$c;$a;$rank{$a}\n";
  }
}
