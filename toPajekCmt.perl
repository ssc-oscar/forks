use strict;
use warnings;

my %rank;
open A, "zcat P2AFullS.nA2Aw.2000.crank.map|"; 
while (<A>){
  chop();
  my ($f, $t)=split(/;/, $_);
  $rank{$f}=$t;
};
  
open A, "zcat P2AFullS.nA2A.2000.names|";
my $i=0;
my %l;
while(<A>){
  chop();
  $l{$i} = $rank{$_} if defined $rank{$_};
  $i++;
}; 
$i=0;
my %cl;
my %cn;
while(<STDIN>){
  chop();
  my ($c,$w) = split(/;/);
  $cl{$i} = $c;
  $cn{$c}{$i}++;
  $i++;
};
my (%n, %e);
$i=0;
open A, "P2AFullS.nA2A.2000.csv"; 
while(<A>){
  chop();
  my ($a, $b, $w) = split(/ /);
  my ($f,$t) = sort ($cl{$a}, $cl{$b});
  $n{$f}++;
  $n{$t}++; 
  $e{"$f,$t"} += $w if $f ne $t;
}
my @ns = keys %n;
my %no2i;
print "*Vertices ".($#ns+1)."\n"; 
for my $no (@ns){
  my @nns = keys %{$cn{$no}};
  my $cnn = "unknown";
  for my $na (@nns){ 
    if (defined $l{$na}) {
      $cnn = $l{$na};
      last;
    }
  }
  $cnn =~ s/["']/ /g;
  print "$i ".'"'."$no:$cnn".'"'." ".$n{$no}."\n";
  $no2i{$no}=$i;
  $i++;
} 
print "*Edges\n"; 
for my $ed (keys %e){
  my ($f, $t) = split (/,/,$ed,-1);
  print "$no2i{$f} $no2i{$t} $e{$ed}\n";
}
