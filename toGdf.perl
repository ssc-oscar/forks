use strict;
use warnings;

my %rank;
open A, "zcat $ARGV[0]|"; 
while (<A>){
  chop();
  my ($f, $t)=split(/;/, $_);
  $rank{$f}=$t;
};
  
open A, "zcat $ARGV[1] |";
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
open A, "$ARGV[2]"; 
while(<A>){
  chop();
  my ($a, $b, $w) = split(/ /);
  my ($f,$t) = sort ($cl{$a}, $cl{$b});
  $n{$f}++;
  $n{$t}++; 
  $e{"$f,$t"} += $w if $f ne $t;
}
my @ns = keys %n;
print "nodedef> name VARCHAR,label VARCHAR,height DOUBLE\n"; 
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
  print "$no,".'"'.$cnn.'"'.",".log($n{$no})."\n";
} 
print "edgedef> node1,node2,weight DOUBLE,directed BOOLEAN\n"; 
for my $ed (keys %e){ 
  print "$ed,$e{$ed},false\n";
}
