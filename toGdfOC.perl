use strict;
use warnings;

my %rank;
my %sz;
open A, "zcat $ARGV[0] |"; 
while (<A>){
  chop();
  my ($f, $t)=split(/;/, $_);
  $rank{$f}=$t;
  $sz{$t}{$f}++;
};

open A, "zcat $ARGV[1]|";
my $i=0;
my (%l, %nn);
while(<A>){
  chop();
  $l{$i} = $rank{$_} if defined $rank{$_};
  $nn{$i} = $_;
  $i++;
}
$i=0;
my %cl;
while (<STDIN>){
  chop();
  next if /^#/;
  my @x = split(/ /);
  for my $j (@x){
    $cl{$nn{$j}} = $i;
  }
  $i++;
}

$i=0;
my (%n, %e, %imap);
open A, "$ARGV[2]"; 
while(<A>){
  chop();
  my ($a, $b, $w) = split(/ /);
  if ($cl{$nn{$a}} ne $cl{$nn{$b}}){
    my $f = defined $rank{$nn{$a}} ? $rank{$nn{$a}} : $nn{$a};
    my $t = defined $rank{$nn{$b}} ? $rank{$nn{$b}} : $nn{$b};
    ($f,$t) = sort ($f, $t);
    if (!defined $imap{$f}) { $imap{$f} = $i; $i++}
    if (!defined $imap{$t}) { $imap{$t} = $i; $i++}
    $n{$f}++; $n{$t}++; 
    $e{"$imap{$f},$imap{$t}"} += $w if ($f ne $t);
    
  }
}
my @ns = keys %n;
print "nodedef> name INT, label VARCHAR,height DOUBLE\n"; 
for my $i (0..$#ns){
  my $no = $ns[$i];
  $no =~ s/["']/ /g;
  print "$imap{$ns[$i]},".'"'.$no.'"'.",".log($n{$ns[$i]})."\n";
} 
print "edgedef> node1,node2,weight DOUBLE,directed BOOLEAN\n"; 
for my $ed (keys %e){ 
  print "$ed,$e{$ed},false\n";
}
