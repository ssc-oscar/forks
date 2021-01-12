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
my @tops = (0, 0, 0, 0, 0);
my @top = ("", "", "", "", "");
while (my ($k, $v) = each %sz){
  my $lg = scalar (keys %$v);
  for my $i (0..$#tops){
    if ($lg > $tops[$i] && ($i ==0 || $lg < $tops[$i-1])){
      $tops[$i] = $lg;
      $top[$i] = $k;
    }
  }
}
print STDERR "@top;@tops\n";

my $which = $ARGV[3]+0;
open A, "zcat $ARGV[1] |";
my $i=0;
my (%l, %nn);
while(<A>){
  chop();
  if (defined $rank{$_}){
    $l{$i} = $rank{$_};
    print STDERR "$i;$_\n" if $_ eq $top[$which];
  }
  $nn{$i} = $_;
  $i++;
}; 
$i=0;
my (%n, %e, %imap, %map);
open A, "$ARGV[2]"; 
while(<A>){
  chop();
  my ($a, $b, $w) = split(/ /);
  #print STDERR "$which;$top[$which];$b;$l{$b}\n" if defined $l{$b};
  if ((defined $l{$a} && $l{$a} eq $top[$which]) || (defined $l{$b} && $l{$b} eq $top[$which])){
    #print STDERR "here $top[$which] $a $b $l{$a} $l{$b}\n"; 
    my ($f,$t) = sort ($nn{$a}, $nn{$b});
    if (!defined $imap{$f}) { $imap{$f} = $i; $map{$i} = $f; $i++}
    if (!defined $imap{$t}) { $imap{$t} = $i; $map{$i} = $t; $i++}
    $n{$f}++; $n{$t}++; 
    $e{"$imap{$f},$imap{$t}"} += $w if $f ne $t;
  }
}
my @ns = keys %n;
print "nodedef> name INT, label VARCHAR,height DOUBLE\n"; 
for my $i (0..$#ns){
  my $no = $ns[$i];
  $no =~ s/["']/ /g;
  print "$imap{$ns[$i]},".'"'.$no.'"'.",".log($n{$ns[$i]})."\n";
} 
print "edgedef> node1,node2,weight DOUBLE,directed BOOLEAN,label VARCHAR\n"; 
for my $ed (keys %e){ 
  my ($a, $b) = split(/,/, $ed);
  print "$ed,$e{$ed},false,".'"'."$map{$a}-$map{$b}".'"'."\n";
}
