use strict;
use warnings;


my (%i2c, %c2i);
my $base = $ARGV[0];
my $sel = $ARGV[2] + 0;
open A, "zcat $base.clones|";
while (<A>){
  chop();
  my ($id, $cl) = split (/;/);
  next if $cl != $sel;
  $i2c{$id} = $cl;
  $c2i{$cl}{$id}++;
}

my (%p, %ip);
my $i = 0;
open A, "zcat $ARGV[1]|";
while (<A>){
  chop();
  $p{$i} = $_; 
  $i++;
}

open A, "$base.csv";
open B, ">$base.$sel.csv";
while (<A>){
  chop();
  my ($a, $b, $w) = split (/ /);
  next if !defined $i2c{$a};
  print B "$a $b $w\n";
  $ip{$p{$a}}{$a}++;
  $ip{$p{$a}}{$b}++;
}

for my $pp (keys %ip){
  my @is = keys %{$ip{$pp}};
  if ($#is > 0){
    print "".(shift @is);
    while ($#is >= 0){
      print " ".(shift @is);
    }
    print "\n";
  }
}



