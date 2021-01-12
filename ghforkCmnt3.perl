use strict;
use warnings;
my (%r2cr, %cr2r);

open A, "zcat ghForkMap.gz|";
while (<A>){
  chop();
  my ($r, $f) = split (/;/);
  $r2cr{$r} = $f;
}

sub mCmp {
  return -1 if length($a) < length($b);
  if (length($a) == length($b)){
    return $a cmp $b;
  }
  return 1;
}

my %r2ccr;
my $pid = "";
my %tmp = ();

#my $i = 1;
#my %i2n;
#open A, "zcat $ARGV[0].names.gz|";
#while (<A>){
#  chop();
#  $i2n{$i} = $_;
#  $i++;
#}
#my %seg;
for my $i ("mmbrshipa"){
  $pid = "";
  %tmp = ();
  open A, "zcat $ARGV[0].$i|";
  while (<A>){
    chop();
    my ($p, $id) = split (/;/);
    #next if $i ne "1" && defined $seg{"1"}{$p};
    #if (!defined $i2n{$p}){
    #  print STDERR "no $p in $ARGV[0].$i\n";
    #  next;
    #}
    #$p=$i2n{$p};
    #$seg{$i}{$p}++;
    if ($pid ne $id && $pid ne ""){
      output();
    }
    $tmp{$p}++ if $p ne ""; 
    $pid = $id;
  }
  output();
}

sub output {
  my @pp = sort mCmp (keys %tmp);
  for my $p (@pp){
    print STDERR "$p;@pp\n" if $p eq "miranagha_js";
    $r2ccr{$p} = $pp[0];
  }
  %tmp = ();
}


while (<STDIN>){
  chop();
  my ($cr, $r) = split (/;/, $_, -1);
  my ($ccr, $cr1, $ccr1, $ccr2) = ("", "", "", "");
  if (defined $r2cr{$r}){
    $cr1 = $r2cr{$r} 
  }
  if (defined $r2ccr{$r}){
    $ccr = $r2ccr{$r}; # not in ghFork
    $ccr1 = $r2cr{$ccr} if defined $r2cr{$ccr}; # get parent if from ghFork
  }
  
  if ($ccr eq ""){
    if ($cr1 ne ""){
      $ccr = $r2ccr{$cr1} if defined $r2ccr{$cr1};
      $ccr2 = $r2cr{$ccr} if defined $r2cr{$ccr};
    }
  }
  print "$r;$cr;$cr1;$ccr;$ccr1;$ccr2\n";
}
  

