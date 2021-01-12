use warnings;
use strict;

my $excludeMr = 1;
$excludeMr = $ARGV[1] if defined $ARGV[1];

open A, "gunzip -c $ARGV[0].names|";
my (@num2f);
while(<A>){
	chop();
  my $n = $_;
	push @num2f, $n;
}
close (A);

my $cn = "";
my %cluster = ();
my $off = 0;

my (@rank, %tr, %best);
while(<STDIN>){
  chop();
  my ($cl, $scr) = split(/;/);
  push @rank, $scr;
  my $cln = $num2f[$cl];
  if (! defined $tr{$cln}){
    $tr{$cln} = -1 
  }
  my $offn = $num2f[$off];
  # print STDERR "$cl $off $cln $offn\n" if $offn =~ /^cl[0-9]+$/;
  if ($offn  !~ /^cl[0-9]+$/ || ! $excludeMr){
    $cluster{$cln}{$offn}++;
    if ($tr{$cln} < $scr || ($tr{$cln} == $scr && defined $best{$cln} && length($best{$cln}) > length($offn))){
      $tr{$cln} = $scr;
      $best{$cln} = $offn;
    }
  }
  $off++;
}
close (A);

if ($off != $#num2f +1){
  print STDERR "wrong length: $off vs ".($#num2f +1)."\n";
  exit (-1);
}

print STDERR "length: $off\n";
while (my ($k, $v) = each %cluster){
  for my $p (keys %{$v}){
    print "$p;$best{$k};$tr{$k}\n";
  }
}
