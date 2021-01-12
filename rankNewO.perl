use warnings;
use strict;

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
open A, "gunzip -c $ARGV[1]|";
while(<A>){
  chop();
  my ($cl, $scr) = split(/;/);
  push @rank, $scr;
}

my %map;
while(<STDIN>){
  chop($_);
  next if $_ =~ /^#/;
  $_ =~ s/\s+$//;
  $_ =~ s/^\s+//;
  $_ =~ s/\s+/ /g;
  my @n = split(/ /, $_, -1);
  if (!defined $rank[$n[0]] || $rank[$n[0]] eq "" || $n[0] > $#rank){
    print STDERR "$_;@n:\n";
    exit();
  }
  my $mrank = $rank[$n[0]];
  my $ctr = $n[0];
  for my $i  (@n){
    if ($rank[$i] >= $mrank && $num2f[$i] !~ /^cl[0-9]+$/){
      $mrank = $rank[$i];
      $ctr = $i;
    }
  }
  my $cln = $num2f[$ctr];
  for my $i (@n){
    my $cln = $num2f[$ctr];
    if ($num2f[$i] !~ /^cl[0-9]+$/){
      if (defined $map{$num2f[$i]}){
        if ($rank[$map{$num2f[$i]}] < $mrank){
          # replace with new center
          $map{$num2f[$i]} = $ctr;
          print "$num2f[$i];$cln\n";
        }else{
          # ignore new center 
          print "$num2f[$i];$num2f[$map{$num2f[$i]}]\n";
        }
        #print STDERR "$num2f[$i];$map{$num2f[$i]};$cln\n";
      }else{
        $map{$num2f[$i]} = $ctr;
        print "$num2f[$i];$cln\n";
      }
    }
  }
}
