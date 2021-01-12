use strict;
use warnings;

my %fp2p;
my %p2fp; 

while(<STDIN>){
  chop();
  s/^ *//;
  s/ /;/;
  my ($nc,$np,$fp,@x) = split(/;/);
  #print "$np;$nc;$fp\n";
  for my $p (@x){
    $fp2p{$fp}{$np}{$p}=$nc;
	 $p2fp{$p}{$fp}=$nc;
  }
}
for my $p (keys %p2fp){
	my @fps = sort keys %{$p2fp{$p}};
	print STDERR "$p;$#fps;@fps\n" if $#fps>1;
}

my %fp2pG;
for my $fp (keys %fp2p){
  for my $np (keys %{$fp2p{$fp}}){
    for my $p (keys %{$fp2p{$fp}{$np}}){
      if (defined $p2fp{$p}){
	     my @fps = sort keys %{$p2fp{$p}};
	     for my $fp1 (@fps){
			  $fp2pG{$fps[0]}{$np}{$p} = $fp2p{$fp}{$np}{$p};
        }
      }
    }
  }
}

my @groups = keys %fp2pG;
#print STDERR "groups=$#groups\n";   
  

for my $fp (keys %fp2pG){
  my @nps = sort {$a <=> $b} (keys %{$fp2pG{$fp}});
  #print STDERR "$fp;$nps[$#nps]\n";
  #print "$fp\;$#nps;@nps\n";
  for my $i (1..($#nps)){
    for my $p (keys %{$fp2pG{$fp}{$nps[$i]}}){
      print "+;$fp;$nps[$i];$nps[$i-1];$p;$fp2pG{$fp}{$nps[$i]}{$p}\n" if !defined $fp2pG{$fp}{$nps[$i-1]}{$p}
    }
  }
  for my $i (0..($#nps-1)){
    for my $p (keys %{$fp2pG{$fp}{$nps[$i+1]}}){
      print "-;$fp;$nps[$i];$nps[$i+1];$p;$fp2pG{$fp}{$nps[$i]}{$p}\n" if !defined $fp2pG{$fp}{$nps[$i+1]}{$p}
    }
  }
}
