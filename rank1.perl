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

my (@rank);
while(<STDIN>){
  chop();
  my ($cl, $scr) = split(/;/);
  push @rank, $scr;
  $cluster{$cl}{$off}++;
  $off++;
}
close (A);

if ($off != $#num2f +1){
  print STDERR "wrong length: $off vs ".($#num2f +1)."\n";
  exit (-1);
}

print STDERR "length: $off\n";
while (my ($k, $v) = each %cluster){
	my %fs = ();
	for my $f (keys %{$v}){
    #print STDERR "$k;$f;$num2f[$f]\n"; 
		$fs{$num2f[$f]} = $rank[$f];
	}
  #print STDERR "$k\n";
	output (\%fs);
}

sub output {
	my $cl = $_[0];
  my %ps;
  my %ps1;
  for my $p (keys %{$cl}){
    if ($p !~ /^PRJ_/){
      $ps{$p}++;
    }
    $ps1{$p}++;
  }
	my @fs = sort { $cl->{$b} <=> $cl->{$a} } (keys %ps);
	my @fs1 = sort { $cl->{$b} <=> $cl->{$a} } (keys %ps1);
	for my $i (0 .. $#fs){
		print "$fs1[$i]\;$fs1[0];$cl->{$fs1[$i]};$cl->{$fs1[0]}\n";
	}
}	
