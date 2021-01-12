use warnings;
use strict;

open A, "gunzip -c $ARGV[0].names|";
my (@num2f);
while(<A>){
	chop();
	push @num2f, $_;
}
close (A);

my $cn = "";
my %cluster = ();
my $ off = 0;
while (<STDIN>){
	chop();
  $cluster{$_}{$off}++;
  $off++;
}

while (my ($k, $v) = each %cluster){
	my %fs = ();
	for my $f (keys %{$v}){
		$fs{$num2f[$f]}++;
	}
	output (\%fs);
}
undef @num2f;

sub output {
	my $cl = $_[0];
	my @fs = sort { length($a) <=> length($b) } (keys %{$cl});
	for my $i (0 .. $#fs){
		print "$fs[$i]\;$fs[0]\n";
	}
}	
