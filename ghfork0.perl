use strict;
use warnings;

my %r2cr;
open A, "zcat projectstr2.gz|";
while (<A>){
  chop();
  my ($r, $cd, $ud, $f, $del, $stars) = split (/;/);
  $r2cr{$r} = $r;
  $r2cr{$r} = $f if $f =~ /_/;  
}

my $pc = "";
my %tmp;
my %noForks;
while(<STDIN>){
  chop();
  my ($c, @x) = split(/;/); 
  if (defined  $r2cr{$x[0]}){
    $x[0] = $r2cr{$x[0]};
  }else{
	 $noForks{$x[0]}++;
  }
  if ($c ne $pc && $pc ne ""){ 
    if (scalar (keys %tmp)>1) { 
      for my $i (keys %tmp){ 
			print "$i\n";
      };
    } 
    %tmp=(); 
  }; 
  $pc = $c;
  my $s = join ';', ($c, @x); 
  $tmp{$s}++;
}

for my $i (sort keys %tmp){ 
  print "$i\n";
}
for my $p (keys %noForks){
  print STDERR "$p\n";
}
