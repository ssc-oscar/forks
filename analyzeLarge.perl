use strict;
use warnings;

my %freq;
open A, "cat fAs|";
while (<A>){
  chop ();
  $freq{$_}++;
}


open A, "zcat /data/basemaps/gz/a2AFullS.s|";
while (<A>){
  chop();
  my ($a, $ra, $bad) = split (/;/);
  if (defined $freq{$ra}){
    print "$a;$ra\n";
  }
}
# zcat /data/basemaps/gz/a2AFullS.s| grep ';0$' |cut -d\; -f2 | uniq -c | sed 's|\s*||;s| |;|' | perl -ane 'chop();($n,$f)=split(/\;/);print "$f\n" if $n > 10' > fAs
