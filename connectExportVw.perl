#!/usr/bin/perl

use warnings;
use strict;
use File::Temp qw/ :POSIX /;

my $badAuthHere =  <<'EOT';
devops <devops@gmail@com>
venudevops <devops@gmail.com>
9tdevops <devops@gmail.com>
devops <devops@gmail.com>
David <David@David-PC>
Tu Nombre <you@example.com>
vagrant <vagrant@scotchbox>
Instant Contiki <user@instant-contiki.(none)>
Mary <mary.example@rypress.com>
Your Name <youremail@domain.com>
Home <Home@Home-PC>
vagrant <vagrant@ironhack>
training_C2D.02.11 <training_C2D.02.11@accenture.com>
Author Name <email@address.com>
name <you@example.com>
Live session user <ubuntu@ubuntu.(none)>
DataWorks-CI <48289519+DataWorks-CI@users.noreply.github.com>
Travis CI <contact@loicortola.com>
user <user@user>
unknown <Daniel@.(none)>
Training_H2A.03.20 <Training_h2a.03.20@accenture.com>
user <user@user.com>
user1 <user1@user1-PC>
Travis-CI <travis@travis>
I <info@remcotolsma.nl>
ubuntu <ubuntu@ubuntu.(none)>
apple <apple@appledeiMac.local>
apple <apple@apples-MacBook-Pro.local>
Logon Aluno <logonrmlocal@fiap.com.br>
VSC <vscavu@microsoft.com>
Demo User <demouser@MacBook-Air.local>
vagrant <vagrant@vagrant-centos65.vagrantup.com>
iyuto <dev@code-check.io>
devops <devops@gmail.com>  
Azure Pages <donotreply@microsoft.com>
jserv <jserv@0xlab.org>
Automated Version Bump <gh-action-bump-version@users.noreply.github.com>
user <user@ubuntu.(none)>
Mac <mac@Macs-MacBook-Pro.local>
Utilisateur <user@debian.arcahe.ovh>
macbook <macbook@macbooks-MacBook-Pro.local>
Demo User <demouser@MacBook-Pro.local>
javascript-ru <47786167+jsru@users.noreply.github.com>
Apprentice <apprentice@devbootcamp.com>
Alex <alex@Alexs-MacBook-Pro.local>
OWASPFoundation <owasp.foundation@owasp.org>
unknown <drdynscript@gmail.com>
Vagrant Default User <vagrant@jessie.vagrantup.com>
Usuario <Usuario@Usuario-PC>
yourusername <your@email.com>
ASUS <ASUS@ASUS-PC>
John Doe <john@doe.org>
Plusieurs textes <email>
pkage <pkage@mit.edu>
Usuario Local <Usuario Local>
source_server <source_server@example.com>
NullDev <NL-Dev@yandex.com>
Alibaba OSS <opensource@alibaba-inc.com>
apple <apple@appledeMacBook-Pro.local>
oscreader <oscreader@OSC>
user <user@debian>
buddy <buddy@buddy.works>
pg <vagrant@packer-debian-7.4-amd64>
User <user@Users-MacBook-Pro.local>
lenovo <lenovo@lenovo-PC>
FIRST_NAME LAST_NAME <MY_NAME@example.com>
pi <pi@raspberrypi>
ubuntu <ubuntu@ubuntu>
mac <mac@macdeMacBook-Pro.local>
Xcode User <xcode@Fasts-MacBook-Pro-2.local>
apple <apple@apples-MacBook-Pro-2.local>
UCloud_Docs <49426641+UCloudDocs@users.noreply.github.com>
student <student@iMac-student.local>
user <user@mail.com>
rails-dev <rails-dev@railsdev-VirtualBox.(none)>
mac <mac@macdeMacBook-Air.local>
MacBook <macbook@MacBook-Pro-MacBook.local>
mac <mac@macs-MacBook-Air.local>
Name <user@example.com>
Server <server@server.com>
info <info@ralfw.de>
student <none@none.com>
unknown <Jason@.(none)>
A Happy User <auser@nothing.none>
codefactor-io <support@codefactor.io>
Travis CI <travis@Traviss-Mac-6.local>
Postman Integration <integration@postman.com>
user <user@DESKTOP-V882PTR>
A U Thor <author@example.com>
unknown <user@.(none)>
acer <acer@boss>
unknown <User@.(none)>
Bj <none@example.com>
- <anybody@ttuwiki.org>
sbx_user1051 <sbx_user1051@169.254.156.1>
Test test (testtest) <noreply@example.com>
bot50 <bot50@users.noreply.github.com>
gituser <git@gituser.com>
 <anonymous@github.com>
EOT

my %bad;
for my $nn (split(/\n/, $badAuthHere)){
  $nn =~ s/^\s*//;
  $nn =~ s/\s*$//;
  $bad{lc($nn)} = 1;
}



open A, "|gzip>$ARGV[0].names";
open B, "|gzip>$ARGV[0].versions";
open C, "|gzip>$ARGV[0].weights";
my (%f2num);
my $n = 0;
my $i = 0;
my @degree;
my @degreeW;

open Z, "zcat /data/basemaps/gz/a2AFullHS.s|";
while(<Z>){
  chop();
  my ($n, $nr, $b, $b1) = split (/;/);
  if ($b + $b1 > 0){
    $n =~ s/\s*$//;
    $nr =~ s/\s*$//;
    $n =~ s/^\s*//;
    $nr =~ s/^\s*//;
    $bad{lc($n)}++;
    $bad{lc($nr)}++;
  }
}

sub isBad {
  my $nn = $_[0];
  $nn =~ s/^\s*//;
  $nn =~ s/\s*$//;
  return 1 if defined $bad{lc($nn)};
  return 0;
}

$i=0;
while(<STDIN>){
  chop();
  my ($w, @vs) = split(/\;/, $_, -1);
  my @vs1 = (); 
  my $cl = "cl$n";# use to identify groups of projects sharing these commits: echo p | getValues -f p2c | cut -d\; -f2 | lsort 1G -u > $p.cs; 
  push @vs, $cl;
  for my $v1 (@vs){
    next if isBad ($v1);
    if (!defined $f2num{$v1}){
      $f2num{$v1} = $i+0;
      print A "$v1\n";
      $i++;
    }
    push @vs1, $v1;
  }
  next if $#vs1 <= 1;
  # my $v0 = shift @vs1a;
  my $v0 = pop @vs1; #to ensure project group is the other node
  for my $v1 (@vs1){
    my $id0 = $f2num{$v0};
    my $id1 = $f2num{$v1};
    print B "$id0 $id1\n";
    $degreeW[$id0]+=$w;
    $degreeW[$id1]+=$w;
    $degree[$id0]++;
    $degree[$id1]++;
    print C "$w\n";
  } 
  $n ++;
  if (!($n%100000000)) { print STDERR "$n lines and $i nodes in $ARGV[0] done\n";}
}
open C, "|gzip>$ARGV[0].degree";
for my $no (0..($i-1)){
  if (defined $degree[$no]){
    print C "$degree[$no];$degreeW[$no]\n";
  }else{
    print C "0;0\n";
  }
}
