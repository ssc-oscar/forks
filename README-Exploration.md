# Detecting git repository forks

We use WoC infrastructure to identify related repositories by commits that 
are shared among the repositories. Git allows arbitrary number of clones to any repository, so it is 
not always clear which sets of repositories represent "independently developed" projects.

It is almost impossible to create two identical commits unintentionally, threfore 
repositories that share a commit must be related. However, if we simply use transitive 
closure of that relationship, a very large fraction og git repositories get into a single cluster.

Why? Several reasons exist. First, some of the repositories contain commits from 
unrelated repositories by intention or accident. Nothing prevents one pushing commits from any repository 
to any other (provided write access). Some repositories appear to serve as a backup for many unrelated projects. 
Other projects intentionally import the full history of the components they use (but not develop) in there project. 
Such nodes connect many seemingly unrelated repositories into a giant megacluster.

One way to operationalize independent (unrelated) repositories would be to use community detection algorithms
on the bi-partite graph where nodes are commits and repositories links represent the fact that a commit belongs to a 
repository.  


The challange is the scale: WoC version Q we use here has almost 2B commits and over 100M projects and has
almost 100B links. Below we describe the process used to accomplish the task. 
Here are main steps:

    - identify GitHub forks from the ghTorrent or via GitHub API
	 - The following steps take time and should be run on 32 different servers in parallel
	 - for each forked repository identify ultimate parent (if a parent is a fork itself, find its parent, until it is not a fork): result ghForks.gz
	 - use WoC databases mapping commits to projects (c2pFullQ{0..31}.s) and replace each project name by the ultimate parent name (if the parent exists): result c2PFullQ{0..31}.s 
    - exclude duplicates from c2PFullQ{0..31}.s: result c2PFullQ{0..31}.ss
	 - convert flat commit;project1 file to one per commit file and ignore commits that do not span projects commit;project1;project2;..: result .p2p files
	 - order projects in each commit (so the same set of projects will always be represented by the same string): result .p2p.gz files
	 - count each set of projects in commits only once (multiple commits that reside in the same set of projects are collapsed into one): result .p2p.s files
	 - (back on one server) merge all 32 .p2p.s files into a single file and remove redundant commits: result c2PFull$ver.p2p.s
	 - encode multiedges edges in c2PFull$ver.p2p.s as integers in a bigraph: lineNumber -> each of the repos on that line
	 - split the edges into ten chunks 139638710 edges each
	 - use iGraph and run cluster_louvain on the resulting graph
	 - output the result 
	 - use the output to assign a cluster name (representing an independently developed project) to all repos in WoC (ultimateMap2.s)
	 - produce various statistics reported in the paper



```
#use ghForks.gz to map each repo to ultimate parent
for j in {0..31}
do zcat c2pFull$ver$j.s | perl ~/lookup/ghfork0.perl 2 | uniq | gzip > c2PFull$ver$j.s &
done
wait
# sort to remove dups
for j in {0..31}
do zcat c2PFull$ver$j.s | lsort 300G -k1,2 -u | gzip > c2PFull$ver$j.ss &
done
wait
# use commits to link repos (commit, list of repos)
for j in {0..31}
do zcat c2PFull$ver$j.ss | perl $HOME/lookup/connectExportPreNoExclude.perl | gzip > c2PFull$ver$j.p2p & 
done
wait
# get rid of commits and sort repos on one line
for j in {0..31}
do zcat c2PFull$ver$j.p2p | perl -ane 'chop();($c,@x) = split(/;/); print "".(join ";", (sort @x))."\n";' | gzip > c2PFull$ver$j.p2p.gz &
done
wait
#remove repeated sets of repos 
for j in {0..31}
do zcat c2PFull$ver$j.p2p.gz | lsort  ${maxM}M -t\| -u | gzip > c2PFull$ver$j.p2p.s &
done
wait
#combine data from all 32 databases
str="lsort $(($maxM*32))M -t\| -u --merge"
for j in {0..31}
do str="$str <(zcat  c2PFull$ver$j.p2p.s)"
done
eval $str | gzip > c2PFull$ver.p2p.s
###############################

#now create names and edges as numbers to reduce mem footprint
zcat c2PFull$ver.p2p.s| perl -e 'my $nn=0;$cl=0;while(<STDIN>){ chop(); (@x)=split(/;/); $n2i{"cl$cl"} = $nn; print STDERR "cl$cl\n"; $nn++; for $n (@x){ if (!defined $n2i{$n}){$n2i{$n}=$nn; print STDERR "$n\n"; $nn++;}} for $i (0..$#x){print $n2i{"cl$cl"}." $n2i{$x[$i]}\n"}; $cl++;}' 2> c2PFull$ver.names | gzip > c2PFull$ver.edges.gz
#split edges into several files
zcat c2PFull$ver.edges.gz | split -l 139638710 -da1 --filter='gzip > $FILE.gz' - ~/c2PFull$ver.edges.
gzip < c2PFull$ver.names > ~/c2PFull$ver.names.gz
##########################
# now run R 
jj="";
nn = scan(paste("c2PFullQ",jj,".names.gz", sep=""), what="");
n=length(nn);
g = make_empty_graph(n, directed = FALSE);
for (i in 0:9){
  x=scan(paste("c2PFullQ.edges.",i,".gz", sep=""), what=0);
  g <- add_edges (g, x+1);
}
g = set_vertex_attr(g, "label", value = nn);
aa1 = V(g)[V(g)$label == "18f5b0ceecb87delgh:e152g2014.11.21002dudadaejannco_addons-server"]   
ga = delete_vertices(g, as.numeric(aa1));
gal = cluster_louvain(ga);
galsz = -sort(-sizes(gal));
mrs = grep ("^cl[0-9]*$",V(ga)$label)
ordr = order(gal$membership[-mrs])
write(paste(V(ga)$label[-mrs][ordr],gal$membership[-mrs][ordr],sep=";"),file=paste("| gzip > c2PFullQ",jj,".mmbrshipa",sep=""),ncol=1);
############################
# use .mmbrshipa to group repos
zcat /data/basemaps/gz/pQnoAdj.map | perl ~/lookup/ghforkCmnt3.perl ~/c2PFullQ | gzip > pQnoAdj.map.a
zcat pQnoAdj.map.a |  perl -ane 'chop();@x=split(/;/); $res=$x[3]; $res = $x[4] if defined  $x[4]; $res = $x[5] if defined  $x[5]; $res=$x[0] if $res eq ""; print "$x[0];$res;$x[1];$x[2];$x[3];$x[4];$x[5]\n";' | lsort 300G -t\; -k1,1 | gzip > ultimateMap2.s

#produce map for WoC
zcat ultimateMap2.s | perl -ane 'chop();@x=split(/;/); print "$x[0];$x[1]\n" if $x[0] ne $x[1]' | gzip > /data/basemaps/gz/p2PQ.s

#report various stats
echo $(zcat ultimateMap2.s|wc -l) "&" $(zcat ultimateMap2.s|cut -d\; -f2 | lsort 10G | uniq -c | lsort 1G -rn > ultimateMap2.cnt; cat ultimateMap2.cnt | wc -l) "&"  $(head -1 ultimateMap2.cnt)


awk '{print $1 "&" $2 "\\\\"}' ultimateMap2.cnt |head -20| sed 's|_|/|'
354920&miranagha/js\\
333645&6101/-\\
241893&ykgm/R\\
211538&jkwonl/test\\
179315&aosp/oz\\
101988&UCF/50\\
94160&rdp/a\\
89602&maiyy/-\\
75231&DT/docs\\
72138&hdl/qfs\\
62686&EEELF/ll\\
62500&kclaw/P6\\
60379&F/hy\\
56126&cuibg/-\\


zcat ultimateMap2.s | perl ~/lookup/grepField.perl test.extra.k 1 | gzip > ultimateMap2.s.test
cat test.extra| cut -d\; -f1,2 | lsort 9G -t\; -k1,1 | uniq | join -t\; - <(zcat ultimateMap2.s.test) > test
cut -d\; -f2,3 test| lsort 4G | uniq -c| sed 's|;| |' | awk  '{print $2";"$3";"$1}' | sort -t\; -k1,1 > test.cnt
cut -d\; -f1  test.cnt | uniq -d | awk '{print "^"$1";"}' > test.cnt.bad

grep -vf test.cnt.bad test.cnt | cut -d\; -f3 | awk '{print i+=$1}' | tail -1
1620790
grep -f test.cnt.bad test.cnt | cut -d\; -f3 | awk '{print i+=$1}' | tail -1
32082

32082/(32082+1620790)
.01940985145855214438

grep '^octocat_Spoon-Knife;' test.cnt | cut -d\; -f3 | awk '{print i+=$1}' | tail -1                                                                                          
9245
da5:/data/play/forks>grep '^rdpeng_ProgrammingAssignment2;' test.cnt | cut -d\; -f3 | awk '{print i+=$1}' | tail -1                                                                                
2717
da5:/data/play/forks>grep '^rdpeng_ExData_Plotting1;' test.cnt | cut -d\; -f3 | awk '{print i+=$1}' | tail -1
1957
da5:/data/play/forks>grep '^spring-projects_spring-boot;' test.cnt | cut -d\; -f3 | awk '{print i+=$1}' | tail -1
1046


zcat diomidis.s | join -t\; - <(zcat ultimateMap2.s|cut -d\; -f1,2) | gzip > merged2.s

#split d
zcat merged2.s|cut -d\; -f2-3| lsort 5G -t\; -k1,1 | uniq | cut -d\; -f1 | uniq -d | wc
100300 
zcat merged2.s|cut -d\; -f2| lsort 5G -u | wc
2036117 
#split my
zcat merged2.s|cut -d\; -f2-3| lsort 5G -t\; -k2,2 | uniq | cut -d\; -f2 | uniq -d | wc
44357
zcat merged2.s|cut -d\; -f3| lsort 5G -u | wc
2124711

#For example, community detection splits kernel into many communities 
zcat merged2.s|cut -d\; -f2-3| lsort 5G -t\; -k1,2 | uniq -c | grep torvalds_linux | wc
629
#At the same time, the d algorithm, splits a lot of related repos:
zcat merged2.s|cut -d\; -f2-3| lsort 5G -t\; -k1,2 | uniq -c | grep aosp_oz | wc  
   1245    2490   57763

e.g:   
   1 04116_linux-artik-custom;aosp_oz
   2 1119553797_sprd-kernel-common;aosp_oz

echo 04116_linux-artik-custom| ~/lookup/getValues -f p2c | cut -d\; -f2 > a
echo 1119553797_sprd-kernel-common| ~/lookup/getValues -f p2c | cut -d\; -f2 > b
join a b | wc
 300390  300390 12315990
wc a b
  568329   568329 23301489 a
  310266   310266 12720906 b

zcat diomidis.s |wc
10649348
zcat merged2.s|wc
8157317
```

# newest attempt after removing 18f5b0ceecb87delgh:e152g2014.11.21002dudadaejannco_addons-server
da5
/data/play/forks
zcat c2pFull$ver.p2p.s1 | perl -e '$pstr="";$nn=0;while(<STDIN>){chop;(@x)=split(/;/);$n=pop @x;$str=join ";", @x; if ($pstr ne $str && $pstr ne ""){ print "$nn;$pstr\n";$pstr=$str;$nn=$n }else{ $pstr=$str;$nn+=$n}};print "$nn;$pstr\n";' | gzip > c2pFull$ver.np2p
zcat c2pFull$ver.np2p | perl connectExportVw.perl c2pFull$ver.np2p1
cd ~/src/networkit
zcat c2pFullR.np2p1.versions| ./cluster 100564768 | gzip > c2pFullR.np2p1.PLMPLP 
modularity=0.960559 nver=100564768 clusters=9424393 largest=334170
modularity=0.960559 nver=100564768 clusters=9424393 largest=334170
zcat c2pFullR.np2p1.PLMPLP|head -100564770 | grep -v '^[BE]' | gzip > c2pFull$ver.np2p1.PLM
zcat c2pFullR.np2p1.PLM | perl rank.perl c2pFullR.np2p1 | gzip > c2pFullR.np2p1.PLM.crank.map
zcat  c2pFullR.np2p1.PLM.crank.map | awk -F\; '{if ($2 != $1)print "$1";"$2 }' | gzip >  c2pFullR.np2p1.PLMmap.forks


zcat c2pFullR.np2p1.PLM | cut -d\; -f1 | lsort 5G | uniq -c | lsort 1G -n | grep -n . | tail
9424384: 105513 812
9424385: 117573 1349
9424386: 138399 817
9424387: 140674 14
9424388: 156518 35
9424389: 171665 39
9424390: 211539 1281
9424391: 223832 17
9424392: 297392 15
9424393: 334170 218

zcat c2pFullR.np2p1.PLM.crank.map|grep -v ';0\.0000' | cut -d\; -f2- | uniq -c | lsort 1g -rn | head                                                                                                                                     
 213254 0000m0000_ProgrammingAssignment2;0.000494
 210946 00-Ling_datasharing;0.000919
 180964 0-6141988_hello-world-ruby-bootcamp-prep-000;0.000150
 128730 0--_Spoon-Knife;0.000220
 128068 00101010_00101010.github.io;0.001743
 103033 0000-bigtree_0000-bigtree.github.com;0.005175
  99796 0-T-0_ps4-linux;1.000020
  78359 0000000111_bootstrap;0.197935
  67199 1nfrag_android_hardware_qcom_audio;0.000103
  65037 0000poiu_github-slideshow;0.000127

#now run updateCFBR.pbs
for i in {0..127}; do scp -p c2PFullR$i.s da3:/data/basemaps/gz/c2PPLMFullR$i.s; done
...

zcat b2PFull$ver.p2p.s | perl -e '$pstr="";$nn=0;while(<STDIN>){chop;(@x)=split(/;/);$n=pop @x;$str=join ";", @x; if ($pstr ne $str && $pstr ne ""){ print "$nn;$pstr\n";$pstr=$str;$nn=$n }else{ $pstr=$str;$nn+=$n}};print "$nn;$pstr\n";' | gzip > b2PFull$ver.np2p
zcat b2PFull$ver.np2p | perl connectExportVw.perl b2PFull$ver.np2p
zcat b2PFull$ver.np2p.names | wc -l
cp -p b2PFullR.np2p.versions ~/src/networkit
cd ~/src/networkit
zcat b2PFullR.np2p.versions| ./cluster 79922413 | gzip > b2PFullR.np2p.PLMPLP 
241404307 edges read
[INFO ]: DynKatz: Sort time: 17895
[INFO ]: parallel coarsening took (49482 ms) 
[WARN ]: move phase aborted after 32 iterations
[INFO ]: parallel coarsening took (18814 ms) 
[INFO ]: parallel coarsening took (2703 ms) 
[INFO ]: parallel coarsening took (566 ms) 
[INFO ]: parallel coarsening took (427 ms) 
[INFO ]: parallel coarsening took (329 ms) 
[INFO ]: parallel coarsening took (366 ms) 
modularity=0.821531 nver=79922413 clusters=1223192 largest=5382218
modularity=0.827331 nver=79922413 clusters=1223192 largest=5548264
#try igraph as well?
zcat b2PFullR.np2p.versions | ../igraph/doAll | gzip > b2PFullR.np2p.commu 
zcat b2PFullR.np2p.commu|cut -c1-100
Modularities:
0.668595 0.800589 0.814301 0.815037 0.815069 0.815071 0.815071
#not as good

zcat ~/src/networkit/b2PFullR.np2p.PLMPLP|head -159844830  | tail -79922415 | grep -v '^[BE]' | gzip > b2PFull$ver.np2p.PLM
zcat b2PFullR.np2p.PLM | perl rank.perl b2PFullR.np2p | gzip > b2PFullR.np2p.PLM.crank.map
zcat b2PFullR.np2p.PLM | cut -d\; -f1 | lsort 5G | uniq -c | lsort 1G -n | grep -n . | tail
1223183:2516886 42
1223184:2556598 103
1223185:2888495 25
1223186:3124864 30
1223187:3365050 11
1223188:4251268 5
1223189:4406530 3
1223190:4502529 17
1223191:4579998 9
1223192:5548264 8

zcat b2PFullR.np2p.PLM.crank.map|grep -v ';0\.0000' | cut -d\; -f2- | uniq -c | lsort 1g -rn | head                                                       
2621382 0000-bigtree_0000-bigtree.github.com;0.864405
2313949 00jy116_spring-boot;0.027938
2016015 0-T-0_ps4-linux;1.005039
1934212 0101adm_cdnjs;0.387243
1660874 Adude11_-tmp-100000-commit-2;0.055651
1540171 04160_react-udemy-demo;0.058654
1405558 07101994_SignalR;0.015550
1303195 00ERNEST00_FFmpeg;0.167898
1290586 1900_mezzanine;0.051218
1166245 0xDEAD_CarND-LaneLines-P1;0.047432



# new attempts

time zcat /home/audris/c2pFullR.np2p.versions| ./cluster 100564796 | gzip > c2pFullR.np2p.PLMPLP                                                                                                                 
10758586476 edges read
[INFO ]: parallel coarsening took (156827 ms) 
[INFO ]: parallel coarsening took (4290 ms) 
[INFO ]: parallel coarsening took (3081 ms) 
[INFO ]: parallel coarsening took (2826 ms) 
[INFO ]: parallel coarsening took (2223 ms) 
modularity=0.960559 nver=100564796 clusters=9396805 largest=466452
modularity=0.960559 nver=100564796 clusters=9396805 largest=468068

zcat /home/audris/src/networkit/c2pFullR.np2p.PLMPLP | head -100564798 | grep -v '^[BE]' | gzip > c2pFullR.np2p.PLMAgain
zcat /home/audris/src/networkit/c2pFullR.np2p.PLMPLP | head -201129596 | tail -$((201129596-100564798)) | grep -v '^[BE]' | gzip > c2pFullR.np2p.PLMPLP

time zcat /home/audris/c2pFullR.np2p.versions  | ./a.out 100564796 |gzip > c2pFullR.np2p.membershipLPD


zcat ~/src/networkit/c2pFullR.np2p.PLMPLP > PLMPLP

zcat c2pFullR.np2p.membershipPLM > PLM
zcat c2pFullR.np2p.membershipPLP > PLP

~/src/igraph/compare 100564796 ~/src/networkit/{PLM,PLP}
VI=0.633809 IVI=0.948729 SJ=13304144.000000 RAND=0.996065 ADJRAND=0.989564
da5:/data/play/forks>~/src/igraph/compare 100564796 ~/src/networkit/{PLM,LPD}
VI=0.353929 IVI=0.972687 SJ=10749880.000000 RAND=0.999964 ADJRAND=0.999905
~/src/igraph/compare 100564796 PLM PLM.doLPI1
VI=0.034496 IVI=0.997266 SJ=1167964.000000 RAND=0.999999 ADJRAND=0.999998
~/src/igraph/compare 100564796 PLM PLMPLP
VI=0.006928 IVI=0.999450 SJ=159805.000000 RAND=0.999999
ADJRAND=0.999996

zcat b2PFullR.np2p.versions | ../igraph/doAll | gzip > b2PFullR.np2p.commu &
zcat b2PFullR.np2p.commu | cut -c1-100
Modularities:
0.635425 0.80517 0.826052 0.826905 0.826941 0.826942 0.826942
...
Leiden found 27088865 clusters using CPM (resolution parameter 0.05 0.01), quality is 0.4075

zcat b2PFullR.np2p.commu|head -9 | tail -1 | perl -ane 's/ /\n/g;print' | gzip > b2PFullR.np2p.PLM

zcat b2PFullR.np2p.PLM | lsort 5G | uniq -c | lsort 1G -n | grep -n . | tail
1231820:2525982 23
1231821:2628933 42
1231822:2785750 9
1231823:3154318 19
1231824:3169093 0
1231825:4287565 7
1231826:4299340 3
1231827:4556907 40
1231828:4713080 6
1231829:4790255 5
zcat b2PFullR.np2p.versions | ./cluster 79777630 > b2PFullR.np2p.PLMPLP                                                                                                                                          
1 edges read
240364900 edges read
[WARN ]: move phase aborted after 32 iterations
[INFO ]: parallel coarsening took (67518 ms) 
[WARN ]: move phase aborted after 32 iterations
[INFO ]: parallel coarsening took (20288 ms) 
[INFO ]: parallel coarsening took (2976 ms) 
[INFO ]: parallel coarsening took (563 ms) 
[INFO ]: parallel coarsening took (365 ms) 
[INFO ]: parallel coarsening took (411 ms) 
[INFO ]: parallel coarsening took (363 ms) 
modularity=0.833420 nver=79777630 clusters=1219263 largest=4915006
modularity=0.839581 nver=79777630 clusters=1219263 largest=5012356

tail -79777631 /home/audris/src/networkit/b2PFullR.np2p.PLMPLP | grep -v '^[BE]' > b2PFullR.np2p.PLM-LP
zcat b2PFullR.np2p.PLM|grep -v '^$' > b2PFullR.np2p.PLM1
~/src/igraph/compare 79777630 b2PFullR.np2p.PLM1 b2PFullR.np2p.PLM-L
VI=1.472765 IVI=0.736541 SJ=22057609.000000 RAND=0.981712 ADJRAND=0.954397


zcat b2PFullR.np2p.versions | ./modularity 79777630 b2PFullR.np2p.PLM1 b2PFullR.np2p.PLM-LP
b2PFullR.np2p.PLM1 modularity=0.826942 nver=79777630 clusters=1231828 largest=4790255
b2PFullR.np2p.PLM-LP modularity=0.839581 nver=79777630 clusters=1219263 largest=5012356


cat b2PFullR.np2p.PLM-LP | perl ~/lookup/connectImportNoCl.perl b2PFullR.np2p | gzip > b2PFullR.np2p.PLM-LP.map

zcat b2PFullR.np2p.PLM.map|cut -d\; -f2 | lsort 50G | uniq -c | lsort 1G -rn > b2PFullR.np2p.PLM.map.cls

zcat b2PFullR.np2p.versions | ./centrality 79777630 b2PFullR.np2p.PLM-LP > b2PFullR.np2p.cranking

zcat /home/audris/c2pFullR.np2p.versions  | ./centrality 100564796 PLM > c2pFullR.np2p.cranking