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
#podman run -it -v /data/play/forks:/data:Z -v /fast/mod:/fast:Z -v /home/audris:/home/audris audris/jupyter-r bash

/data/play/forks
zcat c2pFullS.np2p.s | perl connectExportVw.perl c2pFullS.np2p
zcat c2pFullS.np2p.weights > ~/src/networkit/w
export LD_LIBRARY_PATH=/home/audris/lib64:/home/audris/lib:$(pwd)/build
zcat /data/c2pFullS.np2p.versions | ./clusterw 102087279 10861811185 -1 1| gzip > /data/c2pFullS.np2plw.PLM
#10861811185 edges read
#modularity=0.925592 nver=102087279 clusters=9687490 largest=500235

zcat /data/c2pFullS.np2p.versions | ./clusterw 102087279 10861811185 -1 0| gzip > /data/c2pFullS.np2pw.PLM
#10861811185 edges read
#modularity=0.848492 nver=102087279 clusters=9688257 largest=562499


zcat /data/c2pFullS.np2p.versions | ./cluster 102087279 | gzip > /data/c2pFullS.np2pu.PLM
#10861811185 edges read
modularity=0.958825 nver=102087279 clusters=9686736 largest=314557


zcat c2pFullS.np2pu.PLM | perl rankNew.perl c2pFullS.np2p 1 | gzip > c2pFullS.np2pu.PLM.crank.map
zcat  c2pFullS.np2pu.PLM.crank.map | awk -F\; '{if ($2 != $1)print $1";"$2 }' | gzip >  c2pFullS.np2pu.PLMmap.forks

zcat c2pFullS.np2pw.PLM | perl rankNew.perl c2pFullS.np2p 1 | gzip >  c2pFullS.np2pw.crank.map
zcat  c2pFullS.np2pw.crank.map | awk -F\; '{if ($2 != $1)print $1";"$2 }' | gzip >  c2pFullS.np2pw.PLMmap.forks

zcat c2pFullS.np2plw.PLM | perl rankNew.perl c2pFullS.np2p 1 | gzip >  c2pFullS.np2plw.crank.map
zcat  c2pFullS.np2plw.crank.map | awk -F\; '{if ($2 != $1)print $1";"$2 }' | gzip >  c2pFullS.np2plw.PLMmap.forks


# to identify commit sets responsible for clusters
zcat c2pFullS.np2pw.PLM | perl rankNew.perl c2pFullS.np2p 0 | gzip >  c2pFullS.np2pw.mrcrank.map

#compare weighted and unweighted
zcat /data/basemaps/gz/pS.s| awk -F\; '{print $1";"$1}' | perl ~/lookup/mp.perl 1 c2pFullS.np2pw.PLMmap.forks  | gzip > /data/basemaps/gz/p2PwS.s
zcat /data/basemaps/gz/pS.s| awk -F\; '{print $1";"$1}' | perl ~/lookup/mp.perl 1 c2pFullS.np2plw.PLMmap.forks  | gzip > /data/basemaps/gz/p2PlwS.s

zcat /data/basemaps/gz/pS.s| awk -F\; '{print $1";"$1}' | perl ~/lookup/mp.perl 1 c2pFullS.np2pu.PLMmap.forks  | gzip > /data/basemaps/gz/p2PuS.s
zcat /data/basemaps/gz/p2PuS.s | cut -d\; -f2 | lsort 50G -t\; -u | wc
#79636108

zcat /data/basemaps/gz/p2PwS.s | perl ce1.perl > w
zcat /data/basemaps/gz/p2PuS.s | perl ce1.perl > u

/home/audris/src/OvpNMI/bin/Release/onmi w u
# Average estimated membership in 'w': 1.098 (135162320 / 123154240)
# Average estimated membership in 'u': 1.098 (135162320 / 123153920)
# 0.995246

/home/audris/src/OvpNMI/bin/Release/onmi lw u
# Average estimated membership in 'lw': 1.098 (135162320 / 123154104)
# Average estimated membership in 'u': 1.098 (135162320 / 123153920)
# 0.996826



zcat c2pFullS.np2pu.PLMmap.forks | cut -d\; -f2 |lsort 50G -t\; | uniq -c | lsort 3G -rn | head -20
 213442 0000m0000_ProgrammingAssignment2
 210975 00-Ling_datasharing
 206132 0-Dark-Code-0_github-slideshow
 200054 00101010_00101010.github.io
 169983 0-6141988_ruby-lecture-reading-error-messages-bootcamp-prep-000
 128737 0--_Spoon-Knife
 101793 0-T-0_ps4-linux
  79155 0000000111_bootstrap
  65770 1nfrag_android_hardware_qcom_audio
  62058 0x131_syntastic
  58623 0-admin_tensorflow
  54850 0-deepthought-0_SmartThingsPublic
  53127 0-6141988_javascript-intro-to-functions-lab-bootcamp-prep-000
  51899 00jy116_spring-boot
  51076 001szymon_frontend-nanodegree-resume
  50266 00riddle00_st
  43803 0077cc_gitignore
  43046 0-kaladin_build
  42342 0000-bigtree_0000-bigtree.github.com
  41082 004307ec_caffe

zcat c2pFullS.np2pw.PLMmap.forks | cut -d\; -f2 |lsort 50G -t\; | uniq -c | lsort 3G -rn | head -20
 305253 0-6141988_hello-world-ruby-bootcamp-prep-000
 227117 0000m0000_ProgrammingAssignment2
 215205 0-Dark-Code-0_github-slideshow
 210897 00-Ling_datasharing
 159286 00-Evan_minimal-mistakes
 128062 0--_Spoon-Knife
 105460 0000-bigtree_0000-bigtree.github.com
  97767 0-T-0_ps4-linux
  80343 0000000111_bootstrap
  74090 0x131_syntastic
  58717 0-admin_tensorflow
  54897 0-deepthought-0_SmartThingsPublic
  53127 0-6141988_javascript-intro-to-functions-lab-bootcamp-prep-000
  52790 00jy116_spring-boot
  52363 007sair_gatsby-starter-netlify-cms
  51711 0-T-0_docker
  50266 00riddle00_st
  47919 148R1A0505_proprietary_vendor_motorola
  44584 001szymon_frontend-nanodegree-resume
  43803 0077cc_gitignore

#############################
#do blob sharing where blobs are for programming languages and of at least certain size (450 char long)
zcat P2PFullS.nb2b.2000.s | perl connectExportVw.perl P2PFullS.nb2b.2000 
paste -d\  <(zcat /data/P2PFullS.nb2b.2000.versions) <(zcat /data/P2PFullS.nb2b.2000.weights) | /home/audris/src/networkit/clusterw $(zcat /data/P2PFullS.nb2b.2000.names|wc -l) $(zcat /data/P2PFullS.nb2b.2000.weights|wc -l) | gzip > /data/P2PFullS.nb2bw.2000.PLM
#modularity=0.792458 nver=30128959 clusters=2805378 largest=6652122
zcat P2PFullS.nb2bw.2000.PLM | perl rankNew.perl P2PFullS.nb2b.2000 1 | gzip >  P2PFullS.nb2bw.2000.crank.map  
zcat P2PFullS.nb2b.2000.versions|ssh da3 "bin/connect" | gzip > P2PFullS.nb2b.2000.clones
zcat P2PFullS.nb2b.2000.clones | cut -d\; -f2 | lsort 30G | uniq -c | lsort 1G -rn | head -20
28196096 0
   4181 97898
   1992 78
    940 84
zcat P2PFullS.nb2b.2000.versions | paste  -d\   - <(zcat P2PFullS.nb2b.2000.weights) >  P2PFullS.nb2b.2000.csv 
perl getComponent.perl  P2PFullS.nb2b.2000 P2PFullS.nb2bw.2000.PLM 0 > P2PFullS.nb2bw.2000.0.forOSLO
time /home/audris/src/OSLOM2/oslom_undir -f /data/play/forks/P2PFullS.nb2b.2000.0.csv -hint  /data/play/forks/P2PFullS.nb2bw.2000.0.forOSLO -w -r 1 -hr 1

# do blob knowledge flow where blobs are for programming languages and of at least certain size (450 char long)
# and the next commiter is at least 30 days later
zcat A2AFullS.nfb.30.s | perl connectExportW.perl A2AFullS.nfb.30
zcat A2AFullS.nfb.30.versions | paste  -d\   - <(zcat A2AFullS.nfb.30.weights) >  A2AFullS.nfb.30.csv
cat A2AFullS.nfb.30.csv | /home/audris/src/networkit/clusterw $(zcat /data/A2AFullS.nfb.30.names|wc -l) $(cat /data/A2AFullS.nfb.30.csv|wc -l) | gzip > /data/A2AFullS.nfbw.30.PLM
cat A2AFullS.nfb.30.csv | /home/audris/src/networkit/clusterwd $(zcat /data/A2AFullS.nfb.30.names|wc -l) $(cat /data/A2AFullS.nfb.30.csv|wc -l) | gzip > /data/A2AFullS.nfbwd.30.PLM
#modularity=0.644839 nver=12128613 clusters=86041 largest=4167545
zcat A2AFullS.nfbw.30.PLM | perl rankNew.perl A2AFullS.nfb.30 1 | gzip >  A2AFullS.nfbw.30.crank.map  
zcat A2AFullS.nfb.30.versions|ssh da3 "bin/connect" | gzip > A2AFullS.nfb.30.clones
zcat A2AFullS.nfb.30.clones | cut -d\; -f2 | lsort 30G | uniq -c | lsort 1G -rn | head -20
11942602 0
    199 5677
    131 6422
     76 16373
     48 288
     36 33782
perl getComponent.perl A2AFullS.nfb.30 A2AFullS.nfbw.30.PLM 0 > A2AFullS.nfbw.30.0.forOSLO
time /home/audris/src/OSLOM2/oslom_undir -f /data/play/forks/A2AFullS.nfb.30.0.csv -hint  /data/play/forks/A2AFullS.nfbw.30.0.forOSLO -w -r 1 -hr 1


#############################  
export LD_LIBRARY_PATH=/home/audris/src/networkit/build
c=2000
#  cleanest version in 2000 after manual curation/email join/a2AFullS.s contains a;A;isBad
(zcat P2AFullS.A2A.$c.s | perl -ane 'chop();@x=split(/;/);$w=pop @x;print "$w;".(join ";", @x)."\n";' | perl connectExportVw.perl P2AFullS.nA2A.$c 2> possibleRoots; zcat /data/P2AFullS.nA2A.$c.versions | /home/audris/src/networkit/cluster $(zcat /data/P2AFullS.nA2A.$c.names|wc -l) | gzip > /data/P2AFullS.nA2Au.$c.PLM; paste -d\  <(zcat /data/P2AFullS.nA2A.$c.versions) <(zcat /data/P2AFullS.nA2A.$c.weights) | /home/audris/src/networkit/clusterw $(zcat /data/P2AFullS.nA2A.$c.names|wc -l) $(zcat /data/P2AFullS.nA2A.$c.weights|wc -l) | gzip > /data/P2AFullS.nA2Aw.$c.PLM;zcat P2AFullS.nA2Aw.$c.PLM | perl rankNew.perl P2AFullS.nA2A.$c 1 | gzip >  P2AFullS.nA2Aw.$c.crank.map)
# modularity=0.883691 nver=40196105 clusters=3555419 largest=3149423
# modularity=0.902983 nver=40196105 clusters=3565947 largest=3484825
zcat P2AFullS.nA2A.$c.versions | paste  -d\   - <(zcat P2AFullS.nA2A.$c.weights) >  P2AFullS.nA2A.$c.csv
cd /home/audris/src/verse/python
python3 ./convert.py --format weighted_edgelist /data/P2AFullS.nA2AP.2000.csv /data/P2AFullS.nA2AP.2000.bcsr
cd ../
./verse-weighted -input /data/P2AFullS.nA2AP.2000.bcsr -output /data/P2AFullS.nA2AP.2000.csv.bin
# layout/embed: use ICSE embeddings/create new
python3 ./convert.py --format weighted_edgelist /data/P2AFullS.nA2A.2000.csv /data/P2AFullS.nA2A.2000.bcsr  
./verse-weighted -input /data/P2AFullS.nA2A.2000.bcsr -output /data/P2AFullS.nA2A.2000.csv.bin


#prev version of connectExportVw.perl did not separate project group as the center among devs
zcat P2AFullS.A2A.$c.s | perl -ane 'chop();@x=split(/;/);$w=pop @x;print "$w;".(join ";", @x)."\n";' | perl connectExportVw.perl P2AFullS.nA2AP.$c 2> possibleRootsP
zcat P2AFullS.nA2AP.$c.versions | paste  -d\   - <(zcat P2AFullS.nA2AP.$c.weights) >  P2AFullS.nA2AP.$c.csv 
cat /data/P2AFullS.nA2AP.$c.csv | ./clusterw $(zcat /data/P2AFullS.nA2AP.$c.names|wc -l) $(zcat /data/P2AFullS.nA2AP.$c.weights|wc -l) | gzip > /data/P2AFullS.nA2APw.$c.PLM 
zcat /data/P2AFullS.nA2AP.$c.versions | /home/audris/src/networkit/cluster $(zcat /data/P2AFullS.nA2AP.$c.names|wc -l) | gzip > /data/P2AFullS.nA2APu.$c.PLM
#modularity=0.913354 nver=40187415 clusters=3581556 largest=2454912
#modularity=0.894664 nver=40187415 clusters=3571330 largest=2493313


zcat P2AFullS.nA2APw.$c.PLM | perl rankNew.perl P2AFullS.nA2AP.$c 1 | gzip >  P2AFullS.nA2APw.$c.crank.map
zcat P2AFullS.nA2APw.$c.crank.map | cut -d\; -f2 | /home/audris/bin/lsort 20G | uniq -c | /home/audris/bin/lsort 1G -rn | head -20
1396876 Tomster <tomster@emberjs.com>
1369706 onovy <novy@ondrej.org>
 361039 pascal <pascal@borreli.com>
 328509 Patrick McCarty <pnorcks@gmail.com>
 288493 beckermr <becker.mr@gmail.com>
 260467 mbehling <mb@mariobehling.de>
 257113 Anton <verybigbro@gmail.com>
 219089 Nitesh <nitesh.turaga@gmail.com>
 194444 James Ward <james@jamesward.com>
 182965 cottsay <logans@cottsay.net>
 167255 yours <leezche@gmail.com>
 162328 Steffen Forkmann <sforkmann@gmail.com>
  97857 a <KetanPadegaonkar@gmail.com>
  90945 Wendal <wendal1985@gmail.com>
  90851 mandark <A06438_P5.Training@accenture.com>
  85148 Markus <markus@markus.com>
  82704 drnic <drnicwilliams@gmail.com>
  76520 t-pot <imagire@gmail.com>
  76286 logonrm <logonrm@fiap.com.br>
  74077 alulab14 <alulab14@inf.pucp.edu.pe>

# run oslom
# first connectivity
zcat P2AFullS.nA2AP.$c.versions|ssh da3 "bin/connect" | gzip > P2AFullS.nA2AP.$c.clones

#take largest > 1000 chunks + run oslom on them
zcat P2AFullS.nA2AP.$c.clones | cut -d\; -f2 | lsort 30G | uniq -c | lsort 1G -rn | head -20
25461781 0
  10232 376697
   1268 38514
   1215 5594
   1199 115198
   1197 10362
   1090 104255
   1054 8143
   1011 32555

perl getComponent.perl P2AFullS.nA2AP.$c P2AFullS.nA2APw.$c.PLM 0 > forOSLOP.0
perl getComponent.perl P2AFullS.nA2AP.$c P2AFullS.nA2APw.$c.PLM 376697 > forOSLOP.376697

c=2000  
time /home/audris/src/OSLOM2/oslom_undir -f /data/play/forks/P2AFullS.nA2AP.$c.0.csv -hint forOSLOP.0 -w -r 1 -hr 1
time /home/audris/src/OSLOM2/oslom_undir -f /data/play/forks/P2AFullS.nA2AP.$c.376697.csv -hint forOSLOP.376697 -w -r 1 -hr 1

# layout/embed: use ICSE embeddings/create new
python3 ./convert.py --format weighted_edgelist /data/P2AFullS.nA2A.2000.0.csv /data/P2AFullS.nA2AP.2000.0.bcsr  
./verse-weighted -input /data/P2AFullS.nA2AP.2000.0.bcsr -output /data/P2AFullS.nA2AP.2000.csv.0.bin
#look at 376697?

#try second largest disconnected set
cat P2AFullS.nA2AP.$c.376697.csv_oslo_files/tp | perl rankNewO.perl P2AFullS.nA2AP.${c} P2AFullS.nA2APw.$c.PLM | gzip >  P2AFullS.nA2APOw.${c}.376697.crank.map
perl ./toGdfO.perl P2AFullS.nA2APOw.2000.376697.crank.map P2AFullS.nA2AP.2000.names P2AFullS.nA2AP.$c.376697.csv 0 > 376697.0.gdf

#calculate graphs for the largest set 
cat P2AFullS.nA2AP.$c.0.csv_oslo_files/tp | perl rankNewO.perl P2AFullS.nA2AP.${c} P2AFullS.nA2APw.$c.PLM | gzip >  P2AFullS.nA2APOw.${c}.0.crank.map
perl ./toGdfO.perl P2AFullS.nA2APOw.2000.0.crank.map P2AFullS.nA2AP.2000.names P2AFullS.nA2AP.$c.0.csv 0 > P0.0.gdf
cat P2AFullS.nA2AP.2000.0.csv_oslo_files/tp | perl ./toGdfOC.perl P2AFullS.nA2AOw.2000.0.crank.map P2AFullS.nA2A.2000.names P2AFullS.nA2A.2000.0.csv > Pcmnt.gdf 


#investigate Leuwen first
zcat P2AFullS.nA2APw.2000.PLM|perl clNames.perl > P2AFullS.nA2AP.2000.PLM.clnames
cat P2AFullS.nA2AP.2000.PLM.clnames | cut -d\; -f1|uniq -c|lsort 1G -rn |head -20 > top20
awk '{print $2}' top20 | while read i; do grep '^'$i';' P2AFullS.nA2AP.2000.PLM.clnames| cut -d\; -f2 | ~/lookup/getValues -f A2P |cut -d\; -f2 | lsort 20G -t\; -k1,1 | perl $HOME/bin/grepFieldv.perl PManyA2000.gz 1 |uniq -c |lsort 1G -rn |head -10 > $i.10; done

cat 8.10
   1615 0x38_storybook
   1604 08_fmdb
   1593 01binary_jest
   1591 0204111zhouyi_activeadmin
   1552 04th0809_website
   1488 0-T-0_bootstrap
   1475 0120139_woocommerce
   1438 012Tom_ionic
   1438 0-to-1_redux
   1428 0x73_MagicalRecord

cat 6.10
   1596 00jy116_istio
   1560 10171121_cinder
   1543 06094051_prometheus
   1493 01deyishu_ingress-nginx
   1486 AdamDang_community
   1486 0t3dWCE_boto
   1462 07101040115_celery
   1446 07mehmet0_origin
   1409 007bon_metasploit-framework
   1405 1024122298_openstack-manuals

cat 2.10
   1363 01generator_PrestaShop
   1312 0mars_SonataAdminBundle
   1229 013cbr_FOSUserBundle
   1206 tvdijen_krb5
   1176 3l73_TYPO3.CMS
   1076 00F100_composer
    985 02307_laradock
    981 1001Pharmacies_Sylius
    955 AashishYadavDev_devdocs
    950 0-Yorick-0_ZF2-S-Chazallet
    
cat 7.10
   1825 yarpiin_White-Wolf-SGS9-TW-Pie
   1824 yarpiin_White-Wolf-N9-TW-Pie-Upstream
   1816 linux-mailinglist-archives_linux-arm-kernel.lists.infradead.org.0
   1788 carlitros900_android_kernel_amazon_douglas
   1778 carlitros900_android_kernel_amazon_karnak
   1765 01org_pa-blink
   1751 bq_aquaris-M8
   1742 freeza-inc_bm-galaxy-note10-sd855-pie
   1724 1N4148_hostap
   1717 bq_aquaris-M10-FHD

head top20
1396876 8
1369706 6
 361039 2
 328509 7
 288493 27
 260467 18
 257113 14
 219089 76
 194444 1
 182965 91
   
#inspect the community
zcat P2AFullS.nA2APw.2000.PLM | perl toGdf.perl P2AFullS.nA2APw.2000.crank.map P2AFullS.nA2AP.2000.names P2AFullS.nA2AP.2000.0.csv > 0CMT.gdf
java -Djava.awt.headless=true -Xms10g -Xmx50g -cp forceatlas2.jar:gephi-toolkit-0.9.2-all.jar kco.forceatlas2.Main --input /data/0CMT.gdf --output /data/0CMT2.gdf --2d --format gdf --targetChangePerNode 0.2

for i in {0..3}
do perl ./toGdfO.perl P2AFullS.nA2APw.2000.crank.map P2AFullS.nA2AP.2000.names P2AFullS.nA2AP.$c.0.csv $i > 0.$i.gdf
   java -Djava.awt.headless=true -Xms10g -Xmx50g -cp forceatlas2.jar:gephi-toolkit-0.9.2-all.jar kco.forceatlas2.Main --input /data/0.$i.gdf --output /data/0.$i.20.gdf --2d --format gdf --targetChangePerNode 20
done
   
zcat P2AFullS.nA2APw.2000.PLM | perl toGdf.perl P2AFullS.nA2APw.2000.crank.map P2AFullS.nA2AP.2000.names P2AFullS.nA2AP.2000.0.csv > 0CMT.gdf


perl toPajek.perl P2AFullS.nA2AP.$c.names P2AFullS.nA2AP.2000.376697.csv > 376697.net
#TODO remove: add to findHomonyms.perl from connectExportVw.perl + modify {a2A,A2a}fullHS.s 
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
 <anonymous@github.com>
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

################
zcat P2AFullS.nA2Aw.$c.crank.map | cut -d\; -f2 | /home/audris/bin/lsort 20G | uniq -c | /home/audris/bin/lsort 1G -rn | head -20
2003031  <ash@kambanaria.org>
 455463 --global <tomas.vot@gmail.com>
 374218  <k-okada@jsk.t.u-tokyo.ac.jp>
 362344  <krlmlr+r@mailbox.org>
 315292 * <senirupasan@gmail.com>
 292874  <adam@adamralph.com>
 281212  <gautric@redhat.com>
 246615  <the.sk89q@gmail.com>
 209805  <aneesh.kumar@gmail.com>
 175921  <gmlwjd9405@naver.com>
 157564  <spencerlyon2@gmail.com>
 149159  <kefengong@gmail.com>
 146386 Aluno <Aluno@Digital.net>
 115740 4383 <herveberaud.pro@gmail.com>
 111998  <roberto.tyley@gmail.com>
 103329 aacid <aacid@283d02a7-25f6-0310-bc7c-ecb5cbfe19da>
 101419 0tofu <move@0tofu.com>
  94050 5-5-C <jmarshall.brock@gmail.com>
  89297  <kunal.sumbly@gmail.com>
  84167  <mabedin1@stuy.edu>

#plot layouts

#annotate layouts
# 

# now correlate with API/file extensions/age?
# get data from mongodb for each cluster/use oslom to break largest one first?
- use an algorith similar to AuthSumary.perl, but use group names

  

# run oslom
zcat P2AFullS.nA2A.$c.versions | paste  -d\   - <(zcat P2AFullS.nA2A.$c.weights) >  P2AFullS.nA2A.$c.csv
zcat P2AFullS.nA2Aw.$c.PLM|perl -e '$i=0;while(<STDIN>){chop();($c,$w)=split(/;/);$cl{$c}{$i}++;$i++};while(($k,$v)=each %cl){@a=keys %$v;print "@a\n";}' > forOSLO

cat P2AFullS.nA2A.2000.csv | toPajek.perl > P2AFullS.nA2A.2000.net
cat forOSLO |perl -e '@l=();while(<STDIN>){chop();@x=split(/ /);@l=@x if ($#l < $#x);};for $z (@l){$g{$z}++}; open A, "P2AFullS.nA2A.2000.csv"; while(<A>){chop();($a, $b, $w)=split(/ /);print "$a $b $w\n" if defined  $g{$a} && defined $g{$b};}' > P2AFullS.nA2A.2000.largest.csv                                 
cat forOSLO |perl -e '@l=();while(<STDIN>){chop();@x=split(/ /);@l=@x if ($#l < $#x && $#x<3000000);};print STDERR "$#l\n"; for $z (@l){$g{$z}++}; open A, "P2AFullS.nA2A.2000.csv"; while(<A>){chop();($a, $b, $w)=split(/ /);print "$a $b $w\n" if defined  $g{$a} && defined $g{$b};}' > P2AFullS.nA2A.2000.scndlargest.csv
cat forOSLO |perl -e '@l=();while(<STDIN>){chop();@x=split(/ /);@l=@x if ($#l < $#x && $#x<781387);};print STDERR "$#l\n"; for $z (@l){$g{$z}++}; open A, "P2AFullS.nA2A.2000.csv"; while(<A>){chop();($a, $b, $w)=split(/ /);print "$a $b $w\n" if defined  $g{$a} && defined $g{$b};}' > P2AFullS.nA2A.2000.thrdlargest.csv
611189


#get cluster members (based on Leuwen)
zcat P2AFullS.nA2Aw.2000.PLM|perl clNames.perl > P2AFullS.nA2A.2000.PLM.clnames

#Look at top 20
cat P2AFullS.nA2A.2000.PLM.clnames | cut -d\; -f1|uniq -c|lsort 1G -rn |head -20 > top20
awk '{print $2}' top20 | while read i; do grep '^'$i';' P2AFullS.nA2A.2000.PLM.clnames| cut -d\; -f2 | ~/lookup/getValues -f A2P |cut -d\; -f2 | lsort 20G -t\; -k1,1 | perl $HOME/bin/grepFieldv.perl PManyA2000.gz 1 |uniq -c |lsort 1G -rn |head -10 > $i.10; done

cat 1.10  QT/KDE <ash@kambanaria.org>
   1785 Adikso_nautilus-thumbnail-border
   1703 33905410_qtbase
   1668 08_fmdb
   1631 04th0809_website
   1580 0204111zhouyi_activeadmin
   1573 Dev-Linux_kde-workspace
   1566 01binary_jest
   1541 0dysseas_bedrock
   1504 0x38_storybook
   1498 0x24bin_pyenv

cat 2.10 e-comerce: sonata-project.org/ --global <tomas.vot@gmail.com>
   1405 01generator_PrestaShop
   1318 0mars_SonataAdminBundle
   1263 0ooo_yii
   1252 013cbr_FOSUserBundle
   1192 tvdijen_krb5
   1177 3l73_TYPO3.CMS
   1118 00F100_composer
   1049 115zed81_styles
   1001 2or3_core
    992 1001Pharmacies_Sylius

cat 3.10 robotics <k-okada@jsk.t.u-tokyo.ac.jp>
   1546 130s_rosdistro
   1293 123malin_Paddle
   1198 0123Andrew_gym
   1136 23119841_darknet
   1083 0shine0_OpenNMT-py
   1074 0i0_Mask_RCNN
   1073 0Crap_betaflight
   1000 1002victor_qgroundcontrol
    918 00light00_AirSim
    784 0dysseas_fastai

cat 16.10 bioinformatics conda <krlmlr+r@mailbox.org>
   1592 000CU000_bioconda-recipes
   1276 malcolmbarrett_epibot
   1015 llrs_STATegRa
    999 nturaga_manifest
    963 evince-doc_evince-doc.github.com
    840 4curiosity_dplyr
    739 cderv_cransays
    705 16hoppe_galaxy
    669 0x0L_dask
    532 0000m0000_courses

cat 8.10 starting to program and want to participate in hacktoberfest <senirupasan@gmail.com>
   1391 0limpi0_Hello-Hacktober
   1251 cipta-media_website
   1235 1Nmokhele_Rainbow-Poem
   1182 01samiksha_Data-Structures-And-Algorithms-Hacktoberfest18
   1174 LoserCringe_github-ok
   1139 00x5a4fcf_python-random-quote
   1095 7staff_hacktoberfest2018
   1080 damionx8_Hackoberfest
   1070 1Oozaru1_al-go-rithms
   1027 AHKol_jsonmc

cat 26.10 PowerShell cmdlets for developers and administrators <adam@adamralph.com>
   1467 0xc0re_azure-powershell
   1165 AaronBertrand_sql-docs
    914 AaDake_windows-itpro-docs
    909 0Neji_Umbraco-CMS
    903 07101994_SignalR
    839 0-wiz-0_coreclr
    814 AbhaPatankar_vsts-docs
    789 80er_vsts-tasks
    722 AkJo_microsoft-graph-docs
    712 AaronZhangL_azure-cli

cat 53.10 <mabedin1@stuy.edu>
   1231 0Petya_Cataclysm-DDA Game
    387 ef23_learnddit corse prj /r/IWantToLearn to retrieve relevant comments related to the search query
    367 1Hyena_crawl
    365 MGRiv_euler
    326 mks65_list
    300 Jing907_PCA
    294 nyongtory818_Proj1
    294 JG1990_assignment7
    294 JG1990_assignment5
    293 JG1990_assignment6

# community of clusters induced by developer participation in projects    
zcat P2AFullS.nA2Aw.2000.PLM|perl toGdf.perl P2AFullS.nA2Aw.2000.crank.map P2AFullS.nA2A.2000.names P2AFullS.nA2A.2000.csv > P2AFullS.nA2A.2000CMT.gdf

# Use Lewen as input for much more accurate OSLOM
zcat P2AFullS.nA2Aw.$c.PLM|perl -e '$i=0;while(<STDIN>){chop();($c,$w)=split(/;/);$cl{$c}{$i}++;$i++};while(($k,$v)=each %cl){@a=keys %$v;print "@a\n";}' > forOSLO
time /home/audris/src/OSLOM2/oslom_undir -f /data/play/forks/P2AFullS.nA2A.$c.csv -hint forOSLO -w -r 1 -hr 1
***************************************************************************
CHECK UNIONS AND SIMILAR MODULES DONE
******** module_collection ******** 31197 modules. writing... 
DONE   ****************************
pruning all the modules collected. Partitions found: 1
getting partition from tp-file: /data/play/forks/P2AFullS.nA2A.2000.csv_oslo_files/partitions_level_1
31197 groups found
31197 bss found
checking homeless nodes
writing final solution in file /data/play/forks/P2AFullS.nA2A.2000.csv_oslo_files/short_tp1
******** module_collection ******** 31206 modules. writing... 
DONE   ****************************
hierarchies done ********* 
real    17037m24.576s
user    16539m4.409s
sys     441m58.481s

c=2000
#translate from oslom to author IDs
cat P2AFullS.nA2A.${c}.csv_oslo_files/tp | perl rankNewO.perl P2AFullS.nA2A.${c} P2AFullS.nA2Aw.$c.PLM | gzip >  P2AFullS.nA2AOw.${c}.crank.map
# aggregate individual stats into the group stats
for i in {0..31}; do zcat A2summFullS$i.s; done |perl ~/lookup/AuthToAgg.perl P2AFullS.nA2AOw.2000.crank.map | gzip > P2AFullS.nA2AOw.2000.crank.map.stats
#create csv file to explore via R
zcat P2AFullS.nA2AOw.2000.crank.map.stats|perl -e 'while(<STDIN>){chop();($i,@x)=split(/;/);print "$i";for $k (1..8){$v=shift @x;print ";$v";}for $k (0..$#x){($a, $b)=split(/=/,$x[$i]);$s{$a}=$b;};for $k ("Ruby","Ada","Perl","Clojure","Rust","Go","Kotlin","Erlang","Sql","Julia","OCaml","Java"
,"JavaScript","Scala","Lua","Cobol","TypeScript","Fortran","Python","fml","other","PHP","Dart","R","Basic","Lisp","C/C++","Swift"){print ";$s{$k}";};print "\n"}' > P2AFullS.nA2AOw.2000.crank.map.stats.csv
sed -i 's|^\s*||;s|\r||' P2AFullS.nA2AOw.2000.crank.map.stats.csv
R --no-save
library("data.table")
v=read.table("P2AFullS.nA2AOw.2000.crank.map.stats.csv", sep=";",quote="",comment.char="",header=F)
names(v)=c("a","nA","nc","nf","nfb","np","fc","lc","na","Ruby","Ada","Perl","Clojure","Rust","Go","Kotlin","Erlang","Sql","Julia","OCaml","Java"
,"JavaScript","Scala","Lua","Cobol","TypeScript","Fortran","Python","fml","other","PHP","Dart","R","Basic","Lisp","C/C++","Swift")
sel=v$nA>3&v$nc>200&v$lc-v$fc>3600*24*30;
v1=v[sel,]

quantile(v1$lc/3600/24/365.25+1970)
      0%      25%      50%      75%     100% 
1995.579 2020.409 2020.611 2020.668 2262.271 
 quantile(v1$nA)
   0%   25%   50%   75%  100% 
    4     8    14    27 27474 
> quantile(v1$nc)
      0%      25%      50%      75%     100% 
     201      401      785     1835 12950677 
> quantile(v1$np)
     0%     25%     50%     75%    100% 
      5      39      69     135 1066524 
quantile(v1$fc/3600/24/365.25+1970)
      0%      25%      50%      75%     100% 
1973.169 2012.077 2014.952 2017.200 2020.589 

quantile(v1$lc/3600/24/365.25+1970)
      0%      25%      50%      75%     100% 
1995.579 2020.409 2020.611 2020.668 2262.271 
v1$tot = apply(v1[,c(10:29,31:37)],1,sum, na.rm=T)
#what labguage takes over 30% and 50% of all language-specific files in the community (language-purity)
for (i in c(10:29,31:37)) print (c(names(v1)[i],round(c(sum(v1[,i]/v1[,"tot"]>.3, na.rm=T),sum(v1[,i]/v1[,"tot"]>.5, na.rm=T))/sum(!is.na(v1[,i])),3)));
[1] "Ruby"    "0.157" "0.101"
[1] "Ada"     "0.06"  "0.043"
[1] "Perl"    "0.026" "0.015"
[1] "Clojure" "0.012"   "0.004"  
[1] "Rust"    "0.012" "0.005"
[1] "Go"      "0.027" "0.013"
[1] "Kotlin"  "0.024"  "0.007" 
[1] "Erlang"  "0.014"  "0.004" 
[1] "Sql"     "0.004" "0.002"
[1] "Julia"   "0.001" "0"    
[1] "OCaml"   "0"     "0"    
[1] "Java"    "0.003" "0.001"
[1] "JavaScript" "0.005"      "0.003"     
[1] "Scala"   "0.046" "0.019"
[1] "Lua"     "0.04" "0.02"
[1] "Cobol"   "0.003" "0.001"
[1] "TypeScript" "0.288"      "0.2"       
[1] "Fortran" "0.029"   "0.015"  
[1] "Python"  "0.063"  "0.035" 
[1] "fml"     "0.012" "0.004"
[1] "PHP"     "0.1"   "0.055"
[1] "Dart"    "0"    "0"   
[1] "R" "0"   "0"
[1] "Basic"   "0"     "0"    
[1] "Lisp"    "0.002" "0"    
[1] "C/C++"   "0.003" "0.001"
[1] "Swift"   "0.002" "0.001"
#almost all large clusters have all languages!
#but only two language clusters a pure:  "Ruby" TypeScrip PHP (>10% have > 30% of files in that lang)
for (i in c(10:29,31:37)) print (c(names(v1)[i],round(c(sum(v1[,i]/v1[,"tot"]>.3, na.rm=T),sum(v1[,i]/v1[,"tot"]>.5, na.rm=T))/dim(v1)[1],3)));
[1] "Ruby"    "0.157" "0.101"
[1] "Ada"     "0.06"  "0.043"
[1] "Perl"    "0.026" "0.015"
[1] "Clojure" "0.012" "0.004"  
[1] "Rust"    "0.012" "0.005"
[1] "Go"      "0.027" "0.013"
[1] "Kotlin"  "0.024" "0.007" 
[1] "Erlang"  "0.014" "0.004" 
[1] "Sql"     "0.004" "0.002"
[1] "Julia"   "0.001" "0"    
[1] "OCaml"   "0"     "0"    
[1] "Java"    "0.003" "0.001"
[1] "JavaScript" "0.005"      "0.003"     
[1] "Scala"   "0.046" "0.019"
[1] "Lua"     "0.04"  "0.02"
[1] "Cobol"   "0.003" "0.001"
[1] "TypeScript" "0.288"      "0.2"       
[1] "Fortran" "0.029" "0.015"  
[1] "Python"  "0.063" "0.035" 
[1] "fml"     "0.012" "0.004"
[1] "PHP"     "0.1"   "0.055"
[1] "Dart"    "0"     "0"   
[1] "R"       "0"     "0"
[1] "Basic"   "0"     "0"    
[1] "Lisp"    "0.002" "0"    
[1] "C/C++"   "0.003" "0.001"
[1] "Swift"   "0.002" "0.001"
 		     
#create input for gephi layouts      
perl ./toGdfO.perl P2AFullS.nA2AOw.2000.crank.map P2AFullS.nA2A.2000.names P2AFullS.nA2A.2000.csv 0 > lgst.gdf
perl ./toGdfO.perl P2AFullS.nA2AOw.2000.crank.map P2AFullS.nA2A.2000.names P2AFullS.nA2A.2000.csv 1 > scnd.gdf
perl ./toGdfO.perl P2AFullS.nA2AOw.2000.crank.map P2AFullS.nA2A.2000.names P2AFullS.nA2A.2000.csv 2 > thrd.gdf
perl ./toGdfO.perl P2AFullS.nA2AOw.2000.crank.map P2AFullS.nA2A.2000.names P2AFullS.nA2A.2000.csv 3 > frth.gdf
perl ./toGdfO.perl P2AFullS.nA2AOw.2000.crank.map P2AFullS.nA2A.2000.names P2AFullS.nA2A.2000.csv 4 > ffth.gdf
cat P2AFullS.nA2A.2000.csv_oslo_files/tp | perl ./toGdfOC.perl P2AFullS.nA2AOw.2000.crank.map P2AFullS.nA2A.2000.names P2AFullS.nA2A.2000.csv > cmnt.gdf 

#lets see what are projects most common in each group:
zcat P2AFullS.nA2AOw.${c}.crank.map | lsort 30G -t\; -k2,2 | cut -d\; -f2 | uniq -c |lsort 1G -rn |head -20 > top20O
sed 's|^\s*[0-9]* ||' top20O | while IFS=\;  read i; do grep "$i" <(zcat P2AFullS.nA2AOw.$c.crank.map) | cut -d\; -f1 | ~/lookup/getValues -f A2P |cut -d\; -f2 | lsort 20G -t\; -k1,1 | perl $HOME/bin/grepFieldv.perl PManyA2000.gz 1 |uniq -c |lsort 1G -rn |head -10|awk '{print "'"$i"';"$0}' ; done

 <conrad.irwin@gmail.com>;    577 05BIT008_devise  flexible authentication solution for Rails based on Warden.
 <conrad.irwin@gmail.com>;    550 05BIT008_mongoid
 <conrad.irwin@gmail.com>;    432 01Click_paperclip
 <conrad.irwin@gmail.com>;    423 00kaito_java-docs-samples
 <conrad.irwin@gmail.com>;    402 040mike_rack
 <conrad.irwin@gmail.com>;    399 0bserver07_rails_admin
 <conrad.irwin@gmail.com>;    395 07033320a_select2
 <conrad.irwin@gmail.com>;    392 12spokes_errbit
 <conrad.irwin@gmail.com>;    388 0204111zhouyi_activeadmin
 <conrad.irwin@gmail.com>;    374 0paipai_delayed_job

 <JRadosz@gmail.com>;    656 0101011_atom
 <JRadosz@gmail.com>;    531 01binary_jest
 <JRadosz@gmail.com>;    530 001szymon_coffeescript
 <JRadosz@gmail.com>;    441 216Giorgiy_eslint
 <JRadosz@gmail.com>;    397 0njzy0_flow
 <JRadosz@gmail.com>;    339 0x2206_yargs
 <JRadosz@gmail.com>;    337 01BTC10_npm
 <JRadosz@gmail.com>;    321 00Davo_mocha
 <JRadosz@gmail.com>;    305 00SHUAI00_react-router
 <JRadosz@gmail.com>;    304 0-to-1_redux
 
0xced <cedric.luthi@gmail.com>;    670 000fan000_ReactiveCocoa
0xced <cedric.luthi@gmail.com>;    539 08_fmdb
0xced <cedric.luthi@gmail.com>;    478 0xced_CocoaLumberjack
0xced <cedric.luthi@gmail.com>;    431 007895175_SDWebImage
0xced <cedric.luthi@gmail.com>;    420 002and001_Texture
0xced <cedric.luthi@gmail.com>;    407 0359xiaodong_MBProgressHUD
0xced <cedric.luthi@gmail.com>;    402 0jun0815_SwiftLint
0xced <cedric.luthi@gmail.com>;    368 1103785815_CocoaAsyncSocket
0xced <cedric.luthi@gmail.com>;    354 0x73_MagicalRecord
0xced <cedric.luthi@gmail.com>;    344 0day2010_RestKit

--global <tomas.vot@gmail.com>;    655 0mars_SonataAdminBundle
--global <tomas.vot@gmail.com>;    653 013cbr_FOSUserBundle
--global <tomas.vot@gmail.com>;    463 00577_PHPExcel
--global <tomas.vot@gmail.com>;    435 01viniciusmelo_dolibarr
--global <tomas.vot@gmail.com>;    405 3mg_DoctrineExtensions
--global <tomas.vot@gmail.com>;    405 00F100_composer
--global <tomas.vot@gmail.com>;    389 1001Pharmacies_Sylius
--global <tomas.vot@gmail.com>;    351 0hex_phpunit
--global <tomas.vot@gmail.com>;    346 20uf_http-foundation
--global <tomas.vot@gmail.com>;    312 74Labs_NelmioApiDocBundle

 <jonsimo@aol.com>;   1866 00Delphina00_00Delphina00.github.io
 <jonsimo@aol.com>;   1747 2thumbstech_React-Insta-Clone
 <jonsimo@aol.com>;   1622 00Delphina00_React-UI-Components
 <jonsimo@aol.com>;   1564 120356aa_Sprint-Challenge-React-Smurfs
 <jonsimo@aol.com>;   1556 13JonathanDavid_Tabs-Components
 <jonsimo@aol.com>;   1540 AAsriyan_HTTP-AJAX
 <jonsimo@aol.com>;   1475 AAsriyan_Redux-Todo
 <jonsimo@aol.com>;   1440 13JonathanDavid_Redux-Counter
 <jonsimo@aol.com>;   1281 LambdaSchool_Responsive-Web-Design
 <jonsimo@aol.com>;   1207 13JonathanDavid_Sprint-Challenge-Lambda-Times-React
 
 <the.sk89q@gmail.com>;   1087 0demongamer0_Bukkit
 <the.sk89q@gmail.com>;    802 0uti_MinecraftForge
 <the.sk89q@gmail.com>;    401 1000swoo_SpongeAPI
 <the.sk89q@gmail.com>;    390 0999312_Minecraft-Mod-Language-Package
 <the.sk89q@gmail.com>;    277 1bytee_BungeeCord
 <the.sk89q@gmail.com>;    262 1Rogue_Essentials
 <the.sk89q@gmail.com>;    259 0Laired_WorldEdit-rus-by-DarkFort
 <the.sk89q@gmail.com>;    255 192733213213_BuildCraft
 <the.sk89q@gmail.com>;    220 0Laired_WorldGuard-rus-by-DarkFort
 <the.sk89q@gmail.com>;    215 Agrajagd_Spout
 
 <caniszczyk@gmail.com>;    518 faizann_restcomm.github.io
 <caniszczyk@gmail.com>;    457 Amoshappy_envoy
 <caniszczyk@gmail.com>;    327 0388989777_opensource.guide
 <caniszczyk@gmail.com>;    311 0-T-0_bootstrap
 <caniszczyk@gmail.com>;    308 jain-slee
 <caniszczyk@gmail.com>;    297 Ark-kun_test-infra
 <caniszczyk@gmail.com>;    295 13221325403_Restcomm-Connect
 <caniszczyk@gmail.com>;    291 AdamDang_community
 <caniszczyk@gmail.com>;    283 deev_jain-sip.ha
 <caniszczyk@gmail.com>;    271 Acidburn0zzz_eclipse.platform.ui
 
 <adam@adamralph.com>;    647 FranciscoGileno_haacked.com
 <adam@adamralph.com>;    562 06needhamt_roslyn
 <adam@adamralph.com>;    517 AArnott_visualstudio-docs
 <adam@adamralph.com>;    502 216Giorgiy_cli
 <adam@adamralph.com>;    493 0-wiz-0_coreclr
 <adam@adamralph.com>;    326 0x53A_Paket
 <adam@adamralph.com>;    295 007killyou_ServiceStack
 <adam@adamralph.com>;    281 A-And_corert
 <adam@adamralph.com>;    272 0x53A_FAKE
 <adam@adamralph.com>;    227 AaronBertrand_sql-docs
 
100fue <35096923+100fue@users.noreply.github.com>;   1846 ironhack-labs_lab-bootstrap-cloning-revera
100fue <35096923+100fue@users.noreply.github.com>;   1778 lauphern_debug-test
100fue <35096923+100fue@users.noreply.github.com>;   1768 4ortytwo_lab-express-spotify
100fue <35096923+100fue@users.noreply.github.com>;   1653 3lv27_lab-javascript-memory-game
100fue <35096923+100fue@users.noreply.github.com>;   1616 JorgeMtzCrz_slack
100fue <35096923+100fue@users.noreply.github.com>;   1606 j31_lab-advance-querying-mongo
100fue <35096923+100fue@users.noreply.github.com>;   1589 komi24_ih-contacts
100fue <35096923+100fue@users.noreply.github.com>;   1549 j31_lab-passport-roles
100fue <35096923+100fue@users.noreply.github.com>;   1543 nlevo_project2
100fue <35096923+100fue@users.noreply.github.com>;   1515 elwiiman_theGame

Tomster <tomster@emberjs.com>;    422 0cool321_ember.js
Tomster <tomster@emberjs.com>;    416 0xadada_ember-cli
Tomster <tomster@emberjs.com>;    358 04th0809_website
Tomster <tomster@emberjs.com>;    323 1123456789011_data
Tomster <tomster@emberjs.com>;    208 Acidburn0zzz_guides-1
Tomster <tomster@emberjs.com>;    200 0000marcell_ember-cli-mirage
Tomster <tomster@emberjs.com>;    199 0xadada_ember-simple-auth
Tomster <tomster@emberjs.com>;    171 0000marcell_ember-power-select
Tomster <tomster@emberjs.com>;    170 Alonski_guides-source
Tomster <tomster@emberjs.com>;    160 3tarazona_ember-paper

4383 <herveberaud.pro@gmail.com>;   1017 10171121_cinder
4383 <herveberaud.pro@gmail.com>;    906 1024122298_openstack-manuals
4383 <herveberaud.pro@gmail.com>;    777 21-guns_tempest
4383 <herveberaud.pro@gmail.com>;    766 dmgerman_nova-token: OpenStack Nova fabric controller, supporting a wide variety of virtualization technologies
4383 <herveberaud.pro@gmail.com>;    670 07053220ab_glance
4383 <herveberaud.pro@gmail.com>;    648 08_keystone
4383 <herveberaud.pro@gmail.com>;    636 DingDino_requirements
4383 <herveberaud.pro@gmail.com>;    555 15527370365_aodh
4383 <herveberaud.pro@gmail.com>;    547 249043822_heat
4383 <herveberaud.pro@gmail.com>;    531 02307_phalcon-devtools

 <dloewenherz@gmail.com>;   1502 jerng-org_ruthenium
 <dloewenherz@gmail.com>;    498 0x73_MagicalRecord
 <dloewenherz@gmail.com>;    335 1641731459_scipy
 <dloewenherz@gmail.com>;    261 000fan000_Alamofire
 <dloewenherz@gmail.com>;    259 07cs07_objective-c-style-guide
 <dloewenherz@gmail.com>;    247 002and001_GPUImage
 <dloewenherz@gmail.com>;    222 202works_django-haystack
 <dloewenherz@gmail.com>;    195 007895175_SDWebImage
 <dloewenherz@gmail.com>;    169 07101040115_celery
 <dloewenherz@gmail.com>;    167 007lihegong_IQKeyboardManager
 
  <joan@montane.cat>;    715 01shobitha_vlc
  <joan@montane.cat>;    627 4711_Osmand
  <joan@montane.cat>;    538 1st8_converse.js
  <joan@montane.cat>;    389 Droces_easyssh
  <joan@montane.cat>;    387 0h3r0_jitsi-meet
  <joan@montane.cat>;    381 07101994_KISS
  <joan@montane.cat>;    308 0xf1sh_scummvm
  <joan@montane.cat>;    307 project-draco-hr_Osmand
  <joan@montane.cat>;    296 0359xiaodong_FreeRDP
  <joan@montane.cat>;    269 0x90909_binary-static-www3
  
 <thomas.dubuisson@gmail.com>;    448 0xfinch_Signal-Android
 <thomas.dubuisson@gmail.com>;    412 1930s_ghc
 <thomas.dubuisson@gmail.com>;    354 0x0ece_druid
 <thomas.dubuisson@gmail.com>;    331 2586252322_wesnoth
 <thomas.dubuisson@gmail.com>;    299 0xmohit_cabal
 <thomas.dubuisson@gmail.com>;    293 23Skidoo_stackage
 <thomas.dubuisson@gmail.com>;    238 AndreasPK_containers
 <thomas.dubuisson@gmail.com>;    236 08guye_cas
 <thomas.dubuisson@gmail.com>;    229 1930s_yesod
 <thomas.dubuisson@gmail.com>;    198 0xc0re_Hygieia
 
 <asf@boinkor.net>;    602 0mp_nix
 <asf@boinkor.net>;    472 C4K3_wiki.vg
 <asf@boinkor.net>;    429 0chen0_thrift
 <asf@boinkor.net>;    367 8573_el-get
 <asf@boinkor.net>;    316 4e554c4c_rfcs
 <asf@boinkor.net>;    305 0X1A_cargo
 <asf@boinkor.net>;    261 roblabla_rust-libstd
 <asf@boinkor.net>;    250 2rs2ts_dd-agent
 <asf@boinkor.net>;    195 ANDRON94_sbcl
 <asf@boinkor.net>;    186 1998SYQ_tokio
 <k-okada@jsk.t.u-tokyo.ac.jp>;    802 130s_rosdistro
 <k-okada@jsk.t.u-tokyo.ac.jp>;    522 10eTechnology_librealsense
 <k-okada@jsk.t.u-tokyo.ac.jp>;    423 122689305_ros_comm
 <k-okada@jsk.t.u-tokyo.ac.jp>;    411 6RiverSystems_navigation
 <k-okada@jsk.t.u-tokyo.ac.jp>;    369 130s_moveit-1
 <k-okada@jsk.t.u-tokyo.ac.jp>;    291 130s_gazebo_ros_pkgs
 <k-okada@jsk.t.u-tokyo.ac.jp>;    247 130s_jsk_recognition
 <k-okada@jsk.t.u-tokyo.ac.jp>;    247 0Jiahao_image_pipeline
 <k-okada@jsk.t.u-tokyo.ac.jp>;    207 0table0_rqt_common_plugins
 <k-okada@jsk.t.u-tokyo.ac.jp>;    187 130s_universal_robot
 <Callek@gmail.com>;   1094 0dysseas_bedrock
 <Callek@gmail.com>;    625 Anu-Kannan_https-github.com-dnase-classroom-control-vf
 <Callek@gmail.com>;    623 08_puppet
 <Callek@gmail.com>;    505 AdrianUrsea_www.mozilla.org
 <Callek@gmail.com>;    425 8v060htwyc_build-tools
 <Callek@gmail.com>;    350 TheoChevalier_en-US-central
 <Callek@gmail.com>;    336 TDA_fxa-content-server-l10n
 <Callek@gmail.com>;    329 TheoChevalier_en-US
 <Callek@gmail.com>;    327 TheoChevalier_en-US-test
 <Callek@gmail.com>;    327 0cjs_buildbot
 <wking@tremily.us>;   1151 gentoo_wikiclone
 <wking@tremily.us>;    608 oumpy_oumpy.github.io
 <wking@tremily.us>;    474 Adam-Weiss_ovirt-site
 <wking@tremily.us>;    455 07mehmet0_origin
 <wking@tremily.us>;    414 1235gopalt_openshift-ansible
 <wking@tremily.us>;    412 AdamDang_community
 <wking@tremily.us>;    376 Ark-kun_test-infra
 <wking@tremily.us>;    321 BrandWang_API
 <wking@tremily.us>;    271 AmitRoushan_cloud-provider-openstack
 <wking@tremily.us>;    260 1ForA11_client-go
 <dschlaud@gmail.com>;   1886 A7madXatab_js-ajax-fetch-lab-re-coded-erbil-2018
 <dschlaud@gmail.com>;   1762 Almighty-Mose_css-manifests-lab-v-000
 <dschlaud@gmail.com>;   1727 81Jeremiah_crud-lab-v-000
 <dschlaud@gmail.com>;   1680 Almighty-Mose_js-advanced-functions-intro-to-mocha-readme-v-000
 <dschlaud@gmail.com>;   1672 81Jeremiah_using-to-json-lab-v-000
 <dschlaud@gmail.com>;   1641 81Jeremiah_return-string-data-lab-v-000
 <dschlaud@gmail.com>;   1638 A7madXatab_js-ajax-hitting-apis-lab-re-coded-erbil-2018
 <dschlaud@gmail.com>;   1608 Almighty-Mose_apis-and-faraday-reading-v-000
 <dschlaud@gmail.com>;   1573 AAM77_web-auth-readme-v-000
 <dschlaud@gmail.com>;   1554 81Jeremiah_basic-apis-lab-v-000
 <davidglick@onenw.org>;   1098 collective_imsvdex
 <davidglick@onenw.org>;    816 8bitpanda_intro-to-heroku
 <davidglick@onenw.org>;    489 2silver_Products.CMFPlone
 <davidglick@onenw.org>;    414 28554010_edx-theme
 <davidglick@onenw.org>;    366 25th-floor_plone.app.locales
 <davidglick@onenw.org>;    353 117111302_configuration
 <davidglick@onenw.org>;    344 hoka_zope
 <davidglick@onenw.org>;    320 codeix_z3c.recipe.dev
 <davidglick@onenw.org>;    269 4teamwork_Products.TinyMCE
 <davidglick@onenw.org>;    265 plone_plone.pony


#do unweighted oslom as well
zcat P2AFullS.nA2Au.$c.PLM|perl -e '$i=0;while(<STDIN>){chop();($c,$w)=split(/;/);$cl{$c}{$i}++;$i++};while(($k,$v)=each %cl){@a=keys %$v;print "@a\n";}' > forOSLOu
time /home/audris/src/OSLOM2/oslom_undir -f /data/play/forks/P2AFullS.nA2A.$c.csv -hint forOSLO -uw -r 1 -hr 1
#################################################################################################################################



############# OLD #######################################
for c in 10 100 1000 2000
do zcat P2AFullS.nA2A.$c.s | perl connectExportVw.perl P2AFullS.nA2A.$c
done


export LD_LIBRARY_PATH=/home/audris/src/networkit/build
for c in 10 100 1000 2000
do zcat /data/P2AFullS.nA2A.$c.versions | /home/audris/src/networkit/cluster $(zcat /data/P2AFullS.nA2A.$c.names|wc -l) | gzip > /data/P2AFullS.nA2Au.$c.PLM
   paste -d\  <(zcat /data/P2AFullS.nA2A.$c.versions) <(zcat /data/P2AFullS.nA2A.$c.weights) < | /home/audris/src/networkit/clusterw $(zcat /data/P2AFullS.nA2A.$c.names|wc -l) $(zcat /data/P2AFullS.nA2A.$c.weights|wc -l) | gzip > /data/P2AFullS.nA2Aw.$c.PLM
done
25537691 edges read
modularity=0.999777 nver=28734132 clusters=5175384 largest=34327
50996352 edges read
modularity=0.926016 nver=40769789 clusters=3863554 largest=1683997
59778297 edges read
modularity=0.887575 nver=42388321 clusters=3701295 largest=2355315




#after excluding root admin
c=100

modularity=0.940373 nver=40769789 clusters=4121379 largest=2402178
zcat P2AFullS.nA2Aw.$c.PLM | perl rankNew.perl P2AFullS.nA2A.$c 1 | gzip >  P2AFullS.nA2Aw.$c.crank.map
1362927  <gzlist@googlemail.com>   # is from 3370sohail_gecko-dev
 255915 --get <primetoxinzz@gmail.com> # MC4EP_BetterWithMods
 232945 A-malagon <alejanm09@gmail.com> # dpayo_X-Serv-14.9-CuatroAplis
 225796 Aigars Mahinovs <aigarius@debian.org>
 223459 = <h.eskandari@gmail.com>
 223450 "Enis Afgan ext:(%22) <afgane@gmail.com>
 219292 AJ Alilano <alan.alilano@boom.camp>
 203539  <gerwin.klein@nicta.com.au>
 181806 Andy H. Jung <commaniakr@gmail.com>
 171620  Brian Stansberry <brian.stansberry@redhat.com>

zcat P2AFullS.nA2Aw.$c.PLM | perl rankNew.perl P2AFullS.nA2A.$c 0 | gzip >  P2AFullS.nA2Aw.$c.crankmr.map
zcat P2AFullS.nA2A.$c.s has groups of authors
use P2AFullS.nA2A to reconstruct projects involved

# proceed with forking
#get 100-1K aouthor groups
zcat P2AFullS.nA2Aw.100.crank.map| cut -d\; -f2 | uniq -c | awk '{if ($1 < 1000 && $1 > 100) print $0}' | perl -ane 's|^\s+[0-9]+ ||;print' > P2AFullS.nA2Aw.100.crank.map.100-1000.cl
head -10 P2AFullS.nA2Aw.100.crank.map.100-1000.cl
echo 'Harry king <deon.ruben@ln0hio.com>'  | ~/lookup/getValues -f a2P | cut -d\; -f2  | ~/lookup/getValues -f P2a
#rental page teplate for wp
echo 'Eli Vargas <eevargas@gmail.com>'  | ~/lookup/getValues -f a2P | cut -d\; -f2  | ~/lookup/getValues -f P2a
#looks like real software
echo '126308 <126308@myrp.edu.sg>'  | ~/lookup/getValues -f a2P | cut -d\; -f2  | ~/lookup/getValues -f P2a
#looks like class project
echo 'Benjamin Ozore <65028432+BenOzore@users.noreply.github.com>'  | ~/lookup/getValues -f a2P | cut -d\; -f2  | ~/lookup/getValues -f P2a
#movie database

zcat P2AFullS.nA2Aw.100.crank.map | ~/lookup/s2sBinSorted.perl /fast/a2Aw100Full$ver.tch
zcat P2AFullS.nA2Aw.100.crank.map | perl -ane 'chop();($a,$b)=split(/;/);print "$b;$a\n";' | ~/lookup/s2sBinSorted.perl /fast/A2aw100Full$ver.tch

echo 'Harry king <deon.ruben@ln0hio.com>'  | ~/lookup/getValues -f a2P | cut -d\; -f2  | ~/lookup/getValues -f P2a|wc
echo 'Harry king <deon.ruben@ln0hio.com>'  | ~/lookup/getValues -f A2aw100 | wc

i="712346172@qq.com <712346172@qq.com>"
echo $i | ~/lookup/getValues -f a2P | cut -d\; -f2  | ~/lookup/getValues -f P2a|wc
      3       6     194
echo $i | ~/lookup/getValues -f A2aw100|wc                                         
    511    1557   36669
echo $i | ~/lookup/getValues -f A2aw100| cut -d\; -f2 |~/lookup/getValues -f a2P|cut -d\; -f2 | lsort 1G |uniq -c | lsort 1G -rn | head
     34 L1Z1_2018Android-Study
     10 theoaaa_os-curriculum-design
     10 h2pl_Java-
      9 lazyman-lzy_bvc-scene-api
      9 jianwan_Android-Study
      8 fancn21th_0_85_Big_Screen
head -10 P2AFullS.nA2Aw.100.crank.map.100-1000.cl | while IFS=\; read i; do  echo $i $(echo $i | ~/lookup/getValues -f A2aw100| cut -d\; -f2 |~/lookup/getValues -f a2P|cut -d\; -f2 | lsort 1G |uniq -c | lsort 1G -rn | head -1); done 

Harry king <deon.ruben@ln0hio.com> 73 rental345_english
712346172@qq.com <712346172@qq.com> 34 L1Z1_2018Android-Study
Eli Vargas <eevargas@gmail.com> 16 0000000000000000000009_prework-about-me
Asus <Asus@172.17.82.79> 57 obedottoc_Ex02
126308 <126308@myrp.edu.sg> 86 Dylanwasd_L02-03
abbaswall <67048135+abbaswall@users.noreply.github.com> 21 deotel_tank
Alban <alban.lfbv@gmail.com> 60 hsabouri_42_doom_nukem
L402 <2013091veby@gmail.com> 7 Ikhsandi007_FinalPrakKakas
Athira <aathiraspillai@gmail.com> 27 bitbucket.org_Athira_Pillai_dev-hbpo-server
Benjamin Ozore <65028432+BenOzore@users.noreply.github.com> 23 0therese0_mws-restaurant-stage-1

head -10 P2AFullS.nA2Aw.100.crank.map.100-1000.cl | while IFS=\; read i; do  echo "$i;"$(echo $i | ~/lookup/getValues -f A2aw100| cut -d\; -f2 |~/lookup/getValues -f a2P|cut -d\; -f2 | lsort 1G |uniq -c | lsort 1G -rn | head -3); done 
Harry king <deon.ruben@ln0hio.com>; 73 rental345_english 58 rental345_eng2
712346172@qq.com <712346172@qq.com>; 34 L1Z1_2018Android-Study 10 theoaaa_os-curriculum-design 10 h2pl_Java-
Eli Vargas <eevargas@gmail.com>; 16 0000000000000000000009_prework-about-me 14 BrittaB_Sprout 12 byronferguson_project2
Asus <Asus@172.17.82.79>; 57 obedottoc_Ex02 57 obedottoc_Ex01 55 obedottoc_EX03
126308 <126308@myrp.edu.sg>; 86 Dylanwasd_L02-03 59 Ivannjc_P04-PS 53 Stun4life_Test123
abbaswall <67048135+abbaswall@users.noreply.github.com>; 21 deotel_tank 20 bothdfef_ghjk 16 utakul_lewo
Alban <alban.lfbv@gmail.com>; 60 hsabouri_42_doom_nukem 59 Alblfbv_CoreWar 47 Alblfbv_Lem_in
L402 <2013091veby@gmail.com>; 7 Ikhsandi007_FinalPrakKakas 6 luthfiyaanggraini7_collaborator2 6 jmkl_TestPulsa
Athira <aathiraspillai@gmail.com>; 27 bitbucket.org_Athira_Pillai_dev-hbpo-server 17 leenaNi_estore 9 somoplay_eadating
Benjamin Ozore <65028432+BenOzore@users.noreply.github.com>; 23 0therese0_mws-restaurant-stage-1 16 jniziol_mittflix 13 01ank1998_MyReads-Book-Manager-React

zcat P2AFullS.nA2Aw.100.crank.map| cut -d\; -f2 | uniq -c | awk '{if ($1 > 100000) print $0}' | perl -ane 's|^\s+[0-9]+ ||;print' > P2AFullS.nA2Aw.100.crank.map.100k.cl
cat P2AFullS.nA2Aw.100.crank.map.100k.cl | while IFS=\; read i; do  echo "$i;"$(echo "$i" | ~/lookup/getValues -f A2aw100| cut -d\; -f2 |~/lookup/getValues -f a2P|cut -d\; -f2 | lsort 1G |uniq -c | lsort 1G -rn | head -3); done > P2AFullS.nA2Aw.100.crank.map.100k.cl.top3

Andy H. Jung <commaniakr@gmail.com>; 1526 00101010_00101010.github.io 737 04240_git-training 340 00-Evan_minimal-mistakes
$ <bruce.j.beare@intel.com>; 12115 0-T-0_ps4-linux 7571 1nfrag_android_hardware_qcom_audio 7176 0-kaladin_build
 <graydon@mozilla.com>; 1773 0-admin_bitcoin 677 0DanielTehrani_monero-fork 655 0-0-0-_go-ethereum
 <armenzg@mozilla.com>; 1165 3370sohail_gecko-dev 754 16patsle_gaia 732 0e4ef622_rust
A-malagon <alejanm09@gmail.com>; 1344 00101010_00101010.github.io 811 004307ec_caffe 716 130s_rosdistro
AJ Alilano <alan.alilano@boom.camp>; 2356 098saket_first-contributions 2199 281team29_hacktoberfest 1229 01QueN10_Make-a-Pull-Request
1362927  <gzlist@googlemail.com>;
= <h.eskandari@gmail.com>; 2860 0xc0re_azure-content 908 00101010_00101010.github.io 796 0xced_azure-sdk-for-net
Aigars Mahinovs <aigarius@debian.org>; 3836 0-T-0_ps4-linux 1115 00ERNEST00_FFmpeg 942 0day-ci_xen
"Enis Afgan ext:(%22) <afgane@gmail.com>; 1454 00101010_00101010.github.io 1219 0000m0000_ProgrammingAssignment2 1191 12345678yirgu_2018-06-06-AkU
 Brian Stansberry <brian.stansberry@redhat.com>; 976 0-T-0_ps4-linux 940 0irebRwE_openstack 874 131_puppetlabs-apache
--get <primetoxinzz@gmail.com>; 871 007Axel_spacestation413 656 00101010_00101010.github.io 647 0uti_MinecraftForge
 <gerwin.klein@nicta.com.au>; 950 00101010_00101010.github.io 878 00yk_spacemacs 790 0xbase12_leiningen
Aashish Chaudhary <aashish.chaudhary@kitware.com>; 655 00101010_00101010.github.io 590 0-T-0_julia 519 4gh_METADATA.jl
- - <945852046@qq.com>; 2584 cirosantilli_imagine-all-the-people 706 00101010_00101010.github.io 238 0000wei_element

#fp
python3 embeddingcsv.py /data/fpw.csv
Loading embeddings...
Done loading embeddings (shape: (25181406, 128)).
read (25181406, 128)
UMAP(a=None, angular_rp_forest=True, b=None,
   force_approximation_algorithm=False, init='spectral', learning_rate=1.0,
   local_connectivity=1.0, low_memory=False, metric='cosine',
   metric_kwds=None, min_dist=0.01, n_components=2, n_epochs=11,
   n_neighbors=15, negative_sample_rate=5, output_metric='euclidean',
   output_metric_kwds=None, random_state=None, repulsion_strength=1.0,
   set_op_mix_ratio=1.0, spread=1.0, target_metric='categorical',
   target_metric_kwds=None, target_n_neighbors=-1, target_weight=0.5,
   transform_queue_size=4.0, transform_seed=42, unique=False, verbose=True)
Construct fuzzy simplicial set
Fri Dec 18 15:29:39 2020 Finding Nearest Neighbors
Fri Dec 18 15:29:39 2020 Building RP forest with 256 trees
Sat Dec 19 05:36:06 2020 NN descent for 25 iterations
Segmentation fault (core dumped)

#explore some
echo ":START_ID;:END_ID" > a2P.csv
echo 'Eli Vargas <eevargas@gmail.com>' |  ~/lookup/getValues -f A2aw100 | cut -d\; -f2 |  ~/lookup/getValues -f A2a | cut -d\; -f2 | ~/lookup/getValues -f a2P | lsort 1G -t\; -k1,1 > a2P 
sort -u a2P | ~/lookup/getValues -f a2P 2> /dev/null | awk -F\; '{print $3";"$2}' | sort -u >> a2P.csv

echo "nID:ID;:LABEL" > n.csv
(cut -d\; -f1 a2P.csv | sort -u | perl -ane 'chop();print "$_;A\n"'; cut -d\; -f2 a2P.csv |sort -u |perl -ane 'chop();print "$_;Pr\n"')| grep -v ^: >> n.csv
./bin/neo4j-admin import --separator=\; --nodes=nID=/data/n.csv  --relationships=a2P=/data/a2P.csv

./bin/neo4j-admin import --delimiter=\; --array-delimiter=\| --nodes=nID=/data/n.csv  --relationships=a2P=/data/a2P.csv

for c in 10 100 1000
do zcat P2AFullS.nA2Au.$c.PLM | perl rankNew.perl P2AFullS.nA2A.$c 1 | gzip >  P2AFullS.nA2Au.$c.crank.map
done

for c in 10 100 1000
do zcat P2AFullS.nA2Au.$c.crank.map|cut -d\; -f2 | lsort 2G -u -t\; | wc
   zcat P2AFullS.nA2Au.$c.crank.map|cut -d\; -f2 | lsort 20G -t\; | uniq -c | lsort 3G -rn | head -20
done
5175379
  21047 Charlie <charles.w.marlow@gmail.com>
  18071 Alvaro <rednaxela.007@hotmail.com>
  15851 9614 <46434108+9614@users.noreply.github.com>
  15803 Adrian Lorentzen <Adrian Lorentzen>
  14284 CaiEddie <>
  13130 Darin Buzon <dabuzon@Darins-MacBook-Pro.local>
  12971 Aeranythe Echosong <KatherinaXC@users.noreply.github.com>
  12216 99ro <robindeboer1999@gmail.com>
  11983 Derek Nguyen <48035548+derek-s-nguyen@users.noreply.github.com>
  11249 André H <myersger@googlemail.com>
  11166 Bosslv0 <46979902+Bosslv0@users.noreply.github.com>
  10314 Ainsley Malcolm Pereira <s10186606@connect.np.edu.sg>
  10205 Adam Blank <blank@caltech.edu>
  10185 Chuan Ren (RIT Student) <cr3186@ad.rit.edu>
  10097 Dario Panzuto <dario.panzuto@gmail.com>
   9912 AntonioTonchev <36622262+AntonioTonchev@users.noreply.github.com>
   9109 Bryan Tan <wytan.2017@sis.smu.edu.sg>
   8971 Dhivakar <dhivakar.kanagaraj@solitontech.com>
   8912 DanialVEVO <DanialVEVO@users.noreply.github.com>
   8780 Alejandro Gutierrez <a.guti3@hotmail.com>

(echo ":START_ID,weight,:END_ID"; zcat /da4_data/play/forks/P2AFullS.nA2A.100.versions | paste -d\  - <(zcat /da4_data/play/forks/P2AFullS.nA2A.100.weights)| perl -ane 'chop();($a,$b,$w)=split(/ /);print "$a,$w,$b\n"')  > /data/play/forks/ew.csv
(echo "nID:ID,title"; zcat /da4_data/play/forks/P2AFullS.nA2A.100.names|sed 's|[,\r"]| |g' | sed "s|'| |" | awk '{print i++","$0}') > /data/play/forks/n.csv
(echo ":START_ID,:END_ID"; zcat /da4_data/play/forks/P2AFullS.nA2A.100.versions | sed 's| |,|')  > /data/play/forks/e.csv
bin/neo4j-admin import --nodes=nID=/data/n.csv --relationships=IN=/data/e.csv

./bin/neo4j console
ctrl-C
neo4j-admin set-initial-password  neo4j
./bin/neo4j-admin import --nodes=nID=/data/n.csv  --relationships=IN=/data/e.csv --relationships=INW=/data/ew.csv
./bin/neo4j console

bin/cypher-shell

CALL gds.graph.create( 'grw',  'nID',  { INW: {orientation: 'UNDIRECTED',properties: 'weight'} });
CALL gds.fastRP.write('grw',{ embeddingDimension: 128, iterationWeights:[0,0.5,0.5,1,1,1,1,1,1,1,1], normalizationStrength:-.1, relationshipWeightProperty:'weight', writeProperty:'fastrpw'});
CALL apoc.export.csv.query("MATCH (n) WHERE EXISTS(n.`fastrpw`) RETURN n.nID as id, n.`fastrpw` as fr", "fpw.csv",{});
mv import/fpw.csv data/

CALL gds.graph.create( 'gr',  'nID',  'IN');
CALL gds.fastRP.write('gr',{ embeddingDimension: 128, iterationWeights:[0,1,1,1,1,1,1,1,1,1], normalizationStrength:-.1, writeProperty:'fastrp'});
CALL apoc.export.csv.query("MATCH (n) WHERE EXISTS(n.`fastrp`) RETURN n.nID as id, n.`fastrp` as fr", "fp.csv",{});
mv import/fp.csv data/

scp -p ../fpw.gz da5:/data/play/forks/
grep '^"[0-9]' fp.csv | sed 's|^"||;s|","\[|;|;s|,|;|g;s|\]"||' | gzip > ../fp.gz
scp -p ../fp.gz da5:/data/play/forks/
grep '^"[0-9]' fpw.csv | sed 's|^"||;s|","\[|;|;s|,|;|g;s|\]"||' | gzip > ../fpw.gz
scp -p ../fpw.gz da5:/data/play/forks/

cd /home/audris/src/verse
zcat /data/fp.gz | perl -e 'open A, "zcat /data/P2AFullS.nA2A.100.names|"; $i=0;while (<A>){chop();$b{$i}++ if $_ =~ /^cl[0-9]+$/;$i++}; while(<STDIN>){chop();($i,@x)=split(/;/);print "".(join ";",@x)."\n" if !defined $b{$i};}' | gzip > /data/fp.csv
python3 embeddingcsv.py /data/fp.csv | gzip > /data/fp.csv.embed
zcat /data/fpw.gz | perl -e 'open A, "zcat /data/P2AFullS.nA2A.100.names|"; $i=0;while (<A>){chop();$b{$i}++ if $_ =~ /^cl[0-9]+$/;$i++}; while(<STDIN>){chop();($i,@x)=split(/;/);print "".(join ";",@x)."\n" if !defined $b{$i};}' | gzip > /data/fpw.csv
python3 embeddingcsv.py /data/fpw.csv | gzip > /data/fpw.csv.embed

3882638
1439401 0candy <ngcandy@ca.ibm.com>
 329782 - - <945852046@qq.com>
 276553 --get <primetoxinzz@gmail.com>
 261904 "Nate Koenig ext:(%22) <natekoenig@gmail.com>
 256144 = <linkgeorg@gmail.com>
 234151 "Enis Afgan ext:(%22) <afgane@gmail.com>
 227196 = <h.eskandari@gmail.com>
 189009 Andy H. Jung <commaniakr@gmail.com>
 175868 Aigars Mahinovs <aigarius@debian.org>
 171879  Brian Stansberry <brian.stansberry@redhat.com>
 170811 Aashish Chaudhary <aashish.chaudhary@kitware.com>
 150853  <lucianopf@outlook.com>
 137642 1337GAK <geokonst95@gmail.com>
 117416 $ <bruce.j.beare@intel.com>
 107572 1000TurquoisePogs <30730276+1000TurquoisePogs@users.noreply.github.com>
  97227 Fujita Yuko <fujita.yuko@tis.co.jp>
  94209 ( pancake ) <pancake@youterm.com>
  90562 Albert <al.hu@berkeley.edu>
  85649 Adam <35268982+AdamHawtin@users.noreply.github.com>
  82761 Abimbola <abimbola2010@gmail.com>


3693030 8949710 160095759
1260381  <deller@gmx.de>
1204025  <caniszczyk@gmail.com>
 592828 --global <tomas.vot@gmail.com>
 425531 /orta <orta.therox@gmail.com>
 347293  <adam@adamralph.com>
 328153 "Johannes Schindelin" <Johannes.Schindelin@gmx.de>
 312722  <k-okada@jsk.t.u-tokyo.ac.jp>
 308637  <silviaodwyer119@yahoo.ie>
 288438 "R. Tyler Ballance" <tyler@monkeypox.org>
 285015  <krlmlr+r@mailbox.org>
 276225 3TUSK <urey.s.knowledge@gmail.com>
 222593  <audreyt@audreyt.org>
 205412  <huangdingben_1995@163.com>
 191025 Administrator <Administrator@MSDN-SPECIAL>
 153930 "Luis Henrique Fagundes ext:(%22) <lhfagundes@gmail.com>
 141896  <acyi@ucsd.edu>
 112118  <kunal.sumbly@gmail.com>
 109719 -get <brianhan87usa@gmail.com>
 105699 103066simon <leggsimon@gmail.com>
 103619 207-208 <208_209@outlook.jp>


for c in 10 100 1000
do zcat /data/P2aFullS.na2a.$c.versions | ./cluster $(zcat /data/P2aFullS.na2a.$c.names|wc -l) | gzip > /data/P2aFullS.na2au.$c.PLM
done
37835818 edges read
modularity=0.999531 nver=40512220 clusters=6410130 largest=54034
69057738 edges read
modularity=0.939364 nver=53844628 clusters=4908103 largest=1913913
78315767 edges read
modularity=0.909884 nver=55532050 clusters=4705007 largest=2959284

for c in 10 100 1000
do zcat P2aFullS.na2au.$c.PLM | perl rankNew.perl P2aFullS.na2a.$c 1 | gzip >  P2aFullS.na2au.$c.crank.map
done
c=10
zcat P2aFullS.na2au.$c.crank.map|cut -d\; -f2 | lsort 2G -u -t\; | wc
#6410130
zcat P2aFullS.na2au.$c.crank.map|cut -d\; -f2 | lsort 20G -t\; | uniq -c | lsort 3G -rn | head -20
  34160 Akarsh Hemrajani <ahemrajani3@gatech.edu>
  31059 Aleksandra Sorokina <alexas30s@gmail.com>
  27163 Gnahue <nahuelgalvan@gmail.com>
  25911 CubicNinja64 <abel.carson1@gmail.com>
  24470 Bryan Tan <wytan.2017@sis.smu.edu.sg>
  23723 An Nguyen <An Nguyen>
  22808 --help <taipham@cisco.com>
  22469 Aeranythe Echosong <KatherinaXC@users.noreply.github.com>
  22326 CHOI JOON WON <40337311+Gomyo@users.noreply.github.com>
  19632 Bevin Tang <bevintang@Bevins-MacBook-Air.local>
  18370 Adam Duvin (RIT Student) <axd3010@ad.rit.edu>
  18179 Cao Yizhou <yizhou.cao@finbook.co>
  17387 Marius Brandt <Marius.Brandt@jtl-software.com>
  17352 Alejandro <aledurax@gmail.com>
  17064 Cola Chan <greenzorro@163.com>
  15928 Daniel Yao <ly2328@columbia.edu>
  15755 Alex <alexhors@uw.edu>
  15604 Aleksander Okonski <4145815+GhostOnTheFiber@users.noreply.github.com>
  15386 Alain Milliet <alain.milliet@epfl.ch>
  14776 EMPX <joade361@Student.liu.se>


c=100
zcat P2aFullS.na2au.$c.crank.map|cut -d\; -f2 | lsort 2G -u -t\; | wc
#4908103
zcat P2aFullS.na2au.$c.crank.map|cut -d\; -f2 | lsort 20G -t\; | uniq -c | lsort 3G -rn | head -20
1097040 Aaron Brager <getaaron@gmail.com>
 586792 Aigars Mahinovs <aigarius@debian.org>
 436070 Alexander Opitz <opitz.alexander@googlemail.com>
 325678 Adam Frey <adam@adamfrey.me>
 302861 Adaptivity <verybigbro@gmail.com>
 299496 Agustin Schiavon <v-agschi@microsoft.com>
 292485 Albert Y. Kim <albert.ys.kim@gmail.com>
 277226 Aayush Bisen <aayush23101998@outlook.com>
 239849 Adnan Abdulhussein <adnan@bitnami.com>
 223393 Choi Junwoo <choigo92@nhnent.com>
 213386 Adam Cigánek <adam.ciganek@gmail.com>
 207263 Aaron Blasdel <ablasdel@willowgarage.com>
 192894 Alexander Lutay <1000001@85315e57-5f03-49c3-83f8-201ae2313a75>
 187112 Aashish Chaudhary <aashish.chaudhary@kitware.com>
 185459 Adams <jtyjty99999@126.com>
 163505 52M <52mnee@gmail.com>
 124197 Alex Fisher <alex@linfratech.co.uk>
 120693 00day0 <therandomuser11@gmail.com>
 115273 Allana Caldas <allanacaldas@gmail.com>
 113394 AlexBukach <AlexBukach@1958006.no-reply.drupal.org>

c=1000
zcat P2aFullS.na2au.$c.crank.map|cut -d\; -f2 | lsort 2G -u -t\; | wc
#4705007
zcat P2aFullS.na2au.$c.crank.map|cut -d\; -f2 | lsort 20G -t\; | uniq -c | lsort 3G -rn | head -20
1704216 Aaron Jensen <aaronjensen@gmail.com>
 833349 A S Alam <aalam@users.sf.net>
 660554 Abdul Malik Ikhsan <samsonasik@gmail.com>
 540117 Adam Chainz <adam@adamj.eu>
 521677 98k <18552437190@163.com>
 375495  <darkcoder@darkcoder.localdomain>
 349087 Adam Ralph <adam@adamralph.com>
 340398 Aaron Schumacher <ajschumacher@gmail.com>
 331750 103yiran <1039105206@qq.com>
 318045 Aaron Hill <aa1ronham@gmail.com>
 317288  <@localhost>
 245578 Aaron S. Hawley <aaron.s.hawley@gmail.com>
 215606 1nteger <rhep0820@gmail.com>
 183097 -k <slowdive@me.com>
 175781 Aashish Chaudhary <aashish.chaudhary@kitware.com>
 164751 001u <43672975+001u@users.noreply.github.com>
 164630 18965050 <18965050@qq.com>
 161314 0x00Zhoujialei <35962874+0x00Zhoujialei@users.noreply.github.com>
 159632 Brennen Bearnes <74990+brennen@users.noreply.github.com>
 151746 Aaron D Borden <adborden@a14n.net>



zcat P2aFullS.na2a.s | perl connectExportVw.perl P2aFullS.na2a
zcat /data/P2aFullS.na2a.versions | ./cluster 59388863 | gzip > /data/P2aFullS.na2au.PLM
modularity=0.873406 nver=59388863 clusters=3999512 largest=4108785
zcat P2aFullS.na2au.PLM | perl rankNew.perl P2aFullS.na2a 1 | gzip >  P2aFullS.na2au.crank.map  
zcat P2aFullS.na2au.crank.map|cut -d\; -f2 | lsort 2G -u -t\; | wc
# 3999512
zcat P2aFullS.na2au.crank.map|cut -d\; -f2 | lsort 20G -t\; | uniq -c | lsort 3G -rn | head -20
2201397 GitHub Merge Button <merge-button@github.com>
2051369 dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>
1729189  <>
 699046 Björn Grüning <bjoern@gruenings.eu>
 612134 98k <18552437190@163.com>
 604615 Abdul Malik Ikhsan <samsonasik@gmail.com>
 480856  <Davide@LAPTOP-HEJ8F1OG.localdomain>
 458505  <cys4@DESKTOP-FTIFD5V.localdomain>
 390055 a <
 336546  <annapf@Surface1701.redmond.corp.microsoft.com>
 328224 BuildTools <unconfigured@null.spigotmc.org>
 316221 Google Code Exporter <GoogleCodeExporter@users.noreply.github.com>
 307570 = <=>
 274722 Angular CLI <angular-cli@angular.io>
 267726 Adam Bergmark <adam@bergmark.nl>
 194352 ImgBotApp <ImgBotHelp@gmail.com>
 184917  <jimc@BlackDell.localdomain>
 179311  <JGRAD@JeffGPC.localdomain>
 175287 - <kimhw0820@naver.com>
 170722  <hemanth@PHOENIX.localdomain>



#OSLOM
aa - no result after 31510m41.022s
for c in 10 100 1000
do zcat P2AFullS.nA2A.$c.versions | paste  -d\   - <(zcat P2AFullS.nA2A.$c.weights) > bb$c
done
/home/audris/src/OSLOM2/oslom_undir -f bb10 -w

c=10
cat P2AFullS.nA2A.${c}raw_oslo_files/tp | perl rankNewO.perl P2AFullS.nA2A.${c}raw P2AFullS.nA2Au.$c.PLM |gzip >  P2AFullS.nA2AOw.${c}raw.crank.map
zcat P2AFullS.nA2AOw.${c}raw.crank.map|cut -d\; -f2 | lsort 20G | uniq -c | lsort 1G -rn | head
    628 root <root@074583e83ee3>
    376 Albert Coetzer <albert.coetzer2@gmail.com>
    324 Aaron <Aaron@Aaron-PC>
    298 13518041 Samuel <13518041@std.stei.itb.ac.id>
    264 Angela Hou <angelaho@usc.edu>
    258 Austin Blandford <ABlandford@student.neumont.edu>
    248 * <moumoux001@gmail.com>
    243 Alan Carvalho Galante <alancarvalho@gmail.com>
    222 Achilles Edwin Alfred Saxby <achillessaxby@users.noreply.github.com>
    221 AI-DEV\cap-apreg <annamaria.preg@capgemini.com>

#############################################################################
# OSLOM
#############################################################################

    
c=100
zcat P2AFullS.nA2Au.$c.PLM | perl rankNew.perl P2AFullS.nA2A.$c 1 | gzip >  P2AFullS.nA2Au.$c.crank.map
zcat P2AFullS.nA2Au.${c}.crank.map|cut -d\; -f2 | lsort 20G | uniq -c | lsort 1G -rn | head
1439401 0candy <ngcandy@ca.ibm.com>
 329782 - - <945852046@qq.com>
 276553 --get <primetoxinzz@gmail.com>
 261904 "Nate Koenig ext:(%22) <natekoenig@gmail.com>
 256144 = <linkgeorg@gmail.com>
 234151 "Enis Afgan ext:(%22) <afgane@gmail.com>
 227196 = <h.eskandari@gmail.com>
 189009 Andy H. Jung <commaniakr@gmail.com>
 175868 Aigars Mahinovs <aigarius@debian.org>
 171879  Brian Stansberry <brian.stansberry@redhat.com>

zcat P2AFullS.nA2Aw.$c.PLM | perl rankNew.perl P2AFullS.nA2A.$c 1 | gzip >  P2AFullS.nA2Aw.$c.crank.map
zcat P2AFullS.nA2Aw.${c}.crank.map|cut -d\; -f2 | lsort 20G | uniq -c | lsort 1G -rn | head
1362927  <gzlist@googlemail.com>
 255915 --get <primetoxinzz@gmail.com>
 232945 A-malagon <alejanm09@gmail.com>
 225796 Aigars Mahinovs <aigarius@debian.org>
 223459 = <h.eskandari@gmail.com>
 223450 "Enis Afgan ext:(%22) <afgane@gmail.com>
 219292 AJ Alilano <alan.alilano@boom.camp>
 203539  <gerwin.klein@nicta.com.au>
 181806 Andy H. Jung <commaniakr@gmail.com>
 171620  Brian Stansberry <brian.stansberry@redhat.com>


/home/audris/src/OSLOM2/oslom_undir -f bb100.w -fast -w
***************************************************************************
CHECK UNIONS AND SIMILAR MODULES DONE
******** module_collection ******** 44718 modules. writing...
DONE   ****************************
pruning all the modules collected. Partitions found: 1
getting partition from tp-file: bb100.w_oslo_files/partitions_level_1
44718 groups found
44718 bss found
checking homeless nodes
writing final solution in file bb100.w_oslo_files/short_tp1
******** module_collection ******** 44724 modules. writing...
DONE   ****************************
hierarchies done *********

mv bb$c.w_oslo_files/  P2AFullS.nA2Aw.${c}_oslo_files/
cat P2AFullS.nA2Aw.${c}_oslo_files/tp | perl rankNewO.perl P2AFullS.nA2A.${c} P2AFullS.nA2Aw.$c.PLM |gzip >  P2AFullS.nA2A0w.${c}.crank.map
cat P2AFullS.nA2Aw.${c}_oslo_files/tp1 | perl rankNewO.perl P2AFullS.nA2A.${c} P2AFullS.nA2Aw.$c.PLM |gzip >  P2AFullS.nA2A1w.${c}.crank.map
zcat P2AFullS.nA2A0w.${c}.crank.map|cut -d\; -f2 | lsort 20G | uniq -c | lsort 1G -rn | head
   5140 Aaron Bauman <aaron+github@messageagency.com>
   2652 Bernd Bischl <bernd_bischl@gmx.net>
   2608 'Alexey Khudyakov <alexey.skladnoy@gmail.com>
   2512 18年梦� <getonga2018@gmail.com>
   2413 A K <dev.madand@gmail.com>
   2401 = <mbfrahry@gmail.com>
   2336 --get <primetoxinzz@gmail.com>
   2311  <the.sk89q@gmail.com>
   2136  <philippe.coval@osg.samsung.com>
   2109 00day0 <therandomuser11@gmail.com>

zcat P2AFullS.nA2A1w.${c}.crank.map|cut -d\; -f2 | lsort 20G | uniq -c | lsort 1G -rn | head
 141931 A. M. Lees <Andrew.Lees11@ibm.com>
  59952 Adam Vessey <adam-vessey@users.noreply.github.com>
  39264 Albert <al.hu@berkeley.edu>
  36971 --get <primetoxinzz@gmail.com>
  36463 Adrian F?der <adrian@foeder.de>
  36074  Hryhorii Hrebiniuk <grebinuk@gmail.com>
  31273 (no author) <psiinon@gmail.com>
  29784  <ockham@raz.or.at>
  29116 Adam <adam@novoda.com>
  28503  <the.sk89q@gmail.com>



c=1000
zcat P2AFullS.nA2Au.$c.PLM | perl rankNew.perl P2AFullS.nA2A.$c 1 | gzip >  P2AFullS.nA2Au.${c}.crank.map
zcat P2AFullS.nA2Au.${c}.crank.map|cut -d\; -f2 | lsort 20G | uniq -c | lsort 1G -rn | head
1260381  <deller@gmx.de>
1204025  <caniszczyk@gmail.com>
 592828 --global <tomas.vot@gmail.com>
 425531 /orta <orta.therox@gmail.com>
 347293  <adam@adamralph.com>
 328153 "Johannes Schindelin" <Johannes.Schindelin@gmx.de>
 312722  <k-okada@jsk.t.u-tokyo.ac.jp>
 308637  <silviaodwyer119@yahoo.ie>
 288438 "R. Tyler Ballance" <tyler@monkeypox.org>
 285015  <krlmlr+r@mailbox.org>

zcat P2AFullS.nA2Aw.$c.PLM | perl rankNew.perl P2AFullS.nA2A.$c 1 | gzip >  P2AFullS.nA2Aw.$c.crank.map
zcat P2AFullS.nA2Aw.${c}.crank.map|cut -d\; -f2 | lsort 20G | uniq -c | lsort 1G -rn | head
1292972  <deller@gmx.de>
1160357  <caniszczyk@gmail.com>
 503440 --global <tomas.vot@gmail.com>
 329904  <spencerlyon2@gmail.com>
 326984  <k-okada@jsk.t.u-tokyo.ac.jp>
 326675  <adam@adamralph.com>
 317200 /orta <orta.therox@gmail.com>
 282484  <silviaodwyer119@yahoo.ie>
 258530  <krlmlr+r@mailbox.org>
 251131 3TUSK <urey.s.knowledge@gmail.com>

/home/audris/src/OSLOM2/oslom_undir -f bb1000 -fast -w
***************************************************************************
CHECK UNIONS AND SIMILAR MODULES DONE
******** module_collection ******** 42327 modules. writing... 
DONE   ****************************
pruning all the modules collected. Partitions found: 1
getting partition from tp-file: bb1000_oslo_files/partitions_level_1
42327 groups found
42327 bss found
checking homeless nodes
writing final solution in file bb1000_oslo_files/short_tp1
******** module_collection ******** 42336 modules. writing... 
DONE   ****************************
hierarchies done ********* 

mv bb$c.w_oslo_files/  P2AFullS.nA2Aw.${c}_oslo_files/
cat P2AFullS.nA2Aw.${c}_oslo_files/tp | perl rankNewO.perl P2AFullS.nA2A.${c} P2AFullS.nA2Aw.$c.PLM |gzip >  P2AFullS.nA2A0w.${c}.crank.map
cat P2AFullS.nA2Aw.${c}_oslo_files/tp1 | perl rankNewO.perl P2AFullS.nA2A.${c} P2AFullS.nA2Aw.$c.PLM |gzip >  P2AFullS.nA2A1w.${c}.crank.map
zcat P2AFullS.nA2A0w.${c}.crank.map|cut -d\; -f2 | lsort 20G | uniq -c | lsort 1G -rn | head
  30477 /orta <orta.therox@gmail.com>
  25884  <jdelstrother@gmail.com>
  14368 "Franz Liedke" <franz@develophp.org>
  12864  <aneesh.kumar@gmail.com>
  11226 --global <tomas.vot@gmail.com>
  10869  <adam@adamralph.com>
  10720 #patcon <patrick.c.connolly@gmail.com>
  10510  <ash@kambanaria.org>
   9959  <k-okada@jsk.t.u-tokyo.ac.jp>
   8387 "Amir E. Aharoni" <amir.aharoni@mail.huji.ac.il>

zcat P2AFullS.nA2A1w.${c}.crank.map|cut -d\; -f2 | lsort 20G | uniq -c | lsort 1G -rn | head
 165409  <jdelstrother@gmail.com>
 142025 /orta <orta.therox@gmail.com>
  74593 #patcon <patrick.c.connolly@gmail.com>
  71038  <wangying@google.com>
  65840  <k-okada@jsk.t.u-tokyo.ac.jp>
  65459 --global <tomas.vot@gmail.com>
  60177  <adam@adamralph.com>
  57830  <Callek@gmail.com>
  55077 "Franz Liedke" <franz@develophp.org>
  50408  <reed@reedloden.com>


#############################################################################
#############################################################################
zcat a2PFullS.np2p | perl connectExportVw.perl a2PFullS.np2p
zcat P2aFullS.np2p | perl connectExportVw.perl P2aFullS.np2p
zcat P2aFullS.np2p.versions | paste  -d\ - <(zcat P2aFullS.np2p.weights) > bb
/home/audris/src/OSLOM2/oslom_undir -f bb -w

zcat a2PFullS.np2p.versions | paste  -d\  - <(zcat a2PFullS.np2p.weights) > cc
/home/audris/src/OSLOM2/oslom_undir -f cc -fast -w

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
zcat  c2pFullR.np2p1.PLM.crank.map | awk -F\; '{if ($2 != $1)print $1";"$2 }' | gzip >  c2pFullR.np2p1.PLMmap.forks


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
  
for i in {0..127}; do scp -p c2PFullR$i.s da3:/data/basemaps/gz/c2PPLMFullR$i.s; done


#now run updateCFBR.pbs
cd build
cmake ..
make
make install
cd ..
g++ -O3 -o cluster -I ./include/ -I ./extlibs/tlx/ -L ./build/ cluster.c -lnetworkit
g++ -O3 -o clusterw -I ./include/ -I ./extlibs/tlx/ -L ./build/ clusterw.c -lnetworkit


cat > top1 <<HERE
C;350
Cs;100
F;180
Go;350
JS;1400
Dart;160
Kotlin;400
PY;1700
R;40
Rust;40
Scala;640
TypeScript;1000
ipy;100
java;640
jl;100
pl;40
rb;160
HERE

for LA in Rust jl F Dart ipy pl Kotlin Scala Go; do zcat  /da0_data/play/${LA}thruMaps/b2cPtaPkgR${LA}.*.s | cut -d\; -f3,5 | lsort 100G -t\; -k1,2 | uniq -c | sed 's|^\s*||;s| |;|' |  gzip > P2aR$LA.s; done
zcat P2aRRust.s | grep -v '\[bot\]' | grep -Ev '(Rustu Robot|dependabot|travis-ci|RustyRobot)'| perl connectExportVwP2a.perl P2aRRust
cp -p P2aRRust.versions ~/src/networkit
zcat P2aRRust.weights > ~/src/networkit/w
time zcat P2aRRust.versions | ./clusterw 266641 290857 | gzip > P2aRRust.PLM
modularity=0.897247 nver=266641 clusters=205301 largest=2814
modularity=0.897268 nver=266641 clusters=205301 largest=3140
zcat ~/src/networkit/P2aR$LA.PLM | perl rank1.perl  P2aR$LA | gzip >  P2aR$LA.crank.map

zcat P2aR$LA.crank.map | grep '^[^<;]*;' | lsort 10G -t\; -k1,1 | join -t\; ChrisRust.P2p - > ChrisRust.P2p.joineda
cut -d\; -f3,4 ChrisRust.P2p.joineda | lsort 1G -t\| | uniq -c | sort -n | tail
      8 network-programming;Aceeri_glutin
     10 embedded;Vadim Kaushan <admin@disasm.info>
     10 network-programming;panicbit <panicbit.dev@gmail.com>
     11 embedded;0e4ef622_rust
     11 wasm;11Takanori_cranelift
     14 embedded;0mp_nix
     15 embedded;Borrmann <43264484+2ndTaleStudio@users.noreply.github.com>
     16 embedded;p3we_smoltcp
     36 network-programming;0x1997_hyper
     77 embedded;AssaultCuirass_stdsimd


zcat dist.crank.map | grep '^[^<;]*;' | lsort 10G -t\; -k1,1 | join -t\; ChrisRust.P2p - > ChrisRust.P2p.joinedE
     34 command-line-interface;Kixunil_parse_arg
      6 command-line-interface;alyssais_hyper-stub
      3 command-line-interface;japaric_lm3s6965
     30 command-line-interface;lucazulian_em7180
     17 command-line-interface;myrrlyn_calm_io
     16 command-line-interface;polachok_mio-afpacket
     99 command-line-interface;rv32m1-rust_rv32m1_ri5cy-hal
      6 command-line-interface;vi_rust-stunclient
     13 embedded;Kixunil_parse_arg
      7 embedded;alyssais_hyper-stub
      1 embedded;blm768_wasm-bindgen-console-logger
    219 embedded;eldruin_hdc20xx-rs
      4 embedded;emeraldpay_emerald-vault
    131 embedded;japaric_lm3s6965
     61 embedded;lucazulian_em7180
     19 embedded;myrrlyn_calm_io
      9 embedded;polachok_mio-afpacket
     44 embedded;rv32m1-rust_rv32m1_ri5cy-hal
      7 embedded;vi_rust-stunclient
     53 network-programming;Kixunil_parse_arg
    176 network-programming;alyssais_hyper-stub
      1 network-programming;blm768_wasm-bindgen-console-logger
      1 network-programming;eldruin_hdc20xx-rs
     10 network-programming;japaric_lm3s6965
     56 network-programming;lucazulian_em7180
      8 network-programming;myrrlyn_calm_io
     42 network-programming;polachok_mio-afpacket
     56 network-programming;rv32m1-rust_rv32m1_ri5cy-hal
     32 network-programming;vi_rust-stunclient
     11 wasm;Kixunil_parse_arg
      1 wasm;alyssais_hyper-stub
     30 wasm;blm768_wasm-bindgen-console-logger
     10 wasm;japaric_lm3s6965
      7 wasm;lucazulian_em7180
      5 wasm;myrrlyn_calm_io
      3 wasm;polachok_mio-afpacket
     16 wasm;rv32m1-rust_rv32m1_ri5cy-hal

     


zcat P2mRC  | awk -F\; '{n++;for (i=12; i<= 1000;i*=2){if (NF>i) c[i]++; } }END{for (i=12; i<= 1000;i*=2){print i,c[i]/n}}' 
12 0.421283
24 0.243881
48 0.132938
96 0.061204
192 0.0350636
384 0.0142883
768 0.00877073



zcat P2mRCs  | awk -F\; '{n++;for (i=10; i<= 160;i*=2){if (NF>i) c[i]++; } }END{for (i=10; i<= 160;i*=2){print i,c[i]/n}}' 
10 0.704816
20 0.382408
40 0.127819
80 0.0270784
160 0.00470162

zcat P2mRF  | awk -F\; '{n++;for (i=10; i<= 160;i*=2){if (NF>i) c[i]++; } }END{for (i=10; i<= 160;i*=2){print i,c[i]/n}}' 
10 0.19119
20 0.109169
40 0.055833
80 0.0293892
160 0.0151748

zcat P2mRGo  | awk -F\; '{n++;for (i=10; i<= 640;i*=2){if (NF>i) c[i]++; } }END{for (i=10; i<= 640;i*=2){print i,c[i]/n}}' 
10 0.574852
20 0.29509
40 0.136033
80 0.0656151
160 0.034571
320 0.0173776
640 0.00737377


zcat P2mRJS  | awk -F\; '{n++;for (i=10; i<= 1280;i*=2){if (NF>i) c[i]++; } }END{for (i=10; i<= 1280;i*=2){print i,c[i]/n}}' 
10 0.639941
20 0.527525
40 0.483607
80 0.43637
160 0.398042
640 0.262563
1280 0.0148644
1500;0.00323418

zcat P2mRDart  | awk -F\; '{n++;for (i=10; i<= 640;i*=2){if (NF>i) c[i]++; } }END{for (i=10; i<= 640;i*=2){print i,c[i]/n}}' 
10 0.439696
20 0.241551
40 0.113508
80 0.0325142
160 0.0111281
320 0.00478221
640 0.00215991


zcat P2mRKotlin  | awk -F\; '{n++;for (i=10; i<= 640;i*=2){if (NF>i) c[i]++; } }END{for (i=10; i<= 640;i*=2){print i,c[i]/n}}' 
10 0.697564
20 0.521439
40 0.338474
80 0.175538
160 0.0606591
320 0.0158832
640 0.00379658

zcat P2mRPY  | awk -F\; '{n++;for (i=10; i<= 640;i*=2){if (NF>i) c[i]++; } }END{for (i=10; i<= 640;i*=2){print i,c[i]/n}}' 
10 0.515181
20 0.292578
40 0.148944
80 0.086427
160 0.0475742
320 0.038749
640 0.0243697
zcat P2mRPY  | awk -F\; '{n++;for (i=1000; i<= 2000;i+=200){if (NF>i) c[i]++; } }END{for (i=1000; i<= 2000;i+=200){print i,c[i]/n}}' 
1200 0.0153664
1400 0.0115485
1600 0.00672301
1800 0.00460145


zcat P2mRR  | awk -F\; '{n++;for (i=10; i<= 160;i*=2){if (NF>i) c[i]++; } }END{for (i=10; i<= 160;i*=2){print i,c[i]/n}}' 
10 0.218833
20 0.0653949
40 0.0120985
80 0.00210577
160 0.000443218


zcat P2mRRust  | awk -F\; '{n++;for (i=10; i<= 160;i*=2){if (NF>i) c[i]++; } }END{for (i=10; i<= 160;i*=2){print i,c[i]/n}}' 
10 0.407686
20 0.21676
40 0.0943334
80 0.0358537
160 0.0122366



zcat P2mRScala  | awk -F\; '{n++;for (i=10; i<= 640;i*=2){if (NF>i) c[i]++; } }END{for (i=10; i<= 640;i*=2){print i,c[i]/n}}' 
10 0.602979
20 0.413818
40 0.215903
80 0.102263
160 0.0458566
320 0.0194
640 0.00934548

zcat P2mRTypeScript  | awk -F\; '{n++;for (i=10; i<= 1024;i*=2){if (NF>i) c[i]++; } }END{for (i=10; i<= 1280;i*=2){print i,c[i]/n}}' 
10 0.772294
20 0.64962
40 0.426289
80 0.223612
160 0.12694
320 0.0735866
640 0.0403654
1280 0


zcat P2mRipy  | awk -F\; '{n++;for (i=10; i<= 160;i*=2){if (NF>i) c[i]++; } }END{for (i=10; i<= 160;i*=2){print i,c[i]/n}}' 
10 0.529599
20 0.241883
40 0.073566
80 0.016091
160 0.00271202


zcat P2mRjava  | awk -F\; '{n++;for (i=10; i<= 1280;i*=2){if (NF>i) c[i]++; } }END{for (i=10; i<= 1280;i*=2){print i,c[i]/n}}' 
10 0.757259
20 0.576327
40 0.373395
80 0.19457
160 0.0845198
320 0.0323556
640 0.0115182
1280 0.00385263

zcat P2mRjl  | awk -F\; '{n++;for (i=10; i<= 160;i*=2){if (NF>i) c[i]++; } }END{for (i=10; i<= 160;i*=2){print i,c[i]/n}}' 
10 0.349495
20 0.164012
40 0.057526
80 0.0181676
160 0.00461887

zcat P2mRpl  | awk -F\; '{n++;for (i=10; i<= 160;i*=2){if (NF>i) c[i]++; } }END{for (i=10; i<= 160;i*=2){print i,c[i]/n}}' 
10 0.200477
20 0.0599965
40 0.0175393
80 0.00715389
160 0.00169251

zcat P2mRrb  | awk -F\; '{n++;for (i=10; i<= 160;i*=2){if (NF>i) c[i]++; } }END{for (i=10; i<= 160;i*=2){print i,c[i]/n}}' 
10 0.20421
20 0.0543675
40 0.0249583
80 0.0162167
160 0.0117355

cat top1 | while IFS=\; read LA c; do zcat P2mR$LA | perl connectExportAPI.perl API.$LA.$c $c 2> $LA; done


cat top1 | while IFS=\; read LA c; do
    nn=$(grep ' nodes ' $LA | sed 's|.*Has \([0-9]*\) nodes .*|\1|')
    cp -p API.$LA.$c.versions ~/src/networkit
    echo "time zcat API.$LA.$c.versions | ./cluster $nn | gzip > API.$LA.$c.PLM"
    echo "zcat ~/src/networkit/API.$LA.$c.PLM | perl rank1.perl  API.$LA.$c | gzip >  API.$LA.$c.crank.map &"
	done
	
time zcat API.C.350.versions | ./cluster 6189683 | gzip > API.C.350.PLM
modularity=0.490187 nver=6189683 clusters=2330 largest=1999964
modularity=0.428127 nver=6189683 clusters=2330 largest=2881862


time zcat API.Cs.100.versions | ./cluster 7946538 | gzip > API.Cs.100.PLM
modularity=0.474667 nver=7946538 clusters=6658 largest=1724044
modularity=0.230925 nver=7946538 clusters=6658 largest=6073960

time zcat API.F.180.versions | ./cluster 47078 | gzip > API.F.180.PLM
modularity=0.924190 nver=47078 clusters=1284 largest=3910
modularity=0.924191 nver=47078 clusters=1284 largest=4139

time zcat API.Go.350.versions | ./cluster 2197604 | gzip > API.Go.350.PLM
modularity=0.376564 nver=2197604 clusters=2359 largest=625332
modularity=0.034633 nver=2197604 clusters=2359 largest=1986055

zcat API.JS.1400.versions | ./cluster 7478528 | gzip > API.JS.1400.PLM
modularity=0.205556 nver=7478528 clusters=10199 largest=2880721
modularity=0.000982 nver=7478528 clusters=10199 largest=7249756

time zcat API.Dart.160.versions | ./cluster 1146105 | gzip > API.Dart.160.PLM
modularity=0.652711 nver=1146105 clusters=851 largest=115940
modularity=0.612332 nver=1146105 clusters=851 largest=315411

time zcat API.Kotlin.400.versions | ./cluster 5095621 | gzip > API.Kotlin.400.PLM
modularity=0.516462 nver=5095621 clusters=1859 largest=706239
modularity=0.316415 nver=5095621 clusters=1859 largest=3410469

time zcat API.PY.1400.versions | ./cluster 22790389 | gzip > API.PY.1400.PLM
modularity=0.552178 nver=22790389 clusters=14372 largest=5844473
modularity=0.540139 nver=22790389 clusters=14372 largest=5768252


time zcat API.R.40.versions | ./cluster 586050 | gzip > API.R.40.PLM
modularity=0.380193 nver=586050 clusters=2067 largest=85428
modularity=0.005970 nver=586050 clusters=2067 largest=574391

time zcat API.Rust.200.versions | ./cluster 735936 | gzip > API.Rust.200.PLM
modularity=0.530000 nver=735936 clusters=1719 largest=88391
modularity=0.262738 nver=735936 clusters=1719 largest=490340


time zcat API.Scala.640.versions | ./cluster 2262490 | gzip > API.Scala.640.PLM
modularity=0.650283 nver=2262490 clusters=1519 largest=308120
modularity=0.651281 nver=2262490 clusters=1519 largest=326006

time zcat API.ipy.100.versions | ./cluster 1706676 | gzip > API.ipy.100.PLM
modularity=0.394918 nver=1742647 clusters=4604 largest=321597
modularity=0.035180 nver=1742647 clusters=4604 largest=1604817

zcat API.TypeScript.1000.versions | ./cluster 18128034 | gzip > API.TypeScript.1000.PLM
modularity=0.694366 nver=20478811 clusters=12998 largest=5899462
modularity=0.690437 nver=20478811 clusters=12998 largest=6806967

zcat API.pl.40.versions | ./cluster 582172 | gzip > API.pl.40.PLM
modularity=0.362765 nver=582172 clusters=122 largest=151946
modularity=0.171523 nver=582172 clusters=122 largest=493236

time zcat API.jl.100.versions | ./cluster 133864 | gzip > API.jl.100.PLM
modularity=0.593737 nver=133864 clusters=531 largest=13705
modularity=0.448102 nver=133864 clusters=531 largest=60688

zcat API.rb.160.versions | ./cluster 3762768 | gzip > API.rb.160.PLM
modularity=0.603788 nver=3762768 clusters=22164 largest=935940
modularity=0.601084 nver=3762768 clusters=22164 largest=1028899


modularity=0.577709 nver=3988527 clusters=23096 largest=940544
modularity=0.573104 nver=3988527 clusters=23096 largest=1095478

zcat API.java.640.versions | ./cluster 67614364 | gzip > API.java.640.PLM
modularity=0.567923 nver=67614364 clusters=34122 largest=15148967
modularity=0.557721 nver=67614364 clusters=34122 largest=16488034


####
zcat P2mR* | perl connectExportAPI.perl API.200.100k 200 100000 
#finished at 1592064300 over 5584. Has 103037974 nodes and 403578554 edges
time zcat API.200.100k.versions | ./cluster 103037974 | gzip > API.200.100k.PLMPLP
modularity=0.823989 nver=103037974 clusters=6340393 largest=14049977
modularity=0.828344 nver=103037974 clusters=6340393 largest=14261157
zcat ~/src/networkit/API.200.100k.PLMPLP | grep -nv '^[0-9]' 
1:BEGIN 103037974
103037976:END
103037977:BEGIN 103037974
206075952:END
zcat ~/src/networkit/API.200.100k.PLMPLP | head -206075952  | tail -$((103037974+1)) | grep -v '^[BE]' | gzip >API.200.100k.PLM
zcat API.200.100k.PLM | perl rank1.perl API.200.100k | gzip > API.200.100k.PLM.crank.map&
zcat API.200.100k.PLM | perl rank2.perl API.200.100k | gzip > API.200.100k.PLM.cranka.map&


time zcat API.200.100k.versions | ./cluster 103037974 1000 | gzip > API.200.1k.PLMPLP
121886190 edges removed
modularity=0.844970 nver=103037974 clusters=9132909 largest=10811139
modularity=0.848516 nver=103037974 clusters=9132909 largest=10876603
zcat ~/src/networkit/API.200.1k.PLMPLP | grep -nv '^[0-9]'
1:BEGIN 103037974
103037976:END
103037977:BEGIN 103037974
206075952:END
zcat ~/src/networkit/API.200.1k.PLMPLP | head -206075952  | tail -$((103037974+1)) | grep -v '^[BE]' | gzip >API.200.1k.PLM
zcat API.200.1k.PLM | perl rank1.perl API.200.100k | gzip > API.200.1k.PLM.crank.map&
zcat API.200.1k.PLM | perl rank2.perl API.200.100k | gzip > API.200.1k.PLM.cranka.map&

(time  zcat API.200.100k.versions 10000 | ./cluster 103037974 10000 | gzip > API.200.10k.PLMPLP) &
80495877 edges removed
modularity=0.838909 nver=103037974 clusters=8243244 largest=13451105
modularity=0.842543 nver=103037974 clusters=8243244 largest=13558859
zcat ~/src/networkit/API.200.10k.PLMPLP | grep -nv '^[0-9]'
zcat ~/src/networkit/API.200.10k.PLMPLP | head -206075952  | tail -$((103037974+1)) | grep -v '^[BE]' | gzip >API.200.10k.PLM
zcat API.200.10k.PLM | perl rank1.perl API.200.100k | gzip > API.200.10k.PLM.crank.map&
zcat API.200.10k.PLM | perl rank2.perl API.200.100k | gzip > API.200.10k.PLM.cranka.map&

zcat P2mR* | perl connectExportAPI.perl API.200.1m 200 1000000
#finished at 1592065810 over 5728. Has 103037974 nodes and 684470377 edges
time zcat API.200.1m.versions | ./cluster 103037974 | gzip > API.200.1m.PLMPLP
modularity=0.769493 nver=103037974 clusters=1283446 largest=18068992
modularity=0.769585 nver=103037974 clusters=1283446 largest=18389658
zcat ~/src/networkit/API.200.1m.PLMPLP | grep -nv '^[0-9]'
1:BEGIN 103037974
103037976:END
103037977:BEGIN 103037974
206075952:END
zcat ~/src/networkit/API.200.1k.PLMPLP | head -206075952  | tail -$((103037974+1)) | grep -v '^[BE]' | gzip >API.200.1m.PLM
zcat API.200.1m.PLM | perl rank1.perl API.200.1m | gzip > API.200.1m.PLM.crank.map&
zcat API.200.1m.PLM | perl rank2.perl API.200.1m | gzip > API.200.1m.PLM.cranka.map&




cp -p API.versions ~/src/networkit/
time zcat API.versions | ./cluster 158704920 | gzip > API.PLMPLP
modularity=0.580687 nver=158704920 clusters=51449 largest=89868112
modularity=0.578611 nver=158704920 clusters=51449 largest=89926815

1:BEGIN 158704920
158704922:END
158704923:BEGIN 158704920
317409844:END
zcat ~/src/networkit/API.PLMPLP | head -317409844  | tail -$((158704920+1)) | grep -v '^[BE]' | gzip >API.PLM
zcat API.PLM | perl rank1.perl API | gzip > API.PLM.crank.map


zcat P2mR* | perl connectExportAPI.perl API.200 200
finished at 1591995814 over 25153. Has 92233040 nodes and 31783450 edges
time zcat API.200.versions | ./cluster 92233040 | gzip > API.200.PLMPLP
modularity=0.765188 nver=92233040 clusters=71689 largest=28099907
modularity=0.765183 nver=92233040 clusters=71689 largest=28414623
zcat API.200.PLMPLP | grep -nv '^[0-9]'                                                                          
1:BEGIN 92233040
92233042:END
92233043:BEGIN 92233040
184466084:END
zcat ~/src/networkit/API.200.PLMPLP | head -184466084  | tail -$(( 92233040+1)) | grep -v '^[BE]' | gzip >API.200.PLM

zcat API.200.PLM | perl rank1.perl API.200 | gzip > API.200.PLM.crank.map&
zcat API.200.PLM | perl rank2.perl API.200 | gzip > API.200.PLM.cranka.map&


# PFS is a complete list of blobs (Full) filtered only for relevant languages (Select)
zcat b2PFS$ver.np2p | perl connectExportVw1.perl b2PFS$ver.np2p2
finished at 1592145804 over 10398. Has 47144858 nodes and 166702437 edges
cp -p b2PFS$ver.np2p2.versions ~/src/networkit/
zcat b2PFS$ver.np2p2.weights > ~/src/networkit/w
(i=10; time zcat b2PFS$ver.np2p2.versions | ./clusterw 47144858 166702437 $i 2> b2PFS$ver.np2p2.PLM$i.err | gzip > b2PFS$ver.np2p2.PLM$i ) 
modularity=0.496117 nver=47144858 clusters=30823259 largest=3127060
modularity=0.366850 nver=47144858 clusters=30823259 largest=9777200
zcat ~/src/networkit/b2PFS$ver.np2p2.PLM10 | perl rank.perl b2PFSR.np2p2 | gzip > b2PFSR.np2p2.10.crank.map

zcat b2PFS$ver.np2p2.versions | ./clusterw 47144858 166702437 0 2> b2PFS$ver.np2p2.PLM.err | gzip > b2PFS$ver.np2p2.PLM
modularity=0.528003 nver=47144858 clusters=934377 largest=6036248
modularity=0.411121 nver=47144858 clusters=934377 largest=21467310
zcat ~/src/networkit/b2PFS$ver.np2p2.PLM | perl rank.perl b2PFSR.np2p2 | gzip > b2PFSR.np2p2.crank.map
zcat b2PFSR.np2p2.crank.map | lsort 100G -t\; -rn -k4 | head
zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzia_NetWorkTrans;0--_Spoon-Knife;0.000000;1.000002
zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzia_IntranetTransmission;0--_Spoon-Knife;0.000000;1.000002
zzzzzzzzzzzzz_Printbox;0--_Spoon-Knife;0.000000;1.000002
zzzzzzzzzzzx_MobileWeather;0--_Spoon-Knife;0.000000;1.000002

zcat b2PFSR.np2p2.crank.map | perl -ane 'print if /;cl[0-9]+;/' | lsort 100G -t\; -rn -k4 |head


zcat b2PS$ver.np2p | perl connectExportVw1.perl b2PS$ver.np2p2

for i in 100 50 10 2; do time zcat b2PS$ver.np2p2.versions | ./clusterw  21943048 29593008 $i 2> b2PS$ver.np2p2.PLMPLP$i.err | gzip > b2PS$ver.np2p2.PLMPLP$i; done

zcat ~/src/networkit/b2PSR.np2p2.PLMPLP100 | grep -nv '^[0-9]'                                                                                                                                                        
1:BEGIN 21943048
21943050:END
21943051:BEGIN 21943048
43886100:END
zcat ~/src/networkit/b2PSR.np2p2.PLMPLP100|head -43886100  | tail -$((21943048+1)) | grep -v '^[BE]' | gzip > b2PS$ver.np2p2.PLM100
zcat b2PSR.np2p2.PLM100 | perl rank.perl b2PSR.np2p2 | gzip > b2PSR.np2p2.PLM100.crank.map


zcat b2PSR.np2p2.PLMPLP10 | grep -nv '^[0-9]'
1:BEGIN 20908462
20908464:END
20908465:BEGIN 20908462
41816928:END
zcat ~/src/networkit/b2PSR.np2p2.PLMPLP10|head -41816927  | tail -$((41816928-20908465-1)) | grep -v '^[BE]' | gzip > b2PS$ver.np2p2.PLM10
zcat b2PSR.np2p2.PLM10 | perl rank.perl b2PSR.np2p2 | gzip > b2PSR.np2p2.PLM10.crank.map
zcat b2PSR.np2p2.PLM10 | cut -d\; -f1 | lsort 5G | uniq -c | lsort 1G -n | grep -n . | tail
15078220:  96017 63
15078221: 105027 268
15078222: 107246 47
15078223: 126875 18
15078224: 152734 41
15078225: 166159 55
15078226: 240304 0
15078227: 387769 43
15078228: 775905 39
15078229:2320293 40



zcat b2P$ver.np2p | perl connectExportVw1.perl b2P$ver.np2p2
starting at 1591706950
finished at 1591708132 over 1182. Has 84101124 nodes and 167859513 edges
#if ((g .degree (u) < 200 || g .degree (v) < 200) && w < 5){ remove edge
zcat b2P$ver.np2p2.versions | ./clusterw  84101124 167859513 | gzip > b2P$ver.np2p2.PLMPLP200
#if ((g .degree (u) < 20 || g .degree (v) < 20) && w < 5){ remove edge
zcat b2P$ver.np2p2.versions | ./clusterw  84101124 167859513 | gzip > b2P$ver.np2p2.PLMPLP20

#if ((g .degree (u) < 2 || g .degree (v) < 2) && w < 5){ remove edge
zcat b2P$ver.np2p2.versions | ./clusterw  84101124 167859513 | gzip > b2P$ver.np2p2.PLMPLP
...
removing edge 84061686 5553853 dg0=1 dg1=43 1.000000
removed 55132260 edges
[INFO ]: DynKatz: Sort time: 7886
[INFO ]: DynKatz: Sort time: 11758
[INFO ]: DynKatz: Sort time: 13282
[INFO ]: DynKatz: Sort time: 13427
calculated centrality
[WARN ]: move phase aborted after 32 iterations
[INFO ]: parallel coarsening took (54489 ms) 
[WARN ]: move phase aborted after 32 iterations
[INFO ]: parallel coarsening took (21369 ms) 
[INFO ]: parallel coarsening took (9640 ms) 
[INFO ]: parallel coarsening took (10657 ms) 
[INFO ]: parallel coarsening took (10537 ms) 
[INFO ]: parallel coarsening took (10687 ms) 
modularity=0.858558 nver=84101124 clusters=56405894 largest=4721776
modularity=0.863369 nver=84101124 clusters=56405894 largest=5149650
zcat b2PR.np2p2.PLMPLP | grep -nv '^[0-9]'
1:BEGIN 84101124
84101126:END
84101127:BEGIN 84101124
168202252:END
zcat ~/src/networkit/b2PR.np2p2.PLMPLP|head -168202251  | tail -$((168202252-84101127-1)) | grep -v '^[BE]' | gzip > b2P$ver.np2p2.PLM
zcat b2PR.np2p2.PLM | perl rank.perl b2PR.np2p2 | gzip > b2PR.np2p2.PLM.crank.map
zcat b2PR.np2p2.PLM | cut -d\; -f1 | lsort 5G | uniq -c | lsort 1G -n | grep -n . | tail
56405878: 841194 58
56405879: 876553 0
56405880: 947933 61
56405881:1180262 49
56405882:1222371 46
56405883:1296729 7
56405884:1634271 6
56405885:2232538 45
56405886:3598343 8
56405887:5149650 9

zcat c2pFullR.np2p1.PLM.crank.map| grep ^cran_ | cut -d\; -f2 | lsort 10G -u | gzip > cran
zcat b2PR.np2p2.PLM.crank.map|grep -v ';0\.0000' | cut -d\; -f2- | uniq -c | lsort 1g -rn | head                                                       
zcat b2PR.np2p2.PLM.crank.map | perl ~/bin/grepField.perl cran 1 | lsort 1G -t\; -rn -k3 | head

zcat ~/src/networkit/b2PR.np2p2.PLMPLP20|head -168202251  | tail -$((168202252-84101127-1)) | grep -v '^[BE]' | gzip > b2P$ver.np2p2.PLM20
zcat b2PR.np2p2.PLM20 | perl rank.perl b2PR.np2p2 | gzip > b2PR.np2p2.PLM20.crank.map


zcat b2P$ver.np2p | perl connectExportVw.perl b2P$ver.np2p
zcat b2P$ver.np2p.names | wc -l
84101124
cp -p b2P$ver.np2p.versions ~/src/networkit
cd ~/src/networkit
export LD_LIBRARY_PATH=/home/audris/lib:./build
zcat b2P$ver.np2p.versions| ./cluster  84101124 | gzip > b2P$ver.np2p.PLMPLP 
274076141 edges read
[INFO ]: DynKatz: Sort time: 8657
[INFO ]: DynKatz: Sort time: 19193
[INFO ]: DynKatz: Sort time: 18575
[INFO ]: DynKatz: Sort time: 16438
[INFO ]: DynKatz: Sort time: 16489
calculated centrality
[INFO ]: parallel coarsening took (38915 ms) 
[WARN ]: move phase aborted after 32 iterations
[INFO ]: parallel coarsening took (14473 ms) 
[INFO ]: parallel coarsening took (2307 ms) 
[INFO ]: parallel coarsening took (623 ms) 
[INFO ]: parallel coarsening took (418 ms) 
[INFO ]: parallel coarsening took (326 ms) 
[INFO ]: parallel coarsening took (342 ms) 
modularity=0.820662 nver=84101124 clusters=1278567 largest=5940014
modularity=0.826662 nver=84101124 clusters=1278567 largest=6133073
zcat b2PR.np2p.PLMPLP | grep -nv '^[0-9]'
1:BEGIN 84101124
84101126:END
84101127:BEGIN 84101124
168202252:END

zcat ~/src/networkit/b2PR.np2p.PLMPLP|head -168202251  | tail -$((168202252-84101127-1)) | grep -v '^[BE]' | gzip > b2P$ver.np2p.PLM
zcat b2PR.np2p.PLM | perl rank.perl b2PR.np2p | gzip > b2PR.np2p.PLM.crank.map
zcat b2PR.np2p.PLM | cut -d\; -f1 | lsort 5G | uniq -c | lsort 1G -n | grep -n . | tail
1278558:2251796 4
1278559:2376839 21
1278560:2795153 15
1278561:3058527 23
1278562:3486487 16
1278563:3546301 10
1278564:4288765 18
1278565:4340582 3
1278566:6097687 8
1278567:6133073 5

zcat b2PR.np2p.PLM.crank.map|grep -v ';0\.0000' | cut -d\; -f2- | uniq -c | lsort 1g -rn | head                                                       
2878355 0000-bigtree_0000-bigtree.github.com;0.343017
2816105 0101adm_cdnjs;0.326576
2167448 00jy116_spring-boot;0.013881
1682203 3a-classic_score;0.024765
1652141 Adude11_-tmp-100000-commit-2;0.021895
1476439 07101994_SignalR;0.006226
1271070 00ERNEST00_FFmpeg;0.134429
1241703 1900_mezzanine;0.020577
1174873 09saurabh09_vocal;0.022151
1129123 0-admin_tensorflow;0.021815


zcat c2pFullR.np2p1.PLM.crank.map| grep ^cran_ | cut -d\; -f2 | lsort 10G -u | gzip > cran
zcat b2PR.np2p.PLM.crank.map | perl ~/bin/grepField.perl cran 1 | lsort 1G -t\; -rn -k3 | head

cran_rviewgraph;0-T-0_ps4-linux;1.003315
cran_risaac;0-T-0_ps4-linux;1.003315
cran_gpuR;0-T-0_ps4-linux;1.003315
cran_ggconf;0-T-0_ps4-linux;1.003315
cran_csvread;0-T-0_ps4-linux;1.003315
cran_awspack;0-T-0_ps4-linux;1.003315
cran_XML2R;0-T-0_ps4-linux;1.003315
cran_TAQMNGR;0-T-0_ps4-linux;1.003315
cerebis_HiCseg;0-T-0_ps4-linux;1.003315
Cran_ssh.utils;0-T-0_ps4-linux;1.003315


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

#
cd ~/lookup/
~/src/igraph/compare 116236291 tst.m[12]
ADJRAND=0.484346
./compareSJ 116236291 ~/lookup/tst.m[12] &
./compareNMI 116236291 ~/lookup/tst.m[12]
SJ=15761578.000000
NMI=0.987770

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