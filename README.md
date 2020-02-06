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