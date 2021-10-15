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
#defork
ver=T
zcat c2pFull$ver.np2p.s | perl connectExportVw.perl c2pFull$ver.np2p
zcat c2pFull$ver.np2p.versions | paste  -d\   - <(zcat  c2pFull$ver.np2p.weights) | gzip > c2pFull$ver.np2p.vw
export LD_LIBRARY_PATH=/home/audris/lib64:/home/audris/lib:/home/audris/src/networkit/build
zcat c2pFull$ver.np2p.vw| /home/audris/src/networkit/clusterw $(zcat c2pFull$ver.np2p.names|wc -l) $(zcat /data/c2pFull$ver.np2p.vw|wc -l) -1 0 | gzip > c2pFull$ver.np2pw.PLM
11584040213 edges read
zcat c2pFull$ver.np2pw.PLM | perl rankNew.perl c2pFull$ver.np2p 1 | gzip >  c2pFull$ver.np2pw.crank.map
zcat  c2pFull$ver.np2pw.crank.map | awk -F\; '{if ($2 != $1)print $1";"$2 }' | gzip >  c2pFull$ver.np2pw.PLMmap.forks

zcat c2pFull$ver.np2p.versions | /home/audris/src/networkit/cluster $(zcat c2pFull$ver.np2p.names|wc -l)  | gzip > c2pFull$ver.np2pu.PLM
modularity=0.939244 nver=103612951 clusters=9966309 largest=625309
zcat c2pFull$ver.np2pu.PLM | perl rankNew.perl c2pFull$ver.np2p 1 | gzip >  c2pFull$ver.np2pu.crank.map
zcat  c2pFull$ver.np2pu.crank.map | awk -F\; '{if ($2 != $1)print $1";"$2 }' | gzip >  c2pFull$ver.np2pu.PLMmap.forks

#a2a
c=20 100 3000 30000
zcat P2AFullT.nA2A.$c.s | perl connectExportVw.perl P2AFullT.nA2AP.$c 
zcat P2AFullT.nA2AP.$c.versions | /home/audris/src/networkit/cluster $(zcat P2AFullT.nA2AP.$c.names|wc -l) | gzip > P2AFullT.nA2APu.$c.PLM
paste -d\  <(zcat P2AFullT.nA2AP.$c.versions) <(zcat P2AFullT.nA2AP.$c.weights) | /home/audris/src/networkit/clusterw $(zcat P2AFullT.nA2AP.$c.names|wc -l) $(zcat P2AFullT.nA2AP.$c.weights|wc -l) | gzip > P2AFullT.nA2APw.$c.PLM
zcat P2AFullT.nA2APw.$c.PLM | perl rankNew.perl P2AFullT.nA2AP.$c 1 | gzip >  P2AFullT.nA2APw.$c.crank.map
zcat P2AFullT.nA2APu.$c.PLM | perl rankNew.perl P2AFullT.nA2AP.$c 1 | gzip >  P2AFullT.nA2APu.$c.crank.map
#20u -    modularity=0.985762 nver=30763969 clusters=3914290 largest=151854
#100u -   modularity=0.921918 nver=36981227 clusters=3339331 largest=2458058
#3000u -  modularity=0.875331 nver=38759714 clusters=3201761 largest=3384986
#30000u - modularity=0.868881 nver=39372059 clusters=3123334 largest=2997349
#30000w - modularity=0.893545 nver=39372059 clusters=3112555 largest=2772919
#3000w -  modularity=0.899088 nver=38759714 clusters=3213922 largest=3633416
#100w -   modularity=0.938189 nver=36981227 clusters=3342010 largest=2516262
#20w -    modularity=0.988562 nver=30763969 clusters=3914359 largest=138114

#
c=20
zcat P2AFullT.nA2APw.$c.crank.map | cut -d\; -f2 | lsort 20G | uniq -c | lsort 1G -rn | head -20
  81221 YSR <cumulonimbus96@gmail.com>
  76825 shinara03 <shinara03@gmail.com>
  67693 attarwal <attarwal@UTORONTO.CA>
  62913 Nick <nickinpie@gmail.com>
  51669 weitzman <weitzman>
  49287 tanya1810 <tanya.sood.311@gmail.com>
  48050 cosa65 <mfosc65@gmail.com>
  47841 lassehh <sales@scrapeitout.com>
  46528 Tuan <atuannguyen1101@gmail.com>
  42824 olavoie <olavoie9507@gmail.com>
  42063 tknappek <tknappek@localhost>
  39508 ZackStr <zack.strulovitch@itential.com>
  39275 Unknown <lyndond@uci.edu>
  37938 Ileena Mitra <ileenam@surveymonkey.com>
  36246 Jan Dageförde <jan@dagefor.de>
  34547 tmathern <60901087+tmathern@users.noreply.github.com>
  33485 nlbtdsdso <53069055+nlbtdsdso@users.noreply.github.com>
  32454 phucmh-1356 <49380498+phucmh-1356@users.noreply.github.com>
  31788 Siff <siff.ravn@gmail.com>
  30503 Jura Laakkonen <jura.laakkonen@hiq.fi>

zcat P2AFullT.nA2AP.$c.versions|ssh da3 "bin/connect" | gzip > P2AFullT.nA2AP.$c.clones
zcat P2AFullT.nA2AP.$c.clones | cut -d\; -f2 | lsort 30G | uniq -c | lsort 1G -rn | head -3
13371009 0
    445 383009
    418 676334
    407 91410

c=20    
c=30000
zcat P2AFullT.nA2APw.$c.crank.map | cut -d\; -f2 | lsort 20G | uniq -c | lsort 1G -rn | head -20
1772917 onovy <novy@ondrej.org>
1401109 Tomster <tomster@emberjs.com>
 412298 Force User <lvh3851@g.rit.edu>
 332934 Chris Pitt <chris@silverstripe.com>
 282566 beckermr <becker.mr@gmail.com>
 191827 sc82choi <blissray@gmail.com>
 177113 Anton <verybigbro@gmail.com>
 170220 Nitesh <nitesh.turaga@gmail.com>
 169528 Justin <jflory7@gmail.com>
 159443 a <562826179@qq.com>
 147677 mandark <A06438_P5.Training@accenture.com>
 133136 olprod <olprod@microsoft.com>
 129512 Max Wofford <max@hackedu.us>
 126652 llmel <Aluno@Digital.net>
 120815 Hector <hectorsector@github.com>
 119375 t-pot <imagire@gmail.com>
 109305 MichaelDimmitt <MichaelGDimmitt@gmail.com>
 100329 Federico Aloi <fede@mumuki.org>
  97772 dlockhart <dlockhart@gmail.com>
  92721 Adam Jonas <jonas@chaincode.com>

c=1000  
zcat P2PFullT.nb2b.$c.s | perl connectExportVw.perl P2PFullT.nb2b.$c 
zcat P2PFullT.nb2b.$c.versions | paste  -d\   - <(zcat P2PFullT.nb2b.$c.weights) | /home/audris/src/networkit/clusterw $(zcat P2PFullT.nb2b.$c.names|wc -l) $(zcat P2PFullT.nb2b.$c.weights|wc -l) | gzip > P2PFullT.nb2bw.$c.PLM
#modularity=0.895951 nver=51640553 clusters=1515055 largest=4213396
zcat P2PFullT.nb2bw.$c.PLM | perl rankNew.perl P2PFullT.nb2b.$c 1 | gzip >  P2PFullT.nb2bw.$c.crank.map
zcat P2PFullT.nb2bw.$c.crank.map | cut -d\; -f2 | lsort 20G | uniq -c | lsort 1G -rn | head -20
2097991 cdnjs_cdnjs
1738972 fuel_core
1539576 z-song_laravel-admin
1212948 DasyDong_DasyDong.github.io
1033868 elastic_kibana
 862768 tensorflow_tensorflow
 790681 angular_angular
 563137 TelegramMessenger_Telegram-iOS
 547916 BVLC_caffe
 530790 cwoolcott_cwoolcott.github.io
 504309 buggins_coolreader
 421297 AmartyaU_WordPress
 417589 facebook_react
 415444 rafaelzolt_HugoQuickStart
 398238 cocos2d_cocos2d-x
 389486 openresty_lua-nginx-module
 340091 drupal_drupal
 328303 ARMmbed_mbed-os
 292958 thabataganga_idook
 276958 cvet_fonline

zcat P2PFullT.nb2b.$c.versions|ssh da3 "bin/connect" | gzip > P2PFullT.nb2b.$c.clones
zcat P2PFullT.nb2b.$c.clones | cut -d\; -f2 | lsort 30G | uniq -c | lsort 1G -rn | head -20
zcat P2PFullT.nb2b.$c.versions | paste  -d\   - <(zcat P2PFullT.nb2b.$c.weights) > /data/play/forks/P2PFullT.nb2b.$c.csv
perl getComponent.perl P2PFullT.nb2b.$c P2PFullT.nb2bw.$c.PLM 0 > P2PFullT.nb2bw.$c.0.forOSLO
time /home/audris/src/OSLOM2/oslom_undir -f /data/play/forks/P2PFullT.nb2b.$c.0.csv -hint  /data/play/forks/P2PFullT.nb2bw.$c.0.forOSLO -w -r 1 -hr 1

c=10
modularity=0.975464 nver=35590319 clusters=2401011 largest=1766606
 773837 cdnjs_cdnjs
 507025 opf_openproject
 359456 z-song_laravel-admin
 225837 magicalpanda_MagicalRecord
 223712 AmartyaU_WordPress
 174380 JavascriptBootcamp_group4
 151064 ARMmbed_mbed-os
 149514 mopub_mopub-android-sdk
 106122 kubernetes_kubernetes
 105811 kenmick_WebCrawler
 102712 cndn_intelligent-code-completion
 102279 hrony_GitHub
  99726 angular_code.angularjs.org
  99326 kadirahq_paper-ui
  96991 symfony_symfony
  96429 barryclark_jekyll-now
  94765 laravel_laravel
  91059 BVLC_caffe
  78939 rdpeng_ProgrammingAssignment2
  71513 Sable_mcbench-benchmarks
# too long /home/audris/src/OSLOM2/oslom_undir -f /data/play/forks/P2PFullT.nb2b.$c.csv -w -r 1 -hr 1

base=A2AFullT.nfb.500
zcat $base.s | perl connectExportVw.perl $base
zcat $base.versions | paste  -d\   - <(zcat $base.weights) | /home/audris/src/networkit/clusterw $(zcat $base.names|wc -l) $(zcat $base.weights|wc -l) | gzip > $base.PLM
#modularity=0.733959 nver=452306359 clusters=13681 largest=87041386
zcat $base.PLM | perl rankNew.perl $base 1 | gzip >  $base.crank.map
5181010 grvcodes <gurv0002@gmail.com>
4308454 Taiwo <ayeolakenny@gmail.com>
 970561 Douglas Lovell <dclo@us.ibm.com>
 891289 adedayo2017 <maxistinlove@gmail.com>
 577967 Andreas Amsenius <andreas@amsenius.se>
 471961 Raphael <raphael.rb96@gmail.com>
 413288 Dacio <dacioromero@gmail.com>
 406610 Nathanaël <nathanael.spriet@gmail.com>
 388238 hadley <h.wickham@gmail.com>
 335966 Presto! <presto.core@yandex.com>

base=t2PFullT.np2p
zcat $base.s | perl connectExportVw.perl $base
zcat $base.versions | paste  -d\   - <(zcat $base.weights) | /home/audris/src/networkit/clusterw $(zcat $base.names|wc -l) $(zcat $base.weights|wc -l) | gzip > ${base}w.PLM
#modularity=0.980955 nver=9538003 clusters=2001208 largest=69173
zcat ${base}w.PLM | perl rankNew.perl ${base} 1 | gzip >  ${base}w.crank.map
zcat ${base}w.crank.map | cut -d\; -f2 | lsort 20G | uniq -c | lsort 1G -rn | head -20
  50776 barryclark_jekyll-now
  13369 zuihou_zuihou-admin-cloud
  11861 gitlab.com_claasaug_example-long-term-project
   7030 Jackeagle_kernel_msm-3.18
   6355 TheMattSykes_personal-website
   6208 dcherif_SimpleDevopsProject
   5270 top-think_framework
   5206 rdpeng_ProgrammingAssignment2
   5019 frioux_dotfiles
   4489 laravel_laravel
   3882 forezp_SpringcloudConfig
   3744 gitx_gitx
   3178 git.savannah.gnu.org_git_mediagoblin
   2612 STRYDER-007_aosp_development_aospcaf
   2513 agrosner_DBFlow
   2299 videojs_videojs-contrib-hls
   2283 getlantern_lantern
   2211 thomas-moulard_test_repo
   2065 LhadyChloe_coursera-test
   2038 mad-science_6178

zcat $base.versions | paste  -d\   - <(zcat $base.weights) > $base.csv
time /home/audris/src/OSLOM2/oslom_undir -f $base.csv -w -r 1 -hr 1
real    116m13.715s
cat $base.csv_oslo_files/tp | perl rankNewO.perl $base ${base}w.PLM | gzip > ${base}ow.crank.map
zcat ${base}ow.crank.map | cut -d\; -f2 | lsort 20G | uniq -c | lsort 1G -rn | head -20
  18455 barryclark_jekyll-now
   5720 Jackeagle_kernel_msm-3.18
   5443 TheMattSykes_personal-website
   3435 frioux_dotfiles
   3192 gitlab.com_gitlab-com_www-gitlab-com
   3041 mmistakes_minimal-mistakes
   2750 netlify-templates_gatsby-starter-netlify-cms
   2011 STRYDER-007_aosp_development_aospcaf
   1963 ARMmbed_DAPLink
   1588 vinnymac_vinnymac.github.io
   1587 rdpeng_ProgrammingAssignment2
   1570 thomas-moulard_test_repo
   1506 getlantern_lantern
   1430 JohanSmet_udacity_carnd
   1381 labzero_lunch
   1377 symfony_symfony-standard
   1308 carolmanderson_academic-kickstart
   1205 laravel_laravel
   1183 daattali_beautiful-jekyll
   1144 DevMountain_learn-git

zcat ${base}w.crank.map |grep cran_Epi
cran_Epi;cran_Epi;0.000074
salsa.debian.org_r-pkg-team/r-cran-epi;cran_Epi;0.000074

zcat ${base}ow.crank.map |grep cran_Epi
cran_Epi;cran_Epi


base=T2TFullT.nb2b.20
zcat $base.s | perl connectExportVw.perl $base
zcat $base.versions | paste  -d\   - <(zcat $base.weights) | /home/audris/src/networkit/clusterw $(zcat $base.names|wc -l) $(zcat $base.weights|wc -l) | gzip > ${base}w.PLM
#modularity=0.953600 nver=37981255 clusters=1923369 largest=2158481
zcat ${base}w.PLM | perl rankNew.perl ${base} 1 | gzip >  ${base}w.crank.map
zcat ${base}w.crank.map | cut -d\; -f2 | lsort 20G | uniq -c | lsort 1G -rn | head -20
 952938 JettBurns14_super-super-heroku-bot
 647035 spree_spree
 567755 laravel_laravel
 349635 JulCesMelPin_AngularJS_starters
 338295 frioux_dotfiles
 267115 roguedream_pyc_source
 258407 symfony_symfony-standard
 242865 barryclark_jekyll-now
 192294 EsotericSoftware_spine-runtimes
 190068 klightspeed_Arduino-Libraries
 180382 marmelab_react-admin
 153427 getlantern_lantern
 138478 expo_expo
 131086 hrony_GitHub
 128147 kenmick_WebCrawler
 127551 ant-design_ant-design-pro
 122529 cyrsis_LearnCSharp
 121266 udacity_Sunshine-Version-2
 113158 BVLC_caffe
 111427 git.savannah.gnu.org_git_mediagoblin


zcat $base.versions | paste  -d\   - <(zcat $base.weights) > $base.csv
time /home/audris/src/OSLOM2/oslom_undir -f $base.csv -w -r 1 -hr 1
#too long
base=T2TFullT.nb2b.10
#modularity=0.966172 nver=34231127 clusters=2154971 largest=1494825
 644262 JettBurns14_super-super-heroku-bot
 490972 opf_openproject
 327855 spatie_laravel-feed
 226062 barryclark_jekyll-now
 223753 AFNetworking_AFNetworking
 216180 symfony_symfony-standard
 165154 VagTsop_tutorials-project-collection
 153475 EsotericSoftware_spine-runtimes
 137658 laravel_laravel
 131385 getlantern_lantern
 114050 frioux_dotfiles
 102181 kadirahq_paper-ui
  99238 kenmick_WebCrawler
  87708 micropython_micropython
  83033 hrony_GitHub
  78160 klightspeed_Arduino-Libraries
  78127 imfht_phpsourcecode
  76594 BitBagCommerce_SyliusCmsPlugin
  75357 roguedream_pyc_source
  71051 git.savannah.gnu.org_git_mediagoblin
base=T2TFullT.nb2b.5
#modularity=0.977951 nver=29384072 clusters=2522531 largest=595628
 227979 JettBurns14_super-super-heroku-bot
 217688 SocialiteProviders_Providers
 173142 symfony_symfony-standard
 152170 AFNetworking_AFNetworking
 139473 barryclark_jekyll-now
  96758 Jackeagle_kernel_msm-3.18
  86353 frioux_dotfiles
  84352 ErikMoczi_packages.unity.com
  82328 kaxap_arl
  81981 klightspeed_Arduino-Libraries
  74227 laravel_laravel
  73472 roguedream_pyc_source
  66570 axetroy_greasyfork.py
  64380 git.savannah.gnu.org_git_mediagoblin
  56736 Mattlk13_metacpan-cpan-extracted
  51584 BitBagCommerce_SyliusCmsPlugin
  50501 Wayne-Bai_AST-normalize
  49766 surli_backup-experiments
  49689 cndn_intelligent-code-completion
base=T2TFullT.nb2b.2
#modularity=0.992078 nver=19445607 clusters=3290435 largest=197057
  81888 symfony_symfony-standard
  74347 roguedream_pyc_source
  62663 barryclark_jekyll-now
  57550 JettBurns14_super-super-heroku-bot
  49212 Jackeagle_kernel_msm-3.18
  40875 frioux_dotfiles
  39900 opf_openproject
  39737 laravel_laravel
  37616 getlantern_lantern
  30575 git.savannah.gnu.org_git_mediagoblin
  26581 drupal_drupal
  24517 kaxap_arl
  23962 surli_backup-experiments
  23874 bitcoin_bitcoin
  20261 Homebrew_homebrew-php
  18721 kenmick_WebCrawler
  18714 klightspeed_Arduino-Libraries
  17413 TheMattSykes_personal-website
  15307 Mattlk13_metacpan-cpan-extracted
  15053 zhangkn_Demos4iOS

#determine infrastructure repos
zcat Dl2PfFullT.s |cut -d\; -f1,2|uniq|sed 's/;/|/;s/|Java/|java/' | gzip > Dl2nPT.g
lsort 3G -t\; -k1,1 Infra.100 |cut -d\; -f1 | gzip > Infra.100.n
zcat Dl2PfFullT.s | sed 's/;/|/;s/|Java/|java/' | ~/lookup/grepField.perl  Infra.100.n 1 |gzip > Infra.100.Def
zcat Infra.100.Def|cut -d\; -f1,2|uniq|lsort 3G -t\; -k1,2 -u|cut -d\; -f1|uniq -c|perl -ane 'chop();s/^\s*//;s| |;|;($c,$p)=split(/;/);print "$p;$c\n"'|lsort 3G -t\; -k1,1 | gzip > Infra.100.Def.c
zcat Infra.100.Def.c|join -t\; - <(lsort 3G -t\; -k1,1 Infra.100) |sed 's/|/;/' > Infra.100.DepDef.c



zcat Pkg2lPS{[0-9],[1-3][0-9]}.s| cut -d\; -f1,2 |uniq -c|awk '{if ($1>100) print $0}' |sed 's|^\s*||;s| |;|'| gzip > InfraPackages
zcat InfraPackages | perl -ane '$str=$_;@x=split(/;/,$str);$x[1]=~s/\.\*$//;@y=split(/\./, $x[1]);print "$y[$#y];$str"' | ~/lookup/splitSecCh.perl InfraPackages. 128
for i in {0..127}; do zcat InfraPackages.$i.gz | lsort 1G -t\; -k1,1 | join -t\; - <(zcat export2PFullT$i.s|lsort 2G -t\; -k1,1); done | gzip > InfraPackages.joined
zcat InfraPackages.joined|perl -ane 'chop();@x=split(/;/);@a=split(/\./, $x[2]);pop @a; @b=split(/\./, $x[6]); print "$_\n" if $b[$#b] eq $a[$#a];' |gzip > InfraPackages.clean
#for java
zcat  InfraPackages.clean|awk -F\; '{split($3,a,".");split($7,b,".");if ($2>100 && $4 == $5 && a[length(a)-1]== b[length(b)]){print $0}}'

zcat  InfraPackages.clean|awk -F\; '{if ($2>100 && $4 == $5 && $4 != "java"){print $0}}'
##### Highly shared blobs in: over 1006533 commits
ls /fast/b2cFullT.*.tch.large.*|sed 's|.*tch.large.||'|~/lookup/getValues b2f | awk -F\; '{print $1";"NF";"$NF}' > CommonBlobs
ls -l /fast/b2cFullT.*.tch.large.*|sed 's|.*da \s*||;s| .*tch.large.|;|'|while IFS=\; read s b; do echo $b";"$((($s-20)/20)); done | sort -t\; -k2 -n > CommonBlobs.nc
tail -20 CommonBlobs.nc

dea3013d6710ee273f49ac606a65d5211d480c88;1788867  ISC License
65c5ca88a67c30becee01c5a8816d964b03862f9;1796186  lgplv3
a11777cc471a4344702741ab1c8a588998b1311a;1926624  favicon.ico
6b60c1042f58d9fabb75485aa3624dddcf633b5c;1994195  some svg
dfe0770424b2a19faf507a501ebfc23be8f54e7b;2044669  # Auto detect text files and perform LF normalization\n* text=auto
4d29575de80483b005c29bfcac5061cd2f45313e;2089446  .gitignore
e7b4def49cb53d9aa04228dd3edb14c9e635e003;2123200  include ':app'
7f68460d8b38ac04e3a3224d7c79ef719b1991a9;2374448  RunConfigurationProducerService
94a9ed024d3859793618152ea559a168bbcbb5e2;2415767  gplv3
3c3629e647f5ddf82548912e337bea9826b434af;2430398  node_modules
94fb5490a2ed10b2c69a4a567a4fd2e4f706d841;2546034   glyphicons-halflings-regular.svg
5008ddfcf53c02e82d7eee2e57c38e5672ef89f6;2551838  .DS_Store
b93a4953fff68df523aa7656497ee339d6026d64;2638539  glyphicons-halflings-regular.eot
1413fc609ab6f21774de0cb7e01360095584f65b;2655492  glyphicons-halflings-regular.ttf
9e612858f802245ddcbf59788a0db942224bab35;2672427  glyphicons-halflings-regular.woff
64539b54c3751a6d9adb44c8e3a45ba5a73b77f0;2680048  glyphicons-halflings-regular.woff2
796b96d1c402326528b4ba3c12ee9d92d0e212e9;2891784  "/build"
94a25f7f4cb416c083d265558da75d457237d671;3470163  "xml: VcsDirectoryMappings
8b137891791fe96927ad78e64b0aad7bded08bdc;7290209  "\n"
e69de29bb2d1d6434b8b29ae775ad8c2e48c5391;57126599 ""

https://e.chase.com/T/v6000001799f22108cb0e7bb6e966a3578/df4b90a9c38744eb0000021ef3a0bcc6/df4b90a9-c387-44eb-83fd-26fe8aadead8?__F__=v0fUYvjHMDjRPMSh3tviDHXIoXcPxvDgUUCCPvXMWoX_0JoZLAZABQFwlrGt0pA4PxIRrWOgZzYHN8dssdhH2aKq7iSryUic22
Fwr0QNuBHaGOUH_rsD-9Rzgp--UY3-GIOa0KIJFLod2KlNDMqNjOov-UaS_fzDe-9A459fRGxYSDIOq8-zdqCj9CuFNR5NWSKgPJ4pKiON7GreOoL_aso0pIvl4vIllT2yxNEFTi_gNoX6V1ZrJaj-UtuVJfjRQqac7C7l1fJ1QgzGmAJfb6LlIzKPp-67hQUw1-0dG1L5n8Q1-ZLYRzTAIAt-JzFc6Sc55oxGKhFK4uBS5LWWYTo1lNvg4VNiQ
0neHPiFp2BlbQDqsvzMvdGTIDcvazpTbuasaUlaZ7IAeglJTN3ryBEC4GlBX_rJDN_fin-3hq9B2Vlj9Ue4vOe0sTy6iI92fugSZ3IWs8QGs=


zcat  largeb2PT.s|lsort 100G -t\; -rn -k2 |head -20
e69de29bb2d1d6434b8b29ae775ad8c2e48c5391;22557224 ""
8b137891791fe96927ad78e64b0aad7bded08bdc;3960138  "\n"
94a25f7f4cb416c083d265558da75d457237d671;3053308  
796b96d1c402326528b4ba3c12ee9d92d0e212e9;2333442
7f68460d8b38ac04e3a3224d7c79ef719b1991a9;2042698
dfe0770424b2a19faf507a501ebfc23be8f54e7b;1963691
64539b54c3751a6d9adb44c8e3a45ba5a73b77f0;1921635
9e612858f802245ddcbf59788a0db942224bab35;1916190
1413fc609ab6f21774de0cb7e01360095584f65b;1904289
b93a4953fff68df523aa7656497ee339d6026d64;1893977
94fb5490a2ed10b2c69a4a567a4fd2e4f706d841;1827139
e7b4def49cb53d9aa04228dd3edb14c9e635e003;1814410
4d29575de80483b005c29bfcac5061cd2f45313e;1764745
3c3629e647f5ddf82548912e337bea9826b434af;1761602
5008ddfcf53c02e82d7eee2e57c38e5672ef89f6;1703087
6b60c1042f58d9fabb75485aa3624dddcf633b5c;1652127
cccdd3d517fc5249beaefa600691cf150f2fa3e6;1489139
a11777cc471a4344702741ab1c8a588998b1311a;1452271
59385cdf379bd06a8d2326dcd4de6d5cd5d3f5b0;1444342
b1c56658557b8162aa9f5ba8610ed03a5e558d9d;1365816

#most licenses
zcat ../c2fb/P2LFullT{[0-9],[1-3][0-9]}.s | cut -d\; -f1 | uniq -c | lsort 3G -rn | head
  66247 pombredanne_license-detection-test-data
   7761 barryclark_jekyll-now
   6694 Mattlk13_metacpan-cpan-extracted
   5245 TheMattSykes_personal-website
   5239 gitlab.com_atoomic_CPAN
   3693 caolan_async
   3356 gitlab.com_gitlab-com_www-gitlab-com
   3205 carolmanderson_academic-kickstart
   3106 videojs_videojs-contrib-hls
   3085 drupal_drupal
#Most common licenses
zcat ../c2fb/bL2nPFullT.s|lsort 3G -t\; -k2 -rn |head
e69de29bb2d1d6434b8b29ae775ad8c2e48c5391;22557224 # empty
8b137891791fe96927ad78e64b0aad7bded08bdc;3960138  # \n
796b96d1c402326528b4ba3c12ee9d92d0e212e9;2333442  # /build
dea3013d6710ee273f49ac606a65d5211d480c88;1364470  #ISC
a7ae8ee9b8a30ef2a73ff5a7a80adc3b1a845cae;1248566  #MIT
d4569487a094b9ffb6e42ed157c32f8a5440a07a;1206407  #Nathan LaFreniere
06166077be4d1f620d89b9eb33c76d89e75857da;1178849  #MIT
654d0bfe943437d43242325b1fbcff5f400d84ee;1161826  #MIT
0c068ceecbd48fc4e8279e6451793fec2bf12178;1142651  #MIT
19129e315fe593965a2fdd50ec0d1253bcbd2ece;1130383  #ISC

   

zcat ${base}w.PLM | perl AnnoteCent.perl ${base} | gzip > $base.annotated


perl getType.perl | gzip > P2fP4FullT.s2
zcat P2fP4FullT.s2 | awk -F\; '{if ($1==$2)print $0}'|cut -d\; -f1,3,5,7,9,11,13,14,19,20|gzip > P2fP4FullT.Pstats
zcat P2fP4FullT.Pstats|cut -d\; -f6|lsort 10G |uniq -c
8209126 d
7872040 d70
3829321 d90
60499607 o
9040410 u

zcat P2fP4FullT.s2 | awk -F\; '{if ($1!=$2&&$15/$13 > 0.5 && $15/$17 > .5 && $15>=$16 && $15>1)print $0}' |gzip > P2fP4FullT.StrongReuseLinks05-1

base=P2fP4FullT.StrongReuseLinks05-1
zcat /data/basemaps/gz/$base |  cut -d\; -f1,2  | perl ~/bin/connectBasic.perl $base |gzip > $base.map

zcat $base.versions | /home/audris/src/networkit/cluster $(zcat $base.names|wc -l) | gzip > $base.PLM
#modularity=0.988932 nver=1442864 clusters=387376 largest=70422

zcat $base.PLM | perl rankNew.perl $base 1 | gzip > $base.crank.map

zcat P2fP4FullT.s2 | awk -F\; '{if ($1!=$2 && $15/$13 > 0.1 && $15/$17 > .1 && $15>=$16 && $15>10)print $0}' |gzip > P2fP4FullT.StrongReuseLinks01-10

###########################################################
#time-based blob-sharing


We aim to measure the extent of code sharing at the entire FLOSS
ecosystem level and elicit the reasons for such sharing activity.


We find over 97% of repositories of nontrivial size to share some
blobs and over 10\% percent of the repositories to be 100\% derived.
From among 45M repositories containing over 20 blobs, 9M consist
over 90\%, 12.6M over 70\% and 6M over 50\% of original blobs.

The sharing of blobs among repositories forms a complex network with
over X\% of the repositories sharing blobs with five or more other
repositories and 82% of all repositories formin one giant component that
can be traversed via shared blobs.

The investigation of the nature of code sharing identified several
principal reasons behind wide-scale sharing:

The most widely shared blobs have simple or default content as, for example, and empty
file with just a newline character is shared among 22M repositories.
%Such empty/simple files may not even be intentianally shared but may be created independently.

Once such widely shared artifacts are removed from the blob-project bi-graph
(by considering only the blobs shared among fewer than 100 reopositories)

- The collections of software and artifacts represent
the largest networks of shared source code. For example,
roguedream/pyc_source, gitlab.com/atoomic/CPAN,
cdnjs/cdnjs.

- The next largest network is created by code sharing
among linux kernel repositories (Jackeagle/kernel_msm-3.18).

- The very large collection of SO snippets used to train
DL methods for programming language and bug detection
(aliostad/deep-learning-lang-detection and imanchatterjee/BugDetection) represents another type of
collection where software code is treated as data for DL methods.

- Anoter type of very large network was build bot for chromium that
would get source code from numerous other projects used by chromium
browswer (bloomberg/chromium.bb).

Another type of large-scale sharing involves software collections
curated by well-known individuals and organizations (hadley_cran-packages)

- Reasons for small-scale (few repositories involved) sharing
iiclude code sharing between the primary repository and copy-then-commit
mirrors or forks. For example, the downstream versions of
packages in Debian distribution are not forks of the upstream
repositories. Instead they copy and commit code from the upstream
projects and make minor modifications to ensure that all projects in the
distribution can be copmiled and deployed.

- Artifacts from an author that are related to multiple repos as code
examples from a book in Cran_vars and cran_urca




Only 4.7M of the total of 90M repositories have
such/similar mirrors reducing this 4.7M down to 1.6M of non-duplicate
repositories.


- Yet another drupal_drupal

ver=T
for sel in SL L S
do (time perl getType1.perl $sel |gzip > fP2${sel}PStats$ver.s)
   zcat fP2${sel}PStatsT.s | cut -d\; -f1-6 | uniq | gzip > fP2${sel}PStatsT.Pstats
done


#L
Largest by in-degree
zcat /data/basemaps/gz/$base.s | awk -F\; '{ if ($3 > i){ i=$3; print $0} }'
pombredanne_license-detection-test-data;d;60833;258;64550;480;1;shyler1987_umnenie_0;o70;136;6;4083;3005;0
x0rzkov_dockerfiles-search;o50;156746;497;500244;283834;1;ssweriduk_wordpress-alpine-nginx-scalable-xdebug;o50;6;0;29;18;0

largest by out degree
zcat /data/basemaps/gz/$base.s | awk -F\; '{ if ($4 > o){ o=$4; print $0} }'
TryGhost_Casper;o;976;5034;53403;49105;1;judge_security-drops;o;2;0;150;148;0
2215_he-said-she-said.github.io;o;2093;16089;596119;593274;1;wix_wix-hive-activity-schemas;o;3;1;206;189;0
udacity_frontend-nanodegree-resume;o;8034;17607;415022;385502;1;mmayorivera_angular-sync;u;63;32;155;74;0
frioux_dotfiles;o70;8909;47719;545180;432205;1;zlraymind_vim_config_mac;o50;28;10;107;73;0
gitlab.com_gitlab-com_www-gitlab-com;o;31725;96738;28554613;28209339;2;hamxiaoz_StandOut;o50;1;0;9;7;0
cdnjs_cdnjs;o50;12103;134166;1274921;697481;5;nvbkdw_b3nchai;d70;70;100;2125;510;0
barryclark_jekyll-now;o;91298;209624;13450186;12604996;1;bibhu545_HealthCare;o;21;10;608;584;0

largest by number of artifacts shared out with a single other reopo
zcat /data/basemaps/gz/$base.s | awk -F\; '{ if ($7 > s&& $1 != $8){ s=$7; print $0} }'
eclipse_sumo;o;286;929;465796;423181;112111;trezheur_sumo_svn;d;50;0;150147;1;0
eclipse_jetty.website;o70;40;1;209069;175137;174674;git.eclipse.org_r_www.eclipse.org/jetty;o50;48;137;527649;318069;676
usehotkey_counter;u;3519;1521;997157;450476;450474;one-million-repo_one-million-repo;d;3519;22;997156;2060;2059
freebsd_freebsd;o70;1571;16610;3443707;2857124;810169;mat813_freebsd;d;805;133;1118150;564;86
g0v-data_mirror;o;19;49;6816419;6747858;1024637;g0v-data_mirror-2016;d;8;1;1034766;4698;4698
git-portage_git-portage;o50;1374;4556;3493448;1873632;1675301;gentoo_gentoo;o50;2652;0;18191419;10679708;0
gentoo_gentoo-portage-rsync-mirror;o50;1787;778;10811815;6749143;4225197;gentoo_gentoo;o50;2652;0;18191419;10679708;0

largest by number of artifacts shared in with a single other reopo
zcat /data/basemaps/gz/$base.s | awk -F\; '{ if ($14 > s&& $1 != $8){ s=$14; print $0} }'
bloomberg_chromium.bb;u;5342;60627;9317283;3626734;1394237;chromium_chromium;o50;2053;10999;7370627;4218074;2624148
Oleh-Kravchenko_portage;d90;1699;517;7136663;84107;39732;gentoo_gentoo-portage-rsync-mirror;o50;1787;778;10811815;6749143;3715461

#overall
base=fP2PStatsT
i
zcat /data/basemaps/gz/$base.Pstats| awk -F\; '{ if ($3 > i){ i=$3; print $0} }'
diogom42_d4jdata;o;2;5;303;301
alexey-oblomov_chat_app_backend;o70;15;1;183;160
kylemezzacappa_Gratis;d90;40;1;2145;91
FishAres_Learning-stuff;o70;108;2;679;517
rzfreitas_verde-ghaia;o50;123;9;5616;3750
ashkan18_theReader;d90;184;226;13285;136
nicky132_reactMobx;d90;256;43;26052;1864
bitbucket.org_Carvalier_snapthailand-web;d90;332;8;7643;459
changchiajung_vue-demo;d;431;4;15141;43
MaxYuan85_ChainGamble;d90;433;1;11603;770
swaaminathanm_OTART201802;d90;464;56;15607;1446
vladimireanu_react-e-commerce;d90;502;58;28440;2104
kdv24_Starter_Codes;d;724;31;8036;17
rajnish42413_pizza-front;d70;838;1;25616;5673
gitfer_our-monthly-expenses;d;970;385;15308;22
backpack455_HEARD;d;1626;12;50376;103
JacobcobLee_Fault-reporter;d70;1824;23;71712;14036
rafaelzolt_HugoQuickStart;d;2578;91;41404;71
TryGhost_Casper;o50;3102;18485;95930;49411
Anterotesis_historical-texts;d;32706;227;62340;285
pombredanne_license-detection-test-data;d;67548;258;73723;480
x0rzkov_dockerfiles-search;o50;157323;497;500653;283834


o
zcat /data/basemaps/gz/$base.Pstats| awk -F\; '{ if ($4 > i){ i=$4; print $0} }'
reviewboard_rb-extension-pack;o70;18;50;700;625
UltimateAngular_angular-pro-src;o70;41;134;1276;1106
vyacheslav31_React-Monsters-Tutorial;o50;17;185;58;33
AccordLTN_my-odin-project;u;68;487;271;126
veros7821_web;o70;6;488;36;29
trliner_fake_id;o;0;511;15;15
iominh_mramato.github.io;d70;88;1205;4763;731
camspiers_json-pretty;o;3;2798;46;43
prodev1017_react-bootstrap-master;o50;6;5349;591;381
akbr_prestige;o70;152;13259;4338;3331
gitorious.org_page-six;o70;22;30033;1693;1255
mathiasbynens_jquery-placeholder;o;13;33206;495;479
houdiniproject_houdini;o;130;36022;9462;9181
airbnb_babel-plugin-dynamic-import-node;o;16;36137;520;482
postcss_autoprefixer;o;86;116347;6148;5939
luhur65_web-berita-Rest-API;d90;7;305678;1655;152
sir-haleem_aremuacademy-redesigned;o50;9;549260;56;38
ajbt200128_Attendance-System;o;2;1499507;141;139
adbl_misfire-ui;o70;63;1576626;4343;3684
sprietNathanael_CPE_4_EmbeddedSystems;o50;40;2205859;1596;888
px4_ros-examples;o;4;3118972;474;464
mozilla_gecko-dev;o;4153;18660102;10404110;9467420
Unipisa_test_git_old_dates;o;0;22557223;3;3


gitorious.org_page-six;o70;22;30033;1693;1255;1;sheepcao_ZeroCityShop;d70;616;13;20815;3740;0
mathiasbynens_jquery-placeholder;o;13;33206;495;479;1;second-ob_qiang;d70;85;1;476;48;0
houdiniproject_houdini;o;130;36022;9462;9181;1;osiotestmachine_osio-ci-boost3-BDD-0107-0812-test123;d;21;0;36;1;0 
airbnb_babel-plugin-dynamic-import-node;o;16;36137;520;482;1;chrisyuaners_plan-it-app;d90;489;16;23773;532;0 # .npmignore
postcss_autoprefixer;o;86;116347;6148;5939;3;mojopoly_quotesondev;d90;375;6;14444;1045;0 # "\n\n"
luhur65_web-berita-Rest-API;d90;7;305678;1655;152;5;grmnlrt_secret-santa;d;365;0;12103;34;0 # jquery-3.4.1.min.js
sir-haleem_aremuacademy-redesigned;o50;9;549260;56;38;10;ShafiqHalim007_Shafiq-Halim;u;97;29;423;151;0 #2M glyphicons-halflings-regular.woff
ajbt200128_Attendance-System;o;2;1499507;141;139;7;mainangethe_networking-app;o50;62;11;293;164;0 - initializers/backtrace_silencers.rb 59385cdf379bd06a8d2326dcd4de6d5cd5d3f5b0
adbl_misfire-ui;o70;63;1576626;4343;3684;7;sumitgpl_awsLambda;d70;46;0;385;57;0 # over 20 widely used  with 1.3M The ISC License
sprietNathanael_CPE_4_EmbeddedSystems;o50;40;2205859;1596;888;1;ejdzipi_emotion-jest-serializer-repro;o70;1;0;9;8;0 #.idea\nnode_modules
px4_ros-examples;o;4;3118972;474;464;1;PabloViniegra_ejercicioADEmpleados;o70;4;0;40;36;0 # project file VcsDirectoryMappings 94a25f7f4cb416c083d265558da75d457237d671
mozilla_gecko-dev;o;4153;18660102;10404110;9467420;1;webmproject_vp9-dash;o;3;0;100;97;0 - "\n\n"
Unipisa_test_git_old_dates;o;0;22557223;3;3;1;ruturaaj_jqtoggle;o50;5;0;17;11;0 - empty blob

#most original
zcat /data/basemaps/gz/$base.Pstats| awk -F\; '{ if ($6 > i){ i=$6; print $0} }'
diogom42_d4jdata;o;2;5;303;301
reviewboard_rb-extension-pack;o70;18;50;700;625
rzfreitas_verde-ghaia;o50;123;9;5616;3750
rajnish42413_pizza-front;d70;838;1;25616;5673
arokde_ScratchPadProjects;d70;423;17;19873;5952
d5j6_Geology_BETA;o50;208;33;11624;5983
elymichael_BibliotecaLibros;o;48;7;9340;9177
Mozu_mozu-ios-sdk;o;17;7;12748;12710
Julian-Wyatt_Sudoku;o;7;1;70972;70965
Smithsonian_OpenAccess;o;3;1;138020;138017
CenterForOpenScience_osf.io;o70;963;10852;168668;148335
DoutorGois_BTI;o;11;11;263007;261967
MercifulPotato_mercifulpotato;o;8;2;1787801;1784523
freebsd_freebsd;o70;2337;350634;3665838;2860153
guillermoiglesiashernandez_Dataset;o;12;1;4243527;4243513
grantmakers_profiles;o;43;25;7725253;7725032
elastic_docs;o;160;163;10169275;10070686
whosonfirst-data_whosonfirst-data;o;11;835;15959110;15958156
brbrr_test-dashboard-pages;o;58;2;17216790;17131103
gitlab.com_gitlab-com_www-gitlab-com;o;53298;479814;28924219;28213688
jdtournier_mrtrix3-dev-doc;o;191;162;36998066;36821950
swift-zym_scp-docs;o;18;0;39366489;39120868

no - stats
zcat /data/basemaps/gz/$base.s | awk -F\; '{ if ($6 > s&& $1 != $8){ s=$6; print $0} }'
diogom42_d4jdata;o;2;5;303;301
reviewboard_rb-extension-pack;o70;18;50;700;625
rzfreitas_verde-ghaia;o50;123;9;5616;3750
rajnish42413_pizza-front;d70;838;1;25616;5673
arokde_ScratchPadProjects;d70;423;17;19873;5952
d5j6_Geology_BETA;o50;208;33;11624;5983
elymichael_BibliotecaLibros;o;48;7;9340;9177
Mozu_mozu-ios-sdk;o;17;7;12748;12710
Julian-Wyatt_Sudoku;o;7;1;70972;70965
Smithsonian_OpenAccess;o;3;1;138020;138017
CenterForOpenScience_osf.io;o70;963;10852;168668;148335
DoutorGois_BTI;o;11;11;263007;261967
MercifulPotato_mercifulpotato;o;8;2;1787801;1784523
freebsd_freebsd;o70;2337;350634;3665838;2860153
guillermoiglesiashernandez_Dataset;o;12;1;4243527;4243513
grantmakers_profiles;o;43;25;7725253;7725032
elastic_docs;o;160;163;10169275;10070686
whosonfirst-data_whosonfirst-data;o;11;835;15959110;15958156
brbrr_test-dashboard-pages;o;58;2;17216790;17131103
gitlab.com_gitlab-com_www-gitlab-com;o;53298;479814;28924219;28213688
jdtournier_mrtrix3-dev-doc;o;191;162;36998066;36821950
swift-zym_scp-docs;o;18;0;39366489;39120868


60d7d529a95e:/data/play/forks>zcat /data/basemaps/gz/$base.s | awk -F\; '{ if ($7 > s&& $1 != $8){ s=$7; print $0} }'
diogom42_d4jdata;o;2;5;303;301;1;rwang63_Arbitrage-Finder;o70;13;0;51;38;0
Eiv1nd_MachineLearning;o70;1;1;5;4;3;eriktok_AntcolonyOptimization;d;2;0;4;0;0
DimitraApostolopoulou_PhotoViewer;o50;17;2;49;31;9;dimiap_ImageViewer;u;15;0;39;13;0
reviewboard_rb-extension-pack;o70;18;50;700;625;51;ZHKKKe_se-workshop;d;5;0;60;1;0
UltimateAngular_angular-pro-src;o70;41;134;1276;1106;226;PowerlessCube_Todd-Motto-Angular;u;14;4;474;161;2
UltimateAngular_angular-pro-src;o70;41;134;1276;1106;227;iliassk_UltimateAngular;d70;37;14;601;69;0
UltimateAngular_angular-pro-src;o70;41;134;1276;1106;232;manucho007_UltimateCourses;d70;40;5;660;196;1
gitorious.org_reusable/reusable;o70;62;148;1377;1185;757;cqnu_Reusable;d70;97;1;1481;437;0
WangscGit_learnGit;u;173;6;1698;814;812;WangscGit_gitWarehouse;d;174;0;1698;2;0
ondras_js-like;o;8;120;2289;2133;2130;droidenko_js-like;d70;8;3;2403;265;148
Akai-Kumako_Experiment;o;2;2;22582;22481;4944;kamuiroeru_nitacLT;d;5;3;5086;17;0
Akai-Kumako_Experiment;o;2;2;22582;22481;9809;yamasy1549_exp-nlp;d90;3;1;10023;113;0
sourishbanerjee_APPARELCBIRVGG16;o;1;7;85398;85397;80412;bitbucket.org_udithprabhu_hackathon_dataset;u;2;0;140893;60480;0
eclipse_sumo;o;488;2296;467638;423188;112117;trezheur_sumo_svn;d;83;0;150222;1;0
eclipse_jetty.website;o70;79;1;209232;175137;174674;git.eclipse.org_r_www.eclipse.org/jetty;o50;96;549;527841;318070;677
usehotkey_counter;u;3619;1521;1000003;450476;450474;one-million-repo_one-million-repo;d;3618;22;1000001;2060;2059
ssbattousai_Cuda14;u;5;20;1796523;778157;698925;ssbattousai_Cuda15;u;5;5;1783087;538616;334610
ssbattousai_Cuda14;u;5;20;1796523;778157;754307;ssbattousai_Cuda12;d70;5;20;1832312;322227;297631
freebsd_freebsd;o70;2337;350634;3665838;2860153;812326;mat813_freebsd;d;1323;133;1308058;564;86
g0v-data_mirror;o;21;49;6816423;6747858;967636;g0v-data_mirror-2017;d;7;1;991161;4833;4833
g0v-data_mirror;o;21;49;6816423;6747858;1024637;g0v-data_mirror-2016;d;10;1;1034768;4698;4698
git-portage_git-portage;o50;1534;8178;3494204;1873719;1041325;minaco2_gentoo-portage-mirror;d;1152;2;2269494;22;4
git-portage_git-portage;o50;1534;8178;3494204;1873719;1675388;gentoo_gentoo;o50;2859;18398;18192295;10679971;1415495
gentoo_gentoo-portage-rsync-mirror;o50;1963;778;10812562;6749143;4225197;gentoo_gentoo;o50;2859;18398;18192295;10679971;2856807
9
#Questions what is reused why reused (why not package manager) do therefr where it was copied, does it need to be maintained (license/image), IP/licensing (SO snippets) amount invisible reuse top predict demand (issues)
1) What is the repo for?
2) What artifacts does it contain?
3) Why the artifacts it contains are  useful?
4)   for whom?
5) Who are biggest contributors?
6) What are the most shared artifacts it produced?
7) What are other repos with wich it shares the most artifacts?

Q about blobs
   Usage and Attribution of Stack Overflow Code Snippets in GitHub Projects https://empirical-software.engineering/assets/pdf/emse18-snippets.pdf 
   SOTorrent: Studying the Origin, Evolution, and Usage of Stack Overflow Code Snippets https://arxiv.org/pdf/1809.02814.pdf

   http://www.cs.columbia.edu/~blei/papers/RanganathPerotteElhadadBlei2015.pdf 

#Use histogram in terms of copying
library(data.table)
x=read.table("fP2PStatsT.Pstats.csv", sep=";",quote="",comment.char="")
xL=read.table("fP2LPStatsT.Pstats.csv", sep=";",quote="",comment.char="")
nn = c("prj", "type", "inD", "outD", "nb","nob");
names(x)=nn
names(xL)=nn
quantile(x$inD,0:10/10)                                                                                                                                                                                                                                                                                       
0%    10%    20%    30%    40%    50%    60%    70%    80%    90%   100% 
0      0      1      1      2      3      5      8     13     27 157323 
quantile(x$outD,0:10/10)                                                                                                                                   0%    10%   20%    30%    40%    50%  60%      70%   80%      90%     100% 
0      0     0      0      0      0    0        1     1        4 22557223 
     
quantile(x$nb,0:10/10)
      0%      10%      20%      30%      40%      50%      60%      70% 
       1        2        3        7       12       21       36       61 
     80%      90%     100% 
     113      324 39366489 
quantile(x$nob,0:10/10)
      0%      10%      20%      30%      40%      50%      60%      70% 
       0        1        2        4        7       11       19       32 
     80%      90%     100% 
      61      152 39120868 


> quantile(xL$inD,0:10/10) 
    0%    10%    20%    30%    40%    50%    60%    70%    80%    90%   100% 
     0      0      0      0      0      1      1      2      4      7 156746 
> quantile(xL$outD,0:10/10) 
    0%    10%    20%    30%    40%    50%    60%    70%    80%    90%   100% 
     0      0      0      0      0      0      0      1      1      5 209624 
> quantile(xL$nb,0:10/10) 
      0%      10%      20%      30%      40%      50%      60%      70% 
       1        1        3        6       10       16       25       43 
     80%      90%     100% 
      80      202 39366484 
> quantile(xL$nob,0:10/10) 
      0%      10%      20%      30%      40%      50%      60%      70% 
       0        1        2        4        8       12       20       34 
     80%      90%     100% 
      64      160 39120868 
      
selLL = xL$nb>20
> quantile(xL$outD[selLL],0:10/10)
0%    10%    20%    30%    40%    50%    60%    70%    80%    90%   100% 
     0      0      0      0      0      1      1      2      4     11 209624 
> quantile(xL$inD[selLL],0:10/10) 
     0      0      1      1      2      3      4      5      8     13 156746
selL = selLL & (xL$outD > 4 | xL$inD > 8)

selAL = x$nb>20
quantile(x$outD[selAL],0:10/10)
      0%      10%      20%      30%      40%      50%      60%      70% 
       0        0        0        0        0        0        1        2 
     80%      90%     100% 
       4       10 18660102 
quantile(x$inD[selAL],0:10/10)
    0%    10%    20%    30%    40%    50%    60%    70%    80%    90%   100% 
     0      2      3      5      7     10     13     18     26     50 157323 

selA = selAL &  (x$outD > 10 | x$inD > 50)


#Do first sampling: 40 projects


selLL = xL$nb>20
selL = selLL & (xL$outD > 4 | xL$inD > 8)

selAL = x$nb>20
selA = selAL &  (x$outD > 10 | x$inD > 50)
xSample = sample(as.character(x$prj[selA]),100)
xLSample = sample(as.character(xL$prj[selL]),100)

summary(x[match(xSample,x$prj),])
                                   prj      type         inD        
 2016llf_LLFSwifDemo                 : 1   d  : 3   Min.   :  0.00  
 2040RICARDO_sistemaelectromecanicav1: 1   d70:20   1st Qu.: 19.25  
 ADAMBURTON428_Jedds                 : 1   d90:21   Median : 58.00  
 AnaBrando_SoftBar                   : 1   o  :14   Mean   : 81.87  
 Appefy_Graphico                     : 1   o50: 8   3rd Qu.:102.00  
 Chaitanya-SiteMinder_MEANApp        : 1   o70:22   Max.   :802.00  
 (Other)                             :94   u  :12                   
      outD              nb               nob         
 Min.   :  0.00   Min.   :   24.0   Min.   :   2.00  
 1st Qu.:  2.00   1st Qu.:  183.5   1st Qu.:  59.75  
 Median : 15.00   Median :  668.0   Median : 162.50  
 Mean   : 34.14   Mean   : 1714.3   Mean   : 451.79  
 3rd Qu.: 43.00   3rd Qu.: 2192.2   3rd Qu.: 422.00  
 Max.   :292.00   Max.   :18834.0   Max.   :5841.00  

summary(xL[match(xLSample,xL$prj),])
                                     prj      type         inD       
 1004839068_project                    : 1   d  : 3   Min.   : 0.00  
 AlphaDelete_AguaConsciente            : 1   d70: 9   1st Qu.: 5.00  
 Ariel1990_first_app                   : 1   d90: 3   Median : 9.50  
 B-H-M_Parcel-Sender-MERN              : 1   o  :43   Mean   :14.82  
 BigCatMr_SwiftCat                     : 1   o50: 9   3rd Qu.:20.00  
 Blazingsonic_firecode-work-in-progress: 1   o70:25   Max.   :80.00  
 (Other)                               :94   u  : 8                  
      outD              nb             nob         
 Min.   :  0.00   Min.   :   23   Min.   :    3.0  
 1st Qu.:  1.00   1st Qu.:   95   1st Qu.:   54.0  
 Median :  7.00   Median :  241   Median :  141.5  
 Mean   : 22.66   Mean   : 1031   Mean   :  823.4  
 3rd Qu.: 19.00   3rd Qu.:  525   3rd Qu.:  390.2  
 Max.   :244.00   Max.   :34434   Max.   :31275.0  

 
quantile(x$inD[selAL],90:100/100)
   90%    91%    92%    93%    94%    95%    96%    97%    98%    99%   100% 
    50     55     62     70     81     95    115    146    200    318 157323 
quantile(x$outD[selAL],90:100/100)
     90%      91%      92%      93%      94%      95%      96%      97% 
      10       12       14       17       20       25       33       47 
     98%      99%     100% 
      77      176 18660102 

quantile(xL$inD[selLL],90:100/100)
    90%    91%    92%    93%    94%    95%    96%    97%    98%    99%   100% 
    13     14     15     17     19     21     25     29     38     57 156746 
quantile(xL$outD[selLL],90:100/100)
   90%    91%    92%    93%    94%    95%    96%    97%    98%    99%   100% 
    11     13     15     18     21     26     33     43     61     95 209624 


selLH = selLL & (xL$outD > 95 | xL$inD > 57)
selH = selAL & (x$outD > 176 | x$inD > 318)
hSample = sample(as.character(x$prj[selH]),100)
hLSample = sample(as.character(xL$prj[selLH]),100)
write(hSample,file="hSample",ncol=1)
write(hLSample,file="hLSample",ncol=1)
write(Sample,file="Sample",ncol=1)
write(LSample,file="LSample",ncol=1)

summary(xL[match(hLSample,xL$prj),])
                                prj      type         inD        
 ARCANEDEV_Support                : 1   d  : 4   Min.   :  0.00  
 Andry85_prostositevnua2019       : 1   d70: 9   1st Qu.: 21.00  
 Asad24_me                        : 1   d90: 7   Median : 61.00  
 Brew8it_Twitter-SA-Project       : 1   o  :32   Mean   : 63.25  
 CadiDadi_crypto-portfolio-tracker: 1   o50:13   3rd Qu.: 86.25  
 CarlosjRuiz_guia_practica        : 1   o70:21   Max.   :298.00  
 (Other)                          :94   u  :14                   
      outD             nb               nob         
 Min.   :  0.0   Min.   :   33.0   Min.   :    0.0  
 1st Qu.: 15.5   1st Qu.:  444.5   1st Qu.:  153.8  
 Median : 78.0   Median : 1395.5   Median :  495.0  
 Mean   :100.1   Mean   : 3847.3   Mean   : 2764.3  
 3rd Qu.:135.0   3rd Qu.: 3696.8   3rd Qu.: 1601.5  
 Max.   :751.0   Max.   :81027.0   Max.   :80172.0  

summary(x[match(hSample,x$prj),])
                    prj      type         inD              outD         
 1079107009_MyAndroid : 1   d  : 6   Min.   :   0.0   Min.   :    0.00  
 790396054_zhbj74     : 1   d70:26   1st Qu.:  56.5   1st Qu.:    5.75  
 AlexCout_Turretz     : 1   d90:24   Median : 339.0   Median :  200.50  
 AlexisFinn_linux-home: 1   o  : 4   Mean   : 333.3   Mean   :  670.83  
 Aliis_testpage1K     : 1   o50: 9   3rd Qu.: 488.5   3rd Qu.:  331.00  
 Ant97_MissionSupport : 1   o70:18   Max.   :1332.0   Max.   :14703.00  
 (Other)              :94   u  :13                                      
       nb               nob          
 Min.   :     28   Min.   :     7.0  
 1st Qu.:   1058   1st Qu.:   187.5  
 Median :   7532   Median :   951.0  
 Mean   :  19946   Mean   : 11117.1  
 3rd Qu.:  15433   3rd Qu.:  2080.0  
 Max.   :1042247   Max.   :919854.0  

for i in LSample hLSample; do head -10 $i; done | gzip > selL
for i in Sample hSample; do head -10 $i; done | gzip > sel
(cat fP2PStatsT.Pstats.csv | ~/lookup/grepField.perl sel 1; cat fP2LPStatsT.Pstats.csv|  ~/lookup/grepField.perl selL 1)|sort -R
cat smp1 | python3 query.py

#most projects very short duration
#get projects with 10+ active months
python3 query1.py > ManyMonthsActive
lsort 10G -t\; -k1,1 ManyMonthsActive |join -t\; - <(zcat /data/basemaps/gz/fP2PStatsT.Pstats | lsort 100G -t\; -k1,1) >  ManyMonthsActive.Pstats 
lsort 10G -t\; -k1,1 ManyMonthsActive |join -t\; - <(zcat /data/basemaps/gz/fP2LPStatsT.Pstats | lsort 100G -t\; -k1,1) >  ManyMonthsActiveL.Pstats 

y=read.table("ManyMonthsActive.Pstats", sep=";",quote="",comment.char="")
nn1 = c("prj", "na", "nc", "nCore", "nMnth", "m0", "mn", "gen", "type", "inD", "outD", "nb","nob");
names(y)=nn1
> summary(y)
                      prj                na                  nc          
 0--key_0--key.github.io:      1   Min.   :     1.00   Min.   :      10  
 0--key_lib             :      1   1st Qu.:     2.00   1st Qu.:      66  
 0--key_org-pub         :      1   Median :     3.00   Median :     130  
 0-0-0-_StellarKit      :      1   Mean   :    12.02   Mean   :     651  
 0-0MrLonely_SourceCode :      1   3rd Qu.:     7.00   3rd Qu.:     294  
 0-1-0_twistock         :      1   Max.   :134271.00   Max.   :31610821  
 (Other)                :2107482   NA's   :3064                          
     nCore              nMnth               m0                mn         
 Min.   :    1.00   Min.   :  10.00   2018-03:  28007   2021-02: 277603  
 1st Qu.:    1.00   1st Qu.:  11.00   2018-10:  27658   2021-01: 175150  
 Median :    1.00   Median :  15.00   2018-01:  27399   2021-03: 147578  
 Mean   :    3.73   Mean   :  20.89   2019-03:  27371   2020-12: 118282  
 3rd Qu.:    2.00   3rd Qu.:  23.00   2019-01:  26732   2020-11:  63768  
 Max.   :68872.00   Max.   :1304.00   2017-03:  26636   2020-10:  62138  
 NA's   :3021                         (Other):1943685   (Other):1262969  
     gen           type              inD                outD          
       :1800080   d  :  10967   Min.   :     0.0   Min.   :      0.0  
 female: 307408   d70:  95753   1st Qu.:     4.0   1st Qu.:      1.0  
                  d90:  39317   Median :    12.0   Median :      4.0  
                  o  :1179453   Mean   :    41.5   Mean   :    201.7  
                  o50: 177404   3rd Qu.:    33.0   3rd Qu.:     17.0  
                  o70: 490019   Max.   :133630.0   Max.   :2302297.0  
                  u  : 114575                                         
       nb                nob          
 Min.   :       1   Min.   :       0  
 1st Qu.:     170   1st Qu.:     134  
 Median :     441   Median :     341  
 Mean   :    3128   Mean   :    2247  
 3rd Qu.:    1319   3rd Qu.:     937  
 Max.   :36998066   Max.   :36821950  

yL=read.table("ManyMonthsActiveL.Pstats", sep=";",quote="",comment.char="")
names(yL)=nn1
summary(yL)
                      prj                na                  nc          
 0--key_0--key.github.io:      1   Min.   :     1.00   Min.   :      10  
 0--key_lib             :      1   1st Qu.:     2.00   1st Qu.:      66  
 0--key_org-pub         :      1   Median :     3.00   Median :     130  
 0-0-0-_StellarKit      :      1   Mean   :    12.02   Mean   :     649  
 0-0MrLonely_SourceCode :      1   3rd Qu.:     7.00   3rd Qu.:     294  
 0-1-0_twistock         :      1   Max.   :134271.00   Max.   :31610821  
 (Other)                :2106836   NA's   :3060                          
     nCore              nMnth               m0                mn         
 Min.   :    1.00   Min.   :  10.00   2018-03:  28002   2021-02: 277598  
 1st Qu.:    1.00   1st Qu.:  11.00   2018-10:  27646   2021-01: 175146  
 Median :    1.00   Median :  15.00   2018-01:  27390   2021-03: 147576  
 Mean   :    3.73   Mean   :  20.89   2019-03:  27358   2020-12: 118271  
 3rd Qu.:    2.00   3rd Qu.:  23.00   2019-01:  26725   2020-11:  63762  
 Max.   :68872.00   Max.   :1304.00   2017-03:  26631   2020-10:  62127  
 NA's   :3017                         (Other):1943090   (Other):1262362  
     gen           type              inD                outD          
       :1799454   d  :   7101   Min.   :    0.00   Min.   :     0.00  
 female: 307388   d70:  30199   1st Qu.:    2.00   1st Qu.:     1.00  
                  d90:  12107   Median :    5.00   Median :     4.00  
                  o  :1582586   Mean   :   13.47   Mean   :    28.95  
                  o50:  94710   3rd Qu.:   12.00   3rd Qu.:    16.00  
                  o70: 329764   Max.   :91557.00   Max.   :209627.00  
                  u  :  50375                                         
       nb                nob          
 Min.   :       1   Min.   :       0  
 1st Qu.:     153   1st Qu.:     134  
 Median :     389   Median :     341  
 Mean   :    2601   Mean   :    2245  
 3rd Qu.:    1081   3rd Qu.:     937  
 Max.   :36995595   Max.   :36821950  

yL[yL$nMnth==1304,]
                                prj na    nc nCore nMnth      m0      mn gen
202318 IPDSnelting_velcom-test-repo  3 84861     1  1304 1979-12 2088-07    
       type inD outD nb nob
202318  o70   0    0 28  27

max( yL[yL$nMnth!=1304,"nMnth"])
[1] 1093
yL[yL$nMnth==1093,]
                             prj na     nc nCore nMnth      m0      mn gen type
959326 fcharlie_git-analyze-demo  3 265830     2  1093 2009-12 2100-12      o70
       inD outD  nb nob
959326  13   53 858 746

dim(yL[yL$nMnth>240,]);
[1] 608  13

selyLL = yL$nb>20
selyAL = y$nb>20

quantile(y$inD[selyAL],0:10/10)
    0%    10%    20%    30%    40%    50%    60%    70%    80%    90%   100% 
     0      1      3      5      8     12     18     26     43     86 133630 

quantile(y$outD[selyAL],0:10/10)
      0%     10%     20%     30%     40%     50%     60%     70%     80%     90% 
      0       0       0       1       2       4       7      12      25      78 
   100% 
2302297 

quantile(yL$inD[selyLL],0:10/10)
   0%   10%   20%   30%   40%   50%   60%   70%   80%   90%  100% 
    0     0     1     2     3     5     7    10    15    28 91557 

quantile(yL$outD[selyLL],0:10/10)
    0%    10%    20%    30%    40%    50%    60%    70%    80%    90%   100% 
     0      0      0      1      2      4      7     12     23     58 209627 

ySample = sample(as.character(y$prj[selyAL]),100)
yLSample = sample(as.character(yL$prj[selyLL]),100)

selyL = selyLL & (yL$outD > 15 | yL$inD > 23)
selyA = selyAL &  (y$outD > 78 | y$inD > 86)
hySample = sample(as.character(y$prj[selyA]),100)
hyLSample = sample(as.character(yL$prj[selyL]),100)
write(hySample,file="hySample",ncol=1)
write(hyLSample,file="hyLSample",ncol=1)
write(ySample,file="ySample",ncol=1)
write(yLSample,file="yLSample",ncol=1)
for i in yLSample hyLSample; do head -10 $i; done | gzip > selyL
for i in ySample hySample; do head -10 $i; done | gzip > sely
(cat fP2PStatsT.Pstats.csv | ~/lookup/grepField.perl sely 1; cat fP2LPStatsT.Pstats.csv|  ~/lookup/grepField.perl selyL 1)|sort -R > smpy
cat smpy | python3 query.py


#network of repos that share 20+blobs at a time
base=fP2PStatsT
zcat /data/basemaps/gz/$base.s | awk -F\; '{if ($1!=$8 && $7>20 )print $0}' |gzip > $base.20+
zcat  $base.20+|cut -d\; -f1 |uniq -c |lsort 100G -rn |head
1398003 nakiostudio_Youtube-Downloader
1228434 wujun0919_remote
1075051 QQ986945193_david_hexo_blog
 960731 lijinhai255_nodeBasic
 932881 LokenderSarna_catharsis15

zcat  $base.20+|cut -d\; -f8 |uniq -c |lsort 100G -rn |head
   8115 Mattlk13_metacpan-cpan-extracted
   6483 gitlab.com_atoomic_CPAN
   6271 barryclark_jekyll-now
   5402 drupal_drupal
   5001 AmartyaU_WordPress
   4858 TheMattSykes_personal-website
   4825 Jackeagle_kernel_msm-3.18

#what is being shared
LICENSE


2215_he-said-she-said.github.io;1038;1038;0;5933;d;597072;744;gitlab.com_bmcorser_pars;4
2215_he-said-she-said.github.io;1038;1038;0;5933;d;597072;744;gitlab.com_bmcorser_pars;4
udacity_frontend-nanodegree-resume;0;4037;0;9187;d70;582183;147218;Kottans_frontend-2019-homeworks;46
udacity_frontend-nanodegree-resume;0;4037;0;9187;d70;582183;147218;Kottans_frontend-2019-homeworks;46
frioux_dotfiles;0;1391;0;14007;d90;635860;61600;skeept_dotvim;2944
Adude11_-tmp-100000-commit-2;0;154;0;41556;d70;50000;13570;yangyimeng_test;13570



# means to identify infrastructure (define infrastructure types) - tools to 
# compare from dependencies based infrastructure
# Arun code reuse four types - why
# copy is fine grained: allows excluding unneeded functionality

Add to Prj summaries
+Blob-Linked projects b2B b2b 
+NOriginal blobs P2nfb
LICENSE
zcat fP2PStatsT.StrongReuseLinks011-1.crank.map | awk -F\; '{print $2";"$1}' | lsort 100G | gzip > B2bFullT.s
zcat B2bFullT.s|join -t\; - <(zcat B2bFullT.s)|gzip > B2BFullT.s

Add to Auth summaries 
community?


#Deblobing
time perl getType1.perl | gzip > fP2PStatsT.s
time perl getType1.perl | gzip > fP2LPStatsT.s
time perl getType1.perl | gzip > fP2SLPStatsT.s
time perl getType1.perl | gzip > fP2SPStatsT.s
zcat fP2SLPStatsT.s | cut -d\; -f1-6 | uniq | gzip > fP2SLPStatsT.Pstats
zcat fP2LPStatsT.s | cut -d\; -f1-6 | uniq | gzip > fP2LPStatsT.Pstats
zcat fP2SPStatsT.s | cut -d\; -f1-6 | uniq | gzip > fP2SPStatsT.Pstats
zcat fP2PStatsT.s | cut -d\; -f1-6 | uniq | gzip > fP2PStatsT.Pstats

#by type
for sel in "" S L
do zcat fP2${sel}PStatsT.Pstats | awk -F\; '{ t[$2]++; if($5>20)t20[$2]++; if($5>10)t10[$2]++; if($5>5)t5[$2]++;} END{for (i in t)print i,t[i],t5[i],t10[i],t20[i]; }'&
done

SL
typ all       5+         10+      20+
  o 40711544 24542551 18767869 13831864
o70  9836186  7936200  6645830  5294805
o50  4308962  4308962  3276940  1245493
  u  1966488  1490827  1158141   851447
d70  1528038  1107076   882440   694167
d90   638849   638849   490091   380130
  d  3676685  1392425  1111324   809000

L  
  o 48564010 30232023 24674903 19869258
o70 20103794 16615451 14193532 11492055
o50  7895500 7895500   6294897  2644786
  u  3129304 2528510   2106049  1662689
d70  2154602 1721170   1451796  1189351
d90   892436 892436     722867   571906
d    3429407 1348216   1185651   918449

S
o   18787138 10964502  9030870 7613857
o70 15305134 11208933  9085628 7376845
o50  9081051  9081051  6789467 3157348
u    5703360  4523251  3754207 2959300
d70  5035339  3917862  3340136 2727922
d90  2873749  2873749  2404938 2004529
d   15969710  5363680  4581397 3586242

Full  
o   24488920 13303031 10917881  9334732
o70 22070677 16670002 14431268 12692143
o50 13940010 13940010 11358591  6333512
d70  7872040  6709140  6127959  5458384
d90  3829321  3829321  3482630  3113294
d    8209126  3274916  3032055  2557025
Does not add up to 89M

#blobs not shared
for sel in "" S L
do zcat fP2${sel}PStatsT.Pstats | awk -F\; '{i++;if ($5> 10)i10++; if ($5> 20)i20++;if ($5> 5)i5++;if ($3+$4>0){is++;if ($5> 10)is10++; if ($5> 20)is20++;if ($5> 5)is5++;}}END {print "1+|"i"|"is"|"is/i; print "5+|"i5"|"is5"|"is5/i5; print "10+|"i10"|"is10"|"is10/i10;print "20+|"i20"|"is20"|"is20/i20}'&
done
   
SL
Size|Total    |with sharing|fraction
1+  |62666752 |31550488    |0.503465
5+  |41416890 |25338451    |0.61179
10+ |32332635 |21430142    |0.662802
20+ |23106906 |16524634    |0.715138

L
1+  |86169053 |55978578    |0.649637
5+  |61233306 |47880688    |0.781939
10+ |50629695 |42015385    |0.829857
20+ |38348494 |33618084    |0.876647

S # few blogs account for a lot of connectivity
1+  |72755481 |60784594    |0.835464
5+  |47933028 |43397088    |0.905369
10+ |38986643 |36229075    |0.929269
20+ |29426043 |27951977    |0.949906

Full
1+  |89450504 |74491883    |0.832772
5+  |65542299 |60932347    |0.929664
10+ |56297593 |53780036    |0.955281
20+ |45389449 |44217596    |0.974182


zcat /data/basemaps/gz/$base.s |  cut -d\; -f1,8  | perl ~/bin/connectBasic.perl $base |gzip > $base.map
zcat $base.map|cut -d\; -f2 |lsort 100G |uniq -c |lsort 10G -rn |head
base=fP2PStatsT 73261916/89450504=.82
base=fP2SLPStatsT 26854991/62666752=.43
base=fP2LPStatsT 51945572/86169053=.60
base=fP2SPStatsT 59936632/72755481=.82

export LD_LIBRARY_PATH=/home/audris/lib64:/home/audris/lib:/home/audris/src/networkit/build
zcat $base.versions | /home/audris/src/networkit/cluster $(zcat $base.names|wc -l) | gzip > $base.PLM
#modularity=0.703689 nver=89450504 clusters=15484535 largest=10400230

zcat $base.PLM | perl rankNew.perl $base 1 | gzip > $base.crank.map
base=fP2SLPStatsT
zcat /data/basemaps/gz/$base.s | awk -F\; '{if ($1!=$8 && $7/$5 > 0.1 && $7/$12 > .1 && $7 >= $14 && $7>3)print $0}' |gzip > $base.StrongReuseLinks01-3
base=$base.StrongReuseLinks01-3
zcat $base |  cut -d\; -f1,8  | perl ~/bin/connectBasic.perl $base |gzip > $base.map
zcat $base.versions | /home/audris/src/networkit/cluster $(zcat $base.names|wc -l) | gzip > $base.PLM
# modularity=0.999662 nver=5275622 clusters=1690692 largest=24923
zcat $base.PLM | perl rankNew.perl $base 1 | gzip > $base.crank.map
grep cran_urca <(zcat $base.crank.map)
-frbl_vars;cran_urca;0.029802
-Cran_vars;cran_urca;0.029802
+salsa.debian.org_edd_r-cran-urca;cran_urca;0.029802
-cheaton_vars2;cran_urca;0.029802
-cosname_rfinance;cran_urca;0.029802
+cran_urca;cran_urca;0.029802
+rforge_urca;cran_urca;0.029802
-rforge_vars;cran_urca;0.029802
+bwlewis_urca;cran_urca;0.029802

base=fP2SLPStatsT
zcat /data/basemaps/gz/$base.s | awk -F\; '{if ($1!=$8 && $7/$5 > 0.13 && $7/$12 > .13 && $7 >= $14 && $7>3)print $0}' |gzip > $base.StrongReuseLinks01-10
base=$base.StrongReuseLinks013-3
zcat $base |  cut -d\; -f1,8  | perl ~/bin/connectBasic.perl $base |gzip > $base.map
zcat $base.versions | /home/audris/src/networkit/cluster $(zcat $base.names|wc -l) | gzip > $base.PLM
# modularity=0.999931 nver=4574893 clusters=1540281 largest=9386

salsa.debian.org_edd_r-cran-urca;cran_urca;0.018294
bwlewis_urca;cran_urca;0.018294
rforge_urca;cran_urca;0.018294
cran_urca;cran_urca;0.018294


frbl_vars;Cran_vars;0.018294
cheaton_vars2;Cran_vars;0.018294
Cran_vars;Cran_vars;0.018294
rforge_vars;Cran_vars;0.018294



python3 grabGen.py > cs540vals
(cut -d\, -f2,3 ~/course/cs540-21/bigdata/bv_namsor_ethnicity_2_0_10.csv| sed 's|,| |' | tr '[A-Z]' '[a-z]'| lsort 5G -u -t\; -k1,1  | join -t\; -v1 as.fl.u as.fl.ns > as.fl.u.miss
awk '{if(length($1)>1 && length($2)>1){print $0}}' as.fl.u.miss | lsort 1G -R > as.fl.u.miss.long


base=fP2SLPStatsT
zcat /data/basemaps/gz/$base.s | awk -F\; '{if ($1!=$8 && $7/$5 > 0.12 && $7/$12 > .12 && $7 >= $14 && $7>0)print $0}' |gzip > $base.StrongReuseLinks012-0
#modularity=0.999205 nver=9674252 clusters=2750968 largest=50532

base=fP2SLPStatsT
zcat /data/basemaps/gz/$base.s | awk -F\; '{if ($1!=$8 && $7/$5 > 0.11 && $7/$12 > .11 && $7 >= $14 && $7>0)print $0}' |gzip > $base.StrongReuseLinks011-0
#modularity=0.999034 nver=10351252 clusters=2863192 largest=56649
salsa.debian.org_edd_r-cran-urca;cran_urca;0.013158
cran_urca;cran_urca;0.013158
rforge_urca;cran_urca;0.013158
bwlewis_urca;cran_urca;0.013158

base=fP2PStatsT
zcat /data/basemaps/gz/$base.s | awk -F\; '{if ($1!=$8 && $7/$5 > 0.11 && $7/$12 > .11 && $7 >= $14 && $7>3)print $0}' |gzip > $base.StrongReuseLinks011-3
#modularity=0.976502 nver=12546288 clusters=1935992 largest=762318

grep cran_urca <(zcat $base.crank.map)
salsa.debian.org_edd_r-cran-urca;cran_urca;0.000005
cran_urca;cran_urca;0.000005
bwlewis_urca;cran_urca;0.000005
rforge_urca;cran_urca;0.000005
60d7d529a95e:/data/play/forks>grep Cran_vars <(zcat $base.crank.map)                                                                                                                                                                                                                                                      
Cran_vars;Cran_vars;0.000005
cheaton_vars2;Cran_vars;0.000005
rforge_vars;Cran_vars;0.000005
frbl_vars;Cran_vars;0.000005

base=fP2PStatsT
zcat /data/basemaps/gz/$base.s | awk -F\; '{if ($1!=$8 && $7/$5 > 0.11 && $7/$12 > .11 && $7 >= $14 && $7>1)print $0}' |gzip > $base.StrongReuseLinks011-1
base=$base.StrongReuseLinks011-1
#modularity=0.980753 nver=14313409 clusters=2306614 largest=768584
rforge_urca;cran_urca;0.000005
salsa.debian.org_edd_r-cran-urca;cran_urca;0.000005
bwlewis_urca;cran_urca;0.000005
cran_urca;cran_urca;0.000005


##Woc Update, May 9 ########################################################
Highlights:
- Five hackathon write-ups + one+ MSR paper at MSR 21 
https://2021.msrconf.org/track/msr-2021-technical-papers?track=MSR%20Hackathon

- Version T was collected based on updates/new repositories identified on Feb 12, 2021, and the git objects
have been retrieved by Mar 16. Over 14M new/updated repos were identified and new data exceeded 35TB.

- This required internal changes to the database organization among the servers.

- More attributes in project/author summaries stored in mongodb

- Developer communities based on authors working on the same projects.

- Prototype of additional deforking based on blobs shared among projects.

- License spelled out (forced by the interest from industry)


Details

Almost 300M commits and 1.2B trees were added over these five months. There were over 52M distinct author IDs.
Of them, 1,029,324 were for organizations, used by more than one author, or had insufficient information for
actual author identification. From the remaining 51,326,297 author IDs 35,288,533 distinct authors were identified.
Of the 146M distinct repos 89,840,664 were not clones or forks of others. Using prototype blob-based deforking
this drops the total of independent repos by 10M to 80M. 

In addition to commit-based deforking a prototype blob-based deforking has been implemented. This brings down the
number of independent repositories to 

Number 	type
10457708737 	blob
2595586645 	commit
10687786993 	tree
20216545 	tag
146138675 	projects (distinct repositories)
52355621 	author IDs


Version T - big changes

Fully incorporated da5 (1.3TB RAM 120TB HDD 38TB SSD). Had to re-arranged databases as
blobs no longer fit on a single server.

More attributes in project/author summaries: monthly activity, core team, ...

Developer communities were created by detecting communities of developers who work together on
projects based developer/project bi-graph. 

A prototype deforking based on blob sharing was introduced. Apart from a much larger number of blobs
(than commits or authors) the blob to project graph connects 90% of projects into a single cluster.
The filtering of edges based on the direction and extent of blob reuse in addition to community detection
were needed to get an approximate clustering into independent project clusters - these projects do not share
commits. 
###########################################################


#old wactd_javase-demo;wactd_javase-demo;0;0;2;2;1;1;1;1;o;o;20;19;19;19;20;19;wactd_javase-demo
#new diogom42_d4jdata;o;2;5;303;301;1;rwang63_Arbitrage-Finder;o70;13;0;51;38;0

#this is a more general sharing through commits that create some blobs
zcat PnbFullT.s|wc -l
89450504
zcat PnbFullT.s|awk -F\; '{if ($2> 10)i10++; if ($2> 20)i20++;if ($2> 5)i5++;}END {print i5,i10,i20}'
65542299 56297593 45389449
zcat P2fP4FullT.Pstats | awk -F\; '{ if ($3>1||$5>1)ns++}END{print ns}'
74491883

Total   |with sharing
89450504|74491883
65542299|60932347
56297593|53780036
45389449|44217596

#by type
zcat P2fP4FullT.Pstats | awk -F\; '{ t[$6]++; if($7>20)t20[$6]++; if($7>10)t10[$6]++; if($7>5)t5[$6]++;} END{for (i in t)print i,t[i],t5[i],t10[i],t20[i]; }'
typ all       5+         10+      20+
o   24488920 13303031 10917881  9334732
o70 22070677 16670002 14431268 12692143
o50 13940010 13940010 11358591  6333512
u    9040410  7815879  6947209  5900359
d70  7872040  6709140  6127959  5458384
d90  3829321  3829321  3482630  3113294
d    8209126  3274916  3032055  2557025

zcat P2fP4FullT.Pstats | awk -F\; '{ if (($3>5||$5>5)&&$7>20)ns++}END{print ns}'
34357558

##################################################################
#no blobs removed
###############################################

base=P2fP4FullT.StrongReuseLinks009-0
#modularity=0.960104 nver=26952831 clusters=3080825 largest=2674717
#26952831
cran_urca;cran_urca;0.000001
rforge_urca;cran_urca;0.000001
bwlewi s_urca;cran_urca;0.000001
salsa.debian.org_edd_r-cran-urca;cran_urca;0.000001
#bad
cran_aaSEA;ajbt200128_Attendance-System;0.156510
cran_blrm;Unipisa_test_git_old_dates;1.000001
cran_ggamma;Unipisa_test_git_old_dates;1.000001
cran_mcBFtest;Unipisa_test_git_old_dates;1.000001

#w modularity=0.854355 nver=26952831 clusters=3060858 largest=8468787
cran_urca;cran_urca;0.000001
rforge_urca;cran_urca;0.000001
bwlewis_urca;cran_urca;0.000001
salsa.debian.org_edd_r-cran-urca;cran_urca;0.000001


#bad
base=P2fP4FullT.StrongReuseLinks008-1
#w
#modularity=0.837854 nver=18594745 clusters=2546876 largest=1450086
cheaton_vars2;cran_urca;0.000006
cran_urca;cran_urca;0.000006
bwlewis_urca;cran_urca;0.000006
salsa.debian.org_edd_r-cran-urca;cran_urca;0.000006
rforge_vars;cran_urca;0.000006
Cran_vars;cran_urca;0.000006
rforge_urca;cran_urca;0.000006
frbl_vars;cran_urca;0.000006

#modularity=0.971306 nver=18594745 clusters=2548118 largest=1128437
zcat $base.crank.map|grep cran_urca
cheaton_vars2;cran_urca;0.000006
cran_urca;cran_urca;0.000006
bwlewis_urca;cran_urca;0.000006
salsa.debian.org_edd_r-cran-urca;cran_urca;0.000006
rforge_vars;cran_urca;0.000006
Cran_vars;cran_urca;0.000006
rforge_urca;cran_urca;0.000006
frbl_vars;cran_urca;0.000006

#bad: BUT SHOWS WHERE IT IS CONVERGING TO
base=P2fP4FullT.StrongReuseLinks007-4
#w modularity=0.820646 nver=16345637 clusters=2001156 largest=1682233
zcat ${base}w.crank.map|grep Cran_vars
cheaton_vars2;Cran_vars;0.000008
cran_urca;Cran_vars;0.000008
cosname_rfinance;Cran_vars;0.000008
bwlewis_urca;Cran_vars;0.000008
salsa.debian.org_edd_r-cran-urca;Cran_vars;0.000008
rforge_vars;Cran_vars;0.000008
Cran_vars;Cran_vars;0.000008
rforge_urca;Cran_vars;0.000008
frbl_vars;Cran_vars;0.000008

#modularity=0.956035 nver=16345637 clusters=2002752 largest=1128358
zcat $base.crank.map|grep cran_urca
cran_urca;Cran_vars;0.000008
cheaton_vars2;Cran_vars;0.000008
cran_urca;Cran_vars;0.000008
cosname_rfinance;Cran_vars;0.000008
bwlewis_urca;Cran_vars;0.000008
salsa.debian.org_edd_r-cran-urca;Cran_vars;0.000008
rforge_vars;Cran_vars;0.000008
Cran_vars;Cran_vars;0.000008
rforge_urca;Cran_vars;0.000008
frbl_vars;Cran_vars;0.000008


#bad
base=P2fP4FullT.StrongReuseLinks005-5
#modularity=0.936698 nver=17502650 clusters=1914171 largest=1366398
zcat $base.crank.map|grep cran_urca
cheaton_vars2;cran_urca;0.000011
cran_urca;cran_urca;0.000011
cosname_rfinance;cran_urca;0.000011
bwlewis_urca;cran_urca;0.000011
salsa.debian.org_edd_r-cran-urca;cran_urca;0.000011
rforge_vars;cran_urca;0.000011
Cran_vars;cran_urca;0.000011
githubfun_trading;cran_urca;0.000011
polinas123_vars;cran_urca;0.000011
rforge_urca;cran_urca;0.000011
frbl_vars;cran_urca;0.000011

#bad
base=P2fP4FullT.StrongReuseLinks005-3
#modularity=0.944905 nver=21021539 clusters=2146760 largest=1422192
zcat $base.crank.map|grep cran_urca             
cheaton_vars2;cran_urca;0.000009
cran_urca;cran_urca;0.000009
cosname_rfinance;cran_urca;0.000009
bwlewis_urca;cran_urca;0.000009
salsa.debian.org_edd_r-cran-urca;cran_urca;0.000009
rforge_vars;cran_urca;0.000009
Cran_vars;cran_urca;0.000009
githubfun_trading;cran_urca;0.000009
polinas123_vars;cran_urca;0.000009
rforge_urca;cran_urca;0.000009
frbl_vars;cran_urca;0.000009

#real bad
base=P2fP4FullT.StrongReuseLinks001-1
#modularity=0.881736 nver=37219897 clusters=2055694 largest=2850681
cran_urca;zkan_ml-class-2011;0.010918



base=P2fP4FullT.StrongReuseLinks013-3
#modularity=0.981177 nver=10690866 clusters=1796640 largest=543534
cran_socceR;oddnoc_oddnoc.github.io;0.028821
cran_aaSEA;ajbt200128_Attendance-System;0.446403

base=P2fP4FullT.StrongReuseLinks011-3
#modularity=0.976504 nver=12546288 clusters=1935954 largest=762300

#good
base=P2fP4FullT.StrongReuseLinks01-10
#modularity=0.962836 nver=9946152 clusters=1379003 largest=794267
zcat $base.crank.map|grep cran_urca
cran_urca;cran_urca;0.000006
rforge_urca;cran_urca;0.000006
bwlewis_urca;cran_urca;0.000006
salsa.debian.org_edd_r-cran-urca;cran_urca;0.000006

da5:/data/play/forks>zcat $base.crank.map|grep Cran_vars
rforge_vars;Cran_vars;0.000006
Cran_vars;Cran_vars;0.000006
cheaton_vars2;Cran_vars;0.000006
frbl_vars;Cran_vars;0.000006

#good
base=P2fP4FullT.StrongReuseLinks01-5
#modularity=0.969525 nver=12188240 clusters=1757833 largest=804473
cran_urca;cran_urca;0.000005
rforge_urca;cran_urca;0.000005
bwlewis_urca;cran_urca;0.000005
salsa.debian.org_edd_r-cran-urca;cran_urca;0.000005

#good
base=P2fP4FullT.StrongReuseLinks01-3
#modularity=0.972874 nver=13539286 clusters=1998490 largest=804600
cran_urca;cran_urca;0.000004
rforge_urca;cran_urca;0.000004
bwlewis_urca;cran_urca;0.000004
salsa.debian.org_edd_r-cran-urca;cran_urca;0.000004

rforge_vars;Cran_vars;0.000004
Cran_vars;Cran_vars;0.000004
cheaton_vars2;Cran_vars;0.000004
frbl_vars;Cran_vars;0.000004

#good
base=P2fP4FullT.StrongReuseLinks01-1
zcat /data/basemaps/gz/$base | cut -d\; -f15 |gzip > $base.weights
#modularity=0.977579 nver=15450546 clusters=2383483 largest=805010
cran_urca;cran_urca;0.000004
rforge_urca;cran_urca;0.000004
bwlewis_urca;cran_urca;0.000004
salsa.debian.org_edd_r-cran-urca;cran_urca;0.000004

#w modularity=0.865542 nver=15450546 clusters=2382562 largest=776147
zcat ${base}w.crank.map|grep Cran_vars
rforge_vars;Cran_vars;0.000004
Cran_vars;Cran_vars;0.000004
cheaton_vars2;Cran_vars;0.000004
frbl_vars;Cran_vars;0.000004
da5:/data/play/forks>zcat ${base}w.crank.map|grep cran_urca
cran_urca;cran_urca;0.000004
rforge_urca;cran_urca;0.000004
bwlewis_urca;cran_urca;0.000004
salsa.debian.org_edd_r-cran-urca;cran_urca;0.000004

base=P2fP4FullT.StrongReuseLinks01-0
#modularity=0.963971 nver=24011881 clusters=2975268 largest=2314248
#w modularity=0.865948 nver=24011881 clusters=2957358 largest=7140109
cran_urca;cran_urca;0.000001
rforge_urca;cran_urca;0.000001
bwlewis_urca;cran_urca;0.000001
salsa.debian.org_edd_r-cran-urca;cran_urca;0.000001


#good - use u as reference, but need to separate large:
base=P2fP4FullT.StrongReuseLinks009-1
# how to fix ekstroem_socceR and cran_aaSEA: first copied some fonts
#cran_aaSEA;rohitchandrashekar_cran-package;0;0;1;6;2;0;2;2;o;d70;72;72;66;0;129;25;;0
#modularity=0.974726 nver=17137840 clusters=2475038 largest=1063056
#17137840
cran_urca;cran_urca;0.000004
rforge_urca;cran_urca;0.000004
bwlewis_urca;cran_urca;0.000004
salsa.debian.org_edd_r-cran-urca;cran_urca;0.000004

cran_iq;ndsol_volcano;0.009685#need to separate
cran_socceR;jmandel_jmandel.github.io;0.027642#need to separate
cran_aaSEA;ajbt200128_Attendance-System;0.559323#need to separate


#w modularity=0.851943 nver=17137840 clusters=2474072 largest=992599
#17137840
cran_urca;cran_urca;0.000004
rforge_urca;cran_urca;0.000004
bwlewis_urca;cran_urca;0.000004
salsa.debian.org_edd_r-cran-urca;cran_urca;0.000004



#####################################
# Only code blobs are tracked
base=P2SfP4FullT.StrongReuseLinks013-3
#modularity=0.865994 nver=16468864 clusters=1351953 largest=1575364


base=P2SfP4FullT.StrongReuseLinks005-0
#modularity=0.825861 nver=47677535 clusters=1022173 largest=9372528

base=P2SfP4FullT.StrongReuseLinks009-3
#modularity=0.948164 nver=6533043 clusters=943348 largest=583918
#6533043
zcat ${base}.crank.map| grep cran_urca
#good and only one bad:cran_aaSEA
cran_urca;cran_urca;0.000005
rforge_urca;cran_urca;0.000005
bwlewis_urca;cran_urca;0.000005
salsa.debian.org_edd_r-cran-urca;cran_urca;0.000005

cran_LogicOpt;synthemesc_espresso;0.000036
cran_iq;ndsol_volcano;0.010715 #toshikurauchi_eyeswipe cran_iq share
a bunch of math functions: SparseLUImpl.h MathFunctions.h SparseDenseProduct.h
cran_aaSEA;ajbt200128_Attendance-System;0.683965


base=P2SfP4FullT.StrongReuseLinks013-3
#modularity=0.957997 nver=4178341 clusters=706298 largest=442965
#looses  bwlewis_urca

#wrst ar so-so
cran_iq;reaperhui_CyberSystem_Old;0.007569
cran_aaSEA;ajbt200128_Attendance-System;0.403054
####################
base=P2SfP4FullT.StrongReuseLinks009-3
#modularity=0.948164 nver=6533043 clusters=943353 largest=583918
#6533043
#gains bwlewis_urca
#wrst ar so-so
cran_iq;ndsol_volcano;0.010715
cran_aaSEA;ajbt200128_Attendance-System;0.683965
####################


base=P2SfP4FullT.StrongReuseLinks009-0
#modularity=0.922403 nver=12552042 clusters=1335909 largest=2474981
#12552042
#more wrsr
cran_bayest;ciyanogen_ciyanogen;0.099325
cran_ggpointdensity;ciyanogen_ciyanogen;0.099325
cran_aaSEA;ajbt200128_Attendance-System;0.143641
cran_blrm;Unipisa_test_git_old_dates;1.000000
cran_mcBFtest;Unipisa_test_git_old_dates;1.000000

base=P2SfP4FullT.StrongReuseLinks008-3
#modularity=0.942653 nver=7320733 clusters=1012249 largest=803247
#7320733
#wrst ar so-so
cran_iq;ndsol_volcano;0.010865
cran_aaSEA;ajbt200128_Attendance-System;0.759528


base=P2SfP4FullT.StrongReuseLinks007-3
#modularity=0.939079 nver=8454123 clusters=1081411 largest=843746
#8454123
#bad, mixes in Cran_vars
cran_urca;Cran_vars;0.000010
#wrst get worse
cran_Rcpp;ndsol_volcano;0.010775
cran_iq;ndsol_volcano;0.010775
cran_aaSEA;ajbt200128_Attendance-System;0.827969


#################################################
#################################################
#Most Conservative!!!! ONLY code blobs used in fewer than 100 Ps 
#too restrictive?
base=P2SSfP4FullT.StrongReuseLinks013-3
#modularity=0.999993 nver=1853085 clusters=717811 largest=346
zcat ${base}.crank.map| grep cran_urca 
#looses  bwlewis_urca
cran_urca;cran_urca;0.014185
rforge_urca;cran_urca;0.014185
salsa.debian.org_edd_r-cran-urca;cran_urca;0.014185

#worst pretty good
cran_rgl;dmurdoch_rgl;0.021379
cran_vegan;vegandevs_vegan;0.021379
cran_prodlim;rforge_eventhistory;0.021429
cran_QRM;bpfaff_FRAPO;0.028474
cran_EcoHydRology;mlt_swat;0.035824



base=P2SSfP4FullT.StrongReuseLinks008-2
#modularity=0.999982 nver=3220610 clusters=1180930 largest=3038
#3220610
cran_nametagger;scignscape_kauv-rz;0.058417
cran_EcoHydRology;hawklorry_swat;0.098469
cran_LogicOpt;synthemesc_espresso;0.129174


base=P2SSfP4FullT.StrongReuseLinks008-1 
#modularity=0.999985 nver=3624722 clusters=1318164 largest=3306
#3624722
cran_nametagger;scignscape_kauv-rz;0.058417
cran_EcoHydRology;hawklorry_swat;0.098469
cran_LogicOpt;synthemesc_espresso;0.129174

#choose this
base=P2SSfP4FullT.StrongReuseLinks008-0
#4705944
#modularity=0.999941 nver=4705944 clusters=1602831 largest=4896
zcat ${base}.crank.map| grep cran_urca
cran_urca;cran_urca;0.019356
rforge_urca;cran_urca;0.019356
bwlewis_urca;cran_urca;0.019356
salsa.debian.org_edd_r-cran-urca;cran_urca;0.019356

#wrst great:
cran_nametagger;scignscape_kauv-rz;0.058417
cran_EcoHydRology;hawklorry_swat;0.098469
cran_LogicOpt;synthemesc_espresso;0.129174

#rest bad

base=P2SSfP4FullT.StrongReuseLinks007-3
#modularity=0.999970 nver=3252981 clusters=1177751 largest=4453
#3252981
#bad, mixes in Cran_vars
cran_urca;Cran_vars;0.038305

base=P2SSfP4FullT.StrongReuseLinks005-3
#modularity=0.999906 nver=4056149 clusters=1401330 largest=8254
#4056149
#bad...
cheaton_vars2;cran_urca;0.050728
cran_urca;cran_urca;0.050728
cosname_rfinance;cran_urca;0.050728
bwlewis_urca;cran_urca;0.050728
salsa.debian.org_edd_r-cran-urca;cran_urca;0.050728
rforge_vars;cran_urca;0.050728
Cran_vars;cran_urca;0.050728
githubfun_trading;cran_urca;0.050728
polinas123_vars;cran_urca;0.050728
rforge_urca;cran_urca;0.050728
frbl_vars;cran_urca;0.050728

base=P2SSfP4FullT.StrongReuseLinks005-0
#modularity=0.999749 nver=6972326 clusters=2194036 largest=16288
#not worth checking




base=P2SSfP4FullT.StrongReuseLinks001-1
#modularity=0.995126 nver=10903476 clusters=2615113 largest=129856
#bad
cran_urca;gitlab.com_conradsnicta_armadillo-code;0.421947

base=P2SSfP4FullT.StrongReuseLinks003-1
#modularity=0.999684 nver=7111406 clusters=2203582 largest=21758
bad and bad
cran_urca;cran_urca;0.041208
...

cran_PamBinaries;fotlogo_cv10;0.568654
cran_fastcluster;fotlogo_cv10;0.568654
cran_rTRNG;DeeHants_liveMedia;0.977680




####################################################

zcat /data/basemaps/gz/$base | grep ';[ou];[oud];' | cut -d\; -f1,2 | perl ~/bin/connectBasic.perl ${base}o | gzip > ${base}o.map &
zcat ${base}o.versions | /home/audris/src/networkit/cluster $(zcat ${base}o.names|wc -l) | gzip > ${base}o.PLM
zcat ${base}o.PLM | perl rankNew.perl ${base}o 1 | gzip > ${base}o.crank.map
zcat ${base}o.crank.map|grep '^cran_' | grep -v ';[cC]ran_'|lsort 1G -t\; -k3 -n
#o modularity=0.974226 nver=16697182 clusters=2386354 largest=1034920
#Still an issue
cran_iq;ndsol_volcano;0.009685
cran_socceR;oddnoc_oddnoc.github.io;0.022010
cran_aaSEA;ajbt200128_Attendance-System;0.559323




#perhaps too long?

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
c=10
zcat P2PFullS.nb2b.$c.s | perl connectExportVw.perl P2PFullS.nb2b.$c 
zcat P2PFullS.nb2b.$c.versions | paste  -d\   - <(zcat P2PFullS.nb2b.$c.weights) >  P2PFullS.nb2b.$c.csv 
cat P2PFullS.nb2b.$c.csv | /home/audris/src/networkit/clusterw $(zcat /data/P2PFullS.nb2b.$c.names|wc -l) $(zcat /data/P2PFullS.nb2b.$c.weights|wc -l) | gzip > /data/P2PFullS.nb2bw.$c.PLM
modularity=0.968289 nver=15508231 clusters=1236002 largest=836106
zcat P2PFullS.nb2bw.$c.PLM | perl rankNew.perl P2PFullS.nb2b.$c 1 | gzip >  P2PFullS.nb2bw.$c.crank.map
zcat P2PFullS.nb2b.$c.versions|ssh da3 "bin/connect" | gzip > P2PFullS.nb2b.$c.clones
zcat P2PFullS.nb2b.$c.clones | cut -d\; -f2 | lsort 30G | uniq -c | lsort 1G -rn | head -20
perl getComponent.perl  P2PFullS.nb2b.$c P2PFullS.nb2bw.$c.PLM 0 > P2PFullS.nb2bw.$c.0.forOSLO
time /home/audris/src/OSLOM2/oslom_undir -f /data/play/forks/P2PFullS.nb2b.$c.0.csv -hint  /data/play/forks/P2PFullS.nb2bw.$c.0.forOSLO -w -r 1 -hr 1


c=2000

#do blob sharing where blobs are for programming languages and of at least certain size (450 char long)
zcat P2PFullS.nb2b.2000.s | perl connectExportVw.perl P2PFullS.nb2b.2000 
zcat P2PFullS.nb2b.2000.versions | paste  -d\   - <(zcat P2PFullS.nb2b.2000.weights) >  P2PFullS.nb2b.2000.csv 
cat P2PFullS.nb2b.2000.csv | /home/audris/src/networkit/clusterw $(zcat /data/P2PFullS.nb2b.2000.names|wc -l) $(zcat /data/P2PFullS.nb2b.2000.weights|wc -l) | gzip > /data/P2PFullS.nb2bw.2000.PLM
#modularity=0.792458 nver=30128959 clusters=2805378 largest=6652122
zcat P2PFullS.nb2bw.2000.PLM | perl rankNew.perl P2PFullS.nb2b.2000 1 | gzip >  P2PFullS.nb2bw.2000.crank.map  
zcat P2PFullS.nb2b.2000.versions|ssh da3 "bin/connect" | gzip > P2PFullS.nb2b.2000.clones
zcat P2PFullS.nb2b.2000.clones | cut -d\; -f2 | lsort 30G | uniq -c | lsort 1G -rn | head -20
28196096 0
   4181 97898
   1992 78
    940 84

perl getComponent.perl  P2PFullS.nb2b.2000 P2PFullS.nb2bw.2000.PLM 0 > P2PFullS.nb2bw.2000.0.forOSLO
time /home/audris/src/OSLOM2/oslom_undir -f /data/play/forks/P2PFullS.nb2b.2000.0.csv -hint  /data/play/forks/P2PFullS.nb2bw.2000.0.forOSLO -w -r 1 -hr 1

#crashes, do weight 2+
grep -v ' 1$' P2PFullS.nb2b.2000.csv > P2PFullS.nb2b.2000.2p.csv
awk '{print $1, $2}' P2PFullS.nb2b.2000.2p.csv | ssh da3 "bin/connect" | gzip > P2PFullS.nb2b.2000.2p.clones
zcat P2PFullS.nb2b.2000.2p.clones|cut -d\; -f2 | lsort 30G | uniq -c | lsort 1G -rn | head -20
7232201 5
    578 11720465
    553 16626863
    402 9530741
    
cat P2PFullS.nb2b.2000.2p.csv | /home/audris/src/networkit/clusterw $(zcat /data/P2PFullS.nb2b.2000.names|wc -l) $(zcat /data/P2PFullS.nb2b.2000.weights|wc -l) | gzip > /data/P2PFullS.nb2bw.2000.2p.PLM
#modularity=0.910744 nver=30128959 clusters=21187215 largest=748110
perl getComponent.perl P2PFullS.nb2b.2000.2p P2PFullS.nb2bw.2000.2p.PLM 5 > P2PFullS.nb2bw.2000.2p.5.forOSLO
time /home/audris/src/OSLOM2/oslom_undir -f /data/play/forks/P2PFullS.nb2b.2000.2p.5.csv -hint  /data/play/forks/P2PFullS.nb2bw.2000.2p.5.forOSLO -w -r 1 -hr 1
***************************************************************************
CHECK UNIONS AND SIMILAR MODULES DONE
******** module_collection ******** 34 modules. writing... 
DONE   ****************************
pruning all the modules collected. Partitions found: 1
getting partition from tp-file: /data/play/forks/P2PFullS.nb2b.2000.2p.5.csv_oslo_files/partitions_level_4
34 groups found
34 bss found
checking homeless nodes
writing final solution in file /data/play/forks/P2PFullS.nb2b.2000.2p.5.csv_oslo_files/short_tp4
******** module_collection ******** 34 modules. writing... 
DONE   ****************************
hierarchies done ********* 

real    84519m14.303s
user    83321m34.916s
sys     995m46.387s



#no need to recode, works on original ids
#cat P2PFullS.nb2b.2000.2p.5.csv | perl connectExportW.perl P2PFullS.nb2b.2000.2p.5
#zcat P2PFullS.nb2b.2000.2p.5.versions | paste  -d\   - <(zcat P2PFullS.nb2b.2000.2p.5.weights) >  P2PFullS.nb2b.2000.2p.5recoded.csv 
#time /home/audris/src/OSLOM2/oslom_undir -f /data/play/forks/P2PFullS.nb2b.2000.2p.5recoded.csv -w -r 1 -hr 1

# do blob knowledge flow where blobs are for programming languages and of at least certain size (450 char long)
# and the next commiter is at least 30 days later
zcat A2AFullS.nfb.30.s | perl connectExportW.perl A2AFullS.nfb.30
zcat A2AFullS.nfb.30.versions | paste  -d\  - <(zcat A2AFullS.nfb.30.weights) >  A2AFullS.nfb.30.csv
cat A2AFullS.nfb.30.csv | /home/audris/src/networkit/clusterw $(zcat /data/A2AFullS.nfb.30.names|wc -l) $(cat /data/A2AFullS.nfb.30.csv|wc -l) | gzip > /data/A2AFullS.nfbw.30.PLM
#modularity=0.656848 nver=11980180 clusters=87615 largest=3103190
zcat A2AFullS.nfbw.30.PLM | perl rankNew.perl A2AFullS.nfb.30 1 | gzip >  A2AFullS.nfbw.30.crank.map  
zcat A2AFullS.nfbw.30.crank.map|cut -d\; -f2 | lsort 30G | uniq -c | lsort 1G -rn | head
3103190 ralp-sms <wyred.innovations@gmail.com>
2128103 tea9 <shaomiaom@sina.com>
1046593 Anurag <ajsanuragjain@gmail.com>
 902132 Nathanaël <nathanael.spriet@gmail.com>
 732357 adedayo2017 <maxistinlove@gmail.com>
 642663 zertosh <zertosh@gmail.com>
 479082 JOLY Clement <clement.joly@telecomnancy.eu>
 334570 Oskar Hane <oh@oskarhane.com>
 327431 Unity Technologies <@unity.com>
 317151 WataruSUzuki <wataru0406@gmail.com>

zcat A2AFullS.nfb.30.versions|ssh da3 "bin/connect" | gzip > A2AFullS.nfb.30.clones
zcat A2AFullS.nfb.30.clones | cut -d\; -f2 | lsort 30G | uniq -c | lsort 1G -rn | head -20
11788641 0
    199 5857
    131 6619
     76 16845
perl getComponent.perl A2AFullS.nfb.30 A2AFullS.nfbw.30.PLM 0 > A2AFullS.nfbw.30.0.forOSLO
time /home/audris/src/OSLOM2/oslom_dir -f /data/play/forks/A2AFullS.nfb.30.0.csv -hint  /data/play/forks/A2AFullS.nfbw.30.0.forOSLO -w -r 1 -hr 1


zcat A2AFullS.nfb.0.s | perl connectExportW.perl A2AFullS.nfb.0
zcat A2AFullS.nfb.0.versions | paste  -d\ - <(zcat A2AFullS.nfb.0.weights) >  A2AFullS.nfb.0.csv
cat A2AFullS.nfb.0.csv | /home/audris/src/networkit/clusterw $(zcat /data/A2AFullS.nfb.30.names|wc -l) $(cat /data/A2AFullS.nfb.30.csv|wc -l) | gzip > /data/A2AFullS.nfbw.30.PLM

zcat A2AFullS.nfb.30.s |cut -d\; -f2 |lsort 130G |uniq -c | lsort 1G -rn |head
2695980  <cinnabar@git>
1288126 Nathanaël <nathanael.spriet@gmail.com>
1152753 springsw <web@springsw.com>
1038209 gituser <git@gituser.com>
 987231 a <986945193@QQ.com>
 985800 ralp-sms <wyred.innovations@gmail.com>
 897800 lijinhai255 <18854109500@163.com>
 870039 chuchuyun <happydays0509@gmail.com>
 854262 Anurag <ajsanuragjain@gmail.com>
 804587 Andreas Amsenius <andreas@amsenius.se>

 
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
time /home/audris/src/OSLOM2/oslom_undir -f /data/play/forks/P2AFullS.nA2AP.$c.376697.csv -hint forOSLOP.376697 -w -r 1 -hr 1

time /home/audris/src/OSLOM2/oslom_undir -f /data/play/forks/P2AFullS.nA2AP.$c.0.csv -hint forOSLOP.0 -w -r 1 -hr 1
***************************************************************************
CHECK UNIONS AND SIMILAR MODULES DONE
******** module_collection ******** 203 modules. writing... 
DONE   ****************************
pruning all the modules collected. Partitions found: 1
getting partition from tp-file: /data/play/forks/P2AFullS.nA2AP.2000.0.csv_oslo_files/partitions_level_4
203 groups found
203 bss found
checking homeless nodes
writing final solution in file /data/play/forks/P2AFullS.nA2AP.2000.0.csv_oslo_files/short_tp4
******** module_collection ******** 208 modules. writing... 
DONE   ****************************
hierarchies done ********* 
real    38641m57.425s
user    37641m54.580s
sys     747m20.501s

cat P2AFullS.nA2AP.$c.0.csv_oslo_files/tp | perl rankNewO.perl P2AFullS.nA2AP.${c} P2AFullS.nA2APw.$c.PLM | gzip >  P2AFullS.nA2APOw.${c}.0.crank.map
perl ./toGdfO.perl P2AFullS.nA2APOw.${c}.0.crank.map P2AFullS.nA2AP.${c}.names P2AFullS.nA2AP.$c.0.csv 0 > PP0.0.gdf
cat P2AFullS.nA2AP.${c}.0.csv_oslo_files/tp | perl ./toGdfOC.perl P2AFullS.nA2APOw.${c}.0.crank.map P2AFullS.nA2AP.${c}.names P2AFullS.nA2AP.${c}.0.csv > PPcmnt.gdf 
#investigate the results
zcat P2AFullS.nA2APOw.${c}.0.crank.map|wc -l
15340747
zcat P2AFullS.nA2APOw.${c}.0.crank.map|cut -d\; -f2|lsort 50G -t\; -u | wc -l
883939
zcat P2AFullS.nA2APOw.${c}.0.crank.map|cut -d\; -f2|lsort 50G | uniq -c | lsort 5G -rn | head -20 >top20
head top20
  28659 Tomster <tomster@emberjs.com>
   7934 CarusoGabriel <carusogabriel34@gmail.com>
   7710 Robin Stocker <robin@nibor.org>
   7653 cmr <corey@octayn.net>
   7420 Steve <shade@chemlab.org>
   7322 pascal <pascal@borreli.com>
   7065 onovy <novy@ondrej.org>
   6349 a <Luca@Milanesio.org>
   6326 Fred Zirdung <fred@hackreactor.com>
   5644 edward <edward@4angle.com>

sed 's|^\s*[0-9]* ||'  top20 | while IFS=\; read i; do zcat P2AFullS.nA2APOw.${c}.0.crank.map |grep ";$i\$" | cut -d\; -f1 | ~/lookup/getValues -f A2P |cut -d\; -f2 | lsort 20G -t\; -k1,1 | perl $HOME/bin/grepFieldv.perl PManyA2000.gz 1 |uniq -c |lsort 1G -rn |head -10 > $j.0; j=$((j+1));done
head 0.0
    564 0cool321_ember.js
    480 1123456789011_data
    436 0xadada_ember-cli
    394 04th0809_website
    235 Acidburn0zzz_guides-1
    215 0xadada_ember-simple-auth
    215 0000marcell_ember-cli-mirage
    173 Alonski_guides-source
    169 0000marcell_ember-power-select
    164 22a_ember-template-lint-next-line
head 1.0
    501 php_doc-en
    491 013cbr_FOSUserBundle
    465 php-doc_en
    464 0mars_SonataAdminBundle
    398 00577_PHPExcel
    388 00F100_composer
    341 20uf_http-foundation
    330 3mg_DoctrineExtensions
    262 php_doc-base
    256 0mars_serializer
head 2.0
   1292 bitbucket.org_fusiontestaccount_fusion-renaissance-testing
    989 bitbucket.org_fusiontestaccount_ed-repo
    804 bitbucket.org_fusiontestaccount_rrv2-testing
    660 bitbucket.org_fusiontestaccount_fusiontest2
    603 fusiontestaccount_fusedsl
    337 bitbucket.org_fusiontestaccount_support-readiness
    208 bitbucket.org_mrdon_atlassian-plugins
    196 fusiontestaccount_afouad-repo2
    194 bitbucket.org_200ok_aui
    187 bitbucket.org_atlassian_atlassian-streams
head 3.0
    514 AdamDang_community
    457 000s_https-everywhere
    450 07mehmet0_origin
    422 01deyishu_ingress-nginx
    387 Ark-kun_test-infra
    387 AmitRoushan_cloud-provider-openstack
    386 1235gopalt_openshift-ansible
    338 115100_kops
    313 BrandWang_API
    307 06094051_prometheus


# layout/embed: use ICSE embeddings/create new
python3 ./convert.py --format weighted_edgelist /data/P2AFullS.nA2A.2000.0.csv /data/P2AFullS.nA2AP.2000.0.bcsr  
./verse-weighted -input /data/P2AFullS.nA2AP.2000.0.bcsr -output /data/P2AFullS.nA2AP.2000.csv.0.bin
Calculations took 6359170.50 s to run
real    105986m28.939s
user    96640m48.335s
sys     8954m10.519s
#look at 376697?

#try second largest disconnected set
cat P2AFullS.nA2AP.$c.376697.csv_oslo_files/tp | perl rankNewO.perl P2AFullS.nA2AP.${c} P2AFullS.nA2APw.$c.PLM | gzip >  P2AFullS.nA2APOw.${c}.376697.crank.map
perl ./toGdfO.perl P2AFullS.nA2APOw.2000.376697.crank.map P2AFullS.nA2AP.2000.names P2AFullS.nA2AP.$c.376697.csv 0 > 376697.0.gdf



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

# Use Leuwen as input for much more accurate OSLOM
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