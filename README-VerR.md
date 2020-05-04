# Resolving forks for version R

Unlike for ver Q, no github fork info is used in this case and the 
full graph of project to commit used in conjuction with Louvain modularity
as implemented in networkit (see cluster.c)

```
#first create one record per commit
for j in {0..31}
do zcat c2pFull${ver}$j.s | perl $HOME/lookup/connectExportPreNoExclude.perl | gzip > c2pFull$ver$j.p2p &
done
wait
echo c2pFull$ver$j.p2p

#sort projects within commit
for j in {0..31}
do zcat c2pFull$ver$j.p2p | perl -ane 'chop();($c,@x) = split(/;/); print "".(join ";", (sort @x))."\n";' | gzip > c2pFull$ver$j.p2p.gz &
done
wait
echo c2pFull$ver$j.p2p.gz

# sort lines so we can aggregate commits touching the same projects
for j in {0..31}
do zcat c2pFull$ver$j.p2p.gz | lsort  ${maxM}M -t\| | uniq -c | perl -ane 'chop();s|^\s*([0-9]*) ||;print "$_;$1\n"' | gzip > c2pFull$ver$j.p2p.s &
done
wait
echo c2pFull$ver$j.p2p.s

#merge pieces
str="lsort $(($maxM*32))M -t\| -u --merge"
for j in {0..31}
do str="$str <(zcat  c2pFull$ver$j.p2p.s)"
done
eval $str | gzip > c2pFull$ver.p2p.s

#add weights (commit counts from each piece
zcat  c2pFull$ver.p2p.s | perl -e '$pstr="";$nn=0;while(<STDIN>){chop;(@x)=split(/;/);$n=pop @x;$str=join ";", @x; if ($pstr ne $str && $pstr ne ""){ print "$nn;$pstr\n";$pstr=$str;$nn=$n }else{ $pstr=$str;$nn+=$n}};print "$nn;$pstr\n";' | gzip > c2pFull$ver.np2p
```

The resulting  weighted graph is exported into the format suitable for
PLM, PLM followed by PLP is applied (PLP does not maka any improvements so only 
PLM is used)
```
zcat c2pFull$ver.np2p | perl connectExportVw.perl c2pFull$ver.np2p1
cd ~/src/networkit
zcat c2pFullR.np2p1.versions| ./cluster 100564768 | gzip > c2pFullR.np2p1.PLMPLP 
modularity=0.960559 nver=100564768 clusters=9424393 largest=334170
modularity=0.960559 nver=100564768 clusters=9424393 largest=334170
#extract PLM membership
zcat c2pFullR.np2p1.PLMPLP|head -100564770 | grep -v '^[BE]' | gzip > c2pFull$ver.np2p1.PLM
# label each node by rank (dyn katz centrality)
zcat c2pFullR.np2p1.PLM | perl rank.perl c2pFullR.np2p1 | gzip > c2pFullR.np2p1.PLM.crank.map
#and get the map!
zcat  c2pFullR.np2p1.PLM.crank.map | awk -F\; '{if ($2 != $1)print "$1";"$2 }' | gzip >  c2pFullR.np2p1.PLMmap.forks
```






