#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>

#include <networkit/graph/Graph.hpp>

#include <networkit/structures/Partition.hpp>
#include <networkit/community/CommunityDetectionAlgorithm.hpp>
#include <networkit/community/Modularity.hpp>
#include <networkit/community/PLM.hpp>
#include <networkit/community/PLP.hpp>
#include <networkit/community/LPDegreeOrdered.hpp>
#include <networkit/community/ParallelAgglomerativeClusterer.hpp>
#include <networkit/community/CutClustering.hpp>
#include <networkit/community/LPDegreeOrdered.hpp>
#include <networkit/centrality/DynKatzCentrality.hpp>

int main (int argc, char ** argv){
  long int nnodes = atol (argv[1]);
  NetworKit::Graph g (nnodes, false, false);
  long int fr, to, j = 0;
  while (!feof(stdin)){
    int res = fscanf (stdin, "%li %li\n", &fr, &to);
    if (res != 2) fprintf (stderr, "could not read %s at %li\n", argv[1], j);
    g .addEdge (fr, to);
    //printf ("%li %li\n", fr, to);
    if (!(j++ % 1000000000)) fprintf (stderr, "%li edges read\n", j);
  }
  fprintf (stderr, "%li edges read\n", j);
  long int ne = j;
  NetworKit::DynKatzCentrality c (g, nnodes);
  c. run ();
  fprintf (stderr, "calculated centrality\n");
  NetworKit::PLM m (g); // PLM modularity=0.960559 nver=100564796 clusters=9396335 largest=466522
  //modularity=0.960559 nver=100564796 clusters=9396336 largest=466522
  //
  //NetworKit::PLP m = NetworKit::PLP (g); // does not break up giant cluster; PLP modularity=0.893339 nver=100564796 clusters=9609636 largest=6893550
  //NetworKit::LPDegreeOrdered m = NetworKit::LPDegreeOrdered (g); // modularity=0.956888 nver=100564796 clusters=11349658 largest=291440 runtime 123m
  //NetworKit::ParallelAgglomerativeClusterer m = NetworKit::ParallelAgglomerativeClusterer (g); // this takes 700+GRAM, takes very long 1000; terminate called after throwing an instance of 'std::invalid_argument'
  //  what():  G has self-loops and cannot be processed
  //  try without self-link on the largest cluster c2pFullR.np2p.largesComp.versions.s
  if (0){
    FILE * wg = fopen("w", "r");
    std::vector<NetworKit::edgeweight> w = std::vector<NetworKit::edgeweight>(ne);
    for (j = 0; j < ne; j++){
      double ww;
      int res = fscanf (stdin, "%lf\n", &ww);
      if (res != 1) fprintf (stderr, "could not read w at %li\n", j);
      w [j] = ww;
      if (!(j++ % 1000000000)) fprintf (stderr, "%li weights read\n", j);
    }
    fprintf (stderr, "read weights\n");
    NetworKit::CutClustering m = NetworKit::CutClustering (g, 100); //pick some value
  }
  m .run ();
  NetworKit::Partition p = m .getPartition ();
  NetworKit::Modularity q = NetworKit::Modularity ();
  double quality = q .getQuality (p, g);
  std::map <uint64_t, uint64_t>  ss = p .subsetSizeMap ();
  std::map <uint64_t, uint64_t>::iterator iss = ss .begin ();
  long int mx = 0, which;
  while (iss != ss .end ()){
    if (mx < iss ->second){
      which = iss ->first;
      mx  = iss ->second;
    }
    iss++;
  }
  fprintf (stderr, "modularity=%lf nver=%li clusters=%li largest=%li\n", quality, p .numberOfElements (), p .numberOfSubsets (), mx);
  std::vector<uint64_t> v = p .getVector ();
  printf ("BEGIN %li\n", v .size ());
  for (int i = 0; i < v .size (); i++){
    printf ("%li;%lf\n", v [i], c.score(i));
  }
  printf ("END\n");
  NetworKit::PLP m1 = NetworKit::PLP (g, p);
  m1 .run ();
  NetworKit::Partition p1 = m1 .getPartition ();
  q = NetworKit::Modularity ();
  quality = q .getQuality (p1, g);
  ss = p1 .subsetSizeMap ();
  iss = ss .begin ();
  mx = 0;
  while (iss != ss .end ()){
    if (mx < iss ->second){
      which = iss ->first;
      mx  = iss ->second;
    }
    iss++;
  }
  fprintf (stderr, "modularity=%lf nver=%li clusters=%li largest=%li\n", quality, p .numberOfElements (), p .numberOfSubsets (), mx);
  v = p1 .getVector ();
  printf ("BEGIN %li\n", v .size ());
  for (int i = 0; i < v .size (); i++){
    printf ("%li;%lf\n", v [i], c .score(i));
  }
  printf ("END\n");
  //std::cout << m.toString() << "\n";
}



