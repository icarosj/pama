// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <algorithm>
#include <iostream>
#include <vector>
#include <set>

#include "benchmark.h"
#include "builder.h"
#include "command_line.h"
#include "snapshot.h"
#include "pvector.h"

using namespace std;

int PrintTopology(const MLGraph &g, int num_iterations) {


  std::cout<<"graph has "<<g.num_nodes()<<" nodes and "<<g.num_edges()<<" edges"<<std::endl;
    for (NodeID u=0; u < g.num_nodes(); u++) {
      MLGraph::ml_iterator miter = g.iter();
      std::cout<<"out node("<<g.out_degree(u)<<") "<<u<<" : ";
      for(miter.out_begin(u);!miter.is_end();miter.next())
      {
        std::cout<<"\t"<<*miter;
      }
      std::cout<<std::endl<<std::endl;
    }


    std::cout<<std::endl<<std::endl;

    for (NodeID u=0; u < g.num_nodes(); u++) {

      MLGraph::ml_iterator miter = g.iter();
      std::cout<<"in node("<<g.in_degree(u)<<") "<<u<<" : ";
      for(miter.in_begin(u);!miter.is_end();miter.next())
      {
        std::cout<<"\t"<<*miter;
      }
      std::cout<<std::endl<<std::endl;
    }

  printf("end iterations\n");
  return 0;

}

void PrintTopScores(const MLGraph &g, int nil) {

}

int main(int argc, char* argv[]) {

  install_backtrace_handler();
  CLIterApp cli(argc, argv, "pagerank", 20);
  if (!cli.ParseArgs())
    return -1;
  Builder b(cli);

  MLGraph g = b.MakeSnapshots();
//  printf("|V|:%d  |E|:%d  size index %d  size neigh %d  index[0][0] %d neigh[0] %d %d\n", g.num_nodes(), g.num_edges(),g.size_oi(), g.size_on(), g.dbg_i00(0,9),g.dbg_i00(1,2), g.dbg_n0(9)); 
  auto PRBound = [&cli] (const MLGraph &g) {
    return PrintTopology(g, cli.num_iters());
  };
  BenchmarkKernel(cli, g, PRBound, PrintTopScores);
  return 0;
}
