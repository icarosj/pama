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
typedef float ScoreT;
const float kDamp = 0.85;

pvector<ScoreT> PageRankPull(const MLGraph &g, int num_iterations) {
  const ScoreT init_score = 1.0f / g.num_nodes();
  const ScoreT base_score = (1.0f - kDamp) / g.num_nodes();
  pvector<ScoreT> scores(g.num_nodes(), init_score);
  pvector<ScoreT> contrib(g.num_nodes(), 0);
  printf("SIZE contrib : %ld %ld\n", contrib.size(), contrib.storage());
for (int iter=0; iter < num_iterations; iter++) {
    ScoreT error = 0;
#pragma omp parallel for
    for (NodeID n=0; n < g.num_nodes(); n++)
    {
      //the contrib[] is the only one that gets communicated
      contrib[n] = scores[n] / g.out_degree(n);
    }
/*
    for (NodeID u=0; u < g.num_nodes(); u++) {
      MLGraph::ml_iterator miter = g.iter();
      std::cout<<"out node "<<u<<" : ";
      for(miter.out_begin(u);!miter.is_end();miter.next())
      {
        std::cout<<"\t"<<*miter;
      }
      std::cout<<std::endl<<std::endl;
    }


      std::cout<<std::endl<<std::endl;
*/
#pragma omp parallel for schedule(dynamic,4096) reduction(+ : error)
    for (NodeID u=0; u < g.num_nodes(); u++) {
      ScoreT sum = 0;
      //debug_print("\ninD u is %ld\n",g.in_degree(u));

      MLGraph::ml_iterator miter = g.iter();
 //     std::cout<<"in node "<<u<<" : ";
      for(miter.in_begin(u);!miter.is_end();miter.next())
      {
        NodeID v = *miter;
        sum += contrib[v];
//        std::cout<<"\t"<<*miter;
      }
      ScoreT old_score = scores[u];
      scores[u] = base_score + kDamp * sum;
      error += fabs(scores[u] - old_score);
      /*
      for (NodeID v : g.in_neigh(u))
      {
        debug_print("v is %ld\n",v);
        sum += contrib[v];
        #pragma omp critical 
        {
          acc_cnt[v]++;
          total_cnt++;
        }
      }
      */
    }
    cout << " " << iter << "    " << error << endl;
  }

  printf("end iterations\n");



  return scores;
}

void PrintTopScores(const MLGraph &g, const pvector<ScoreT> &scores) {
  vector<pair<NodeID, ScoreT>> score_pairs(g.num_nodes());
  for (NodeID n=0; n < g.num_nodes(); n++) {
    score_pairs[n] = make_pair(n, scores[n]);
  }
  int k = 5;
  vector<pair<ScoreT, NodeID>> top_k = TopK(score_pairs, k);
  k = min(k, static_cast<int>(top_k.size()));
  for (auto kvp : top_k)
    cout << kvp.second << ":" << kvp.first << endl;
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
    return PageRankPull(g, cli.num_iters());
  };
  BenchmarkKernel(cli, g, PRBound, PrintTopScores);
  return 0;
}
