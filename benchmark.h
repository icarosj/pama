// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef BENCHMARK_H_
#define BENCHMARK_H_

#include <algorithm>
#include <cinttypes>
#include <functional>
#include <random>
#include <utility>
#include <vector>

#include "builder.h"
#include "graph.h"
#include "print_util.h"
#include "timer.h"
#include "writer.h"


/*
GAP Benchmark Suite
File:   Benchmark
Author: Scott Beamer

Various helper functions to ease writing of kernels
*/


// Default type signatures for commonly used types
typedef int64_t NodeID;
typedef NodeID WeightT;
typedef NodeWeight<NodeID, WeightT> WNode;

typedef Snapshots<NodeID> MLGraph;
typedef CSRGraph<NodeID> Graph;
typedef CSRGraph<NodeID, WNode> WGraph;

typedef BuilderBase<NodeID, NodeID, NodeID> Builder;
typedef BuilderBase<NodeID, WNode, WeightT> WeightedBuilder;

typedef WriterBase<NodeID, NodeID> Writer;
typedef WriterBase<NodeID, WNode> WeightedWriter;

class pairCMP
{
  public:
  bool operator()(std::pair<NodeID, int> a, std::pair<NodeID, int> b)
  {
    return a.second >= b.second;
  }
};



// Used to pick random non-zero degree starting points for search algorithms
template<typename GraphT_>
class SourcePicker {
 public:
  explicit SourcePicker(const GraphT_ &g, NodeID given_source = -1)
      : given_source(given_source), rng(8), udist(0, g.num_nodes()-1), g_(g) {}

  NodeID PickNext() {
    if (given_source != -1)
      return given_source;
    NodeID source;
    do {
      source = udist(rng);
    } while (g_.out_degree(source) == 0);
    return source;
  }

 private:
  NodeID given_source;
  std::mt19937 rng;
  std::uniform_int_distribution<NodeID> udist;
  const GraphT_ &g_;
};


// Returns k pairs with largest values from list of key-value pairs
template<typename KeyT, typename ValT>
std::vector<std::pair<ValT, KeyT>> TopK(
    const std::vector<std::pair<KeyT, ValT>> &to_sort, size_t k) {
  std::vector<std::pair<ValT, KeyT>> top_k;
  ValT min_so_far = 0;
  for (auto kvp : to_sort) {
    if ((top_k.size() < k) || (kvp.second > min_so_far)) {
      top_k.push_back(std::make_pair(kvp.second, kvp.first));
      std::sort(top_k.begin(), top_k.end(),
                std::greater<std::pair<ValT, KeyT>>());
      if (top_k.size() > k)
        top_k.resize(k);
      min_so_far = top_k.back().first;
    }
  }
  return top_k;
}


// Calls (and times) kernel according to command line arguments
template<typename GraphT_, typename GraphFunc, typename AnalysisFunc>
void BenchmarkKernel(const CLApp &cli, const GraphT_ &g,
                     GraphFunc kernel, AnalysisFunc stats) {
  g.PrintStats();
  double total_seconds = 0;
  Timer trial_timer;
  for (int iter=0; iter < cli.num_trials(); iter++) {
    trial_timer.Start();
    auto result = kernel(g);
    trial_timer.Stop();
    PrintTime("Trial Time", trial_timer.Seconds());
    if (cli.do_analysis() && (iter == (cli.num_trials()-1))) {
      stats(g, result);
    }
    total_seconds += trial_timer.Seconds();
  }
  PrintTime("Average Time", total_seconds / cli.num_trials());
  /*
  for (int i=0;i<10;i++)
  {
    for (auto j : g.out_neigh(i)) {
      std::cout << j << " ";
    }
    std::cout<<std::endl<<std::endl;
  }
*/
}

#endif  // BENCHMARK_H_
