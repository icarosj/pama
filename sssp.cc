// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <cinttypes>
#include <limits>
#include <iostream>
#include <vector>

#include "benchmark.h"
#include "bucket.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "platform_atomics.h"
#include "pvector.h"
#include "timer.h"

using namespace std;

const WeightT kDistInf = numeric_limits<WeightT>::max()/2;


// reduces barriers down to 2
pvector<WeightT> DeltaStep(const WGraph &g, NodeID source, WeightT delta) {
  Timer t;
  //pvector seems to be a priority queue 
  //implemented with a vector
  //no, just a vector
  pvector<WeightT> dist(g.num_nodes(), kDistInf);
  dist[source] = 0;
  // two element arrays for double buffering curr=iter%2, next=(iter+1)%2
  //THIS seems to be a PQ
  Bucket<NodeID> shared_bins[2];
  size_t shared_indexes[2] = {0, kDistInf};
  shared_bins[0].push_back(source);
  int64_t num_checks = 0;
  t.Start();
  #pragma omp parallel reduction(+ : num_checks)
  {
    vector<vector<NodeID>> local_bins;
    size_t iter = 0;
    while (shared_indexes[iter%2] != kDistInf) {
      size_t curr_bin_index = shared_indexes[iter%2];
      #pragma omp for nowait schedule(dynamic, 64)
      //iterate over the checked bucket?
      for (auto it = shared_bins[iter%2].begin();
           it < shared_bins[iter%2].end(); ++it) {
        NodeID u = *it;
        //if not converged?
        //the default value of delta is 1
        if (dist[u] >= delta*static_cast<WeightT>(curr_bin_index)) {
          num_checks += g.out_degree(u);
          for (WNode wn : g.out_neigh(u)) {
            WeightT old_dist = dist[wn.v];
            WeightT new_dist = dist[u] + wn.w;
            if (new_dist < old_dist) {
              bool changed_dist = true;
              while (!compare_and_swap(dist[wn.v], old_dist, new_dist)) {
                old_dist = dist[wn.v];
                if (old_dist <= new_dist) {
                  //a smaller dist was set by another thr
                  changed_dist = false;
                  break;
                }
              }
              if (changed_dist) {
                //by doing this, we choose smallest dist node
                size_t dest_bin = new_dist/delta;
                if (dest_bin >= local_bins.size()) {
                  local_bins.resize(dest_bin+1);
                }
                //put same-distanced ones in the same bin
                local_bins[dest_bin].push_back(wn.v);
              }
            }
          }
        }
      }
      for (size_t i=curr_bin_index; i < local_bins.size(); i++) {
        if (!local_bins[i].empty()) {
          #pragma omp critical
            shared_indexes[(iter+1)%2] = min(shared_indexes[(iter+1)%2], i);
            break;
        }
      }
      #pragma omp barrier
      #pragma omp single nowait
      {
        t.Stop();
        PrintStep(curr_bin_index, t.Millisecs(), shared_bins[iter%2].size());
        t.Start();
        shared_bins[iter%2].clear();
        shared_indexes[iter%2] = kDistInf;
      }
      //does it mean "anything exists" in the index?
      //yeah it seems so
      if (shared_indexes[(iter+1)%2] < local_bins.size())
        //is this swap thread-safe??
        //yeah it is..
        shared_bins[(iter+1)%2].swap_vector_in(
            local_bins[shared_indexes[(iter+1)%2]]);
      iter++;
      #pragma omp barrier
    }
  }
  return dist;
}


void PrintSSSPStats(const WGraph &g, const pvector<WeightT> &dist) {
  auto NotInf = [](WeightT d) { return d != kDistInf; };
  int64_t num_reached = count_if(dist.begin(), dist.end(), NotInf);
  cout << "SSSP Tree reaches " << num_reached << " nodes" << endl;
}


int main(int argc, char* argv[]) {
  CLDelta cli(argc, argv, "single-source shortest-path");
  if (!cli.ParseArgs())
    return -1;
  WeightedBuilder b(cli);
  WGraph g = b.MakeGraph();
  SourcePicker<WGraph> sp(g, cli.start_vertex());
  auto SSSPBound = [&sp, &cli] (const WGraph &g) {
    return DeltaStep(g, sp.PickNext(), cli.delta());
  };
  BenchmarkKernel(cli, g, SSSPBound, PrintSSSPStats);
  return 0;
}
