// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef BUILDER_H_
#define BUILDER_H_

#include <algorithm>
#include <cinttypes>
#include <fstream>
#include <functional>
#include <type_traits>
#include <utility>
#include <fcntl.h>
#include <sys/mman.h>

#include "command_line.h"
#include "generator.h"
#include "snapshot.h"
#include "platform_atomics.h"
#include "print_util.h"
#include "pvector.h"
#include "reader.h"
#include "timer.h"


/*
GAP Benchmark Suite
Class:  BuilderBase
Author: Scott Beamer

Given arguements from the command line (cli), returns a built graph
 - MakeGraph() will parse cli and obtain edgelist and call
   MakeGraphFromEL(edgelist) to perform actual graph construction
 - edgelist can be from file (reader) or synthetically generated (generator)
 - Common case: BuilderBase typedef'd (w/ params) to be Builder (benchmark.h)
*/


template <typename NodeID_, typename DestID_ = NodeID_,
          typename WeightT_ = NodeID_, bool invert = true>
class BuilderBase {
  typedef EdgePair<NodeID_, DestID_> Edge;
  typedef pvector<Edge> EdgeList;

  const CLBase &cli_;
  bool symmetrize_;
  bool needs_weights_;
  int64_t num_nodes_ = -1;

 public:
  explicit BuilderBase(const CLBase &cli) : cli_(cli) {
    symmetrize_ = cli_.symmetrize();
    needs_weights_ = !std::is_same<NodeID_, DestID_>::value;
  }

  DestID_ GetSource(EdgePair<NodeID_, NodeID_> e) {
    return e.u;
  }

  DestID_ GetSource(EdgePair<NodeID_, NodeWeight<NodeID_, WeightT_>> e) {
    return NodeWeight<NodeID_, WeightT_>(e.u, e.v.w);
  }

  NodeID_ FindMaxNodeID(const EdgeList &el) {
    NodeID_ max_seen = 0;
    #pragma omp parallel for reduction(max : max_seen)
    for (auto it = el.begin(); it < el.end(); it++) {
      Edge e = *it;
      max_seen = std::max(max_seen, e.u);
      max_seen = std::max(max_seen, (NodeID_) e.v);
    }
    return max_seen;
  }

  pvector<NodeID_> CountDegrees(const EdgeList &el, bool transpose) {
    pvector<NodeID_> degrees(num_nodes_, 0);
    #pragma omp parallel for
    for (auto it = el.begin(); it < el.end(); it++) {
      Edge e = *it;
      if (symmetrize_ || (!symmetrize_ && !transpose))
        fetch_and_add(degrees[e.u], 1);
      if (symmetrize_ || (!symmetrize_ && transpose))
        fetch_and_add(degrees[(NodeID_) e.v], 1);
    }
    return degrees;
  }

  static
  pvector<SGOffset> PrefixSum(const pvector<NodeID_> &degrees) {
    pvector<SGOffset> sums(degrees.size() + 1);
    SGOffset total = 0;
    for (size_t n=0; n < degrees.size(); n++) {
      sums[n] = total;
      total += degrees[n];
    }
    sums[degrees.size()] = total;
    return sums;
  }

  static
  pvector<SGOffset> ParallelPrefixSum(const pvector<NodeID_> &degrees) {
    const size_t block_size = 1<<20;
    const size_t num_blocks = (degrees.size() + block_size - 1) / block_size;
    pvector<SGOffset> local_sums(num_blocks);
    #pragma omp parallel for
    for (size_t block=0; block < num_blocks; block++) {
      SGOffset lsum = 0;
      size_t block_end = std::min((block + 1) * block_size, degrees.size());
      for (size_t i=block * block_size; i < block_end; i++)
        lsum += degrees[i];
      local_sums[block] = lsum;
    }
    pvector<SGOffset> bulk_prefix(num_blocks+1);
    SGOffset total = 0;
    for (size_t block=0; block < num_blocks; block++) {
      bulk_prefix[block] = total;
      total += local_sums[block];
    }
    bulk_prefix[num_blocks] = total;
    pvector<SGOffset> prefix(degrees.size() + 1);
    #pragma omp parallel for
    for (size_t block=0; block < num_blocks; block++) {
      SGOffset local_total = bulk_prefix[block];
      size_t block_end = std::min((block + 1) * block_size, degrees.size());
      for (size_t i=block * block_size; i < block_end; i++) {
        prefix[i] = local_total;
        local_total += degrees[i];
      }
    }
    prefix[degrees.size()] = bulk_prefix[num_blocks];
    return prefix;
  }

  // Removes self-loops and redundant edges
  // Side effect: neighbor IDs will be sorted
 /* 
  void SquishCSR(const CSRGraph<NodeID_, DestID_, invert> &g, bool transpose,
                 DestID_*** sq_index, DestID_** sq_neighs) {
    pvector<NodeID_> diffs(g.num_nodes());
    DestID_ *n_start, *n_end;
    #pragma omp parallel for private(n_start, n_end)
    for (NodeID_ n=0; n < g.num_nodes(); n++) {
      if (transpose) {
        n_start = g.in_neigh(n).begin();
        n_end = g.in_neigh(n).end();
      } else {
        n_start = g.out_neigh(n).begin();
        n_end = g.out_neigh(n).end();
      }
      std::sort(n_start, n_end);
      DestID_ *new_end = std::unique(n_start, n_end);
      new_end = std::remove(n_start, new_end, n);
      diffs[n] = new_end - n_start;
    }
    pvector<SGOffset> sq_offsets = ParallelPrefixSum(diffs);
    E_storage = sizeof(DestID_)*sq_offsets[g.num_nodes()]/2;
    printf("gen neigh length %d size %lld account for %lld\n",sq_offsets[g.num_nodes()],E_storage, E_storage/2); 
    *sq_neighs = new DestID_[sq_offsets[g.num_nodes()]];
    *sq_index = CSRGraph<NodeID_, DestID_>::GenIndex(sq_offsets, *sq_neighs);
    #pragma omp parallel for private(n_start)
    for (NodeID_ n=0; n < g.num_nodes(); n++) {
      if (transpose)
        n_start = g.in_neigh(n).begin();
      else
        n_start = g.out_neigh(n).begin();
      std::copy(n_start, n_start+diffs[n], (*sq_index)[n]);
    }
  }

  CSRGraph<NodeID_, DestID_, invert> SquishGraph(
      const CSRGraph<NodeID_, DestID_, invert> &g) {

    vt_elem* out_vertex_table,  DestID_* out_edge_table ;
    vt_elem* in_vertex_table, DestID_* in_edge_table ;


    SquishCSR(g, false, &out_vertex_table, &out_edge_table);
    if (g.directed()) {
      if (invert)
        SquishCSR(g, true, &in_vertex_table, &in_edge_table);
      return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), out_index,
                                                out_neighs, in_index,
                                                in_neighs);
    } else {
      return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), out_index,
                                                out_neighs);
    }
  }
*/
  /*
  Graph Bulding Steps (for CSR):
    - Read edgelist once to determine vertex degrees (CountDegrees)
    - Determine vertex offsets by a prefix sum (ParallelPrefixSum)
    - Allocate storage and set points according to offsets (GenIndex)
    - Copy edges into storage
  */
  void MakeCSR(const EdgeList &el, bool transpose, 
      vt_elem** p_vertex_table, DestID_** p_edge_table ) {
    pvector<NodeID_> degrees = CountDegrees(el, transpose);
    pvector<SGOffset> offsets = ParallelPrefixSum(degrees);
    std::cout<<"alloc edgelist "<<offsets[num_nodes_]<<std::endl;



/*
     int f = open((cli_.dbpath()+"/base_e_"+std::to_string(transpose)).c_str(),
          O_RDWR);
    if(f < 0 )
    {
      f = open((cli_.dbpath()+"/base_e_"+std::to_string(transpose)).c_str(),
          O_CREAT | O_EXCL | O_RDWR, 0777);
    }
    ASSERT(f, > , 0); 

		int r = fallocate(f, 0, 0, sizeof(DestID_)*offsets[num_nodes_]);
		ASSERT(r, == ,0) ;
    

    *p_edge_table = (DestID_*)mmap(0, sizeof(DestID_)*offsets[num_nodes_], 
        PROT_READ | PROT_WRITE, MAP_SHARED, f, 0);

*/
    *p_edge_table = (DestID_*)mmap_alloc("base_e_"+std::to_string(transpose), 
        cli_.dbpath(), 0, sizeof(DestID_)*offsets[num_nodes_]);



//    *p_edge_table = new DestID_[offsets[num_nodes_]];
    *p_vertex_table = CSRGraph<NodeID_, DestID_>::
      GenVertexTable(offsets, *p_edge_table,transpose, cli_.dbpath());
    #pragma omp parallel for
    for (auto it = el.begin(); it < el.end(); it++) {
      Edge e = *it;

      if (symmetrize_ || (!symmetrize_ && !transpose))
        (*p_edge_table)[fetch_and_add(offsets[e.u], 1)] = e.v;
      if (symmetrize_ || (!symmetrize_ && transpose))
        (*p_edge_table)[fetch_and_add(offsets[static_cast<NodeID_>(e.v)], 1)] =
            GetSource(e);
    }
  }

  CSRGraph<NodeID_, DestID_, invert> MakeGraphFromEL(EdgeList &el) {
    vt_elem* vertex_table, *inv_vertex_table;
    DestID_* edge_table, *inv_edge_table;

    Timer t;
    t.Start();
    if (num_nodes_ == -1)
      num_nodes_ = FindMaxNodeID(el)+1;
    //if (needs_weights_)
    //  Generator<NodeID_, DestID_, WeightT_>::InsertWeights(el);

    std::cout<<"EL size : "<<el.size()<<std::endl;
    MakeCSR(el, false, &vertex_table, &edge_table);
    if (!symmetrize_ && invert)
      MakeCSR(el, true, &inv_vertex_table, &inv_edge_table);
    t.Stop();
    PrintTime("Build Time", t.Seconds());
    if (symmetrize_)
      return CSRGraph<NodeID_, DestID_, invert>(num_nodes_, vertex_table, edge_table);
    else
      return CSRGraph<NodeID_, DestID_, invert>(num_nodes_, vertex_table, edge_table,
          inv_vertex_table, inv_edge_table);
  }



 CSRGraph<NodeID_, DestID_, invert> MakeGraph() {
    CSRGraph<NodeID_, DestID_, invert> g;
    {  // extra scope to trigger earlier deletion of el (save memory)
      EdgeList el;
      if (cli_.filename() != "") {
        Reader<NodeID_, DestID_, WeightT_, invert> r(cli_.filename());
        el = r.ReadFile(needs_weights_);
      } else if (cli_.scale() != -1) {
        assert(0);
      }
      g = MakeGraphFromEL(el);
    }
    return g;
    //return SquishGraph(g);
  }



 Snapshots<NodeID_, DestID_, invert> MakeSnapshots() {
   Snapshots<NodeID_, DestID_, invert> mlg; //multi-level graph
   {  // extra scope to trigger earlier deletion of el (save memory)
     EdgeList el;
     //TODO : start from the DB level
     int level = 0;
     if(cli_.loadDB())
     {
       level = mlg.read_database(cli_.dbpath());
     }
     for (auto filename : cli_.files())
     { 
       Reader<NodeID_, DestID_, WeightT_, invert> r(filename);
       el = r.ReadFile(needs_weights_);
       if(level == 0)
       {
         //g = new CSRGraph<NodeID_, DestID_, invert>(MakeGraphFromEL(el)); ;
         // mlg.add_graph(g);
         mlg.create_base_from_EL(el, cli_.dbpath());
         //         mlg.add_graph(CSRGraph<NodeID_, DestID_, invert>(MakeGraphFromEL(el)));
       }
       else 
       {
         std::cout<<"make delta"<<std::endl;
         mlg.add_delta_from_EL(el, cli_.dbpath());
         std::cout<<"make delta done"<<std::endl;
       }
       printf("acccnt is sized as %ld\n", mlg.num_nodes());
       level++;
     } //graphs for snapshots have to be input from a file. no generation supported.

     mlg.write_meta_file(cli_.dbpath());
     mlg.write_indirection_file(cli_.dbpath());
     std::cout<<level<<" levels created "<<std::endl;

   }
   return mlg;
   //return SquishGraph(g);
 }




 // Relabels (and rebuilds) graph by order of decreasing degree
 /*
    static
    CSRGraph<NodeID_, DestID_, invert> RelabelByDegree(
    const CSRGraph<NodeID_, DestID_, invert> &g) {
    if (g.directed()) {
    std::cout << "Cannot relabel directed graph" << std::endl;
    std::exit(-11);
    }
    Timer t;
    t.Start();
    typedef std::pair<int64_t, NodeID_> degree_node_p;
    pvector<degree_node_p> degree_id_pairs(g.num_nodes());
    #pragma omp parallel for
    for (NodeID_ n=0; n < g.num_nodes(); n++)
      degree_id_pairs[n] = std::make_pair(g.out_degree(n), n);
    std::sort(degree_id_pairs.begin(), degree_id_pairs.end(),
              std::greater<degree_node_p>());
    pvector<NodeID_> degrees(g.num_nodes());
    pvector<NodeID_> new_ids(g.num_nodes());
    #pragma omp parallel for
    for (NodeID_ n=0; n < g.num_nodes(); n++) {
      degrees[n] = degree_id_pairs[n].first;
      new_ids[degree_id_pairs[n].second] = n;
    }
    pvector<SGOffset> offsets = ParallelPrefixSum(degrees);
    DestID_* neighs = new DestID_[offsets[g.num_nodes()]];
    DestID_** index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, neighs);
    #pragma omp parallel for
    for (NodeID_ u=0; u < g.num_nodes(); u++) {
      for (NodeID_ v : g.out_neigh(u))
        neighs[offsets[new_ids[u]]++] = new_ids[v];
      std::sort(index[new_ids[u]], index[new_ids[u]+1]);
    }
    t.Stop();
    PrintTime("Relabel", t.Seconds());
    return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), index, neighs);
  }
  */
};

#endif  // BUILDER_H_
