// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <iostream>
#include <vector>

#include "benchmark.h"
#include "bitmap.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "platform_atomics.h"
#include "pvector.h"
#include "sliding_queue.h"
#include "timer.h"

using namespace std;


/*

Using optimization of precomputing degrees in bulk in beginning and storing
them in parent array as negative numbers. Thus the encoding of parent is:
parent[x] < 0 implies it is -out_degree(x)
parent[x] >= 0 implies it is parent(x)

*/


vector<int> acc_cnt; //(g.num_nodes(), 0);
int total_cnt;

int64_t BUStep(const Graph &g, pvector<NodeID> &parent, Bitmap &front,
               Bitmap &next) {
  int64_t awake_count = 0;
  next.reset();
  #pragma omp parallel for reduction(+ : awake_count) schedule(dynamic, 1024)
  for (NodeID u=0; u < g.num_nodes(); u++) {
    if (parent[u] < 0) {
      for (NodeID v : g.in_neigh(u)) {
        if (front.get_bit(v)) {
          parent[u] = v;
          awake_count++;
          next.set_bit(u);
          break;
        }
      }
    }
  }
  return awake_count;
}


int64_t TDStep(const Graph &g, pvector<NodeID> &parent,
               SlidingQueue<NodeID> &queue) {
               //SlidingQueue<NodeID> &queue, vector<int>& acc_cnt) {
  int64_t scout_count = 0;
  #pragma omp parallel
  {
    QueueBuffer<NodeID> lqueue(queue);
    #pragma omp for reduction(+ : scout_count)
    for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
      NodeID u = *q_iter;
      for (NodeID v : g.out_neigh(u)) {
        NodeID curr_val = parent[v];
        #pragma omp critical 
        {
          acc_cnt[v]++;
          total_cnt++;
        }
        if (curr_val < 0) {
          if (compare_and_swap(parent[v], curr_val, u)) {
            lqueue.push_back(v);
            scout_count += -curr_val;
          }
        }
      }
    }
    lqueue.flush();
  }
  return scout_count;
}


void QueueToBitmap(const SlidingQueue<NodeID> &queue, Bitmap &bm) {
  #pragma omp parallel for
  for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
    NodeID u = *q_iter;
    bm.set_bit_atomic(u);
  }
}

void BitmapToQueue(const Graph &g, const Bitmap &bm,
                   SlidingQueue<NodeID> &queue) {
  #pragma omp parallel
  {
    QueueBuffer<NodeID> lqueue(queue);
    #pragma omp for
    for (NodeID n=0; n < g.num_nodes(); n++)
      if (bm.get_bit(n))
        lqueue.push_back(n);
    lqueue.flush();
  }
  queue.slide_window();
}

pvector<NodeID> InitParent(const Graph &g) {
  pvector<NodeID> parent(g.num_nodes());
  #pragma omp parallel for
  for (NodeID n=0; n < g.num_nodes(); n++)
    parent[n] = g.out_degree(n) != 0 ? -g.out_degree(n) : -1;
  return parent;
}


pvector<NodeID> DOBFS(const Graph &g, NodeID source, int alpha = 26,
                      int beta = 72) {
  cout << "source: " << source << endl;
  Timer t;
  t.Start();
  pvector<NodeID> parent = InitParent(g);
  printf("SIZE parent: %d %d\n", parent.size(), parent.storage());
  data_storage = parent.storage();
  t.Stop();
  PrintStep("i", t.Seconds());
  parent[source] = source;
  SlidingQueue<NodeID> queue(g.num_nodes());
  queue.push_back(source);
  queue.slide_window();
  Bitmap curr(g.num_nodes());
  curr.reset();
  Bitmap front(g.num_nodes());
  front.reset();
  int64_t edges_to_check = g.num_edges_directed();
  int64_t scout_count = g.out_degree(source);


  total_cnt=0;
  acc_cnt = vector<int>(g.num_nodes(), 0);

  while (!queue.empty()) {
    //if (scout_count > edges_to_check / alpha) {
    if (0) {
      int64_t awake_count;
      TIME_OP(t, QueueToBitmap(queue, front));
      PrintStep("e", t.Seconds());
      queue.slide_window();
      do {
        t.Start();
        awake_count = BUStep(g, parent, front, curr);
        front.swap(curr);
        t.Stop();
        PrintStep("bu", t.Seconds(), awake_count);
      } while (awake_count > g.num_nodes() / beta);
      TIME_OP(t, BitmapToQueue(g, front, queue));
      PrintStep("c", t.Seconds());
      scout_count = 1;
    } else {
      t.Start();
      edges_to_check -= scout_count;
      scout_count = TDStep(g, parent, queue);
      queue.slide_window();
      t.Stop();
      PrintStep("td", t.Seconds(), queue.size());
    }
  }

  printf("end iterations\n");
  std::set<pair<NodeID, int>, pairCMP> accesses;
  int acc=0;
  for (NodeID u=0; u < g.num_nodes(); u++) 
  {
    pair<NodeID, int> tmp(u, acc_cnt[u]);
    accesses.insert(tmp);
//    printf("adding (%d, %d)\n", tmp.first, tmp.second);
    acc+=acc_cnt[u];
  }
  int total_graph_size = E_storage+V_storage+data_storage;
printf("  total graph size :%.6fMB vertex %2.2f%% edge %2.2f%% data %2.2f%%\n"
    , (float)total_graph_size/1000000, (float)V_storage/total_graph_size*100,  (float)E_storage/total_graph_size*100,  (float)data_storage/total_graph_size*100);     

  printf("total cnt : %d , %d\n", total_cnt, acc);
  acc=0;
  printf("acc,NodeID, deg, acc/deg, access, self\%, node ratio\%, accumulate\%, ref to whole\n");
  int cnt=0;
  for(auto it : accesses)
  {
    cnt++;
    acc+=it.second;
    printf("acc,%d,  %d, %d, %2.3f, %2.3f\%, %2.3f\%, %2.3f\%, %2.3f\%\n", it.first, it.second, g.out_degree(it.first) , (float)(it.second)/(g.out_degree(it.first)), (float)(it.second)/total_cnt*100, (float)(cnt)/g.num_nodes()*100,(float)(acc)/total_cnt*100,(float)(acc)/total_cnt*data_storage/total_graph_size*100);
  }

  printf("acc again : %d\n", acc);








  return parent;
}


void PrintBFSStats(const Graph &g, const pvector<NodeID> &bfs_tree) {
  int64_t tree_size = 0;
  int64_t n_edges = 0;
  for (NodeID n=0; n < g.num_nodes(); n++) {
    if (bfs_tree[n] >= 0) {
      n_edges += g.out_degree(n);
      tree_size++;
    }
  }
  cout << "BFS Tree has " << tree_size << " nodes and ";
  cout << n_edges << " edges" << endl;
}


int main(int argc, char* argv[]) {
  CLApp cli(argc, argv, "breadth-first search");
  if (!cli.ParseArgs())
    return -1;
  Builder b(cli);
  Graph g = b.MakeGraph();
  SourcePicker<Graph> sp(g, cli.start_vertex());
  auto BFSBound = [&sp] (const Graph &g) { return DOBFS(g, sp.PickNext()); };
  BenchmarkKernel(cli, g, BFSBound, PrintBFSStats);
  return 0;
}
