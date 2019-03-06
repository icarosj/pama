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

#define _ll_cas_asm(ptr, oldval, newval) \
    ({ \
      __typeof__(ptr) ___p = (ptr); \
      __typeof__(*___p) ___oldval = (oldval); \
      __typeof__(*___p) ___newval = (newval); \
      register unsigned char ___result; \
      register __typeof__(*___p) ___readval; \
      if (sizeof(*___p) == 4) { \
        __asm__ __volatile__ ("lock; cmpxchgl %3,%1; sete %0" \
                              : "=q"(___result), "=m"(*___p), "=a"(___readval) \
                              : "r"(___newval), "m"(*___p), "2"(___oldval) \
                              : "memory"); \
      } else if (sizeof(*___p) == 8) { \
        __asm__ __volatile__ ("lock; cmpxchgq %3,%1; sete %0" \
                              : "=q"(___result), "=m"(*___p), "=a"(___readval) \
                              : "r"(___newval), "m"(*___p), "2"(___oldval) \
                              : "memory"); \
      } else { \
        abort(); \
      } \
      (___result==1); \
    })


static inline bool _ll_atomic_compare_and_swap(float *dest, float old_val,
		float new_val) {
    return _ll_cas_asm(dest, old_val, new_val);
}



template<typename T>
inline void ATOMIC_ADD(T* target, T value) {

    if (value == 0) return;

    T oldValue, newValue;
    do {
        oldValue = *target;
        newValue = oldValue + value;
    } while (_ll_atomic_compare_and_swap((T*) target, oldValue, newValue) == false);
}


using namespace std;
typedef float ScoreT;
const float kDamp = 0.85;

ScoreT* PageRankPull(const MLGraph &g, int num_iterations) {
  const ScoreT init_score = 1.0f / g.num_nodes();
  const ScoreT base_score = (1.0f - kDamp) / g.num_nodes();
//  pvector<ScoreT> scores(g.num_nodes(), init_score);
  ScoreT* G_pg_rank = new ScoreT[g.num_nodes()];
  ScoreT* G_pg_rank_nxt = new ScoreT[g.num_nodes()];
  ScoreT diff = 0.0 ;
  int32_t cnt = 0 ;
  ScoreT N = (ScoreT)(g.num_nodes()) ;
//  pvector<ScoreT> contrib(g.num_nodes(), 0);
#pragma omp parallel for
    for (NodeID n=0; n < g.num_nodes(); n++)
    {
      G_pg_rank[n] = init_score;
    }


    printf("start iteration\n");

		do
		{
			diff = 0.000000 ;
#pragma omp parallel
			{
				ScoreT diff_prv = 0.0 ;

				diff_prv = 0.000000 ;

#pragma omp for nowait schedule(dynamic,4096)
				for (NodeID t = 0; t < g.num_nodes(); t ++) 
				{
					ScoreT val = 0.0 ;
					ScoreT __S1 = 0.0 ;

					__S1 = 0.000000 ;

          MLGraph::ml_iterator miter = g.iter();
          for(miter.in_begin(t);!miter.is_end();miter.next())
          {
            NodeID w = *miter;
						__S1 = __S1 + G_pg_rank[w] / ((ScoreT)((g.out_degree(w)))) ;
          }


					val = (1 - kDamp) / N + kDamp * __S1 ;
					diff_prv = diff_prv +  std::abs((val - G_pg_rank[t]))  ;
					G_pg_rank_nxt[t] =  val;
				}
				ATOMIC_ADD<ScoreT>(&diff, diff_prv);
			}

#pragma omp parallel for
			for (NodeID i3 = 0; i3 < g.num_nodes(); i3 ++) 
				G_pg_rank[i3] =  G_pg_rank_nxt[i3];

			cnt = cnt + 1 ;
      cout << " " << cnt << "    " << diff<< endl;
    }
		while (cnt < num_iterations);
		

  printf("end iterations\n");



  return G_pg_rank;
}

void PrintTopScores(const MLGraph &g, ScoreT* scores) {
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
