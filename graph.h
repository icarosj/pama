// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef GRAPH_H_
#define GRAPH_H_

#include <cinttypes>
#include <iostream>
#include <type_traits>
#include <set>
#include <assert.h>

#include "pvector.h"
#include "error.h"
#ifdef DEBUG 
#define debug_print(fmt, ...) \
              do { fprintf(stdout, fmt, __VA_ARGS__); } while (0)
#else
#define debug_print(fmt, ...) ;
#endif

/*#define DEBUG 0
#define debug_print(fmt, ...) \
              do { if (DEBUG) fprintf(stdout, fmt, __VA_ARGS__); } while (0)*/
/*
GAP Benchmark Suite
Class:  CSRGraph
Author: Scott Beamer

Simple container for graph in CSR format
 - Intended to be constructed by a Builder
 - To make weighted, set DestID_ template type to NodeWeight
 - MakeInverse parameter controls whether graph stores its inverse
*/


//level_n_index structure
#define VT_IDX_BITS            64
#define VT_IDX_LEVEL_BITS      10
#define VT_IDX_INDEX_BITS      (VT_IDX_BITS-VT_IDX_LEVEL_BITS) 
#define INDEX_MASK             ((1ull<<VT_IDX_INDEX_BITS)-1)

#define CONT_NULL -1ll


//indirection structure
#define IND_VERTICE_BITS       9 //this make a page 2^12B for base, 2^13B for delta
//#define IND_VERTICE_BITS       6 //this make a page 2^12B for base, 2^13B for delta
#define IND_VERTICE_MASK       ((1ull << IND_VERTICE_BITS)-1) 
#define IND_VERTICES_PER_PAGE  (1ull<<IND_VERTICE_BITS)  
#define BLOCK_SIZE             (((IND_VERTICES_PER_PAGE)) * (16))

#define IND_V_OFFSET(x)        ((x) & (IND_VERTICE_MASK)) 
#define IND_V_PG(x)            ((x) >> (IND_VERTICE_BITS))         
#define IND_V_PG_FILTER(x)     (((x) >> (IND_VERTICE_BITS)) <<(IND_VERTICE_BITS))        




//Make it page-aligned!!!
void* mmap_alloc(std::string filename, std::string dbpath, uint32_t level,  uint64_t size)
{
  std::string filepath = dbpath + "/" + filename + std::to_string(level);
  int f = open(filepath.c_str(),   O_RDWR);
  if(f < 0 )
  {
    f = open(filepath.c_str(),
        O_CREAT | O_EXCL | O_RDWR, 0777);
  }
  ASSERT(f, > , 0); 

  size_t s = size;
  std::cout<<"alloc size "<<s<<std::endl;
  if ((s & (BLOCK_SIZE-1)) != 0) {
    size_t x = BLOCK_SIZE - (s & (BLOCK_SIZE-1));
    s += x;
    std::cout<<"alloc size changed by "<<x<<" and now "<<s<<std::endl;
  }

  int r = fallocate(f, 0, 0, s);
//  int r = ftruncate(f, size);
  if(r!=0)
    std::cout<<"f "<<f<<"size "<<size<<" errno "<<errno<<std::endl;
  ASSERT(r, == ,0) ;


  return  mmap(0, s, PROT_READ | PROT_WRITE, MAP_SHARED, f, 0);

}



void* mmap_open(std::string filename, std::string dbpath, uint32_t level,  uint64_t size)
{
  std::string filepath = dbpath + "/" + filename + std::to_string(level);
  int f = open(filepath.c_str(),   O_RDWR);
  ASSERT(f, > , 0); 


  size_t s = size;
  std::cout<<"open size "<<s<<std::endl;
  if ((s & (BLOCK_SIZE-1)) != 0) {
    size_t x = BLOCK_SIZE - (s & (BLOCK_SIZE-1));
    s += x;
    std::cout<<"open size changed by "<<x<<" and now "<<s<<std::endl;
  }

  //return malloc(size);
  std::cout<<"mmap open size "<<size<<std::endl;
  return  mmap(0, s, PROT_READ | PROT_WRITE, MAP_SHARED, f, 0);

}





// Used to hold node & weight, with another node it makes a weighted edge
template <typename NodeID_, typename WeightT_>
struct NodeWeight {
  NodeID_ v;
  WeightT_ w;
  NodeWeight() {}
  NodeWeight(NodeID_ v) : v(v), w(1) {}
  NodeWeight(NodeID_ v, WeightT_ w) : v(v), w(w) {}

  bool operator< (const NodeWeight& rhs) const {
    return v == rhs.v ? w < rhs.w : v < rhs.v;
  }

  // doesn't check WeightT_s, needed to remove duplicate edges
  bool operator== (const NodeWeight& rhs) const {
    return v == rhs.v;
  }

  // doesn't check WeightT_s, needed to remove self edges
  bool operator== (const NodeID_& rhs) const {
    return v == rhs;
  }

  operator NodeID_() {
    return v;
  }
};

template <typename NodeID_, typename WeightT_>
std::ostream& operator<<(std::ostream& os,
                         const NodeWeight<NodeID_, WeightT_>& nw) {
  os << nw.v << " " << nw.w;
  return os;
}

template <typename NodeID_, typename WeightT_>
std::istream& operator>>(std::istream& is, NodeWeight<NodeID_, WeightT_>& nw) {
  is >> nw.v >> nw.w;
  return is;
}



// Syntatic sugar for an edge
template <typename SrcT, typename DstT = SrcT>
struct EdgePair {
  SrcT u;
  DstT v;

  EdgePair() {}

  EdgePair(SrcT u, DstT v) : u(u), v(v) {}
};

// SG = serialized graph, these types are for writing graph to file
typedef int32_t SGID;
typedef EdgePair<SGID> SGEdge;
typedef int64_t SGOffset;

typedef uint64_t vt_elem;


template <class NodeID_, class DestID_ = NodeID_, bool MakeInverse = true>
class Snapshots;
//class Snapshots::ml_iterator;

template <class NodeID_, class DestID_ = NodeID_, bool MakeInverse = true>
class CSRGraph {
  friend class Snapshots<NodeID_>;
  // Used to access neighbors of vertex, basically sugar for iterators
  protected:
    bool directed_;
    int64_t num_nodes_;
    int64_t num_edges_;

    //points to middle of vertex_table
    //  vt_elem* indirection[VERTEX_PER_PAGE];

    //stores indices of edge_table
    vt_elem* out_vertex_table_;
    vt_elem* in_vertex_table_;

    //stores destination node ids
    DestID_* out_edge_table_;
    DestID_* in_edge_table_;



    class Neighborhood {
      NodeID_ n_;
      vt_elem* vertex_table;
      DestID_* edge_table;
      public:
      Neighborhood(){}
      Neighborhood(NodeID_ n, vt_elem* vertex_table, DestID_* edge_table) : n_(n), vertex_table(vertex_table), edge_table(edge_table){}
      typedef DestID_* iterator;
      iterator begin() 
      { 
        //debug_print("begin n_ %ld, g[n_] %ld * %ld\n",
//            n_, vertex_table[n_], edge_table[vertex_table[n_]]); 
        return &edge_table[vertex_table[n_]]; 
      }

      iterator end()   
      { 
       // debug_print("end n_ %ld, g[n_+1] %ld * %ld\n",
//            n_+1, vertex_table[n_+1], edge_table[vertex_table[n_+1]-1]);
        return &edge_table[vertex_table[n_+1]]; 
      }
    };

    void ReleaseResources() {
      if (out_vertex_table_ != nullptr)
        delete[] out_vertex_table_;
      if (out_edge_table_ != nullptr)
        delete[] out_edge_table_;
      if(directed_)
      {
        if (in_vertex_table_ != nullptr)
          delete[] in_vertex_table_;
        if (in_edge_table_ != nullptr)
          delete[] in_edge_table_;
      }
    }


  public:
    CSRGraph() : directed_(false), num_nodes_(-1), num_edges_(-1),
    out_vertex_table_(nullptr),in_vertex_table_(nullptr), 
    out_edge_table_(nullptr) , in_edge_table_(nullptr){}



    //Undirected edges
    CSRGraph(int64_t num_nodes, vt_elem* vertex_table, DestID_* edge_table ) :
      directed_(false), num_nodes_(num_nodes),
      out_vertex_table_(vertex_table), in_vertex_table_(vertex_table), 
      out_edge_table_(edge_table),     in_edge_table_(edge_table) {
        std::cout<<"undirected constructor"<<std::endl;
        num_edges_ = (vertex_table[num_nodes_] - vertex_table[0])/2;
      }

    //Directed graphs
    CSRGraph(int64_t num_nodes, 
        vt_elem* out_vertex_table,  DestID_* out_edge_table ,
        vt_elem* in_vertex_table, DestID_* in_edge_table ) :
      directed_(false), num_nodes_(num_nodes),
      out_vertex_table_(out_vertex_table), in_vertex_table_(in_vertex_table), 
      out_edge_table_(out_edge_table),     in_edge_table_(in_edge_table) {
        std::cout<<"directed constructor"<<std::endl;
        num_edges_ = (out_vertex_table_[num_nodes_] - out_vertex_table_[0]);
      }



    //Directed graphs, explicit num_edges
    CSRGraph(int64_t num_nodes, int64_t num_edges, 
        vt_elem* out_vertex_table,  DestID_* out_edge_table ,
        vt_elem* in_vertex_table, DestID_* in_edge_table ) :
      directed_(false), num_nodes_(num_nodes),
      out_vertex_table_(out_vertex_table), in_vertex_table_(in_vertex_table), 
      out_edge_table_(out_edge_table),     in_edge_table_(in_edge_table) {
        std::cout<<"directed constructor"<<std::endl;
        num_edges_ = num_edges;
      }


    CSRGraph(const CSRGraph& other) : directed_(other.directed_),
    num_nodes_(other.num_nodes_), num_edges_(other.num_edges_),
    out_vertex_table_(other.out_vertex_table_),in_vertex_table_(other.in_vertex_table_),
    out_edge_table_(other.out_edge_table_) ,   in_edge_table_(other.in_edge_table_) {
    }

/*
    CSRGraph(CSRGraph&& other) : directed_(other.directed_),
    num_nodes_(other.num_nodes_), num_edges_(other.num_edges_),
    in_vertex_table_(other.in_vertex_table_),out_vertex_table_(other.out_vertex_table_),
    in_edge_table_(other.in_edge_table_) ,   out_edge_table_(other.out_edge_table_) {
      other.num_edges_ = -1;
      other.num_nodes_ = -1;
      other.in_vertex_table_= nullptr;
      other.out_vertex_table_= nullptr;
      other.in_edge_table_= nullptr;
      other.out_edge_table_= nullptr;
    }
*/
    ~CSRGraph() {
      //ReleaseResources();
    }

    CSRGraph& operator=(CSRGraph&& other) {
      if (this != &other) {
        ReleaseResources();
        directed_ = other.directed_;
        num_edges_ = other.num_edges_;
        num_nodes_ = other.num_nodes_;
        in_vertex_table_ = other.in_vertex_table_;
        in_edge_table_ =   other.in_edge_table_;
        out_vertex_table_ = other.out_vertex_table_;
        out_edge_table_ =   other.out_edge_table_;
        other.num_edges_ = -1;
        other.num_nodes_ = -1;
        other.in_vertex_table_= nullptr;
        other.in_edge_table_= nullptr;
        other.out_vertex_table_= nullptr;
        other.out_edge_table_= nullptr;
      }
      return *this;
    }

    virtual bool directed() const {
      return directed_;
    }

    virtual int64_t num_nodes() const {
      return num_nodes_;
    }

    virtual int64_t num_edges() const {
      return num_edges_;
    }

    virtual int64_t num_edges_directed() const {
      return directed_ ? num_edges_ : 2*num_edges_;
    }

    virtual int64_t out_degree(NodeID_ v) const {
      return out_vertex_table_[v+1] - out_vertex_table_[v];
    }


    virtual int64_t in_degree(NodeID_ v) const {
      return in_vertex_table_[v+1] - in_vertex_table_[v];
    }

    /*
       int64_t in_degree(NodeID_ v) const {
       static_assert(MakeInverse, "Graph inversion disabled but reading inverse");
       debug_print("getting degree : u %d  in_indexv+1 %d  in_indexv %d res %d\n", v, in_index_[v+1], in_index_[v],in_index_[v+1] - in_index_[v]);
       return in_index_[v+1] - in_index_[v];
       }*/

    Neighborhood out_neigh(NodeID_ n) const {
      return Neighborhood(n, out_vertex_table_, out_edge_table_);
    }

    Neighborhood in_neigh(NodeID_ n) const {
      return Neighborhood(n, in_vertex_table_, in_edge_table_);
      /*debug_print("looking for in neighbor of %d\n", n);
        static_assert(MakeInverse, "Graph inversion disabled but reading inverse");
        return Neighborhood(n, in_index_);*/
    }

    vt_elem* out_vt(uint64_t offset)
    {
      return out_vertex_table_+offset; 
    }

    vt_elem* in_vt(uint64_t offset)
    {
      return in_vertex_table_+offset; 
    }

    virtual void PrintStats() const {
      std::cout << "Graph has " << num_nodes_ << " nodes and "
        << num_edges_ << " ";
      if (!directed_)
        std::cout << "un";
      std::cout << "directed edges for degree: ";
      std::cout << num_edges_/num_nodes_ << std::endl;
    }

    virtual void PrintTopology() const {
      for (NodeID_ i=0; i < num_nodes_; i++) {
        std::cout << i << ": ";
        for (DestID_ j : out_neigh(i)) {
          std::cout << j << " ";
        }
        std::cout << std::endl;
      }
    }


    static vt_elem* GenVertexTable(const pvector<SGOffset> &offsets, 
        DestID_* edge_table, bool transpose, std::string dbpath) {
      NodeID_ length = offsets.size();

      /*
         int f = open((dbpath+"/base_v_"+std::to_string(transpose)).c_str(),
         O_RDWR);
         if(f < 0 )
         {
         f = open((dbpath+"/base_v_"+std::to_string(transpose)).c_str(),
         O_CREAT | O_EXCL | O_RDWR, 0777);
         }
         ASSERT(f, > , 0); 

         int r = fallocate(f, 0, 0, sizeof(DestID_)*length);
         ASSERT(r, == ,0) ;


         vt_elem* vertex_table = (vt_elem*)mmap(0, sizeof(vt_elem)*length, 
         PROT_READ | PROT_WRITE, MAP_SHARED, f, 0);

*/
      vt_elem* vertex_table = (vt_elem*)mmap_alloc("base_v_"+std::to_string(transpose), dbpath, 0, sizeof(vt_elem)*length);


#pragma omp parallel for
      for (NodeID_ n=0; n < length; n++)
          vertex_table[n] = offsets[n];
      return vertex_table;
    }

    pvector<SGOffset> VertexOffsets(bool in_graph = false) const {
      assert(0);//not used. not thati know of
      pvector<SGOffset> offsets(num_nodes_+1);
      for (NodeID_ n=0; n < num_nodes_+1; n++)
        if (in_graph)
          offsets[n] = in_vertex_table_[n] - in_vertex_table_[0];
        else
          offsets[n] = out_vertex_table_[n] - out_vertex_table_[0];
      return offsets;
    }

    virtual size_t size_on()
    {
      return sizeof(out_vertex_table_);
    }
    virtual size_t size_oi()
    {
      return sizeof(out_edge_table_);
    }

    virtual DestID_ dbg_i00(int x, int y)
    {
      return out_edge_table_[out_vertex_table_[x]+y];
    }
    virtual DestID_ dbg_n0(int x)
    {
      return out_vertex_table_[x];
    }

    //  int* acc_cnt;
    // private:
  protected:
};

#endif  // GRAPH_H_
