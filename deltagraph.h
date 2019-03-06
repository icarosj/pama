
#ifndef DELTAGRAPH_H_
#define DELTAGRAPH_H_

#include "graph.h"


//deltagraph has : 
//length annotated for vertex tables
//vertex table entries:
//[snapshot# : offset] [delta langth]





//a 64-bit storage. upper bits store the levels
//and lower bits store the indices
class level_n_index
{
  private:
    //signed int to represent null continuation recore (-1)
    int64_t level_n_index_;
  public:
    level_n_index(){}
    level_n_index(int64_t content):level_n_index_(content) {}
    level_n_index(uint32_t level, uint64_t index)
    {
      level_n_index_ = (((uint64_t)level << VT_IDX_INDEX_BITS) | index);
    }

    int64_t content() const
    {
      return level_n_index_;
    }

    vt_elem index()
    {
      return (level_n_index_ & INDEX_MASK) ;
    }
    uint32_t level()
    {
      return (level_n_index_ >> VT_IDX_INDEX_BITS);
    }
    bool cont_record_null()
    {
      if (level_n_index_ == CONT_NULL) return 1;//end
      else return 0;//not end
    }
};



//now it has store level, index and length (degree for optional)
//the idea is to reduce space for the base graph.
//only the index for base. level_n_index and length, degree for deltas 
//64bit for base elem, 128bit for delta elem
class delta_vt_elem
{
  private:

    //I could also have level + offset for cont.record 
    //and save an indirection step
    //but have to pay 128bit for the continuation 
    //(instead of 64)
    //----i changed it to have 128 bit later.
    level_n_index level_n_index_;

    uint32_t length_;
    uint32_t degree_;

  public:
    delta_vt_elem() {}
    delta_vt_elem(int64_t level_n_index):level_n_index_(level_n_index) {}
 //   delta_vt_elem(int64_t level_n_index, uint32_t length):
 //     level_n_index_(level_n_index), length_(length) {}
    delta_vt_elem(int64_t level,uint32_t index, 
        uint32_t length, uint32_t degree):
      level_n_index_(level,index), 
      length_(length),degree_(degree) {}

    vt_elem index()
    {
      return (level_n_index_.index() ) ;
    }
    uint32_t level()
    {
      return (level_n_index_.level());
    }
    bool cont_record_null()
    {
      return level_n_index_.cont_record_null();
    }

    uint32_t length()
    {
      return length_;
    }

    uint32_t degree()
    {
      return degree_;
    }
    
    
};


//graph structure for more than 1 levels
template <class NodeID_, class DestID_ = NodeID_, bool MakeInverse = true>
class DeltaGraph: public CSRGraph<NodeID_> //should i inherit?
{
  friend Snapshots<NodeID_>;


  private:
    delta_vt_elem* out_delta_vt_;
    delta_vt_elem* in_delta_vt_;

    uint64_t out_delta_size_; 
    uint64_t  in_delta_size_; 
    uint64_t out_et_size_; 
    uint64_t  in_et_size_; 

    //inherited from parent
    //DestID_* out_edge_table_;
    //DestID_* in_edge_table_;

    class Neighborhood{
      NodeID_ n_;
      delta_vt_elem vt_entry_;
      DestID_* edge_table_;
      public:
      Neighborhood(){}
      Neighborhood(delta_vt_elem vt,  DestID_* edge_table) : 
       vt_entry_(vt), edge_table_(edge_table){}
      typedef DestID_* iterator;
      iterator begin() 
      { 
        //debug_print("begin n_ %lu, g[n_] %ld \n",
//            vt_entry_.index(), edge_table_[vt_entry_.index()]); 

        return &edge_table_[vt_entry_.index()]; 
      }

      //at the end, it stores delta_vt_elem for the continuation recored
      //make sure delta_vt_elem is twice the size of DestID_ 
      iterator end()   
      { 
        //debug_print("end n_ %lu, g[n_+1] %u \n",
//            vt_entry_.index(), vt_entry_.length());
        return &edge_table_[vt_entry_.index() + vt_entry_.length()]; 
      }
    };






  public:

    DeltaGraph() : CSRGraph<NodeID_>(), out_delta_vt_(nullptr), in_delta_vt_(nullptr)
  {}
 

    //Directed graphs
    DeltaGraph(int64_t num_nodes, uint64_t num_edges,
        delta_vt_elem* out_delta_vt, 
        delta_vt_elem* in_delta_vt,   
        DestID_* out_et,  
        DestID_* in_et,
        uint64_t out_delta_size,  uint64_t in_delta_size,   
        uint64_t out_et_size,     uint64_t in_et_size 
        ) :
      CSRGraph<NodeID_>(num_nodes,num_edges,
          nullptr, out_et, nullptr, in_et), 
      out_delta_vt_(out_delta_vt),    
      in_delta_vt_  (in_delta_vt),
      out_delta_size_(out_delta_size),
      in_delta_size_(in_delta_size),
      out_et_size_(out_et_size),
      in_et_size_(in_et_size)
      {
        std::cout<<"directed delta constructor"<<std::endl;
        //      num_edges_ = (out_vertex_table_[num_nodes_] - out_vertex_table_[0]);
      }


    DeltaGraph(const DeltaGraph& other) : CSRGraph<NodeID_>(other),
    out_delta_vt_(other.out_delta_vt_), in_delta_vt_(other.in_delta_vt_),
      out_delta_size_(other.out_delta_size_),
      in_delta_size_(other.in_delta_size_),
      out_et_size_(other.out_et_size_),
      in_et_size_(other.in_et_size_)
  {
  }

/*
    DeltaGraph(DeltaGraph&& other) : CSRGraph<NodeID_>(other),
    out_delta_vt_(other.out_delta_vt_), in_delta_vt_(other.in_delta_vt_)
  {
    other.in_delta_vt_= nullptr;
    other.out_delta_vt_ = nullptr;
  }
*/
    ~DeltaGraph() {
      //ReleaseResources();
    }

    DeltaGraph& operator=(DeltaGraph&& other) {
      if (this != &other) {
        ReleaseResources();
        this->directed_ = other.directed_;
        this->num_edges_ = other.num_edges_;
        this->num_nodes_ = other.num_nodes_;
        this->in_vertex_table_ = other.in_vertex_table_;
        this->in_edge_table_ =   other.in_edge_table_;
        this->out_vertex_table_ = other.out_vertex_table_;
        this->out_edge_table_ =   other.out_edge_table_;
        other.num_edges_ = -1;
        other.num_nodes_ = -1;
        other.in_vertex_table_= nullptr;
        other.in_edge_table_= nullptr;
        other.out_vertex_table_= nullptr;
        other.out_edge_table_= nullptr;


        out_delta_vt_=other.out_delta_vt_;
        in_delta_vt_=other.in_delta_vt_;
        other.in_delta_vt_= nullptr;
        other.out_delta_vt_ = nullptr;
      }
      return *this;
    }


    uint64_t out_et_size()
    {
      return out_et_size_;
    }
    uint64_t in_et_size()
    {
      return in_et_size_;
    }

    uint64_t out_delta_size()
    {
      return out_delta_size_;
    }
    uint64_t in_delta_size()
    {
      return in_delta_size_;
    }

    int64_t out_degree(NodeID_ v) const {
      return out_delta_vt_[v].degree();
    }


    int64_t in_degree(NodeID_ v) const {
      return in_delta_vt_[v].degree();
    }

    //from the data from the passed vt entry (it could be from other snapshot), follow its neighbors. 
    Neighborhood out_delta_neigh(delta_vt_elem vt) const {
      return Neighborhood(vt, this->out_edge_table_);
    }

    Neighborhood in_delta_neigh(delta_vt_elem vt) const {
      return Neighborhood(vt, this->in_edge_table_);
     }

    delta_vt_elem* out_delta_vt(uint64_t offset)
    {
      return out_delta_vt_+offset;
    }
    delta_vt_elem* in_delta_vt(uint64_t offset)
    {
      return in_delta_vt_+offset;
    }


    void PrintStats() const {
      std::cout << "DeltaGraph has " << this->num_nodes_ << " nodes and "
        << this->num_edges_ << " ";
      if (!this->directed_)
        std::cout << "un";
      std::cout << "directed edges for degree: ";
      std::cout << this->num_edges_/this->num_nodes_ << std::endl;
    }

    void PrintTopology() const {
      for (NodeID_ i=0; i < this->num_nodes_; i++) {
        std::cout << i << ": ";
        for (DestID_ j : out_delta_neigh(i)) {
          std::cout << j << " ";
        }
        std::cout << std::endl;
      }
    }

    size_t size_on()
    {
      return sizeof(out_delta_vt_);
    }
    size_t size_oi()
    {
      return sizeof(this->out_edge_table_);
    }

    DestID_ dbg_i00(int y)
    {
      return this->out_edge_table_[y];
    }
    DestID_ dbg_n0(int x)
    {
      return this->out_delta_vt_[x].index();
    }


  private:
    //is it okay to have the same name?
    //i think so :)
    //it's only 'delta's of the neighborhood. 

    void ReleaseResources() {
      if(out_delta_vt_ != nullptr)
        delete[] out_delta_vt_;
      if(in_delta_vt_ != nullptr)
        delete[] in_delta_vt_;
    }
};



#endif //DELTAGRAPH_H_



